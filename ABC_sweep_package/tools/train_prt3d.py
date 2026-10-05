#!/usr/bin/env python3
"""Train, validate and test a PRT-DeepONet on the ABC dataset, and give a verdict.

    python3 tools/train_prt3d.py --dataset dataset/dataset.h5 --species C

Standalone: NumPy, h5py and PyTorch, nothing else. It reads the same .h5 the
PRT-DeepONet-GUI reads, honours the train/validation/test split recorded in it,
and finishes with a number that says whether the whole chain, from CompLaB3D
through the collector to the network, produced something a model can learn.

    python3 3D/model/train.py --data dataset/dataset.h5 --species C --distance gdf

is the other way to do this, and the better one once you are past the first
dataset: it is the code the published results came from, it has the three feature
switches, and it is what the GUI drives. Use it. This file exists so the package
can answer "does the whole chain work" without needing that repository present,
and so the answer is a number rather than an impression.

THE MODEL
---------
The published 3D architecture, at this dataset's resolution:

    branch 1   the pore space, (1, 32, 32, 32)
               five conv blocks, 16 32 64 128 256 channels, SiLU, average pooling
               32 halved five times is 1, so the flatten is 1*1*1*256 = 256
    branch 2   the parameters, (Pe, Da), three linear layers of 128
    trunk      (x, y, z, t_norm, gdf) at sampled pore voxels, eight linear of 128
    output     the dot product of branch and trunk, one scalar per point

One model per species, as in the 2D release: the network has a single output
field. --species picks which.

WHY THE TRUNK IS SAMPLED
------------------------
Evaluating it at all 32768 voxels of every sample is most of the cost and almost
none of the information: neighbouring voxels of a smooth field are nearly the
same question asked twice. Drawing n_points random PORE voxels per sample per
step costs a fraction of that and converges to the same thing. Solid voxels are
excluded because there is nothing there to predict.

WHAT THE VERDICT MEANS
----------------------
The test score is a relative L2 error over held-out GEOMETRIES, against two
baselines:

    mean        predict the dataset mean everywhere. A model that does not beat
                this has learned nothing at all.
    per-sample  predict each sample's own mean. This one is harder: beating it
                means the model has learned the SHAPE of the field and not just
                how much of it there is.

A short run on a small dataset will not give a good number, and it is not meant
to. What it is meant to give is a number that is clearly better than both
baselines, which is what says the inputs carry the information the outputs need,
which is what says the CompLaB3D sweep did its job.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

PORE = 2


# ------------------------------------------------------------------- data ---
class ABCDataset(object):
    """Loaded into memory: the whole sweep is about 1 GB as float32 and being
    able to index it freely is worth more here than streaming it."""

    def __init__(self, path, species, split=None, n_points=2048, seed=0):
        import h5py
        with h5py.File(path, "r") as h:
            names = [s.decode() if isinstance(s, bytes) else str(s)
                     for s in h.attrs["species"]]
            if species not in names:
                raise SystemExit("no species %r in the dataset; it has %s"
                                 % (species, names))
            ci = names.index(species)

            self.material = h["geom/material"][:].astype(np.float32)
            self.gdf = np.nan_to_num(h["geom/gdf"][:]).astype(np.float32)
            gidx = h["samples/geom_index"][:]
            params = h["samples/params"][:].astype(np.float32)
            t_norm = h["samples/t_norm"][:].astype(np.float32)
            scale = np.asarray(h["samples"].attrs["conc_scale"], np.float32)

            if split is not None and "split" in h["samples"]:
                sp = np.array([x.decode() if isinstance(x, bytes) else str(x)
                               for x in h["samples/split"][:]])
                sel = np.where(sp == split)[0]
            else:
                sel = np.arange(len(gidx))
            if sel.size == 0:
                raise SystemExit("no samples in split %r" % split)

            conc = np.empty((sel.size,) + h["samples/conc"].shape[1:2]
                            + h["samples/conc"].shape[3:], np.float32)
            for k, s in enumerate(sel):
                conc[k] = h["samples/conc"][s, :, ci].astype(np.float32)

        self.conc = conc / max(float(scale[ci]), 1e-30)
        self.geom_index = gidx[sel]
        self.params = params[sel]
        self.t_norm = t_norm[sel]
        self.species = species
        self.conc_scale = float(scale[ci])
        self.n_points = n_points
        self.rng = np.random.default_rng(seed)

        # normalise the trunk's geometry column once, not per item
        self.gdf_n = self.gdf / max(float(self.gdf.max()), 1.0)
        self.pore = [np.argwhere(self.material[g] == PORE)
                     for g in range(self.material.shape[0])]
        self.S, self.T = self.conc.shape[0], self.conc.shape[1]
        self.shape = self.conc.shape[2:]
        self.index = [(s, t) for s in range(self.S) for t in range(self.T)]

    def __len__(self):
        return len(self.index)

    def item(self, i, full=False):
        s, t = self.index[i]
        g = int(self.geom_index[s])
        pts = self.pore[g]
        if not full and self.n_points < len(pts):
            pick = self.rng.choice(len(pts), self.n_points, replace=False)
            pts = pts[pick]
        x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
        nx, ny, nz = self.shape
        trunk = np.stack([
            x / (nx - 1.0), y / (ny - 1.0), z / (nz - 1.0),
            np.full(len(pts), self.t_norm[s, t], np.float32),
            self.gdf_n[g][x, y, z],
        ], 1).astype(np.float32)
        branch1 = (self.material[g] == PORE).astype(np.float32)[None]
        target = self.conc[s, t][x, y, z].astype(np.float32)
        return branch1, self.params[s], trunk, target


def batches(ds, batch_size, shuffle=True):
    order = np.arange(len(ds))
    if shuffle:
        ds.rng.shuffle(order)
    for k in range(0, len(order), batch_size):
        chunk = order[k:k + batch_size]
        b1, b2, tr, ta = zip(*[ds.item(int(i)) for i in chunk])
        yield (np.stack(b1), np.stack(b2), np.stack(tr), np.stack(ta))


# ------------------------------------------------------------------ model ---
def build_model(n_params, n_trunk, width=128):
    import torch.nn as nn

    class Branch1(nn.Module):
        def __init__(self):
            super().__init__()
            chans, layers, cin = (16, 32, 64, 128, 256), [], 1
            for c in chans:
                layers += [nn.Conv3d(cin, c, 3, padding=1), nn.SiLU(),
                           nn.AvgPool3d(2)]
                cin = c
            self.conv = nn.Sequential(*layers)
            self.head = nn.Sequential(nn.Flatten(), nn.LazyLinear(width))

        def forward(self, x):
            return self.head(self.conv(x))

    def mlp(nin, depth):
        import torch.nn as nn
        seq, prev = [], nin
        for _ in range(depth):
            seq += [nn.Linear(prev, width), nn.SiLU()]
            prev = width
        return nn.Sequential(*seq)

    class DeepONet(nn.Module):
        def __init__(self):
            super().__init__()
            self.b1 = Branch1()
            self.b2 = mlp(n_params, 3)
            self.trunk = mlp(n_trunk, 8)
            self.bias = nn.Parameter(__import__("torch").zeros(1))

        def forward(self, g, p, t):
            # (B, W) * (B, W) -> (B, W), then dotted with (B, N, W)
            b = self.b1(g) * self.b2(p)
            return (self.trunk(t) * b[:, None, :]).sum(-1) + self.bias

    return DeepONet()


# THE METRIC IS AGGREGATED OVER A SPLIT, NOT AVERAGED OVER SAMPLES.
#
#     rel L2 = sqrt( sum over every point (pred - true)^2 / sum over every point true^2 )
#
# The obvious alternative, the mean of each sample's own relative error, is
# unusable here and was tried first. Early snapshots hold almost no product: C at
# the second frame of a low-Damkohler run is a field of near-zeros, its own norm
# is near zero, and its relative error is enormous however good the prediction
# is. One such sample then dominates the average and the reported score is four
# figures. Aggregating first and dividing once gives the error of the FIELD,
# which is the quantity anyone means by "how wrong is it", and it cannot be
# driven by a sample that has nothing in it.
def _sq_sums(pred, true):
    return float(((pred - true) ** 2).sum()), float((true ** 2).sum())


def evaluate(model, ds, batch_size, device):
    """Aggregate relative L2 over the whole split."""
    import torch
    model.eval()
    num = den = 0.0
    with torch.no_grad():
        for b1, b2, tr, ta in batches(ds, batch_size, shuffle=False):
            b1 = torch.as_tensor(b1, device=device)
            b2 = torch.as_tensor(b2, device=device)
            tr = torch.as_tensor(tr, device=device)
            ta = torch.as_tensor(ta, device=device)
            n, d = _sq_sums(model(b1, b2, tr), ta)
            num += n
            den += d
    return float(np.sqrt(num / max(den, 1e-30)))


def baselines(ds, batch_size):
    """Two things the model has to beat to have learned anything.

    The dataset mean is the weaker one: beating it only means the model knows
    roughly how much product there is. The per-sample mean is the one that
    matters, because beating THAT means the model has learned the SHAPE of the
    field rather than its magnitude.
    """
    import torch
    tot = cnt = 0.0
    for _, _, _, ta in batches(ds, batch_size, shuffle=False):
        tot += float(ta.sum())
        cnt += ta.size
    gmean = tot / max(cnt, 1)
    ng = dg = ns = 0.0
    for _, _, _, ta in batches(ds, batch_size, shuffle=False):
        t = torch.as_tensor(ta)
        a, d = _sq_sums(torch.full_like(t, gmean), t)
        b, _ = _sq_sums(t.mean(-1, keepdim=True).expand_as(t), t)
        ng += a
        ns += b
        dg += d
    return (float(np.sqrt(ng / max(dg, 1e-30))),
            float(np.sqrt(ns / max(dg, 1e-30))))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default="dataset/dataset.h5")
    ap.add_argument("--species", default="C")
    ap.add_argument("--out", default="model")
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--n-points", type=int, default=2048)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--patience", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    try:
        import torch
    except ImportError:
        sys.exit("PyTorch is not installed. pip install torch, or train with "
                 "3D/model/train.py in the PRT-DeepONet-GUI instead.")
    torch.manual_seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    os.makedirs(args.out, exist_ok=True)

    tr = ABCDataset(args.dataset, args.species, "train", args.n_points, args.seed)
    va = ABCDataset(args.dataset, args.species, "val", args.n_points, args.seed + 1)
    te = ABCDataset(args.dataset, args.species, "test", args.n_points, args.seed + 2)
    print("species %s   train %d / val %d / test %d samples   on %s"
          % (args.species, len(tr), len(va), len(te), device))
    print("  geometries: train %s  val %s  test %s"
          % (sorted(set(tr.geom_index.tolist())),
             sorted(set(va.geom_index.tolist())),
             sorted(set(te.geom_index.tolist()))))

    model = build_model(tr.params.shape[1], 5).to(device)
    # LazyLinear needs one forward pass before the optimiser can see its weights
    b1, b2, t4, _ = tr.item(0)
    model(torch.as_tensor(b1[None], device=device),
          torch.as_tensor(b2[None], device=device),
          torch.as_tensor(t4[None], device=device))
    nparam = sum(p.numel() for p in model.parameters())
    print("  %d parameters" % nparam)

    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    best, best_ep, since = float("inf"), -1, 0
    ckpt = os.path.join(args.out, "best_%s.pt" % args.species)
    history = []

    for ep in range(args.epochs):
        model.train()
        t0, tot, n = time.time(), 0.0, 0
        for bb1, bb2, btr, bta in batches(tr, args.batch_size):
            bb1 = torch.as_tensor(bb1, device=device)
            bb2 = torch.as_tensor(bb2, device=device)
            btr = torch.as_tensor(btr, device=device)
            bta = torch.as_tensor(bta, device=device)
            opt.zero_grad()
            # Plain mean squared error to train on: it is smooth, it is finite for
            # a target of all zeros, and it does not let a nearly empty early
            # snapshot dominate the gradient. The relative L2 above is the metric,
            # not the objective.
            loss = torch.nn.functional.mse_loss(model(bb1, bb2, btr), bta)
            loss.backward()
            opt.step()
            tot += float(loss.detach()) * len(bta)
            n += len(bta)
        v = evaluate(model, va, args.batch_size, device)
        history.append({"epoch": ep, "train_mse": tot / max(n, 1), "val_rel_l2": v})
        mark = ""
        if v < best:
            best, best_ep, since = v, ep, 0
            torch.save(model.state_dict(), ckpt)
            mark = "  <- best"
        else:
            since += 1
        print("  epoch %3d   train mse %.5g   val relL2 %.4f   %.1fs%s"
              % (ep, tot / max(n, 1), v, time.time() - t0, mark))
        if since >= args.patience:
            print("  stopping: no improvement in %d epochs" % args.patience)
            break

    model.load_state_dict(torch.load(ckpt, map_location=device))
    test = evaluate(model, te, args.batch_size, device)
    bg, bs = baselines(te, args.batch_size)

    print()
    print("=" * 72)
    print("held-out TEST relative L2, species %s" % args.species)
    print("  model                      %.4f" % test)
    print("  baseline, dataset mean     %.4f" % bg)
    print("  baseline, per-sample mean  %.4f" % bs)
    print()
    if test < 0.5 * min(bg, bs):
        verdict = ("PASS. The model is comfortably better than both baselines on "
                   "geometries it never saw, so the sweep's inputs carry the "
                   "information its outputs need.")
        rc = 0
    elif test < min(bg, bs):
        verdict = ("MARGINAL. Better than both baselines but not by much. Train "
                   "longer, or look at whether the sweep spans enough of the "
                   "parameter space to be worth generalising over.")
        rc = 0
    else:
        verdict = ("FAIL. The model does not beat predicting a mean. Either the "
                   "training run was far too short, or the dataset does not relate "
                   "its inputs to its outputs -- run tools/verify_dataset.py, which "
                   "tests that directly.")
        rc = 1
    print(verdict)
    print("=" * 72)

    json.dump({"species": args.species, "test_rel_l2": test,
               "baseline_dataset_mean": bg, "baseline_sample_mean": bs,
               "best_val": best, "best_epoch": best_ep,
               "epochs_run": len(history), "parameters": nparam,
               "conc_scale": tr.conc_scale, "verdict": verdict.split(".")[0],
               "history": history},
              open(os.path.join(args.out, "result_%s.json" % args.species), "w"),
              indent=2)
    print("wrote %s" % os.path.join(args.out, "result_%s.json" % args.species))
    return rc


if __name__ == "__main__":
    sys.exit(main())
