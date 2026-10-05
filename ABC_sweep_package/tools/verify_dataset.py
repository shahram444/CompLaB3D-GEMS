#!/usr/bin/env python3
"""Decide whether the collected dataset is worth training on.

    python3 tools/verify_dataset.py --dataset dataset/dataset.h5

This is the verdict on the CompLaB3D half of the chain. It opens the .h5 and asks
the questions that a dataset can answer about itself, in rough order of how badly
a No would hurt:

  1  DOES DAMKOHLER DO ANYTHING?
     For every geometry and Peclet, each neighbouring pair of Damkohler cases must differ. If they
     do not, the solver was built against a header that ignores PRT_KABIO, every
     case ran identical chemistry while recording a different Damkohler, and the
     dataset teaches the network that one of its two inputs is noise. Every run
     succeeded, every file is present, and nothing else anywhere says so. This is
     the check the whole package is arranged around.

  2  DOES PECLET DO ANYTHING?
     The same question for the other input. A dataset where it does not is one
     where <delta_P> was zero, or the flow solver never ran.

  3  DOES THE GEOMETRY DO ANYTHING?
     At fixed Peclet and Damkohler, two pore spaces must give different fields.
     If they do not, the geometries are not reaching the solver.

  4  ARE THE BOUNDARY CONDITIONS THE ONES THAT WERE ASKED FOR?
     A should sit at the feed value on the inlet plane and at zero on the outlet
     plane; B the mirror image. This catches a swapped boundary, a Neumann that
     should have been Dirichlet, and a geometry written with the x axis reversed.

  5  IS ANYTHING NEGATIVE?
     The v1.3.2 clamp should have made this impossible. Anything above a trace is
     a sign the clamp is not in the binary.

  6  IS THE PRODUCT ACTUALLY BEING MADE?
     C must be non-zero in every sample past the first, and it must be where A
     and B overlap rather than at a boundary.

  7  IS float16 HOLDING THE FIELD?
     Reports the smallest positive stored value and whether it fell into the
     subnormal range, where relative error is several per cent.

  8  IS THE SPLIT HONEST?
     Train, validation and test must share no geometry. A network tested on a
     pore space it trained on is being scored on memory rather than transfer.

Exit status 0 means the dataset is worth training on.
"""
import argparse
import sys

import numpy as np

PORE = 2


class Verdict(object):
    def __init__(self):
        self.fail, self.warn, self.ok = [], [], []

    def check(self, condition, name, detail_ok="", detail_bad="", warn_only=False):
        if condition:
            self.ok.append((name, detail_ok))
        elif warn_only:
            self.warn.append((name, detail_bad))
        else:
            self.fail.append((name, detail_bad))
        return bool(condition)

    def show(self):
        for n, d in self.ok:
            print("  ok    %-34s %s" % (n, d))
        for n, d in self.warn:
            print("  WARN  %-34s %s" % (n, d))
        for n, d in self.fail:
            print("  FAIL  %-34s %s" % (n, d))
        print()
        if self.fail:
            print("VERDICT: NOT usable. %d check(s) failed." % len(self.fail))
            return 1
        if self.warn:
            print("VERDICT: usable, with %d thing(s) to read first." % len(self.warn))
            return 0
        print("VERDICT: the CompLaB3D sweep produced a sound dataset.")
        return 0


def rel_diff(a, b, mask):
    """Relative L2 difference between two fields, over the pore space only."""
    a, b = a[mask].astype(np.float64), b[mask].astype(np.float64)
    den = np.sqrt((a * a).sum()) + np.sqrt((b * b).sum())
    if den == 0:
        return 0.0
    return float(2.0 * np.sqrt(((a - b) ** 2).sum()) / den)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default="dataset/dataset.h5")
    ap.add_argument("--feed", type=float, default=1.0,
                    help="the Dirichlet feed the campaign used")
    ap.add_argument("--sensitivity", type=float, default=0.02,
                    help="relative difference below which two cases count as "
                         "the same field")
    args = ap.parse_args(argv)

    import h5py
    V = Verdict()
    h = h5py.File(args.dataset, "r")

    species = [s.decode() if isinstance(s, bytes) else str(s)
               for s in h.attrs["species"]]
    mat = h["geom/material"][:]
    gidx = h["samples/geom_index"][:]
    gid = h["geom/gid"][:]
    par = h["samples/params"][:]
    tn = h["samples/t_norm"][:]
    conc = h["samples/conc"]
    S, T, C = conc.shape[0], conc.shape[1], conc.shape[2]
    shape = conc.shape[3:]
    bnd = h.attrs.get("boundaries", b"unknown")
    bnd = bnd.decode() if isinstance(bnd, bytes) else str(bnd)

    def probe(n):
        """n sample indices spread evenly across the campaign, not the first n.

        Samples are written in run order and a run id is
        (cases per geometry) * gid + ..., so the FIRST n rows are the first two or
        three geometries only: one porosity level, all of it in the train split.
        A check that reads conc[0:24] would then never see a high-porosity pack, a
        held-out geometry, or the largest Peclet in anything but the first pack.
        A stride spans every gid, both ends of each parameter ladder and all three
        splits for the same number of reads.
        """
        if S <= n:
            return list(range(S))
        step = S / float(n)
        return sorted({int(i * step) for i in range(n)})

    print("%s" % args.dataset)
    print("  %d samples, %d geometries, %d snapshots, %d species %s, shape %s"
          % (S, len(gid), T, C, species, tuple(shape)))
    print("  params %s" % [s.decode() if isinstance(s, bytes) else str(s)
                           for s in h.attrs["param_names"]])
    print()

    # ---- structure --------------------------------------------------------
    V.check(par.shape == (S, 2), "parameter table shape",
            "(%d, 2)" % S, "expected (%d, 2), got %s" % (S, par.shape))
    V.check(tn.shape == (S, T), "time axis shape",
            "(%d, %d)" % (S, T), "expected (%d, %d), got %s" % (S, T, tn.shape))
    V.check(np.all(np.diff(tn, axis=1) > 0), "time axis increases",
            "strictly, in every sample", "some sample's t_norm does not increase")
    V.check(np.allclose(tn, tn[0]), "time axis is common",
            "every sample is on the same t_norm grid",
            "samples are on different time grids, so slot k does not mean one "
            "thing across the dataset", warn_only=True)
    V.check("conc_scale" in h["samples"].attrs, "conc_scale recorded",
            str(np.asarray(h["samples"].attrs.get("conc_scale", [])).round(5).tolist()),
            "dataset_reader.py divides the target by conc_scale and will refuse "
            "a file without it")
    V.check("velocity" in h["samples"], "flow field present",
            "/samples/velocity", "no /samples/velocity, so the flow-proxy and "
            "dimension-free switches cannot be used with this file",
            warn_only=True)

    pore = (mat == PORE)                              # (G, nx, ny, nz)

    # ---- 5. negativity ----------------------------------------------------
    worst_neg, worst_s = 0.0, -1
    for s in range(S):
        m = float(conc[s].min())
        if m < worst_neg:
            worst_neg, worst_s = m, s
    V.check(worst_neg > -1e-6 * args.feed, "nothing went negative",
            "min over the whole file %.3g" % worst_neg,
            "sample %d reaches %.3g. The v1.3.2 positivity clamp should make this "
            "impossible; check the [CLAMP] line in that run's log."
            % (worst_s, worst_neg))

    # ---- 4. boundary conditions ------------------------------------------
    # x = 0 is the inlet plane and x = nx-1 the outlet. Measured on pore voxels
    # only: a solid voxel carries no dynamics and whatever its density reports is
    # not a concentration.
    #
    # Only the two FEED faces are held in every boundary variant, so only those
    # two can be tested against a known number. The rest is tested by direction:
    # A must fall from its feed face to the far one and B must rise, whether the
    # far face is a held zero or an open outflow.
    iA, iB = species.index("A"), species.index("B")
    a_in, a_out, b_in, b_out = [], [], [], []
    for s in probe(24):
        g = gidx[s]
        pin, pout = pore[g][0], pore[g][-1]
        last = conc[s, -1]
        if pin.any():
            a_in.append(float(last[iA][0][pin].mean()))
            b_in.append(float(last[iB][0][pin].mean()))
        if pout.any():
            a_out.append(float(last[iA][-1][pout].mean()))
            b_out.append(float(last[iB][-1][pout].mean()))
    f = args.feed
    V.check(abs(np.mean(a_in) - f) < 0.02 * f, "A is held at the inlet",
            "mean %.4g on x = 0 against a feed of %g" % (np.mean(a_in), f),
            "mean %.4g on x = 0, not the feed of %g. The inlet Dirichlet value is "
            "not what the XML asked for, or the x axis is reversed."
            % (np.mean(a_in), f))
    V.check(abs(np.mean(b_out) - f) < 0.02 * f, "B is held at the outlet",
            "mean %.4g on x = nx-1 against a feed of %g" % (np.mean(b_out), f),
            "mean %.4g on x = nx-1, not the feed of %g" % (np.mean(b_out), f))
    V.check(np.mean(a_out) < np.mean(a_in), "A falls across the domain",
            "%.4g at the inlet down to %.4g at the outlet"
            % (np.mean(a_in), np.mean(a_out)),
            "A is %.4g at the inlet and %.4g at the outlet. It should be consumed "
            "on the way; check the flow direction and the reaction."
            % (np.mean(a_in), np.mean(a_out)))
    V.check(np.mean(b_in) < np.mean(b_out), "B falls the other way",
            "%.4g at the outlet down to %.4g at the inlet"
            % (np.mean(b_out), np.mean(b_in)),
            "B is %.4g at the outlet and %.4g at the inlet, so it is not entering "
            "from the outlet face" % (np.mean(b_out), np.mean(b_in)))

    # Under the all-Dirichlet set the FAR face of each reactant is held at zero
    # too, and a held zero is an exact number rather than a direction. Checking it
    # catches a swapped or dropped plane that the two directional tests above
    # would pass: if A's far face were left open, A would still fall across the
    # domain and the test would say nothing. This check only exists when the file
    # records that boundary set, which the collector writes from the campaign.
    if bnd == "dirichlet":
        V.check(np.mean(a_out) < 0.02 * f, "A's far face is held at zero",
                "mean %.3g on x = nx-1, against a feed of %g" % (np.mean(a_out), f),
                "mean %.3g on x = nx-1. The campaign records boundaries = dirichlet, "
                "so that plane should be held at 0.0. It is not, which means the "
                "rendered XML and the recorded boundary set disagree."
                % np.mean(a_out))
        V.check(np.mean(b_in) < 0.02 * f, "B's far face is held at zero",
                "mean %.3g on x = 0, against a feed of %g" % (np.mean(b_in), f),
                "mean %.3g on x = 0, which should be held at 0.0 under "
                "boundaries = dirichlet" % np.mean(b_in))
        iCf = species.index("C")
        cf = []
        for s in probe(24):
            g = gidx[s]
            for plane, msk in ((0, pore[g][0]), (-1, pore[g][-1])):
                if msk.any():
                    cf.append(float(conc[s, -1, iCf][plane][msk].mean()))
        V.check(max(cf) < 0.02 * f, "C is held at zero on both x faces",
                "largest face mean %.3g" % max(cf),
                "C reaches %.3g on an x face that boundaries = dirichlet holds at "
                "0.0" % max(cf))

    # ---- 5b. is a boundary manufacturing mass? ---------------------------
    # The open-outflow variant uses Neumann planes, and
    # FlatAdiabaticBoundaryFunctional3D re-defines such a plane from the layer
    # inside it every step, which tops a field rising from zero up from nothing
    # and streams the surplus back in. Nothing here is fed above the feed
    # concentration and nothing can concentrate, so anything above it came from
    # a boundary rather than from the chemistry.
    # Every run, not a sample of them. This is the check campaign01 would have
    # failed, and the run that reached 67 x the feed was run_0010: the 11th, which
    # a first-24 window happens to include and a first-6 window would not. One
    # frame set per run is about 4 MB, so the whole pass is a single streamed read.
    peak_of = [float(np.asarray(conc[s], np.float32).max()) for s in range(S)]
    peak = max(peak_of)
    over = sum(1 for x in peak_of if x > 1.05 * f)
    # Above a quarter over the feed this stops being something to note and
    # becomes a reason not to train on the file. A campaign of the mixed boundary
    # set reached 67 x the feed at Peclet 2, with half the pore voxels of one run
    # above 1.05; a field like that is not the reaction it claims to be anywhere,
    # not just near the plane that made it.
    V.check(peak <= 1.05 * f, "nothing exceeds the feed",
            "the largest value anywhere is %.4g against a feed of %g, over all %d "
            "runs" % (peak, f, S),
            "%d of %d runs exceed the feed of %g and the worst reaches %.4g. Nothing "
            "in this problem can concentrate, so it was manufactured at a boundary. "
            "Rebuild with --boundaries dirichlet, which is the default."
            % (over, S, f, peak),
            warn_only=(peak <= 1.25 * f))

    # ---- 6. is the product being made? -----------------------------------
    iC = species.index("C")
    # every run, since this claims every run: one frame each, about 17 MB in all
    made = [float(conc[s, -1, iC][pore[gidx[s]]].max()) for s in range(S)]
    V.check(min(made) > 1e-8 * f, "C is produced in every run",
            "peak C between %.3g and %.3g" % (min(made), max(made)),
            "at least one run ends with essentially no C (peak %.3g). A and B "
            "never met, or the reaction never fired." % min(made))
    # and it is made INSIDE, not at a face
    interior_fraction = []
    for s in probe(12):
        g = gidx[s]
        c = conc[s, -1, iC]
        tot = float(c[pore[g]].sum())
        if tot <= 0:
            continue
        faces = float(c[0][pore[g][0]].sum() + c[-1][pore[g][-1]].sum())
        interior_fraction.append(1.0 - faces / tot)
    V.check(np.mean(interior_fraction) > 0.9, "C is made inside the domain",
            "%.1f%% of it is away from the two x faces"
            % (100 * np.mean(interior_fraction)),
            "only %.1f%% of C is away from the x faces, so it is being made at a "
            "boundary rather than where A and B overlap"
            % (100 * np.mean(interior_fraction)), warn_only=True)

    # ---- 1, 2, 3. does each input do anything? ---------------------------
    key = {}
    for s in range(S):
        key[(int(gidx[s]), round(float(par[s, 0]), 6), round(float(par[s, 1]), 6))] = s

    def spread(pairs, label, what):
        """Largest relative difference of the final C field over the given pairs."""
        diffs = []
        for s1, s2, g in pairs:
            diffs.append(rel_diff(conc[s1, -1, iC], conc[s2, -1, iC], pore[g]))
        if not diffs:
            V.check(False, label, "", "no pair of runs to compare, so %s could not "
                    "be tested at all" % what, warn_only=True)
            return
        lo, hi = min(diffs), max(diffs)
        V.check(lo > args.sensitivity, label,
                "relative difference %.3f to %.3f over %d pair(s)"
                % (lo, hi, len(diffs)),
                "two runs differing only in %s give fields that agree to %.4f. "
                "%s" % (what, lo, DETAIL[label]))

    gids_present = sorted({int(g) for g in gidx})
    pes = sorted({round(float(p), 6) for p in par[:, 0]})
    das = sorted({round(float(p), 6) for p in par[:, 1]})

    # CONSECUTIVE pairs, not first against last.
    #
    # With two values on an axis those are the same thing. With three they are not:
    # comparing only 0.1 against 10 leaves every Da = 1.0 run, a third of the
    # campaign, never compared against anything. The failure this whole package
    # exists to catch is a binary that runs one chemistry while recording three,
    # and a version of it that spared the middle rung would pass a first-against-
    # last test. Each neighbouring pair has to differ on its own.
    pairs = []
    for g in gids_present:
        for pe in pes:
            for d1, d2 in zip(das, das[1:]):
                s1, s2 = key.get((g, pe, d1)), key.get((g, pe, d2))
                if s1 is not None and s2 is not None:
                    pairs.append((s1, s2, g))
    spread(pairs, "Damkohler changes the answer", "Damkohler")

    pairs = []
    for g in gids_present:
        for da in das:
            for p1, p2 in zip(pes, pes[1:]):
                s1, s2 = key.get((g, p1, da)), key.get((g, p2, da))
                if s1 is not None and s2 is not None:
                    pairs.append((s1, s2, g))
    spread(pairs, "Peclet changes the answer", "Peclet")

    # geometry: compare two different pore spaces, masked by the pore space they
    # share, so the difference is in the field and not merely in where pore is
    pairs_g = []
    if len(gids_present) >= 2:
        for pe in pes:
            for da in das:
                s1 = key.get((gids_present[0], pe, da))
                s2 = key.get((gids_present[-1], pe, da))
                if s1 is not None and s2 is not None:
                    shared = pore[gidx[s1]] & pore[gidx[s2]]
                    if shared.any():
                        pairs_g.append((s1, s2, shared))
    diffs = [rel_diff(conc[a, -1, iC], conc[b, -1, iC], m) for a, b, m in pairs_g]
    V.check(bool(diffs) and min(diffs) > args.sensitivity,
            "the pore space changes the answer",
            "relative difference %.3f to %.3f over %d pair(s)"
            % (min(diffs) if diffs else 0, max(diffs) if diffs else 0, len(diffs)),
            "two different pore spaces at the same Peclet and Damkohler agree to "
            "%.4f (over %d pair(s), threshold %.3f). Either the geometry is not "
            "reaching the solver -- check that each case's input/geometry.dat is "
            "its own -- or the two pore spaces are too alike to be worth "
            "distinguishing."
            % (min(diffs) if diffs else 0.0, len(diffs), args.sensitivity))

    # ---- 7. float16 headroom ---------------------------------------------
    # The question is not whether ANY value is subnormal: the tail of a front is
    # arbitrarily small and always will be, and a value at 1e-6 of the feed does
    # not matter however badly it is stored. The question is whether the BULK of
    # the field sits in the subnormal range, which is what happens if the feed
    # concentration is chosen at 1e-3 or below.
    vals = []
    for s in range(min(S, 8)):
        a = np.asarray(conc[s], np.float32)
        p = a[a > 0]
        if p.size:
            vals.append(p if p.size < 400000 else p[::max(1, p.size // 400000)])
    if vals:
        v = np.concatenate(vals)
        p90 = float(np.percentile(v, 90))
        frac_sub = float((v < 6.1e-5).mean())
        V.check(p90 > 6.1e-5, "float16 headroom",
                "90th percentile of the positive field is %.3g, well above float16's "
                "smallest normal; %.1f%% of values are in the subnormal tail"
                % (p90, 100 * frac_sub),
                "90 per cent of the positive field is below float16's smallest "
                "normal (6.1e-5), where the relative error runs to several per "
                "cent. Raise the feed concentration: A0 = 1.0 puts five decades of "
                "the field in the normal range.", warn_only=True)

    # ---- 8. is the split honest? -----------------------------------------
    if "split" in h["samples"]:
        sp = [x.decode() if isinstance(x, bytes) else str(x)
              for x in h["samples/split"][:]]
        by = {}
        for s, name in enumerate(sp):
            by.setdefault(name, set()).add(int(gidx[s]))
        overlap = set()
        names = list(by)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                overlap |= by[names[i]] & by[names[j]]
        V.check(not overlap, "the split shares no geometry",
                "  ".join("%s %d" % (k, len(by[k])) for k in sorted(by)),
                "geometries %s are in more than one split, so the score measures "
                "memory rather than transfer" % sorted(overlap))
        V.check(all(by.get(k) for k in ("train", "val", "test")),
                "all three splits are populated",
                "  ".join("%s %d geometries" % (k, len(by[k])) for k in sorted(by)),
                "one of train, validation or test has no geometry in it")
    else:
        V.check(False, "the split is recorded", "",
                "no /samples/split, so the training script has to invent one",
                warn_only=True)

    h.close()
    print()
    return V.show()


DETAIL = {
    "Damkohler changes the answer":
        "That is the silent failure this package exists to catch: the solver was "
        "built against config/kinetics/defineAbioticKinetics.default.hh, which has "
        "its rate compiled in and ignores PRT_KABIO. Every case ran the same "
        "chemistry. Rebuild with kinetics/defineAbioticKinetics.hh from this "
        "package and run tools/check_campaign.py check --complab first.",
    "Peclet changes the answer":
        "Either <delta_P> was zero, in which case complab_functions.hh silently "
        "sets Pe = 0 and every case was pure diffusion, or the flow solver did not "
        "converge. Check the [NS] Pe achieved line in the run logs.",
}

if __name__ == "__main__":
    sys.exit(main())
