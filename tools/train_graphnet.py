#!/usr/bin/env python3
"""
train_graphnet.py -- fit a bipartite graph network to a reaction network and write a .gnn file.

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
Meile Lab, University of Georgia.

WHY A GRAPH AND NOT A PLAIN NETWORK
    A dense network from concentrations to rates has to learn, from data, which species affect
    which.  That information is already known: it is the stoichiometric matrix.  Putting it in as
    structure rather than making the network infer it means fewer weights, less training data, and
    a model that cannot invent a coupling the chemistry does not have.

THE GRAPH
    Bipartite.  One node per species, one node per reaction, an edge wherever the stoichiometric
    coefficient s[i][r] is non-zero, carrying that coefficient as its weight.  This is the
    stoichiometric matrix read as a graph, so nothing is approximated.

ONE ROUND OF MESSAGE PASSING
    reactions   hR_r = tanh( Wsr . sum_i s[i][r] * hS_i  +  Wda * da_r  +  bR )
    species     hS_i = tanh( Wss . hS_i  +  Wrs . sum_r s[i][r] * hR_r  +  bS )

    The sums are weighted by the stoichiometry, so a species that consumes two moles per reaction
    speaks twice as loudly into that reaction as one that consumes one.  da_r is the per-reaction
    Damkohler number, an edge-level parameter that lets one trained network cover a range of
    transport regimes instead of one.

INPUT AND OUTPUT
    In:  one concentration per species.
    Out: one rate per species, plus optionally a growth rate.
    Both are scaled the way the .srg files scale theirs, so the same reasoning applies: a network
    is fitted over a box and is meaningless outside it, so the box travels in the file.

    python train_graphnet.py --stoich stoich.csv --data samples.csv --out network.gnn
"""

import argparse
import datetime
import sys

import numpy as np


# ------------------------------------------------------------------------------------------------
# the model
# ------------------------------------------------------------------------------------------------
class GraphNet(object):
    def __init__(self, S, width, rounds, rng, n_out):
        """S: (n_species, n_reactions) stoichiometric matrix."""
        self.S = np.asarray(S, dtype=np.float64)
        self.nS, self.nR = self.S.shape
        self.W = width
        self.L = rounds
        self.n_out = n_out

        def he(shape, fan_in):
            return rng.normal(0.0, np.sqrt(1.0 / max(fan_in, 1)), size=shape)

        # one scalar concentration per species, lifted to width W
        self.Wenc = he((width, 1), 1)
        self.benc = np.zeros(width)

        self.Wsr, self.Wda, self.bR = [], [], []
        self.Wss, self.Wrs, self.bS = [], [], []
        for _ in range(rounds):
            self.Wsr.append(he((width, width), width))
            self.Wda.append(he((width, 1), 1))
            self.bR.append(np.zeros(width))
            self.Wss.append(he((width, width), width))
            self.Wrs.append(he((width, width), width))
            self.bS.append(np.zeros(width))

        # The readout has ONE row shared by every species node, plus a second row for growth when
        # there is one.  Not one row per species: the species are told apart by their own hidden
        # state and by the per-species output scaling, and giving each its own row would be nS-1
        # rows the gradient never reaches.
        self.n_head = 2 if n_out > 1 else 1
        self.Wout = he((self.n_head, width), width)
        self.bout = np.zeros(self.n_head)

    def params(self):
        p = [self.Wenc, self.benc]
        for l in range(self.L):
            p += [self.Wsr[l], self.Wda[l], self.bR[l], self.Wss[l], self.Wrs[l], self.bS[l]]
        p += [self.Wout, self.bout]
        return p

    def forward(self, x, da):
        """x: (batch, nS) scaled concentrations.  da: (nR,).  Returns (batch, n_out)."""
        B = x.shape[0]
        # encode: every species node starts from its own scaled concentration
        hS = np.tanh(x[:, :, None] * self.Wenc[:, 0][None, None, :] + self.benc[None, None, :])
        for l in range(self.L):
            # species -> reactions, weighted by stoichiometry
            mR = np.einsum('ir,bih->brh', self.S, hS)
            hR = np.tanh(mR @ self.Wsr[l].T
                         + da[None, :, None] * self.Wda[l][:, 0][None, None, :]
                         + self.bR[l][None, None, :])
            # reactions -> species, same weights
            mS = np.einsum('ir,brh->bih', self.S, hR)
            hS = np.tanh(hS @ self.Wss[l].T + mS @ self.Wrs[l].T + self.bS[l][None, None, :])
        # read out: one value per species node, then the growth slot if there is one
        y = (hS @ self.Wout[0]) + self.bout[0]
        if self.n_out > 1:
            g = (hS.mean(axis=1) @ self.Wout[1]) + self.bout[1]
            return np.concatenate([y, g[:, None]], axis=1)
        return y


# ------------------------------------------------------------------------------------------------
# fitting
# ------------------------------------------------------------------------------------------------
def numerical_fit(net, X, Y, da, epochs, lr, rng, verbose=True):
    """Adam on the squared error.  Gradients by central differences on a random parameter subset
    each step: slow per step but exact to the model, and it keeps this script dependency-free."""
    ps = net.params()
    m = [np.zeros_like(p) for p in ps]
    v = [np.zeros_like(p) for p in ps]
    b1, b2, eps = 0.9, 0.999, 1e-8

    def loss():
        return float(np.mean((net.forward(X, da) - Y) ** 2))

    best, bestp = loss(), [p.copy() for p in ps]
    stale = 0
    for t in range(1, epochs + 1):
        for k, p in enumerate(ps):
            flat = p.reshape(-1)
            n = flat.size
            idx = rng.choice(n, size=min(n, 24), replace=False)
            g = np.zeros(n)
            h = 1e-5
            for i in idx:
                o = flat[i]
                flat[i] = o + h; lp = loss()
                flat[i] = o - h; lm = loss()
                flat[i] = o
                g[i] = (lp - lm) / (2 * h)
            g = g.reshape(p.shape)
            m[k] = b1 * m[k] + (1 - b1) * g
            v[k] = b2 * v[k] + (1 - b2) * g * g
            mh = m[k] / (1 - b1 ** t)
            vh = v[k] / (1 - b2 ** t)
            p -= lr * mh / (np.sqrt(vh) + eps)
        cur = loss()
        if cur < best:
            best, bestp = cur, [p.copy() for p in ps]
            stale = 0
        else:
            # Finite-difference gradients on a random subset are noisy, so a step size that was
            # fine at the start will eventually overshoot and the loss climbs away.  Rather than
            # leave that for the user to notice in the printout, back the step off and go back to
            # the best parameters seen.  Without this the run still ENDS at the best fit -- that
            # is restored below -- but it wastes every epoch after the turn.
            stale += 1
            if stale >= 5:
                lr *= 0.5
                for p, bp in zip(ps, bestp):
                    p[...] = bp
                m = [np.zeros_like(p) for p in ps]
                v = [np.zeros_like(p) for p in ps]
                stale = 0
                if verbose:
                    print("  epoch %5d   stalled; step size halved to %.4g" % (t, lr))
        if verbose and (t % max(1, epochs // 10) == 0):
            print("  epoch %5d   mse %.6e   best %.6e" % (t, cur, best))
    for p, bp in zip(ps, bestp):
        p[...] = bp
    return best


# ------------------------------------------------------------------------------------------------
# the file
# ------------------------------------------------------------------------------------------------
def write_gnn(path, net, species, reactions, da, xoff, xgain, yoff, ygain,
              trainmin, trainmax, units, growth, provenance):
    def row(f, name, a):
        f.write(name + " " + " ".join("%.17g" % v for v in np.asarray(a).reshape(-1)) + "\n")

    with open(path, "w") as f:
        f.write("# CompLaB3D graph network\n")
        f.write("# THE NETWORK IS ONLY VALID INSIDE trainmin..trainmax. Outside that box a\n")
        f.write("# neural network does not fail, it returns confident nonsense.\n")
        # '#' lines are comments and are dropped on load. 'provenance' lines are kept and
        # printed in the run log, so a result carries the network that produced it.
        for line in provenance:
            f.write("provenance " + line + "\n")
        f.write("version 1\n")
        f.write("units %s\n" % units)
        f.write("species %s\n" % " ".join(species))
        f.write("reactions %s\n" % " ".join(reactions))
        f.write("rounds %d\n" % net.L)
        f.write("width %d\n" % net.W)
        f.write("growth %d\n" % (1 if growth else 0))
        f.write("stoich %d %d\n" % (net.nS, net.nR))
        for i in range(net.nS):
            f.write(" ".join("%.17g" % v for v in net.S[i]) + "\n")
        row(f, "da", da)
        row(f, "xoffset", xoff)
        row(f, "xgain", xgain)
        f.write("xymin -1\n")
        row(f, "yoffset", yoff)
        row(f, "ygain", ygain)
        f.write("yymin -1\n")
        row(f, "trainmin", trainmin)
        row(f, "trainmax", trainmax)
        row(f, "Wenc", net.Wenc)
        row(f, "benc", net.benc)
        for l in range(net.L):
            f.write("layer %d\n" % l)
            row(f, "Wsr", net.Wsr[l]); row(f, "Wda", net.Wda[l]); row(f, "bR", net.bR[l])
            row(f, "Wss", net.Wss[l]); row(f, "Wrs", net.Wrs[l]); row(f, "bS", net.bS[l])
        row(f, "Wout", net.Wout)
        row(f, "bout", net.bout)


def read_gnn(path):
    """Read a .gnn back.  Returns (net, meta) where meta carries the scaling, da and units.

    This exists so the file is round-trippable on the Python side too: it is what
    tests/xval_gnn.py uses to check the C++ forward pass against this one, and it is how you
    would load a fitted network to plot it or to compare two fits."""
    toks, i = [], 0
    for line in open(path):
        if line.startswith("#"):
            continue
        if line.startswith("provenance"):
            continue
        toks += line.split()

    M = {"units": "per_second"}
    S, layers = None, []

    def take(n):
        nonlocal i
        v = np.array([float(x) for x in toks[i:i + n]], dtype=np.float64)
        i += n
        return v

    while i < len(toks):
        k = toks[i]; i += 1
        if k == "version":
            i += 1
        elif k == "units":
            M["units"] = toks[i]; i += 1
        elif k in ("species", "reactions"):
            names = []
            while i < len(toks) and toks[i] not in ("reactions", "rounds", "width",
                                                    "growth", "stoich"):
                names.append(toks[i]); i += 1
            M[k] = names
        elif k in ("rounds", "width", "growth"):
            M[k] = int(toks[i]); i += 1
        elif k == "stoich":
            M["nS"], M["nR"] = int(toks[i]), int(toks[i + 1]); i += 2
            S = take(M["nS"] * M["nR"]).reshape(M["nS"], M["nR"])
        elif k == "da":       M["da"] = take(M["nR"])
        elif k == "xoffset":  M["xoff"] = take(M["nS"])
        elif k == "xgain":    M["xgain"] = take(M["nS"])
        elif k == "xymin":    M["xymin"] = float(toks[i]); i += 1
        elif k == "yoffset":  M["yoff"] = take(M["nS"] + M["growth"])
        elif k == "ygain":    M["ygain"] = take(M["nS"] + M["growth"])
        elif k == "yymin":    M["yymin"] = float(toks[i]); i += 1
        elif k == "trainmin": M["trainmin"] = take(M["nS"])
        elif k == "trainmax": M["trainmax"] = take(M["nS"])
        elif k == "Wenc":     M["Wenc"] = take(M["width"]).reshape(M["width"], 1)
        elif k == "benc":     M["benc"] = take(M["width"])
        elif k == "layer":    layers.append({}); i += 1
        elif k in ("Wsr", "Wss", "Wrs"):
            layers[-1][k] = take(M["width"] ** 2).reshape(M["width"], M["width"])
        elif k == "Wda":      layers[-1][k] = take(M["width"]).reshape(M["width"], 1)
        elif k in ("bR", "bS"): layers[-1][k] = take(M["width"])
        elif k == "Wout":
            nhead = 2 if M["growth"] else 1
            M["Wout"] = take(nhead * M["width"]).reshape(nhead, M["width"])
        elif k == "bout":
            M["bout"] = take(2 if M["growth"] else 1)
        else:
            raise ValueError("unknown keyword '%s' in %s" % (k, path))

    net = GraphNet(S, M["width"], M["rounds"], np.random.default_rng(0),
                   M["nS"] + M["growth"])
    net.Wenc, net.benc = M["Wenc"], M["benc"]
    for l, L in enumerate(layers):
        net.Wsr[l], net.Wda[l], net.bR[l] = L["Wsr"], L["Wda"], L["bR"]
        net.Wss[l], net.Wrs[l], net.bS[l] = L["Wss"], L["Wrs"], L["bS"]
    net.Wout, net.bout = M["Wout"], M["bout"]
    M["unitscale"] = 1.0 / 3600.0 if M["units"] == "per_hour" else 1.0
    return net, M


def predict(net, M, C):
    """Concentrations in, rates out, in real units -- the same path the C++ takes: clamp to the
    training box, scale in, forward, scale out, apply the unit scale."""
    C = np.clip(np.atleast_2d(C), M["trainmin"], M["trainmax"])
    Xs = (C - M["xoff"]) * M["xgain"] + M["xymin"]
    Ys = net.forward(Xs, M["da"])
    return ((Ys - M["yymin"]) / M["ygain"] + M["yoff"]) * M["unitscale"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stoich", required=True,
                    help="CSV, one row per species, one column per reaction")
    ap.add_argument("--data", required=True,
                    help="CSV with a header: one column per species concentration, then one "
                         "rate column per species named <species>_rate, then growth if used")
    ap.add_argument("--out", default="network.gnn")
    ap.add_argument("--width", type=int, default=8)
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--lr", type=float, default=0.02)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--units", default="per_hour", choices=["per_second", "per_hour"])
    ap.add_argument("--da", default="", help="comma-separated Damkohler number per reaction")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    S = np.loadtxt(args.stoich, delimiter=",", ndmin=2)
    nS, nR = S.shape

    with open(args.data) as f:
        header = f.readline().strip().split(",")
    D = np.loadtxt(args.data, delimiter=",", skiprows=1, ndmin=2)

    species = [h for h in header if not h.endswith("_rate") and h != "growth"]
    if len(species) != nS:
        sys.exit("the data has %d concentration columns but the stoichiometry has %d species"
                 % (len(species), nS))
    growth = "growth" in header
    X = D[:, [header.index(s) for s in species]]
    Ycols = [header.index(s + "_rate") for s in species]
    Y = D[:, Ycols]
    if growth:
        Y = np.concatenate([Y, D[:, [header.index("growth")]]], axis=1)

    da = (np.array([float(v) for v in args.da.split(",")]) if args.da
          else np.ones(nR))
    if da.size != nR:
        sys.exit("--da has %d values but there are %d reactions" % (da.size, nR))

    lo, hi = X.min(axis=0), X.max(axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    xoff, xgain = lo, 2.0 / span
    ylo, yhi = Y.min(axis=0), Y.max(axis=0)
    yspan = np.where(yhi > ylo, yhi - ylo, 1.0)
    yoff, ygain = ylo, 2.0 / yspan

    Xs = (X - xoff) * xgain - 1.0
    Ys = (Y - yoff) * ygain - 1.0

    net = GraphNet(S, args.width, args.rounds, rng, Y.shape[1])
    print("fitting %d species, %d reactions, width %d, %d rounds, %d samples"
          % (nS, nR, args.width, args.rounds, X.shape[0]))
    mse = numerical_fit(net, Xs, Ys, da, args.epochs, args.lr, rng)

    pred = net.forward(Xs, da)
    r = np.corrcoef(pred.reshape(-1), Ys.reshape(-1))[0, 1]
    print("final scaled mse %.6e   R %.6f" % (mse, r))

    write_gnn(args.out, net, species,
              ["R%d" % (i + 1) for i in range(nR)], da,
              xoff, xgain, yoff, ygain, lo, hi, args.units, growth,
              ["fitted %s from %s" % (datetime.date.today().isoformat(), args.data),
               "samples %d, scaled mse %.4e, R %.5f" % (X.shape[0], mse, r)])
    print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
