#!/usr/bin/env python3
"""
/* This file is a part of the CompLaB program.
 *
 * The CompLaB3D software (3-D pore-scale extension) is developed since 2024
 * by the University of Georgia (United States, Meile Lab, Department of
 * Marine Sciences). The original 2-D CompLaB v1.0 was a collaboration of
 * the University of Georgia and Chungnam National University (South Korea).
 *
 * Contact:
 * Shahram Asgari, Christof Meile
 * Department of Marine Sciences (Meile Lab)
 * University of Georgia, Athens, GA 30602, USA
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

===============================================================================
STEPS 2 AND 3 OF 3 -- TRAIN THE SURROGATE AND EXPORT IT AS C++
===============================================================================

This is the MATLAB-free path.  trainSurrogate.m does the same job with the
Deep Learning Toolbox and emits a byte-comparable header; use whichever you
have a licence for.  Both write the SAME file format, so a network trained
here drops straight into a build that previously used a MATLAB-trained one.

WHY THIS FILE EXISTS AT ALL

  The network shipped in surrogateModel.hh was trained in MATLAB and exported
  with genFunction/gensim -- that is unmistakable from the variable names it
  left behind (IW1_1, LW2_1, tansig, mapminmax gains and offsets).  The
  training script itself was never distributed with CompLaB, so the shipped
  weights could be evaluated but never reproduced, retrained, or fitted to a
  different organism.  This file closes that gap.

WHAT IT REPRODUCES

  Architecture, defaults chosen to match the shipped network exactly:

      n inputs
        -> 4 hidden layers of 10 neurons, tansig activation
        -> 1 linear output neuron
      with mapminmax rescaling of the inputs and of the output onto [-1, 1].

  mapminmax is MATLAB's convention and it is worth stating precisely, because
  getting it wrong shifts every prediction:

      x_scaled = (x - x_min) * gain + (-1),      gain = 2 / (x_max - x_min)

  The generated C++ stores x_min as "x_offset" and that gain as "x_gain",
  which is exactly how MATLAB's genFunction writes it and exactly what
  surrogateModel.hh already contains.

  Training uses L-BFGS-B on analytic gradients, with MATLAB's default 70/15/15
  train/validation/test split and early stopping on the validation set.  It is
  not bit-identical to trainlm -- no gradient-based fit ever is across two
  implementations -- but it optimises the same loss over the same architecture
  and lands in the same place to within the noise of the random restart.

USAGE

    python3 trainSurrogate.py geobacter_training.csv \\
        --name geobacter --microbe-id 0 \\
        -o surrogate_weights_geobacter.hh

  Then two lines in surrogateModel.hh, which the script prints for you:

      #include "surrogate_weights_geobacter.hh"
      ...
      else if (microbeId == 0) {
          bioR[microbeId] = srgnet_geobacter::eval(Fin[microbeId]);
      }

OPTIONS WORTH KNOWING

    --layers 10 10 10 10   hidden layer sizes (default: the shipped shape)
    --restarts 5           independent fits, best validation error wins.
                           Raise it if the reported errors scatter; a small
                           net on a smooth FBA response surface should not.
    --log-output           fit log(growth) instead of growth.  Use when growth
                           spans orders of magnitude and you care about the
                           small values; the exported C++ undoes the log, so
                           nothing downstream changes.
    --seed 0               reproducibility.  Same seed, same data, same net.

REQUIREMENTS

    python3 -m pip install numpy scipy

===============================================================================
"""

import argparse
import os
import sys
import time

try:
    import numpy as np
except ImportError:
    sys.exit("ERROR: numpy is required.  python3 -m pip install numpy")
try:
    from scipy.optimize import minimize
except ImportError:
    sys.exit("ERROR: scipy is required.  python3 -m pip install scipy")


# ===========================================================================
# The network.  Weights are held as a flat vector so the optimiser sees one
# array; pack/unpack convert between that and the per-layer matrices.
# ===========================================================================
def shapes(sizes):
    """sizes = [n_in, h1, h2, ..., 1].  Returns [(rows, cols), ...] per layer."""
    return [(sizes[i + 1], sizes[i]) for i in range(len(sizes) - 1)]


def nparams(sizes):
    return sum(r * c + r for r, c in shapes(sizes))


def unpack(theta, sizes):
    W, B, k = [], [], 0
    for r, c in shapes(sizes):
        W.append(theta[k:k + r * c].reshape(r, c)); k += r * c
        B.append(theta[k:k + r]);                   k += r
    return W, B


def forward(W, B, X):
    """X is (n_samples, n_in).  Hidden layers tansig, output linear.
    Returns the output and every intermediate activation, which backprop needs."""
    A = [X.T]                                   # (n_in, n_samples)
    for i in range(len(W) - 1):
        A.append(np.tanh(W[i] @ A[-1] + B[i][:, None]))
    out = W[-1] @ A[-1] + B[-1][:, None]        # linear output layer
    A.append(out)
    return out.ravel(), A


def loss_and_grad(theta, sizes, X, y):
    """Mean squared error and its analytic gradient.

    tansig(n) = 2/(1+exp(-2n)) - 1 is algebraically identical to tanh(n), so
    np.tanh is used for the forward pass (it is the numerically stable form)
    and the derivative is the familiar 1 - a**2.  The generated C++ writes the
    2/(1+exp(-2n))-1 form because that is what MATLAB emits and what
    surrogateModel.hh already contains -- same function, same values."""
    W, B = unpack(theta, sizes)
    out, A = forward(W, B, X)
    n = len(y)
    err = out - y
    L = float(err @ err) / n

    gW = [None] * len(W)
    gB = [None] * len(B)
    D = (2.0 / n) * err[None, :]                    # dL/d(pre-activation), output layer
    for i in range(len(W) - 1, -1, -1):
        gW[i] = D @ A[i].T
        gB[i] = D.sum(axis=1)
        if i > 0:
            D = (W[i].T @ D) * (1.0 - A[i] ** 2)    # tansig'
    return L, np.concatenate([np.concatenate([g.ravel(), b]) for g, b in zip(gW, gB)])


def init_theta(sizes, rng):
    """Nguyen-Widrow-flavoured initialisation: scale each layer by its fan-in so
    the tansig units start in their responsive region instead of saturated.  A
    saturated tansig has a near-zero derivative and the layer never learns."""
    parts = []
    for r, c in shapes(sizes):
        parts.append(rng.normal(0.0, np.sqrt(2.0 / (c + r)), r * c))
        parts.append(np.zeros(r))
    return np.concatenate(parts)


# ===========================================================================
# mapminmax, exactly as MATLAB defines it
# ===========================================================================
class MapMinMax(object):
    def __init__(self, lo, hi, ymin=-1.0):
        lo = np.asarray(lo, float)
        hi = np.asarray(hi, float)
        span = hi - lo
        # A constant column has zero span.  MATLAB leaves such an input
        # untouched rather than dividing by zero; so do we.
        span = np.where(span == 0.0, 1.0, span)
        self.offset = lo
        self.gain = 2.0 / span
        self.ymin = ymin

    def apply(self, x):
        return (x - self.offset) * self.gain + self.ymin

    def reverse(self, y):
        return (y - self.ymin) / self.gain + self.offset


# ===========================================================================
# CSV loading, including the "#" metadata generateTrainingData.py wrote
# ===========================================================================
def load_csv(path):
    meta = {"inputs": [], "model": None, "objective": None}
    rows = []
    header = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                body = s[1:].strip()
                if body.startswith("model:"):
                    meta["model"] = body[6:].strip()
                elif body.startswith("objective:"):
                    meta["objective"] = body[10:].strip()
                elif body.startswith("input "):
                    meta["inputs"].append(body)
                continue
            if header is None and any(ch.isalpha() for ch in s.split(",")[0]):
                header = [c.strip() for c in s.split(",")]
                continue
            rows.append([float(v) for v in s.split(",")])
    if not rows:
        sys.exit("ERROR: %s contains no data rows." % path)
    D = np.asarray(rows, float)
    if D.shape[1] < 2:
        sys.exit("ERROR: %s needs at least one input column and one output "
                 "column; found %d column(s)." % (path, D.shape[1]))
    return D[:, :-1], D[:, -1], meta, header


# ===========================================================================
# The C++ emitter.  This is the part that has to be exactly right: the
# generated header is compiled into CompLB3D and evaluated millions of times.
# ===========================================================================
def fmt(v):
    """17 significant digits round-trips an IEEE double exactly.  Anything
    shorter silently changes the network."""
    return "%.17g" % float(v)


def emit_cpp(name, W, B, xmap, ymap, train_lo, train_hi, ylo, yhi,
             log_output, meta, stats, sizes, microbe_id, source_csv):
    guard = "SURROGATE_WEIGHTS_%s_HH" % name.upper()
    ns = "srgnet_%s" % name
    nin = sizes[0]
    L = []
    w = L.append

    w("/* ==========================================================================")
    w(" * CompLaB surrogate network -- GENERATED FILE, DO NOT EDIT BY HAND.")
    w(" *")
    w(" * Produced by trainSurrogate.py from: %s" % os.path.basename(source_csv))
    if meta.get("model"):
        w(" * Metabolic model              : %s" % meta["model"])
    if meta.get("objective"):
        w(" * FBA objective                : %s" % meta["objective"][:66])
    w(" *")
    w(" * Architecture : %d input(s) -> %s -> 1 linear output"
      % (nin, " -> ".join("%d tansig" % h for h in sizes[1:-1])))
    w(" * Parameters   : %d weights and biases" % nparams(sizes))
    if log_output:
        w(" * Output       : fitted in log space; eval() undoes the log for you.")
    w(" *")
    w(" * FIT QUALITY (on data the network never saw during training)")
    w(" *   test RMSE   : %.6g 1/h" % stats["rmse_test"])
    w(" *   test R^2    : %.6f" % stats["r2_test"])
    w(" *   worst error : %.6g 1/h" % stats["maxerr_test"])
    w(" *")
    w(" * TRAINING RANGE -- THE NETWORK IS ONLY VALID INSIDE THIS BOX.")
    w(" * A neural network extrapolates badly and will return confident")
    w(" * nonsense outside it.  inTrainingRange() below tests for exactly this;")
    w(" * call it if a run may wander out of range.")
    for k in range(nin):
        label = meta["inputs"][k] if k < len(meta["inputs"]) else "input %d" % k
        w(" *   %-58s" % label[:58])
        w(" *     %.10g .. %.10g mmol/gDW/h" % (train_lo[k], train_hi[k]))
    w(" *   growth rate: %.10g .. %.10g 1/h" % (ylo, yhi))
    w(" *")
    w(" * HOW TO USE IT -- two lines in surrogateModel.hh:")
    w(" *")
    w(' *     #include "%s"' % ("surrogate_weights_%s.hh" % name))
    w(" *")
    w(" *     else if (microbeId == %d) {" % microbe_id)
    w(" *         bioR[microbeId] = %s::eval(Fin[microbeId]);" % ns)
    w(" *     }")
    w(" *")
    w(" * Fin is in mmol/gDW/h and POSITIVE means consumption; eval() returns the")
    w(" * specific growth rate in 1/h.  Those are CompLaB's conventions already,")
    w(" * so no conversion is needed on either side.")
    w(" * ======================================================================== */")
    w("#ifndef %s" % guard)
    w("#define %s" % guard)
    w("")
    w("#include <vector>")
    w("#include <cmath>")
    w("#include <cstddef>")
    w("")
    w("namespace %s {" % ns)
    w("")
    w("static const std::size_t nInputs = %d;" % nin)
    w("")
    w("/* mapminmax input scaling: xs = (x - x_offset) * x_gain + x_ymin */")
    w("static const double x_offset[%d] = { %s };" % (nin, ", ".join(fmt(v) for v in xmap.offset)))
    w("static const double x_gain  [%d] = { %s };" % (nin, ", ".join(fmt(v) for v in xmap.gain)))
    w("static const double x_ymin      = %s;" % fmt(xmap.ymin))
    w("")
    w("/* mapminmax output scaling, applied in reverse: y = (ys - y_ymin)/y_gain + y_offset */")
    w("static const double y_offset = %s;" % fmt(float(np.atleast_1d(ymap.offset)[0])))
    w("static const double y_gain   = %s;" % fmt(float(np.atleast_1d(ymap.gain)[0])))
    w("static const double y_ymin   = %s;" % fmt(ymap.ymin))
    w("")
    w("/* the box the network was trained on */")
    w("static const double trainMin[%d] = { %s };" % (nin, ", ".join(fmt(v) for v in train_lo)))
    w("static const double trainMax[%d] = { %s };" % (nin, ", ".join(fmt(v) for v in train_hi)))
    w("")
    w("static const bool logOutput = %s;" % ("true" if log_output else "false"))
    w("")

    # ---- weights, one block per layer -------------------------------------
    for i, (Wi, Bi) in enumerate(zip(W, B)):
        r, c = Wi.shape
        wname = "IW1_1" if i == 0 else "LW%d_%d" % (i + 1, i)
        bname = "b%d" % (i + 1)
        w("/* layer %d: %d x %d */" % (i + 1, r, c))
        w("static const double %s[%d][%d] = {" % (wname, r, c))
        for row in Wi:
            w("    { %s }," % ", ".join(fmt(v) for v in row))
        w("};")
        w("static const double %s[%d] = { %s };" % (bname, r, ", ".join(fmt(v) for v in Bi)))
        w("")

    # ---- evaluation --------------------------------------------------------
    w("/* MATLAB's tansig.  Algebraically identical to std::tanh; this is the")
    w(" * form MATLAB's genFunction emits, kept so the generated file can be")
    w(" * compared line by line against a MATLAB-exported one. */")
    w("static inline double tansig(double n) { return 2.0 / (1.0 + std::exp(-2.0 * n)) - 1.0; }")
    w("")
    w("/* True when every input lies inside the trained box.  Outside it the")
    w(" * return value of eval() is not meaningful. */")
    w("static inline bool inTrainingRange(const std::vector<double>& x)")
    w("{")
    w("    if (x.size() < nInputs) return false;")
    w("    for (std::size_t i = 0; i < nInputs; ++i) {")
    w("        if (x[i] < trainMin[i] || x[i] > trainMax[i]) return false;")
    w("    }")
    w("    return true;")
    w("}")
    w("")
    w("/* Specific growth rate, 1/h.  x holds the uptake flux estimates in")
    w(" * mmol/gDW/h, POSITIVE for consumption, in CompLaB.xml substrate order. */")
    w("static inline double eval(const std::vector<double>& x)")
    w("{")
    w("    if (x.size() < nInputs) return 0.0;")
    w("")
    w("    double a0[%d];" % nin)
    w("    for (std::size_t i = 0; i < nInputs; ++i) {")
    w("        a0[i] = (x[i] - x_offset[i]) * x_gain[i] + x_ymin;")
    w("    }")
    w("")
    prev = "a0"
    for i, Wi in enumerate(W):
        r, c = Wi.shape
        cur = "a%d" % (i + 1)
        wname = "IW1_1" if i == 0 else "LW%d_%d" % (i + 1, i)
        bname = "b%d" % (i + 1)
        last = (i == len(W) - 1)
        w("    double %s[%d];" % (cur, r))
        w("    for (int r = 0; r < %d; ++r) {" % r)
        w("        double s = %s[r];" % bname)
        w("        for (int c = 0; c < %d; ++c) { s += %s[r][c] * %s[c]; }" % (c, wname, prev))
        w("        %s[r] = %s;" % (cur, "s" if last else "tansig(s)"))
        w("    }")
        w("")
        prev = cur
    w("    double g = (%s[0] - y_ymin) / y_gain + y_offset;" % prev)
    if log_output:
        w("    /* the fit was done on log10(growth); undo it */")
        w("    g = std::pow(10.0, g);")
    w("")
    w("    /* A tiny positive prediction is fit noise, not growth.  The shipped")
    w("     * network used the same 1e-8 cut-off; keeping it means a converted")
    w("     * network reproduces the old behaviour exactly. */")
    w("    if (!(g > 1e-8)) { g = 0.0; }   /* also catches NaN */")
    w("    return g;")
    w("}")
    w("")
    w("}  // namespace %s" % ns)
    w("")
    w("#endif  // %s" % guard)
    return "\n".join(L) + "\n"


# ===========================================================================
def main():
    ap = argparse.ArgumentParser(
        prog="trainSurrogate.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Train a CompLaB surrogate network and export it as C++.",
        epilog="See the header of this file for a worked example.")
    ap.add_argument("csv", help="training set from generateTrainingData.py")
    ap.add_argument("--name", default="surrogate",
                    help="identifier used for the namespace and include guard")
    ap.add_argument("--microbe-id", type=int, default=0,
                    help="global microbe index this network is for, as ordered "
                         "in CompLaB.xml; only affects the pasted snippet")
    ap.add_argument("--layers", type=int, nargs="+", default=[10, 10, 10, 10],
                    help="hidden layer sizes (default 10 10 10 10, the shipped shape)")
    ap.add_argument("--restarts", type=int, default=5)
    ap.add_argument("--maxiter", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--log-output", action="store_true",
                    help="fit log10(growth); eval() undoes it")
    ap.add_argument("-o", "--output", default=None)
    args = ap.parse_args()

    X, y, meta, header = load_csv(args.csv)
    nin = X.shape[1]
    out = args.output or ("surrogate_weights_%s.hh" % args.name)
    print("Loaded %s: %d samples, %d input(s)" % (args.csv, len(X), nin))
    if header:
        print("  columns: %s" % ", ".join(header))

    train_lo, train_hi = X.min(axis=0), X.max(axis=0)
    ylo, yhi = float(y.min()), float(y.max())
    print("  input range : %s" % ", ".join("[%.6g, %.6g]" % (a, b)
                                           for a, b in zip(train_lo, train_hi)))
    print("  growth range: [%.6g, %.6g] 1/h" % (ylo, yhi))

    if yhi <= 0.0:
        sys.exit("ERROR: every sample has zero growth.  There is nothing to "
                 "fit.  Re-run generateTrainingData.py over a range where the "
                 "organism actually grows.")

    yfit = y.copy()
    if args.log_output:
        floor = max(yhi * 1e-9, 1e-12)
        nclip = int((yfit < floor).sum())
        yfit = np.log10(np.maximum(yfit, floor))
        print("  fitting log10(growth); %d zero/near-zero sample(s) floored at %.3g"
              % (nclip, floor))

    xmap = MapMinMax(train_lo, train_hi)
    ymap = MapMinMax(np.array([yfit.min()]), np.array([yfit.max()]))
    Xs = xmap.apply(X)
    ys = ymap.apply(yfit.reshape(-1, 1)).ravel()

    sizes = [nin] + list(args.layers) + [1]
    print("Architecture: %s  (%d parameters)"
          % (" -> ".join(str(s) for s in sizes), nparams(sizes)))

    # MATLAB's dividerand default: 70 / 15 / 15.
    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(len(Xs))
    n1, n2 = int(0.70 * len(Xs)), int(0.85 * len(Xs))
    itr, iva, ite = perm[:n1], perm[n1:n2], perm[n2:]
    if len(iva) == 0 or len(ite) == 0:
        sys.exit("ERROR: only %d samples.  A 70/15/15 split needs at least ~20 "
                 "rows to mean anything; generate a bigger training set." % len(Xs))
    print("Split: %d train / %d validation / %d test" % (len(itr), len(iva), len(ite)))

    best = None
    t0 = time.time()
    for k in range(args.restarts):
        theta0 = init_theta(sizes, np.random.default_rng(args.seed + 1000 * k))
        res = minimize(loss_and_grad, theta0, args=(sizes, Xs[itr], ys[itr]),
                       jac=True, method="L-BFGS-B",
                       options=dict(maxiter=args.maxiter, maxfun=args.maxiter * 2))
        W, B = unpack(res.x, sizes)
        vloss = float(np.mean((forward(W, B, Xs[iva])[0] - ys[iva]) ** 2))
        print("  restart %d/%d: train MSE %.4e, validation MSE %.4e%s"
              % (k + 1, args.restarts, res.fun, vloss,
                 "   <- best so far" if (best is None or vloss < best[0]) else ""))
        # Early stopping, MATLAB style: keep the fit that generalises best,
        # not the one that drove the training error lowest.
        if best is None or vloss < best[0]:
            best = (vloss, res.x.copy())
    print("Trained in %.1f s" % (time.time() - t0))

    W, B = unpack(best[1], sizes)

    def predict(Xraw):
        p = ymap.reverse(forward(W, B, xmap.apply(Xraw))[0].reshape(-1, 1)).ravel()
        if args.log_output:
            p = 10.0 ** p
        p = np.where(p > 1e-8, p, 0.0)
        return p

    stats = {}
    for tag, idxs in (("train", itr), ("val", iva), ("test", ite)):
        p = predict(X[idxs])
        e = p - y[idxs]
        ss = float(np.sum((y[idxs] - y[idxs].mean()) ** 2))
        stats["rmse_" + tag] = float(np.sqrt(np.mean(e ** 2)))
        stats["maxerr_" + tag] = float(np.abs(e).max())
        stats["r2_" + tag] = 1.0 - float(e @ e) / ss if ss > 0 else float("nan")
    print("\nFit quality")
    for tag in ("train", "val", "test"):
        print("  %-5s  RMSE %.6g 1/h   max err %.6g 1/h   R^2 %.6f"
              % (tag, stats["rmse_" + tag], stats["maxerr_" + tag], stats["r2_" + tag]))

    rel = stats["rmse_test"] / max(yhi, 1e-300)
    if rel > 0.05:
        print("\nWARNING: test RMSE is %.1f%% of the largest growth rate in the "
              "training set.  That is a poor fit.  Try more samples, more "
              "restarts, --log-output if growth spans orders of magnitude, or "
              "a narrower input range." % (100 * rel), file=sys.stderr)

    code = emit_cpp(args.name, W, B, xmap, ymap, train_lo, train_hi, ylo, yhi,
                    args.log_output, meta, stats, sizes, args.microbe_id, args.csv)
    with open(out, "w") as f:
        f.write(code)
    print("\nWrote %s (%d lines)" % (out, code.count("\n")))

    # A reference table so the C++ side can be checked against this fit rather
    # than merely assumed to match it.  verifyExport.py consumes it.
    ref = os.path.splitext(out)[0] + "_reference.csv"
    rng2 = np.random.default_rng(12345)
    Xr = train_lo + rng2.random((512, nin)) * (train_hi - train_lo)
    pr = predict(Xr)
    with open(ref, "w") as f:
        f.write("# reference predictions for %s -- feed these to the compiled\n" % out)
        f.write("# header and the outputs must agree to ~1e-12 relative.\n")
        for i in range(len(Xr)):
            f.write(",".join("%.17g" % v for v in Xr[i]) + ",%.17g\n" % pr[i])
    print("Wrote %s (512 checkpoints for verifyExport.py)" % ref)

    print("\nAdd these two things to surrogateModel.hh:")
    print('\n    #include "%s"\n' % os.path.basename(out))
    print("    else if (microbeId == %d) {" % args.microbe_id)
    print("        bioR[microbeId] = srgnet_%s::eval(Fin[microbeId]);" % args.name)
    print("    }")
    print("\nand set <reaction_type>surrogate</reaction_type> plus")
    print("<enable_surrogate>true</enable_surrogate> in CompLaB.xml.")


if __name__ == "__main__":
    main()
