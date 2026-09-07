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
    out = W[-1] @ A[-1] + B[-1][:, None]        # linear output layer, (n_out, n_samples)
    A.append(out)
    return out, A


def loss_and_grad(theta, sizes, X, Y):
    """Mean squared error and its analytic gradient.

    tansig(n) = 2/(1+exp(-2n)) - 1 is algebraically identical to tanh(n), so
    np.tanh is used for the forward pass (it is the numerically stable form)
    and the derivative is the familiar 1 - a**2.  The generated C++ writes the
    2/(1+exp(-2n))-1 form because that is what MATLAB emits and what
    surrogateModel.hh already contains -- same function, same values."""
    W, B = unpack(theta, sizes)
    out, A = forward(W, B, X)                       # (n_out, n_samples)
    n_out, n = out.shape
    err = out - Y                                   # Y is (n_out, n_samples)
    L = float((err * err).sum()) / (n * n_out)

    gW = [None] * len(W)
    gB = [None] * len(B)
    # Every output contributes to the same shared hidden layers, so the
    # gradient is simply summed over outputs -- which is what the matrix
    # product below does. No per-output weighting: the targets have already
    # been mapped onto -1..1 individually, so a flux in mmol/gDW/h and a growth
    # rate in 1/h carry equal weight in the loss. That is deliberate. Weighting
    # them by magnitude instead would let the largest flux dominate the fit.
    D = (2.0 / (n * n_out)) * err                   # dL/d(pre-activation), output layer
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
    meta = {"inputs": [], "outputs": [], "model": None, "objective": None,
            "n_inputs": None}
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
                elif body.startswith("n_inputs:"):
                    meta["n_inputs"] = int(body.split(":", 1)[1])
                elif body.startswith("input "):
                    meta["inputs"].append(body)
                elif body.startswith("output "):
                    meta["outputs"].append(body)
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

    # HOW MANY COLUMNS ARE INPUTS.  generateTrainingData.py writes "# n_inputs:"
    # so this is exact.  A CSV without it is assumed to be the old single-output
    # form -- everything but the last column is an input -- which keeps every
    # training set written before this change loading unchanged.
    nin = meta["n_inputs"]
    if nin is None:
        nin = D.shape[1] - 1
        meta["n_inputs"] = nin
    if nin >= D.shape[1]:
        sys.exit("ERROR: %s says n_inputs %d but has only %d column(s)."
                 % (path, nin, D.shape[1]))

    # Name the outputs, from the header row if there is one, else positionally.
    nout = D.shape[1] - nin
    if header is not None and len(header) == D.shape[1]:
        out_names = [h.strip() for h in header[nin:]]
    else:
        out_names = ["growth"] + ["out_%d" % k for k in range(1, nout)]
    meta["out_names"] = out_names

    return D[:, :nin], D[:, nin:], meta, header


# ===========================================================================
# The C++ emitter.  This is the part that has to be exactly right: the
# generated header is compiled into CompLB3D and evaluated millions of times.
# ===========================================================================
def fmt(v):
    """17 significant digits round-trips an IEEE double exactly.  Anything
    shorter silently changes the network."""
    return "%.17g" % float(v)


def emit_cpp(name, W, B, xmap, ymap, train_lo, train_hi, ylo, yhi,
             log_output, meta, stats, sizes, microbe_id, source_csv,
             out_names=None):
    guard = "SURROGATE_WEIGHTS_%s_HH" % name.upper()
    ns = "srgnet_%s" % name
    nin = sizes[0]
    nout = sizes[-1]
    if not out_names:
        out_names = ["growth"] + ["out_%d" % k for k in range(1, nout)]
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
    w(" * Architecture : %d input(s) -> %s -> %d linear output(s)"
      % (nin, " -> ".join("%d tansig" % h for h in sizes[1:-1]), nout))
    w(" * Outputs      : %s" % ", ".join(out_names))
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
    if nout > 1 and stats.get("per_out"):
        w(" *")
        w(" * PER-OUTPUT TEST ERROR. Read this before the headline number above:")
        w(" * the flux columns are harder to fit than growth, because they are")
        w(" * piecewise linear with kinks where the binding constraint changes.")
        for nm in out_names:
            d = stats["per_out"].get(nm)
            if d:
                w(" *   %-20s RMSE %-12.6g span %-12.6g R^2 %.5f"
                  % (nm, d["rmse"], d["span"], d["r2"]))
    w(" *")
    w(" * HOW TO USE IT -- in surrogateModel.hh:")
    w(" *")
    w(' *     #include "%s"' % ("surrogate_weights_%s.hh" % name))
    w(" *")
    if nout == 1:
        w(" *     else if (microbeId == %d) {" % microbe_id)
        w(" *         bioR[microbeId] = %s::eval(Fin[microbeId]);" % ns)
        w(" *     }")
    else:
        w(" *     else if (microbeId == %d) {" % microbe_id)
        w(" *         double o[%s::nOutputs];" % ns)
        w(" *         %s::evalAll(Fin[microbeId], o);" % ns)
        w(" *         bioR[microbeId] = o[0];               // %s" % out_names[0])
        for k, nm in enumerate(out_names[1:], start=1):
            w(" *         Fout[microbeId][%d] = o[%d];           // %s" % (k - 1, k, nm))
        w(" *     }")
        w(" *")
        w(" *     THE INDEX MAPPING IS YOURS TO GET RIGHT. Output k above is the")
        w(" *     flux of %s," % ", ".join(out_names[1:]))
        w(" *     in that order. Fout[microbeId][j] is substrate j as ordered in")
        w(" *     CompLaB.xml. Those two orders are only the same if you swept the")
        w(" *     exchanges in substrate order. Check it once, in the XML, rather")
        w(" *     than assuming it: a transposed pair here is silent and fatal.")
    w(" *")
    w(" * Fin is in mmol/gDW/h and POSITIVE means consumption; output 0 is the")
    w(" * specific growth rate in 1/h and the rest are fluxes in mmol/gDW/h,")
    w(" * positive = consumed. Those are CompLaB's conventions already, so no")
    w(" * conversion is needed on either side.")
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
    w("static const std::size_t nInputs  = %d;" % nin)
    w("static const std::size_t nOutputs = %d;" % nout)
    w("")
    w("/* what each output IS, in order. Kept as text so the file is readable")
    w(" * on its own and a mis-wired index is visible in a diff. */")
    w("static const char* const outputNames[%d] = { %s };"
      % (nout, ", ".join('"%s"' % n for n in out_names)))
    w("")
    w("/* mapminmax input scaling: xs = (x - x_offset) * x_gain + x_ymin */")
    w("static const double x_offset[%d] = { %s };" % (nin, ", ".join(fmt(v) for v in xmap.offset)))
    w("static const double x_gain  [%d] = { %s };" % (nin, ", ".join(fmt(v) for v in xmap.gain)))
    w("static const double x_ymin      = %s;" % fmt(xmap.ymin))
    w("")
    w("/* mapminmax output scaling, applied in reverse: y = (ys - y_ymin)/y_gain + y_offset")
    w(" * ONE ENTRY PER OUTPUT. Each output was mapped onto -1..1 separately, so")
    w(" * a growth rate in 1/h and a flux in mmol/gDW/h carry equal weight in the")
    w(" * loss instead of the larger one dominating. */")
    w("static const double y_offset[%d] = { %s };"
      % (nout, ", ".join(fmt(v) for v in np.atleast_1d(ymap.offset))))
    w("static const double y_gain  [%d] = { %s };"
      % (nout, ", ".join(fmt(v) for v in np.atleast_1d(ymap.gain))))
    w("static const double y_ymin      = %s;" % fmt(ymap.ymin))
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
    # ---- the forward pass, written once and shared by eval/evalAll ---------
    w("/* The forward pass. Hidden layers are tansig, the output layer linear.")
    w(" * out[] receives nOutputs values in PHYSICAL units. */")
    w("static inline void evalAll(const std::vector<double>& x, double* out)")
    w("{")
    w("    for (std::size_t k = 0; k < nOutputs; ++k) out[k] = 0.0;")
    w("    if (x.size() < nInputs) return;")
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
    w("    for (std::size_t k = 0; k < nOutputs; ++k) {")
    w("        out[k] = (%s[k] - y_ymin) / y_gain[k] + y_offset[k];" % prev)
    w("    }")
    if log_output:
        w("")
        w("    /* the fit was done on log10(growth); undo it. ONLY output 0: a flux")
        w("     * may legitimately be zero or negative and was never logged. */")
        w("    out[0] = std::pow(10.0, out[0]);")
    w("")
    w("    /* A tiny positive prediction is fit noise, not growth. The shipped")
    w("     * network used the same 1e-8 cut-off; keeping it means a converted")
    w("     * network reproduces the old behaviour exactly. The cut-off applies")
    w("     * to GROWTH ONLY -- clamping a flux would corrupt it. */")
    w("    if (!(out[0] > 1e-8)) { out[0] = 0.0; }   /* also catches NaN */")
    if nout > 1:
        w("    for (std::size_t k = 1; k < nOutputs; ++k) {")
        w("        if (out[k] != out[k]) { out[k] = 0.0; }   /* NaN only */")
        w("    }")
    w("}")
    w("")
    w("/* Growth alone. Unchanged in meaning from the single-output networks, so")
    w(" * a caller that only wants growth needs no edit. */")
    w("static inline double eval(const std::vector<double>& x)")
    w("{")
    w("    double o[nOutputs];")
    w("    evalAll(x, o);")
    w("    return o[0];")
    w("}")
    w("")
    w("}  // namespace %s" % ns)
    w("")
    w("#endif  // %s" % guard)
    return "\n".join(L) + "\n"


def emit_srg(path, name, W, B, xmap, ymap, train_lo, train_hi, Y,
             log_output, meta, out_names, in_names):
    """The same fitted network in the run-time format.

    WHY BOTH.  The .hh is compiled into the solver and is the fastest thing
    there is; the .srg is read at start-up, so swapping network is an edit to
    CompLaB.xml rather than a rebuild.  They hold identical numbers -- this
    writes them from the same fit in the same pass, precisely so they cannot
    drift -- and tests/test_surrogate_parity.cpp checks that the two evaluators
    agree on them.

    Without this the run-time path could only ever hold a growth-only network,
    because the in-run trainer inside the solver fits growth alone.  A
    multi-output .srg has to come from here.
    """
    nin = len(np.atleast_1d(xmap.offset))
    nout = len(np.atleast_1d(ymap.offset))
    sizes = [nin] + [w.shape[0] for w in W]
    L = []
    a = L.append
    a("# CompLB3D surrogate network")
    a("# fitted by trainSurrogate.py from %s%s"
      % (os.path.basename(meta.get("source", "")),
         (", model %s" % meta["model"]) if meta.get("model") else ""))
    a("# THE NETWORK IS ONLY VALID INSIDE trainMin..trainMax. Outside that box a")
    a("# neural network does not fail, it returns confident nonsense.")
    a("version 1")
    a("inputs %d" % nin)
    a("layers " + " ".join(str(v) for v in sizes))
    a("logoutput %d" % (1 if log_output else 0))
    if in_names:
        a("inputnames " + " ".join(in_names))
    a("xoffset " + " ".join("%.17g" % v for v in np.atleast_1d(xmap.offset)))
    a("xgain "   + " ".join("%.17g" % v for v in np.atleast_1d(xmap.gain)))
    a("xymin %.17g" % xmap.ymin)
    a("outputs %d" % nout)
    a("outputnames " + " ".join(out_names))
    a("yoffset " + " ".join("%.17g" % v for v in np.atleast_1d(ymap.offset)))
    a("ygain "   + " ".join("%.17g" % v for v in np.atleast_1d(ymap.gain)))
    a("yymin %.17g" % ymap.ymin)
    a("trainmin " + " ".join("%.17g" % v for v in train_lo))
    a("trainmax " + " ".join("%.17g" % v for v in train_hi))
    a("outrange " + " ".join("%.17g %.17g" % (Y[:, k].min(), Y[:, k].max())
                             for k in range(nout)))
    for i, (Wi, Bi) in enumerate(zip(W, B)):
        a("W %d" % i)
        for r in range(Wi.shape[0]):
            a(" ".join("%.17g" % v for v in Wi[r]))
        a("B %d" % i)
        a(" ".join("%.17g" % v for v in Bi))
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    return path


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
    ap.add_argument("--srg", default=None, metavar="PATH",
                    help="ALSO write the network in the run-time .srg format. "
                         "Point <weights_file> at it and the solver reads it at "
                         "start-up, so switching network needs no rebuild. Both "
                         "files come from the same fit and hold the same numbers.")
    ap.add_argument("-o", "--output", default=None)
    args = ap.parse_args()

    X, Y, meta, header = load_csv(args.csv)
    nin = X.shape[1]
    nout = Y.shape[1]
    out_names = meta.get("out_names") or ["growth"]
    y = Y[:, 0]                                   # column 0 is always growth
    out = args.output or ("surrogate_weights_%s.hh" % args.name)
    print("Loaded %s: %d samples, %d input(s), %d output(s)"
          % (args.csv, len(X), nin, nout))
    if header:
        print("  columns: %s" % ", ".join(header))

    train_lo, train_hi = X.min(axis=0), X.max(axis=0)
    ylo, yhi = float(y.min()), float(y.max())
    print("  input range : %s" % ", ".join("[%.6g, %.6g]" % (a, b)
                                           for a, b in zip(train_lo, train_hi)))
    for k, nm in enumerate(out_names):
        print("  output %-16s [%+.6g, %+.6g]" % (nm, Y[:, k].min(), Y[:, k].max()))

    if yhi <= 0.0:
        sys.exit("ERROR: every sample has zero growth.  There is nothing to "
                 "fit.  Re-run generateTrainingData.py over a range where the "
                 "organism actually grows.")

    # A column that never moves cannot be fitted and cannot be useful: mapminmax
    # would divide by a zero span, and the network would spend capacity learning
    # a constant. Say so and stop, rather than producing a header whose extra
    # outputs are silently meaningless.
    flat = [nm for k, nm in enumerate(out_names)
            if k > 0 and (Y[:, k].max() - Y[:, k].min()) <= 0.0]
    if flat:
        sys.exit("ERROR: output column(s) %s are constant across every sample.\n"
                 "  Nothing can be learned from them. Re-run the sweep without "
                 "them, or over a range where they vary." % ", ".join(flat))

    Yfit = Y.astype(float).copy()
    if args.log_output:
        floor = max(yhi * 1e-9, 1e-12)
        nclip = int((Yfit[:, 0] < floor).sum())
        # ONLY the growth column is logged. A flux column can legitimately be
        # zero or negative, and log10 of either is not a number.
        Yfit[:, 0] = np.log10(np.maximum(Yfit[:, 0], floor))
        print("  fitting log10(growth); %d zero/near-zero sample(s) floored at %.3g"
              % (nclip, floor))

    xmap = MapMinMax(train_lo, train_hi)
    ymap = MapMinMax(Yfit.min(axis=0), Yfit.max(axis=0))
    Xs = xmap.apply(X)
    ys = ymap.apply(Yfit).T                       # (n_out, n_samples)

    sizes = [nin] + list(args.layers) + [nout]
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
        res = minimize(loss_and_grad, theta0, args=(sizes, Xs[itr], ys[:, itr]),
                       jac=True, method="L-BFGS-B",
                       options=dict(maxiter=args.maxiter, maxfun=args.maxiter * 2))
        W, B = unpack(res.x, sizes)
        vloss = float(np.mean((forward(W, B, Xs[iva])[0] - ys[:, iva]) ** 2))
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
        """Returns (n_samples, n_out) in physical units."""
        p = ymap.reverse(forward(W, B, xmap.apply(Xraw))[0].T)
        if args.log_output:
            p[:, 0] = 10.0 ** p[:, 0]
        # The growth floor applies to growth only. A flux may legitimately be
        # small, zero, or negative, and clamping one would corrupt it.
        p[:, 0] = np.where(p[:, 0] > 1e-8, p[:, 0], 0.0)
        return p

    stats = {}
    per_out = {}
    for tag, idxs in (("train", itr), ("val", iva), ("test", ite)):
        P = predict(X[idxs])
        E = P - Y[idxs]
        # headline numbers stay the growth column, so an old header and a new
        # one report the same thing in the same place
        e = E[:, 0]
        ss = float(np.sum((y[idxs] - y[idxs].mean()) ** 2))
        stats["rmse_" + tag] = float(np.sqrt(np.mean(e ** 2)))
        stats["maxerr_" + tag] = float(np.abs(e).max())
        stats["r2_" + tag] = 1.0 - float(e @ e) / ss if ss > 0 else float("nan")
        if tag == "test":
            for k, nm in enumerate(out_names):
                ek = E[:, k]
                ssk = float(np.sum((Y[idxs, k] - Y[idxs, k].mean()) ** 2))
                per_out[nm] = dict(
                    rmse=float(np.sqrt(np.mean(ek ** 2))),
                    maxerr=float(np.abs(ek).max()),
                    r2=(1.0 - float(ek @ ek) / ssk) if ssk > 0 else float("nan"),
                    span=float(Y[:, k].max() - Y[:, k].min()))
    stats["per_out"] = per_out

    print("\nFit quality (growth column)")
    for tag in ("train", "val", "test"):
        print("  %-5s  RMSE %.6g 1/h   max err %.6g 1/h   R^2 %.6f"
              % (tag, stats["rmse_" + tag], stats["maxerr_" + tag], stats["r2_" + tag]))

    if nout > 1:
        # PER OUTPUT, ON THE TEST SET. A single overall error hides the one
        # column that did not fit, and the flux columns are much harder to fit
        # than growth: they are piecewise-linear with kinks where the binding
        # constraint changes, and a network smooths exactly those.
        print("\nPer-output test error")
        print("  %-20s %12s %12s %9s" % ("output", "RMSE", "span", "R^2"))
        for nm in out_names:
            d = per_out[nm]
            print("  %-20s %12.6g %12.6g %9.5f" % (nm, d["rmse"], d["span"], d["r2"]))
        bad = [nm for nm in out_names
               if per_out[nm]["span"] > 0 and
                  per_out[nm]["rmse"] / per_out[nm]["span"] > 0.05]
        if bad:
            print("\nWARNING: %s fitted to worse than 5%% of their own span. Those "
                  "columns are the weak link, not the growth column. More nodes "
                  "and more samples help; so does dropping a column you do not "
                  "actually need." % ", ".join(bad), file=sys.stderr)

    rel = stats["rmse_test"] / max(yhi, 1e-300)
    if rel > 0.05:
        print("\nWARNING: test RMSE is %.1f%% of the largest growth rate in the "
              "training set.  That is a poor fit.  Try more samples, more "
              "restarts, --log-output if growth spans orders of magnitude, or "
              "a narrower input range." % (100 * rel), file=sys.stderr)

    code = emit_cpp(args.name, W, B, xmap, ymap, train_lo, train_hi, ylo, yhi,
                    args.log_output, meta, stats, sizes, args.microbe_id, args.csv,
                    out_names)
    with open(out, "w") as f:
        f.write(code)
    print("\nWrote %s (%d lines)" % (out, code.count("\n")))

    if args.srg:
        meta["source"] = args.csv
        in_names = []
        for line in meta.get("inputs", []):
            # "input 0: EX_glc__D_e (index 27) range ..." -> the reaction id
            parts = line.split(":", 1)
            if len(parts) == 2:
                in_names.append(parts[1].strip().split()[0])
        if len(in_names) != nin:
            in_names = []
        emit_srg(args.srg, args.name, W, B, xmap, ymap, train_lo, train_hi, Y,
                 args.log_output, meta, out_names, in_names)
        print("Wrote %s (the same network, for <weights_file>)" % args.srg)

    # A reference table so the C++ side can be checked against this fit rather
    # than merely assumed to match it.  verifyExport.py consumes it.
    ref = os.path.splitext(out)[0] + "_reference.csv"
    rng2 = np.random.default_rng(12345)
    Xr = train_lo + rng2.random((512, nin)) * (train_hi - train_lo)
    Pr = predict(Xr)
    with open(ref, "w") as f:
        f.write("# reference predictions for %s -- feed these to the compiled\n" % out)
        f.write("# header and the outputs must agree to ~1e-12 relative.\n")
        f.write("# n_inputs: %d\n" % nin)
        f.write("# n_outputs: %d\n" % nout)
        f.write("# outputs: %s\n" % " ".join(out_names))
        for i in range(len(Xr)):
            f.write(",".join("%.17g" % v for v in Xr[i]) + ","
                    + ",".join("%.17g" % v for v in Pr[i]) + "\n")
    print("Wrote %s (512 checkpoints for verifyExport.py)" % ref)

    ns = "srgnet_%s" % args.name
    print("\nAdd these two things to surrogateModel.hh:")
    print('\n    #include "%s"\n' % os.path.basename(out))
    print("    else if (microbeId == %d) {" % args.microbe_id)
    if nout == 1:
        print("        bioR[microbeId] = %s::eval(Fin[microbeId]);" % ns)
        print("    }")
    else:
        print("        double o[%s::nOutputs];" % ns)
        print("        %s::evalAll(Fin[microbeId], o);" % ns)
        print("        bioR[microbeId] = o[0];               // %s" % out_names[0])
        for k, nm in enumerate(out_names[1:], start=1):
            print("        Fout[microbeId][%d] = o[%d];           // %s" % (k - 1, k, nm))
        print("    }")
        print("\nCHECK THE INDEX MAPPING BEFORE YOU RUN IT. Fout[microbeId][j] is")
        print("substrate j in <name_of_substrates> order. The outputs above are in")
        print("the order the exchanges were swept:")
        for k, nm in enumerate(out_names[1:], start=1):
            print("    o[%d] -> %s" % (k, nm))
        print("Those two orders agree only if you swept in substrate order. A")
        print("transposed pair is silent: both numbers stay plausible and every")
        print("result is wrong.")
    print("\nand set <reaction_type>surrogate</reaction_type> plus")
    print("<enable_surrogate>true</enable_surrogate> in CompLaB.xml.")
    if nout > 1:
        print("\nWhat changes at run time: the solver now takes the substrate")
        print("consumption from the network instead of falling back to a Monod")
        print("guess. Where the swept bound was not the binding constraint -- most")
        print("of the domain -- those two differ, and the Monod guess was the one")
        print("that was wrong.")


if __name__ == "__main__":
    main()
