#!/usr/bin/env python3
"""
===============================================================================
INSPECT AN EXISTING SURROGATE NETWORK
===============================================================================

Recovers, from a surrogateModel.hh (or a generated surrogate_weights_*.hh),
the facts about a network that nobody wrote down:

  * its architecture,
  * the RANGE OF INPUTS IT WAS TRAINED ON, and
  * the range of growth rates it can produce.

The training range is not stored anywhere as such -- but mapminmax IS
invertible, so it can be recovered exactly.  MATLAB's mapminmax stores

      gain = 2 / (x_max - x_min),      offset = x_min

so the box the network saw is

      x_min = offset,      x_max = offset + 2 / gain

That matters because a network outside its training box does not fail, it
lies.  Before trusting a network somebody else fitted -- including the one
CompLaB has shipped since v1.0 -- run this and check that your simulation
stays inside the box it prints.

    python3 inspectSurrogate.py ../surrogateModel.hh
    python3 inspectSurrogate.py surrogate_weights_geobacter.hh --eval 2.0 0.1

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
===============================================================================
"""

import argparse
import re
import sys

try:
    import numpy as np
except ImportError:
    sys.exit("ERROR: numpy is required.  python3 -m pip install numpy")

NUM = r'-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?'


def parse(path):
    t = open(path).read()

    def vec(name):
        # The "=" must follow the name directly (allowing an array bound), or a
        # mention of the same name inside a COMMENT matches first and silently
        # returns the wrong block. That is not hypothetical: the generated
        # header documents the scaling as "(x - x_offset) * x_gain + x_ymin"
        # a few lines above the real declarations.
        m = re.search(r'\b' + name + r'\b\s*(?:\[\s*\d+\s*\])?\s*=\s*\{(.*?)\};', t, re.S)
        if not m:
            raise KeyError(name)
        return np.array([float(x) for x in re.findall(NUM, m.group(1))])

    def mat(name):
        decl = r'\b' + name + r'\b\s*(?:\[\s*\d+\s*\]){0,2}\s*=\s*'
        m = re.search(decl + r'\{(.*?)\n\s*\};', t, re.S)
        if not m:
            m = re.search(decl + r'\{(.*?)\};', t, re.S)
        if not m:
            raise KeyError(name)
        rows = re.findall(r'\{([^{}]*)\}', m.group(1))
        if not rows:
            # A single-row "matrix" written as a plain vector, which is how the
            # shipped surrogateModel.hh stores its output layer (LW5_4).
            return np.array([[float(x) for x in re.findall(NUM, m.group(1))]])
        return np.array([[float(x) for x in re.findall(NUM, r)] for r in rows])

    def scal(name, default=None):
        m = re.search(r'\b' + name + r'\b\s*=\s*(' + NUM + r')', t)
        if not m:
            if default is None:
                raise KeyError(name)
            return default
        return float(m.group(1))

    W, B = [mat('IW1_1')], [vec('b1')]
    i = 2
    while True:
        try:
            Wi = mat('LW%d_%d' % (i, i - 1))
        except KeyError:
            break
        # The shipped surrogateModel.hh writes the OUTPUT layer as a plain
        # vector LW5_4 with a SCALAR bias b5, not as a 1 x n matrix with a
        # 1-element bias array.  Both spellings have to be accepted, or this
        # tool cannot read the very file it exists to inspect.
        try:
            Bi = vec('b%d' % i)
        except KeyError:
            Bi = np.array([scal('b%d' % i)])
        W.append(np.atleast_2d(Wi))
        B.append(np.atleast_1d(Bi))
        i += 1

    nin = W[0].shape[1]
    try:
        xoff = vec('x_offset')
        xgain = vec('x_gain')
    except KeyError:
        xoff = np.array([scal('x_offset%d' % k) for k in range(nin)])
        xgain = np.array([scal('x_gain%d' % k) for k in range(nin)])

    return dict(W=W, B=B, xoff=xoff, xgain=xgain,
                xymin=scal('x_ymin', -1.0), yymin=scal('y_ymin', -1.0),
                ygain=scal('y_gain'), yoff=scal('y_offset', 0.0),
                logout=('logOutput = true' in t))


def evaluate(p, x):
    a = (np.asarray(x, float) - p['xoff']) * p['xgain'] + p['xymin']
    for W, b in zip(p['W'][:-1], p['B'][:-1]):
        a = 2.0 / (1.0 + np.exp(-2.0 * (W @ a + b))) - 1.0
    a = float((p['W'][-1] @ a + p['B'][-1]).ravel()[0])
    g = (a - p['yymin']) / p['ygain'] + p['yoff']
    if p['logout']:
        g = 10.0 ** g
    return 0.0 if not (g > 1e-8) else g


def main():
    ap = argparse.ArgumentParser(prog="inspectSurrogate.py")
    ap.add_argument("header")
    ap.add_argument("--eval", nargs="+", type=float, metavar="V",
                    help="evaluate at these uptake rates, mmol/gDW/h")
    args = ap.parse_args()

    p = parse(args.header)
    nin = p['W'][0].shape[1]
    sizes = [nin] + [w.shape[0] for w in p['W']]
    npar = sum(w.size + b.size for w, b in zip(p['W'], p['B']))

    print("File         : %s" % args.header)
    print("Architecture : %s   (%d parameters)"
          % (" -> ".join(str(s) for s in sizes), npar))
    print("Output       : %s" % ("log10(growth), undone on evaluation"
                                 if p['logout'] else "growth rate, 1/h"))
    print("\nRECOVERED TRAINING RANGE  (mapminmax inverted; the network is only")
    print("valid inside this box)")
    lo = np.empty(nin)
    hi = np.empty(nin)
    for k in range(nin):
        lo[k] = p['xoff'][k]
        hi[k] = p['xoff'][k] + 2.0 / p['xgain'][k]
        print("  input %d : %.10g .. %.10g mmol/gDW/h" % (k, lo[k], hi[k]))
    ylo = p['yoff']
    yhi = p['yoff'] + 2.0 / p['ygain']
    if p['logout']:
        print("  output  : %.10g .. %.10g 1/h  (10^ of the fitted range)"
              % (10.0 ** ylo, 10.0 ** yhi))
    else:
        print("  output  : %.10g .. %.10g 1/h" % (ylo, yhi))

    if args.eval:
        if len(args.eval) != nin:
            sys.exit("ERROR: this network takes %d input(s); %d given."
                     % (nin, len(args.eval)))
        g = evaluate(p, args.eval)
        inbox = all(lo[k] <= args.eval[k] <= hi[k] for k in range(nin))
        print("\nAt %s:" % ", ".join("%g" % v for v in args.eval))
        print("  growth rate = %.10g 1/h" % g)
        if not inbox:
            print("  WARNING: OUTSIDE the training box. This number is an "
                  "extrapolation and should not be trusted.")
    else:
        # A quick sanity sweep along each axis, others held at their midpoint.
        print("\nResponse along each axis (others at their midpoint)")
        mid = 0.5 * (lo + hi)
        for k in range(nin):
            xs = np.linspace(lo[k], hi[k], 7)
            vals = []
            for v in xs:
                x = mid.copy()
                x[k] = v
                vals.append(evaluate(p, x))
            print("  input %d: %s" % (k, "  ".join("%.4g" % v for v in vals)))
        print("  (input values: %s)"
              % "  ".join("%.3g" % v for v in np.linspace(lo[0], hi[0], 7)))


if __name__ == "__main__":
    main()
