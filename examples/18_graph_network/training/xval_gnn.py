#!/usr/bin/env python3
"""
xval_gnn.py -- check that complab3d_graphnet.hh and train_graphnet.py compute the same thing.

Two implementations of the same forward pass, written at different times in different languages,
are two chances to be wrong.  This runs both over the same network and the same inputs and reports
the largest disagreement.  Inputs are drawn half inside the training box and half outside it, so
the clamp is compared as well as the arithmetic.

    python tests/xval_gnn.py

Needs numpy and a C++ compiler.  Nothing else -- no Palabos.
"""
# ------------------------------------------------------------------------------------------------
#  DETERMINISM.  This block runs before numpy is imported, which is why it sits above the other
#  imports rather than among them.
#
#  A seed alone does not make a fitted result reproducible.  numpy and scipy hand their linear
#  algebra to a threaded BLAS, and a threaded reduction adds its terms in whatever order the
#  threads happen to finish in.  The last bits of a fitted constant therefore differ between two
#  runs of the same command, and in a search that ranks candidates against one another a
#  difference in the last bits decides which candidate survives.  From there the two runs diverge
#  completely: the same command with the same --seed returns a different answer, which is the
#  failure this block exists to remove.
#
#  One thread per pool fixes the reduction order and makes the run reproducible.  It is slower on
#  a large problem, so COMPLAB_THREADS is the way out, and anything other than 1 is announced on
#  stderr rather than applied quietly.
#
#  Reproducible means: same command, same seed, same machine, same library versions.  The version
#  numbers are written into every fitted file this repository produces, because a different scipy
#  can converge to a different local minimum and no environment variable can prevent that.
# ------------------------------------------------------------------------------------------------
import os as _os

COMPLAB_THREADS = _os.environ.get("COMPLAB_THREADS", "1")
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    _os.environ[_var] = COMPLAB_THREADS
if COMPLAB_THREADS != "1":
    import sys as _sys
    _sys.stderr.write(
        "note: COMPLAB_THREADS=%s. This run is faster and NOT reproducible; unset it to get\n"
        "      the same answer every time from the same seed.\n" % COMPLAB_THREADS)


def complab_versions():
    """The library versions a fitted result depends on, as one short string for a file header."""
    import platform
    out = ["python " + platform.python_version()]
    for _m in ("numpy", "scipy", "torch", "sklearn"):
        try:
            out.append("%s %s" % (_m, __import__(_m).__version__))
        except Exception:
            pass
    return ", ".join(out)


import os, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC  = os.path.join(ROOT, "src")
TOOLS = os.path.join(ROOT, "tools")
sys.path.insert(0, TOOLS)
import train_graphnet as tg

DRIVER = r'''
#include "complab3d_graphnet.hh"
#include <cstdio>
int main(int, char** argv) {
    complab_gnn::Network N; std::string err;
    if (!complab_gnn::load(N, argv[1], &err)) { std::printf("LOAD FAIL: %s\n", err.c_str()); return 1; }
    std::FILE* f = std::fopen(argv[2], "r");
    std::vector<double> c((size_t)N.nS), out;
    for (;;) {
        bool ok = true;
        for (int i = 0; i < N.nS; ++i)
            if (std::fscanf(f, "%lf", &c[(size_t)i]) != 1) { ok = false; break; }
        if (!ok) break;
        N.eval(c, out, true, 0);
        for (size_t i = 0; i < out.size(); ++i)
            std::printf("%.17g%c", out[i], i + 1 == out.size() ? '\n' : ' ');
    }
    std::fclose(f);
    return 0;
}
'''


def main():
    tmp = tempfile.mkdtemp(prefix="gnnxval_")
    rng = np.random.default_rng(7)

    # a stoichiometry with real zeros, so the "no edge" shortcut is exercised
    nS, nR = 5, 4
    S = np.zeros((nS, nR))
    S[0, 0], S[1, 0], S[2, 0] = -1, -2, 1
    S[2, 1], S[3, 1] = -1, 1.5
    S[0, 2], S[4, 2] = -0.5, 1
    S[3, 3], S[4, 3] = -1, 2
    np.savetxt(os.path.join(tmp, "stoich.csv"), S, delimiter=",", fmt="%.17g")

    sp = ["A", "B", "C", "D", "E"]
    N = 40
    X = rng.uniform(0.1, 2.0, size=(N, nS))
    Y = rng.normal(0.0, 1.0, size=(N, nS + 1))
    hdr = ",".join(sp + [s + "_rate" for s in sp] + ["growth"])
    np.savetxt(os.path.join(tmp, "samples.csv"), np.hstack([X, Y]),
               delimiter=",", header=hdr, comments="", fmt="%.17g")

    # BOTH readouts are checked, not only the default.
    #
    # The species readout is what every .gnn written before the extent readout existed means, and
    # those files must keep evaluating to exactly what they always did. A cross-check that only
    # covers the new mode would let the old one drift silently, and the drift would surface as
    # someone's year-old network quietly returning different numbers.
    gnns = {}
    for mode in ("extent", "species"):
        g = os.path.join(tmp, "net_%s.gnn" % mode)
        rc = subprocess.run([sys.executable, os.path.join(TOOLS, "train_graphnet.py"),
                             "--stoich", os.path.join(tmp, "stoich.csv"),
                             "--data", os.path.join(tmp, "samples.csv"),
                             "--out", g, "--width", "4", "--rounds", "2",
                             "--epochs", "0", "--units", "per_hour",
                             "--readout", mode,
                             "--da", "1.0,0.35,2.5,0.7"],
                            capture_output=True, text=True)
        if rc.returncode:
            print(rc.stdout, rc.stderr); return 1
        gnns[mode] = g

    src = os.path.join(tmp, "driver.cpp")
    open(src, "w").write(DRIVER)
    exe = os.path.join(tmp, "driver")
    cc = subprocess.run(["g++", "-O2", "-Wall", "-Wextra", "-std=c++11",
                         "-I", SRC, "-o", exe, src], capture_output=True, text=True)
    if cc.returncode:
        print(cc.stderr); return 1
    if cc.stderr.strip():
        print("compiler had something to say:\n" + cc.stderr)

    ok = True
    for mode in ("extent", "species"):
        gnn = gnns[mode]
        net, M = tg.read_gnn(gnn)
        # half inside the training box, half well outside it
        C = np.vstack([rng.uniform(M["trainmin"], M["trainmax"], size=(30, nS)),
                       rng.uniform(-1.0, 4.0, size=(20, nS))])
        inp = os.path.join(tmp, "inputs_%s.txt" % mode)
        np.savetxt(inp, C, fmt="%.17g")

        Ypy = tg.predict(net, M, C)
        run = subprocess.run([exe, gnn, inp], capture_output=True, text=True)
        if run.returncode or not run.stdout.strip():
            print(run.stdout, run.stderr); return 1
        Ycc = np.array([[float(v) for v in ln.split()] for ln in run.stdout.strip().split("\n")])

        if Ycc.shape != Ypy.shape:
            print("shape mismatch: C++ %s, python %s" % (Ycc.shape, Ypy.shape)); return 1

        d = np.abs(Ycc - Ypy)
        scale = max(np.abs(Ypy).max(), 1e-300)
        print("  readout %s: %d inputs x %d outputs" % ((mode,) + Ycc.shape))
        print("    typical output magnitude   %.6e" % np.abs(Ypy).mean())
        print("    largest absolute disagreement %.6e" % d.max())
        print("    as a fraction of the largest output %.3e" % (d.max() / scale))

        if mode == "extent":
            # The claim of this readout, checked on the C++ side rather than the Python one: the
            # rates it returns must lie in the column space of S, at every input including the
            # ones outside the training box.
            Sp = np.linalg.pinv(S)
            R = Ycc[:, :nS]
            off = R - (R @ Sp.T) @ S.T
            den = float(np.sum(R ** 2))
            dev = float(np.sqrt(np.sum(off ** 2) / den)) if den > 0 else 0.0
            print("    off the stoichiometric subspace: %.3e of the rate norm" % dev)
            if dev > 1e-12:
                print("    ** extent mode is not stoichiometrically exact in the C++ evaluator **")
                ok = False

        if d.max() / scale >= 1e-12:
            ok = False

    print("\n%s" % ("the two implementations agree to machine precision, in both readouts"
                    if ok else "** THE TWO IMPLEMENTATIONS DISAGREE **"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
