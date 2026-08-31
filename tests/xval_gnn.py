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

    gnn = os.path.join(tmp, "net.gnn")
    rc = subprocess.run([sys.executable, os.path.join(TOOLS, "train_graphnet.py"),
                         "--stoich", os.path.join(tmp, "stoich.csv"),
                         "--data", os.path.join(tmp, "samples.csv"),
                         "--out", gnn, "--width", "4", "--rounds", "2",
                         "--epochs", "0", "--units", "per_hour",
                         "--da", "1.0,0.35,2.5,0.7"],
                        capture_output=True, text=True)
    if rc.returncode:
        print(rc.stdout, rc.stderr); return 1

    src = os.path.join(tmp, "driver.cpp")
    open(src, "w").write(DRIVER)
    exe = os.path.join(tmp, "driver")
    cc = subprocess.run(["g++", "-O2", "-Wall", "-Wextra", "-std=c++11",
                         "-I", SRC, "-o", exe, src], capture_output=True, text=True)
    if cc.returncode:
        print(cc.stderr); return 1
    if cc.stderr.strip():
        print("compiler had something to say:\n" + cc.stderr)

    net, M = tg.read_gnn(gnn)
    # half inside the training box, half well outside it
    C = np.vstack([rng.uniform(M["trainmin"], M["trainmax"], size=(30, nS)),
                   rng.uniform(-1.0, 4.0, size=(20, nS))])
    inp = os.path.join(tmp, "inputs.txt")
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
    print("  %d inputs x %d outputs" % Ycc.shape)
    print("  typical output magnitude   %.6e" % np.abs(Ypy).mean())
    print("  largest absolute disagreement %.6e" % d.max())
    print("  as a fraction of the largest output %.3e" % (d.max() / scale))

    ok = d.max() / scale < 1e-12
    print("\n%s" % ("the two implementations agree to machine precision"
                    if ok else "** THE TWO IMPLEMENTATIONS DISAGREE **"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
