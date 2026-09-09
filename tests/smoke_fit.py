#!/usr/bin/env python3
"""
smoke_fit.py -- check that fit_symbolic.py runs and writes a .sym the C++ side will accept.

Deliberately tiny: a real search takes minutes, and what is being checked here is the plumbing,
not the quality of the fit.  The quality claim is in the README and rests on pipelines/B_offline_models/B2_surrogate_network/expected/ecoli_sweep.csv,
where the search recovers the dual-Monod law it was never told about to within 0.03%.
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
tmp = tempfile.mkdtemp(prefix="fitsmoke_")
rng = np.random.default_rng(0)

a = rng.uniform(0.1, 2.0, 60)
b = rng.uniform(0.1, 2.0, 60)
np.savetxt(os.path.join(tmp, "d.csv"), np.column_stack([a, b, 3.0 * a * b]),
           delimiter=",", header="A,B,growth", comments="", fmt="%.17g")

out = os.path.join(tmp, "r.sym")
r = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "fit_symbolic.py"),
                    "--data", os.path.join(tmp, "d.csv"), "--target", "growth",
                    "--out", out, "--pop", "60", "--gens", "6", "--depth", "3"],
                   capture_output=True, text=True)
if r.returncode or not os.path.exists(out):
    print(r.stdout, r.stderr); sys.exit(1)

drv = os.path.join(tmp, "d.cpp")
open(drv, "w").write('''
#include "complab3d_symbolic.hh"
#include <cstdio>
int main(int, char**argv){
    complab_sym::Program P; std::string e;
    if(!complab_sym::load(P,argv[1],&e)){ std::printf("REFUSED: %s\\n",e.c_str()); return 1; }
    std::printf("loaded, %d variable(s), %d rate(s)\\n",(int)P.vars.size(),(int)P.rates.size());
    return 0;
}''')
exe = os.path.join(tmp, "d")
c = subprocess.run(["g++", "-O1", "-std=c++11", "-I", os.path.join(ROOT, "src"),
                    "-o", exe, drv], capture_output=True, text=True)
if c.returncode:
    print(c.stderr); sys.exit(1)
v = subprocess.run([exe, out], capture_output=True, text=True)
print("  fit_symbolic.py wrote a .sym and the C++ loader accepted it:", v.stdout.strip())
sys.exit(0 if v.returncode == 0 else 1)
