#!/usr/bin/env python3
"""
smoke_fit.py -- check that fit_symbolic.py runs and writes a .sym the C++ side will accept.

Deliberately tiny: a real search takes minutes, and what is being checked here is the plumbing,
not the quality of the fit.  The quality claim is in the README and rests on examples/ecoli_sweep.csv,
where the search recovers the dual-Monod law it was never told about to within 0.03%.
"""
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
