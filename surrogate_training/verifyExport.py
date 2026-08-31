#!/usr/bin/env python3
"""
===============================================================================
VERIFY THAT THE GENERATED C++ AGREES WITH THE TRAINER
===============================================================================

An exporter that silently drops a digit, transposes a weight matrix, or gets
the mapminmax direction backwards produces a header that compiles, runs, and
is wrong.  Nothing downstream would catch it: the growth rates would simply be
different numbers, and there is no independent reference to compare them to.

So the trainer writes 512 reference predictions alongside the header, and this
script compiles the header, evaluates it at those same points, and compares.
Run it after every export.  It needs a C++ compiler and nothing else.

    python3 verifyExport.py surrogate_weights_geobacter.hh

It exits non-zero if the two disagree by more than 1e-10 relative, which is
far above the ~1e-16 that differing floating-point summation orders explain
and far below any error that would matter biologically.
===============================================================================
"""

import os
import subprocess
import sys
import tempfile


DRIVER = r"""
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "%(header)s"
int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: driver <points.csv>\n"); return 2; }
    std::FILE* f = std::fopen(argv[1], "r");
    if (!f) { std::fprintf(stderr, "cannot open %%s\n", argv[1]); return 2; }
    char line[65536];
    while (std::fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        std::vector<double> v;
        char* p = line;
        while (*p) {
            char* end = 0;
            double x = std::strtod(p, &end);
            if (end == p) break;
            v.push_back(x);
            p = end;
            while (*p == ',' || *p == ' ') ++p;
        }
        if (v.size() < %(ns)s::nInputs + 1) continue;
        v.pop_back();                       /* drop the trailer, the expected value */
        std::printf("%%.17g\n", %(ns)s::eval(v));
    }
    std::fclose(f);
    return 0;
}
"""


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: python3 verifyExport.py <surrogate_weights_*.hh>")
    header = os.path.abspath(sys.argv[1])
    if not os.path.exists(header):
        sys.exit("ERROR: %s not found." % header)
    ref = os.path.splitext(header)[0] + "_reference.csv"
    if not os.path.exists(ref):
        sys.exit("ERROR: %s not found.  It is written by trainSurrogate.py "
                 "next to the header; re-export to regenerate it." % ref)

    ns = None
    for line in open(header):
        if line.startswith("namespace "):
            ns = line.split()[1].strip("{ ").strip()
            break
    if ns is None:
        sys.exit("ERROR: no namespace found in %s -- is this a generated "
                 "surrogate header?" % header)

    tmp = tempfile.mkdtemp(prefix="srgverify_")
    src = os.path.join(tmp, "driver.cpp")
    exe = os.path.join(tmp, "driver")
    with open(src, "w") as f:
        f.write(DRIVER % {"header": header, "ns": ns})

    cc = os.environ.get("CXX", "g++")
    cmd = [cc, "-O2", "-std=c++11", "-Wall", "-Wextra", "-o", exe, src]
    print("Compiling: %s" % " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.stdout.strip() or r.stderr.strip():
        print(r.stdout + r.stderr)
    if r.returncode != 0:
        sys.exit("FAIL: the generated header does not compile.")
    if r.stderr.strip():
        print("NOTE: the compiler emitted warnings above at -Wall -Wextra.")

    out = subprocess.run([exe, ref], capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit("FAIL: the driver crashed.\n" + out.stderr)

    got = [float(x) for x in out.stdout.split()]
    want = []
    for line in open(ref):
        if line.startswith("#") or not line.strip():
            continue
        want.append(float(line.rsplit(",", 1)[1]))

    if len(got) != len(want):
        sys.exit("FAIL: the C++ produced %d values, the reference has %d."
                 % (len(got), len(want)))

    worst = 0.0
    worst_i = -1
    for i, (a, b) in enumerate(zip(got, want)):
        d = abs(a - b) / max(abs(b), 1e-30)
        if d > worst:
            worst, worst_i = d, i
    print("Checked %d points." % len(got))
    print("Worst relative difference: %.3e (point %d: C++ %.17g, trainer %.17g)"
          % (worst, worst_i, got[worst_i], want[worst_i]))

    if worst > 1e-10:
        sys.exit("FAIL: the compiled network does not reproduce the trainer.")
    print("PASS: the compiled network reproduces the trainer to %.1e relative."
          % worst)


if __name__ == "__main__":
    main()
