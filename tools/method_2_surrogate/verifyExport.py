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
#include <cstddef>
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
        if (v.size() < %(ns)s::nInputs + %(ns)s::nOutputs) continue;
        /* the row is nInputs values then nOutputs expected values; keep the
           inputs, drop the expectations, and print what the header computes */
        v.resize(%(ns)s::nInputs);
        double o[%(ns)s::nOutputs];
        %(ns)s::evalAll(v, o);
        for (std::size_t k = 0; k < %(ns)s::nOutputs; ++k)
            std::printf("%%.17g%%s", o[k], (k + 1 == %(ns)s::nOutputs) ? "\n" : " ");
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

    # How many inputs and outputs the reference carries. Files written before
    # multi-output have neither line and are single-output by definition.
    nin = nout = None
    names = None
    for line in open(ref):
        if not line.startswith("#"):
            break
        b = line[1:].strip()
        if b.startswith("n_inputs:"):
            nin = int(b.split(":", 1)[1])
        elif b.startswith("n_outputs:"):
            nout = int(b.split(":", 1)[1])
        elif b.startswith("outputs:"):
            names = b.split(":", 1)[1].split()
    if nout is None:
        nout = 1

    got = [[float(v) for v in ln.split()]
           for ln in out.stdout.splitlines() if ln.strip()]
    want = []
    for line in open(ref):
        if line.startswith("#") or not line.strip():
            continue
        vals = [float(v) for v in line.split(",")]
        want.append(vals[-nout:])

    if len(got) != len(want):
        sys.exit("FAIL: the C++ produced %d rows, the reference has %d."
                 % (len(got), len(want)))
    if got and len(got[0]) != nout:
        sys.exit("FAIL: the header returns %d output(s), the reference has %d."
                 % (len(got[0]), nout))

    if names is None:
        names = ["output %d" % k for k in range(nout)]

    # PER OUTPUT, not pooled. A single worst-case number over all outputs would
    # let a correct growth column hide a transposed flux column.
    print("Checked %d points, %d output(s)." % (len(got), nout))
    overall = 0.0
    for k in range(nout):
        worst, worst_i = 0.0, -1
        for i in range(len(got)):
            a, b = got[i][k], want[i][k]
            d = abs(a - b) / max(abs(b), 1e-30)
            if d > worst:
                worst, worst_i = d, i
        overall = max(overall, worst)
        print("  %-20s worst %.3e  (point %d: C++ %.17g, trainer %.17g)"
              % (names[k] if k < len(names) else "output %d" % k,
                 worst, worst_i, got[worst_i][k], want[worst_i][k]))

    if overall > 1e-10:
        sys.exit("FAIL: the compiled network does not reproduce the trainer.")
    print("PASS: the compiled network reproduces the trainer to %.1e relative."
          % overall)


if __name__ == "__main__":
    main()
