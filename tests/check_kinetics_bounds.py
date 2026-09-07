#!/usr/bin/env python3
"""Every shipped rate law must survive a substrate list shorter than it expects.

defineKinetics.hh and defineAbioticKinetics.hh are the two files a user is meant
to edit, and they are the only place in the program where a plain C++ subscript
is written against a vector whose length comes from CompLaB.xml. Write
`subsR[3] = ...` in a case that declares three substrates and the solver reads
one double past the end of a heap block: no compiler warning, no crash most of
the time, and a rate that is whatever happened to be in the next allocation.

That is not a hypothetical. The comment at the top of the shipped default header
records why it exists at all -- so an example could not accidentally compile
against a network indexing C[0] to C[94].

Every shipped header currently guards its own lengths on entry. Nothing checked
that, which is exactly how the next one comes to be written without a guard. So
this compiles each header on its own under AddressSanitizer and calls it with
every substrate count from zero up to one more than it needs, and every microbe
count from zero to three. A read or write past the end aborts the process and is
reported here with the header that did it.

    python3 tests/check_kinetics_bounds.py

Silence and a zero exit mean every shipped rate law is safe on a short list.
Run from the repository root; called by tests/run_tests.sh.
"""
import glob
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# The two entry points, and the driver body that exercises each. Both are called
# once per voxel per iteration by the solver, with vectors sized from the XML.
BIOTIC = """
    std::vector<double> B(nb, 0.7), C(nc, 0.5);
    std::vector<double> subsR(nc, 0.0), bioR(nb, 0.0);
    defineRxnKinetics(B, C, subsR, bioR, mask);
"""
ABIOTIC = """
    std::vector<double> C(nc, 0.5);
    std::vector<double> subsR(nc, 0.0);
    defineAbioticRxnKinetics(C, subsR, mask);
    double mineral = 1.0, mineralR = 0.0;
    defineDissolutionRate(0, C, mineral, subsR, mineralR, mask);
    (void) nb;
"""

DRIVER = """
#include <vector>
#include <cstdio>
#include "palabos_stub.hh"
namespace plb {{ PlbOut pcout; }}

#include "{header}"

static void one(std::size_t nc, std::size_t nb, plb::plint mask)
{{
{body}
}}

int main()
{{
    /* mask 0 and 1 are solid and wall; 2 upwards is open pore or biofilm. A rate
     * law is called in all of them and must be safe in all of them. */
    for (plb::plint mask = 0; mask <= 3; ++mask)
        for (std::size_t nc = 0; nc <= 5; ++nc)
            for (std::size_t nb = 0; nb <= 3; ++nb)
                one(nc, nb, mask);
    std::printf("ok\\n");
    return 0;
}}
"""


def headers():
    found = []
    for pat in ("config/kinetics/*.hh", "examples/*/kinetics/*.hh"):
        found += sorted(glob.glob(os.path.join(ROOT, pat)))
    return found


def main():
    bad = []
    checked = 0
    with tempfile.TemporaryDirectory() as tmp:
        for h in headers():
            rel = os.path.relpath(h, ROOT)
            body = ABIOTIC if "Abiotic" in os.path.basename(h) else BIOTIC
            src = os.path.join(tmp, "drv.cpp")
            exe = os.path.join(tmp, "drv")
            with open(src, "w") as f:
                f.write(DRIVER.format(header=h, body=body))
            cc = subprocess.run(
                ["g++", "-O1", "-g", "-std=c++11", "-fsanitize=address,undefined",
                 "-fno-omit-frame-pointer", "-I", HERE, "-I", os.path.join(ROOT, "src"),
                 "-o", exe, src],
                capture_output=True, text=True)
            if cc.returncode != 0:
                bad.append((rel, "does not compile", cc.stderr.strip()[-1500:]))
                continue
            run = subprocess.run([exe], capture_output=True, text=True,
                                 env=dict(os.environ, ASAN_OPTIONS="detect_leaks=0"))
            checked += 1
            if run.returncode != 0 or "ok" not in run.stdout:
                bad.append((rel, "reads or writes out of bounds",
                            (run.stderr or run.stdout).strip()[-2500:]))

    for rel, what, detail in bad:
        print("  %-64s %s" % (rel, what))
        for line in detail.splitlines()[:20]:
            print("      " + line)
    if bad:
        print("\n%d of %d shipped rate law(s) are not safe on a short substrate list"
              % (len(bad), checked + len(bad)))
        return 1
    print("all %d shipped rate laws are safe on a short substrate list" % checked)
    return 0


if __name__ == "__main__":
    sys.exit(main())
