#!/bin/sh
# Every check in one command.  Needs g++ and python3 with numpy.  Does NOT need Palabos.
set -e
cd "$(dirname "$0")"
rm -f t13 t14 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 2>/dev/null || true
echo "### the expression language"
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t1 test_sym.cpp                  && ./t1
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t2 test_file.cpp                 && ./t2
echo "### the .sym reaction block: the solver writes the substrate lines"
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t14 test_sym_stoich.cpp          && ./t14
echo "### the surrogate: more than growth out of one network"
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t8 test_surrogate_multi.cpp && ./t8
echo "### the surrogate: the run-time and compiled paths agree"
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t9 test_surrogate_parity.cpp && ./t9
echo "### the surrogate: the compiled path holds its training box too"
g++ -O2 -Wall -Wextra -std=c++11 -I. -I.. -I../src -o t12 test_surrogate_range.cpp && ./t12
echo "### the graph network"
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t3 test_gnn.cpp                  && ./t3
echo "### the two processors, run on a stub lattice"
g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src    -o t4 test_symbolic_processor.cpp   && ./t4
g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src    -o t5 test_graphnet_processor.cpp   && ./t5
echo "### the thermodynamic gate"
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t11 test_thermo.cpp             && ./t11

echo
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t13 test_upscale.cpp            && ./t13
echo "### the wiring: every rate path is reachable, not just implemented"
g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src    -o t10 test_wiring.cpp             && ./t10
echo "### the abiotic sweep: reacts with no biomass present"
g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src    -o t6 test_abiotic.cpp               && ./t6
echo "### C++ against Python, same network same inputs"
python3 xval_gnn.py
echo "### the shipped example against the law it was fitted from"
python3 check_aom.py
echo "### the symbolic fitter writes a file the solver accepts"
python3 smoke_fit.py
echo "### every shipped rate law survives a substrate list shorter than it wants"
python3 check_kinetics_bounds.py
echo "### the shipped example rate laws load"
g++ -O2 -Wall -Wextra -std=c++11 -I../src -o t7 check_examples.cpp && ./t7
rm -f t13 t14 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 sym_*.sym b*.sym *.thm *_test.srg 2>/dev/null || true          # test_file.cpp leaves its malformed samples behind
echo
echo "everything passed"
