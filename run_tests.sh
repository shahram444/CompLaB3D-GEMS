#!/bin/sh
# Every check in one command.  Needs g++ and python3 with numpy.  Does NOT need Palabos.
set -e
cd "$(dirname "$0")/tests"
rm -f t1 t2 t3 t4 t5 t6 t7 2>/dev/null || true
echo "### the expression language"
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t1 test_sym.cpp                  && ./t1
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t2 test_file.cpp                 && ./t2
echo "### the graph network"
g++ -O2 -Wall -Wextra -std=c++11 -I../src        -o t3 test_gnn.cpp                  && ./t3
echo "### the two processors, run on a stub lattice"
g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src    -o t4 test_symbolic_processor.cpp   && ./t4
g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src    -o t5 test_graphnet_processor.cpp   && ./t5
echo "### the abiotic sweep: reacts with no biomass present"
g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src    -o t6 test_abiotic.cpp               && ./t6
echo "### dissolution: the mineral inventory actually falls"
g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src    -o t7 test_dissolution.cpp             && ./t7
echo "### C++ against Python, same network same inputs"
python3 xval_gnn.py
echo "### the shipped example against the law it was fitted from"
python3 "../examples/16_learned_rate_laws/fitting/check_aom.py"
echo "### the symbolic fitter writes a file the solver accepts"
python3 smoke_fit.py
rm -f t1 t2 t3 t4 t5 t6 t7 b*.sym 2>/dev/null || true          # test_file.cpp leaves its malformed samples behind
echo
echo "everything passed"
