# Installing

## What you need

| | Needed for | Notes |
|---|---|---|
| A C++17 compiler | everything | GCC 9+, Clang 10+, Intel 2021+ |
| CMake 3.16+ | everything | |
| [Palabos](https://palabos.unige.ch) 2.3.0 | everything | the lattice Boltzmann library underneath |
| MPI | parallel runs | OpenMPI or MPICH; serial runs need none |
| [GLPK](https://www.gnu.org/software/glpk/) | flux balance, GLPK path | `libglpk-dev` |
| CPython 3.8+ with [COBRApy](https://opencobra.github.io/cobrapy/) | flux balance, COBRApy path | only if you use that back end |
| Python 3.8+ with NumPy | all the offline tools | |

Nothing beyond the compiler and Palabos is needed to run precipitation,
dissolution or compiled kinetics.

## Palabos

```bash
wget https://palabos.unige.ch/downloads/palabos-v2.3.0.zip
unzip palabos-v2.3.0.zip
export PALABOS_ROOT=$PWD/palabos-v2.3.0
```

## Build

```bash
git clone https://github.com/shahram444/CompLaB3D-GEMS.git
cd CompLaB3D-GEMS

cmake -B build -S . \
      -DPALABOS_ROOT=$PALABOS_ROOT \
      -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Optional back ends are off by default, so the code builds on a machine that has
neither GLPK nor Python. Turn on only what you have installed and intend to use:

```bash
cmake -B build -S . -DPALABOS_ROOT=$PALABOS_ROOT \
      -DENABLE_GLPK=ON \
      -DENABLE_COBRAPY=ON
```

| Option | Default | What it does |
|---|---|---|
| `ENABLE_GLPK` | OFF | Flux balance analysis through GLPK, in process. Needs `libglpk-dev`. |
| `ENABLE_COBRAPY` | OFF | Flux balance analysis through an embedded Python interpreter. Needs the Python development headers and COBRApy. |
| `ENABLE_MPI` | ON | Parallel runs. Turn off for a serial-only build. |
| `FBA_BULK_ONLY` | OFF | Restrict the FBA processors to the bulk domain, skipping envelope cells. Typically 20–40% faster; verify against a default run on your own case before trusting it. |

The surrogate, symbolic and graph-network paths need no build option — they are
always available.

Turning an option ON only makes the solver **available**. You still switch it on
for a given run in `CompLaB.xml`, under `<simulation_mode>`. If the XML asks for
a solver that was not built in, the program stops at start-up with a message
naming the option to enable, rather than ignoring you silently.

`comp.sh` is a thin wrapper over the same commands if you would rather not think
about CMake.

## Check it works

```bash
./tests/check_repo.sh     # the tree, in one second, no dependencies
./tests/run_tests.sh      # the solver
```

Seven suites, 154 checks, no Palabos needed — the rate processors compile
against a minimal stub in `tests/palabos_stub.hh` and run on lattices the tests
fill directly.

Then run something real:

```bash
./scripts/setup_case.sh 01_flow_only run/first
cd run/first && cmake -B build -S . && cmake --build build -j
./complab CompLaB.xml
```

## Common problems

**`Palabos headers not found`** — `PALABOS_ROOT` is not set, or points at the
directory above the unzipped one. It should contain `src/palabos3D.h`.

**`glp_create_prob was not declared`** — GLPK is missing while `ENABLE_GLPK=ON`.
Install `libglpk-dev`, or build without that back end.

**COBRApy path aborts at start-up** — the embedded interpreter could not import
`cobra`. Check that the Python it found is the one with COBRApy installed:
the solver prints the interpreter path in its first few lines.

**The run does nothing and finishes fast** — usually a geometry with no
percolating path. Check it before blaming the chemistry:

```bash
python tools/geometry.py --inspect input/geometry.dat
```
