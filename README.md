# CompLB3D

**CompLB3D: 3D pore-scale reactive transport with microbial metabolism**

A 3D pore-scale reactive transport model: lattice-Boltzmann flow and solute
transport coupled to microbial metabolism, aqueous speciation, and mineral
precipitation and dissolution that change the pore geometry as the simulation
runs.

CompLB3D extends [CompLaB](https://bitbucket.org/MeileLab/complab) from 2D to
3D and adds an equilibrium chemistry solver, an abiotic kinetics path, and
mineral dissolution.

---

## Authors

| | |
|---|---|
| **Shahram Asgari** | Department of Marine Sciences, University of Georgia, Athens, GA, USA — <shahram.asgari@uga.edu> |
| **Christof Meile** | Department of Marine Sciences, University of Georgia, Athens, GA, USA |

Developed in the Meile Lab, University of Georgia. CompLB3D builds on the
two-dimensional CompLaB, developed by Heewon Jung (Chungnam National
University) and co-workers with the University of Georgia.

---

## What it does

| | |
|---|---|
| **Flow** | Navier–Stokes by lattice Boltzmann (D3Q19), in an arbitrary pore geometry |
| **Transport** | advection–diffusion for any number of solutes (D3Q7) |
| **Biomass** | three interchangeable solvers: cellular automaton, finite difference, lattice Boltzmann |
| **Metabolism** | hand-written kinetics, flux balance analysis through **GLPK** or **COBRApy**, or a trained **surrogate** network — chosen per organism |
| **Speciation** | a components-and-species equilibrium tableau, solved in every pore voxel |
| **Geometry change** | mineral precipitation clogs pores; dissolution reopens them; the flow is re-solved either way |

Every one of those is optional and switched from `CompLaB.xml`. A run with all
of them off is a plain flow-and-transport solver.

---

## Requirements

| | |
|---|---|
| Palabos | **v2.3.0**, downloaded separately (see below) |
| CMake | ≥ 3.10 |
| A C++11 compiler | GCC or Clang |
| MPI | optional but on by default |
| GLPK | optional, only for `-DENABLE_GLPK=ON` |
| Python 3 + cobrapy | optional, only for `-DENABLE_COBRAPY=ON` |

**Palabos is not in this repository.** It is a third-party library of about
100 MB under its own licence; vendoring it would make this repo unclonable on a
slow link and muddy the licence question. Download it from
<https://palabos.unige.ch/> and point CMake at it.

## Building

```bash
git clone https://github.com/shahram444/CompLB3D-lattice-Boltzmann-FBA-surrogate-COBRApy-GLPK.git
cd CompLB3D-lattice-Boltzmann-FBA-surrogate-COBRApy-GLPK

# Edit lines 85 to 93 of CMakeLists.txt to point at your Palabos v2.3.0 tree.

mkdir build && cd build
cmake ..                          # plain build
cmake -DENABLE_GLPK=ON ..         # with flux balance analysis through GLPK
cmake -DENABLE_COBRAPY=ON ..      # with flux balance analysis through COBRApy
make -j
```

`CompLaB.xml` is read from the working directory, by that exact name.

```bash
cd ..
./complab              # or:  srun ./complab
```

`comp.sh` is a Slurm template; `examples/runExamples.slurm` submits the whole
example suite.

---

## Start here

```bash
python3 examples/makeGeometry.py
python3 examples/makeExamples.py
./examples/runAllExamples.sh .
```

Fifteen cases, one per capability, each running in seconds on a laptop. Each is
a self-contained folder with its own `CompLaB.xml`, geometry, kinetics headers
and a README saying what to expect. `examples/README.md` is the map.

Four of them have answers that do not come from this code — an analytic
diffusion profile, mass balances, and a flux balance case whose optimum is
`min(Vs, 2·Vo)` by hand. Those are the ones to trust.

---

## Configuring a run

Everything is in `CompLaB.xml`. `CompLaB_reference_template.xml` documents every
tag with its units, defaults and known limits; it is the reference to keep open.

Two choices are per organism and independent of each other:

```xml
<microbe0>
    <solver_type>CA</solver_type>          <!-- how biomass moves: CA, FD, LBM -->
    <reaction_type>glpk</reaction_type>    <!-- how it metabolises -->
</microbe0>
```

`<reaction_type>` accepts `none`, `kinetics`, `glpk`, `cobrapy`, `surrogate`,
and the combinations `glpk_and_kinetics`, `cobrapy_and_kinetics`,
`surrogate_and_kinetics`, which run both paths and add the rates. Each also
needs its global switch on (`<enable_fba_glpk>` and friends), and a mismatch
stops the run with a message rather than doing something quietly wrong.

**Rate laws are compiled in.** `defineKinetics.hh` and
`defineAbioticKinetics.hh` are `#include`d, not read at run time, so changing
your chemistry means a rebuild. Every example carries its own pair.

---

## Fitting your own surrogate

`surrogate_training/` has the offline FBA sweep, a MATLAB trainer and a
licence-free Python one that emit byte-identical C++, and a verifier that
compiles the generated header and checks it against the fit.

`inspectSurrogate.py` reports any network's **valid input range**, recovered by
inverting the `mapminmax` scaling. Worth running before trusting a network
someone else fitted, including the one shipped in `surrogateModel.hh`: outside
its training box a network does not fail, it returns confident nonsense.

---

## Licence and third-party code

CompLB3D is **AGPL-3.0-or-later**, the same as CompLaB and Palabos. See
`LICENSE`.

`src/complab3d_glpkcpp.hh` derives from GLPKMEX / the Octave-Forge `glpkcc.cpp`
GLPK interface, which is GPL-licensed. It is included here on that basis and
its provenance is recorded in the file header. Anyone redistributing should
satisfy themselves that this is correct for their situation.

Palabos is a separate work under its own licence and is not distributed here.

---

## Citing

Use the **"Cite this repository"** button in the sidebar, which reads
`CITATION.cff`, or cite directly:

> Asgari, S. and Meile, C. (2026). *CompLB3D: a 3D pore-scale reactive
> transport model with microbial metabolism* (version 1.0.0).
> https://github.com/shahram444/CompLB3D-lattice-Boltzmann-FBA-surrogate-COBRApy-GLPK

```bibtex
@software{asgari_complb3d_2026,
  author  = {Asgari, Shahram and Meile, Christof},
  title   = {{CompLB3D: a 3D pore-scale reactive transport model with
             microbial metabolism}},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/shahram444/CompLB3D-lattice-Boltzmann-FBA-surrogate-COBRApy-GLPK}
}
```

Please cite the 2D CompLaB paper alongside this one: CompLB3D extends that work
rather than replacing it.

## Contributing

See `CONTRIBUTING.md`. The short version: run the example suite before opening
a pull request.
