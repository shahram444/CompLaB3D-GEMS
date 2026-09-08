# 10 - FBA through COBRApy: the same linear program, a different solver

## 1. The scenario

This is example 09's case again, the same pore space, the same three
solutes, the same four-reaction metabolic model, run through a different
linear-program solver. Flux balance analysis, FBA for short, asks a linear
program in every voxel and at every step how fast an organism's metabolic
network can run given what diffusion has delivered there; example 09 calls
GLPK from inside the C++ solver, this one hands the stoichiometry to
COBRApy, a Python package widely used in the metabolic-modelling community
for exactly this kind of constraint-based analysis, running in an embedded
Python interpreter, and reads the fluxes back.

A reactive-transport modeller can read this pair as a verification exercise:
two independent implementations of the same coupling, checked against each
other rather than against a hand-derived answer. A computational biologist
will recognise COBRApy itself as the tool the field already reads models
with, so this is also the path that reaches everything COBRApy can do that a
bare linear-program call cannot.

The two must agree. They share the geometry, the transport step and the
uptake-bound construction, and they are handed the same linear program; the
only thing that differs is who solves it. So a disagreement between 09 and
10 is not a modelling choice, it is a bug in one of the two couplings, and
that is the strongest check available on the flux-balance path in this code.
Why keep both: GLPK is roughly fifty times faster and needs no interpreter to
configure, so it is what you run day to day; COBRApy is what you check
against, because it is the tool the field trusts.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um
   (the pore space itself is 24 x 24 x 6; preprocess.py adds a one voxel
   wall on the y and z faces the solver conditions with nothing, so NY
   and NZ come out two larger than the pore space it declares)

   x=0                                                          x=23
    |                                                             |
S,O |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> |
1.0 |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> | Neumann,
    |>>>  ====                                          >>>>>>>> | open outlet
    |               ####    ####    ####    ####    ####         |
    |               ====                                         |
    +-------------------------------------------------------------+
      grain blocks alternate between the y=0 and y=23 walls every
      four voxels along x, so diffusion has to weave to get through

      #### = solid grain (wall)     ==== = initial Toybug patch
      porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry boundary conditions. The
other four faces are an inert wall, drawn by `preprocess.py`, because the
solver gives a face with no declared condition nothing at all, and a lattice
Boltzmann field with nothing conditioning it streams off the edge of the
domain and reads back values that nothing updates. This is the identical
geometry file example 09 uses, so the two back ends can be compared voxel by
voxel.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `S`, `O`, `P` |
| **Diffusivity** | 5e-10 m2/s for all three species, in pore and in biofilm alike |
| **Chemistry** | none written by hand. `<enable_kinetics>false</enable_kinetics>`, `<enable_abiotic_kinetics>false</enable_abiotic_kinetics>` |
| **Biology** | `Toybug`, one organism, `<reaction_type>cobrapy</reaction_type>`, attached biofilm seeded on its own material number, spreading by finite-difference diffusion |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 500 advection-diffusion steps, output every 100 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **S** (donor) | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **O** (acceptor) | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **P** (product) | closed | closed | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, which is a reservoir. *Neumann with a zero
gradient* means nothing flows across that face, an open outlet the species
never actually reaches in this run. *Closed* means the same thing on both
ends: P is made inside the domain by the reaction and stays there.

## 4. The reaction

`input/toy_model.xml`, the same file example 09 uses, small enough to check
by hand: three metabolites, four reactions. A metabolic model is a
stoichiometric matrix, the table of how many molecules of each metabolite
each reaction consumes or produces, together with an objective, the one
combination of reaction rates (fluxes) the linear program is asked to
maximise.

```
    metabolites      S_c    O_c    P_c

    reactions
      0  EX_S        S_c ->                    exchange reaction; negative
                                                flux imports the donor
      1  EX_O        O_c ->                    negative flux imports the
                                                acceptor
      2  EX_P        P_c ->                    positive flux exports the
                                                product
      3  BIO         S_c + 0.5 O_c -> P_c      the objective: this flux IS
                                                growth
```

The linear program, in words: maximise the objective flux, subject to every
internal metabolite being at steady state, and subject to each flux staying
between a lower and an upper bound set by the local uptake bounds. Because
the exchange reactions are tied to the objective by the steady-state
constraint, the whole system reduces to one line worth checking by hand:

```
    maximise   c' v          the objective flux, reaction 3
    subject to S v  = 0      every internal metabolite at steady state
               lb <= v <= ub  uptake bounded by the LOCAL concentrations

    growth  =  min( Vs , 2 * Vo )
```

where Vs and Vo are the uptake bounds rebuilt every step from the local
concentration through a Michaelis-Menten term, `V = Vmax * C / (Kc + C)`. That
number comes from the stoichiometry, not from either solver, which is
exactly what makes this pair a test rather than a demonstration.

How COBRApy is reached from C++:

```
    C++ side                                Python side
    --------                                -----------
    read input/toy_model.xml once     -->   build a cobra.Model from the arrays
    per voxel, per step:
      build lb/ub from local C        -->   set Reaction.bounds
      call run_FBA()                  -->   model.optimize()
      read growth + exchange fluxes   <--   solution.fluxes
```

The model file is opened by the C++ side, not by COBRApy: the stoichiometric
matrix, the bound vectors and the objective position cross the boundary as
plain arrays, and the Python half assembles a `cobra.Model` from them. That
is what makes the two back ends comparable, they are handed the same
numbers, in the same order, from the same parser. The contract for that
argument vector is written out in full at the top of
`training/complab3d_cobrapy.py`, and duplicated verbatim in
`src/complab3d_pythonAPI.hh`; if you change one, change the other.

## 5. What happens each step

The solver repeats this 500 times:

1. **Stream and collide** each species on its own D3Q7 lattice, advancing
   diffusion by one step, and diffuse the biomass field the same way.
2. **Apply the boundary conditions**, so S and O are re-pinned to 1.0 at the
   left face.
3. **Build the uptake bounds**, per voxel, per substrate: a Michaelis-Menten
   term on the local concentration, capped by `<fba_maximum_uptake_flux>`.
4. **Hand the bounds to the embedded Python interpreter**, which sets them as
   `Reaction.bounds` on a pre-built `cobra.Model` and calls
   `model.optimize()`. This is the expensive step and it is more expensive
   than example 09's: a Python call and a general-purpose solver, in every
   voxel that holds biomass, every iteration.
5. **Read the growth rate and the exchange fluxes back** from
   `solution.fluxes`, scale by the local biomass and the time step through
   `<biomass_molar_mass>`, and write the increments into the change
   lattices.
6. **Add the increments**, clamped so a step cannot draw a voxel's substrate
   below what it actually holds.

## 6. What comes out

```
output/
  S_0000100.vti  ...  S_0000500.vti        concentration of the donor
  O_*.vti  P_*.vti                          the other two species
  Toybug_*.vti                              the biomass field
  rate_S_*.vti  rate_O_*.vti  rate_P_*.vti  the reaction rate as a field, mol/L/s
  summary.csv                               one row every 50 steps
  run.log                                   the whole run, including the checks
```

| Quantity | What it should do |
|---|---|
| Total S, total O | Both **fall** relative to what a Dirichlet supply alone would give, because the organism is consuming them |
| Total P | Rises from 0 and keeps rising; P is the only closed species |
| Total biomass | Rises where both S and O are present, flattens or falls where either runs out |
| S : O consumption ratio | Near 1 : 0.5, the stoichiometry of the objective reaction, wherever the program is unconstrained |
| Every field, compared to example 09 | Agrees to solver tolerance. Anything larger is a bug in one coupling, not a property of either solver |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Slower than example 09's GLPK run, on the order of tens of times slower, because each solve now crosses into Python and a general-purpose solver |

## 7. What to check

1. **09 and 10 agree.** Run both and compare `output/summary.csv` field by
   field; they should agree to solver tolerance.
2. **The growth rate approaches the stoichiometric bound**, `min(Vs, 2*Vo)`,
   the same relation example 09 checks.
3. **P accumulates at the rate biomass is made.** The objective reaction
   exports exactly one P per unit growth, so the two curves are the same
   curve.
4. **S and O are drawn in a 1 : 0.5 ratio** wherever the program is
   unconstrained.
5. **No `[NEG!]` warnings.** Uptake is bounded by what is locally present and
   the increment is clamped again before it is applied.

`postprocess.py` reads `output/summary.csv` and `run.log`, reports every
field's change, flags any negative value, and prints the log lines that
mention the linear program, so checks 3, 4 and 5 are largely done for you.
Check 1 needs both runs side by side, and check 2 needs a look at the fields;
both are yours to do.

## 8. What this case demonstrates

**The same coupling on a second solver, as a cross-check rather than a second
demonstration.** Two things cost people time on this path specifically.

The interpreter: COBRApy is imported by name at run time from `<src_path>`,
so `src/complab3d_cobrapy.py` has to be there, and the embedded interpreter
has to be the same Python you installed `cobra` into. The solver prints the
interpreter path in its first few lines and `offline.sh` prints the shell's,
so a mismatch is visible before the run rather than during it.

Positional exchange indices: `<exchange_reaction_indices>` is positional,
`0 1 2` is correct only for this exact model file. Insert one reaction into
the model and every index after it points somewhere else, silently.
`<exchange_reaction_names>` names the reactions instead and resolves them
against the model at start-up, so a wrong name stops the run and prints the
near matches; it needs an SBML model, which this flat matrix file is not, so
this case still uses indices. Example 16 shows the named form on a
genome-scale model.

**What it deliberately leaves out:** flow, an abiotic reaction, and any
change to the pore space, same as example 09. Example 11 replaces the linear
program with a network fitted to its answers. Example 12 runs one organism
through both a flux-balance path and a hand-written kinetics path in the
same voxels.

## 9. How to build and run

Two files do this, and they are deliberately separate because compiling and
running belong in different places on a cluster.

| File | What it is | Where you run it |
|---|---|---|
| `COMPILE.txt` | The interactive session that builds `./complab`, step by step | Once, by hand, on an interactive node |
| | `run.sh` | The SLURM batch job: modules, pre-processing, the solver, the checks. Submit with `sbatch run.sh`. |
| | `COMPILE.txt` | The interactive session that builds `./complab`, one step at a time. |

**First time in this case folder,** open `COMPILE.txt` and follow it. It asks
for an interactive node, loads the modules, runs cmake with the flags this case
actually needs, and compiles. It also says exactly when you have to come back
and recompile, which for most cases is almost never.

**Every run after that** is one command:

```bash
sbatch run.sh
squeue -u $USER                 # watch it
tail -f output/run.log          # read it while it runs
```

`run.sh` carries a comment on every line: the SLURM header, the module loads
that must match what you compiled with, the pre-processing, the solver call and
the post-processing checks. It runs on one rank on purpose, and the note at the
bottom of the file says why.

To work in a scratch copy instead of dirtying this folder:

```bash
./scripts/setup_case.sh 10_fba_cobrapy run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
Palabos v2.3.0, Python development headers (`python3-dev` on Debian/Ubuntu),
and a Python interpreter with `cobra` importable. This case's
`<reaction_type>` is `cobrapy`, so the executable must be built with the flag
that links in the embedded Python interpreter:

```bash
cmake -B build -S . \
      -DCMAKE_BUILD_TYPE=Release \
      -DPALABOS_ROOT=/path/to/palabos-v2.3.0 \
      -DENABLE_COBRAPY=ON
cmake --build build -j
```

**When you must recompile.** `-DENABLE_COBRAPY=ON` only makes the embedded
interpreter available in the executable; `<enable_fba_cobrapy>true</enable_fba_cobrapy>`
in `CompLaB.xml` switches it on for a given run, and that line can be edited
without a rebuild. But the two have to agree: a build without the flag
compiles cleanly and then refuses at start-up, because turning the XML path
on is not enough if the flag was never passed to CMake. Separately from the
build flag, `src/complab3d_cobrapy.py` must exist at run time (`setup_case.sh`
copies it there from `tools/`) and the embedded interpreter must be able to
`import cobra`, or the run fails at start-up with "failed to import module
complab3d_cobrapy" after the build has already succeeded. Editing
`input/toy_model.xml` or `input/geometry.dat` needs no rebuild; only a change
to `CMakeLists.txt` or the solver sources does.

On an EasyBuild cluster, loading a Python module is not enough for this
path: the module system's `PYTHONPATH` hook applies to your interactive
shell, not to the embedded interpreter the compiled executable starts, so
`cobra` can be importable from your terminal and still missing to the
running solver. `run.sh` shows the explicit `PYTHONPATH` export this needs.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Prints which Python `cobra` imports from, and checks the model file parses. A mismatch between that Python and the embedded interpreter's is the usual failure on this path. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| | `input/toy_model.xml` | The metabolic model of section 4, in the flat `<Metabolic_Model>` format. Read by the C++ side, not by COBRApy. |
| | `src/complab3d_cobrapy.py` | The Python half of the bridge, imported by name at run time. `setup_case.sh` copies it from `tools/` into `src/`. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
