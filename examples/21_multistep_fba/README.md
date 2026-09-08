# 21 - Multi-step FBA: pinning the by-products that plain FBA leaves free

## 1. The scenario

Flux balance analysis, FBA for short, predicts a growth rate by solving a
linear program: an optimisation problem in which both the objective, the
quantity to maximise, and every constraint are straight lines rather than
curves. Given a metabolic model, a stoichiometric matrix listing how many
molecules of each internal metabolite every reaction consumes or produces,
and a bound on how much of each substrate a cell can take up, the linear
program asks for the one combination of internal reaction rates, called
fluxes, that maximises growth while keeping every internal metabolite exactly
balanced. Case 09 solves exactly one such program per voxel, per step, and
uses whatever the solver returns.

For the growth rate itself that is enough, because in a typical model the
optimal growth rate is unique: every combination of fluxes that achieves the
maximum growth rate achieves the same value of the objective, so there is
nothing left to choose between them. The by-products are a different story.
Offered glucose and oxygen together at generous bounds, the E. coli core model
used in these cases can hit the same maximum growth rate while exporting
anywhere along a whole range of acetate and formate combinations: the two
by-products absorb the same spare electrons, the objective does not care which
one carries them, and the plain linear program returns whichever combination
its internal pivoting happens to land on first. Nothing about that answer is
wrong, it genuinely achieves the optimal growth rate, but it is arbitrary. It
need not survive a change of solver version, and there is no way to compare it
against a measured yield, because the number reported is not determined by
the model at all.

This case removes that arbitrariness with lexicographic optimisation: instead
of one linear program, a short chain of them, each optimising the next
quantity of interest while holding every earlier one at, or near, its own
optimum. After the last stage in the chain, nothing in the model is free to
vary any more, and the same model at the same bounds returns the same numbers
on any solver, every time. This is the multi-step FBA of Song et al. (2025);
the identical device, under the name lexicographic optimisation, appears in
DFBAlab (Gomez, Hoffner and Barton 2014).

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um
   (the pore space itself is 24 x 24 x 6; preprocess.py adds a one voxel
   wall on the y and z faces the solver conditions with nothing, so NY
   and NZ come out two larger than the pore space it declares)

   x=0                                                          x=23
    |                                                             |
glc |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> |
o2  |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> |
1.0 |>>>  ====                                          >>>>>>>> | Neumann,
    |               ####    ####    ####    ####    ####         | open outlet
    |               ====                                         |
    +-------------------------------------------------------------+
      grain blocks alternate between the y=0 and y=25 walls every
      four voxels along x, so diffusion has to weave to get through

      #### = solid grain (wall)     ==== = initial Ecoli patch
      porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry boundary conditions. The
other four faces are an inert wall, drawn by `preprocess.py`, for the same
reason as every other case in this set: the solver gives a face with no
declared condition nothing at all, and a lattice Boltzmann field left
unconditioned there streams off the edge of the block and reads back values
nothing updates.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `glc`, `o2`, `ac`, `for` |
| **Chemistry** | none written by hand. Growth and every exchange flux come from the metabolic model |
| **Biology** | `Ecoli`, one organism, `<reaction_type>glpk</reaction_type>`, `e_coli_core.xml` read natively as SBML, attached biofilm spreading by finite-difference diffusion, three-stage lexicographic chain per voxel |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 500 advection-diffusion steps, output every 100 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **glc** (donor) | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **o2** (acceptor) | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **ac** (product) | closed | Neumann, zero gradient | 0 |
| **for** (product) | closed | Neumann, zero gradient | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, which is a reservoir. *Neumann with a zero
gradient* means nothing is forced to flow across that face, an open outlet.
*Closed* on the left means acetate and formate have no source there; both
are made only by the organism and can leave only through the right face,
which is exactly what lets the two stages that pin them constrain a quantity
that actually reaches the transport step.

## 4. The reaction

`input/e_coli_core.xml` is the E. coli core metabolic model, read directly in
SBML format: 72 metabolites, 95 reactions, including the four this case
watches: `EX_glc__D_e`, `EX_o2_e`, `EX_ac_e`, `EX_for_e` for glucose, oxygen,
acetate and formate exchange, and `Biomass_Ecoli_core` for growth. A plain
single-stage linear program on this model would read:

```
    maximise   biomass                    the objective, one flux out of 95
    subject to S v  = 0                   every internal metabolite at steady state
               lb <= v <= ub              uptake and export bounded by local transport
```

which leaves acetate and formate export free within the range described in
section 1. The chain this case runs instead is three such programs, solved in
order, each adding one more constraint fixing what the previous stage already
achieved:

```
    1.  max  biomass                                  ->  mu*
    2.  max  acetate export,   biomass >= a1 x mu*
    3.  max  formate export,   biomass >= a1 x mu*
                                acetate >= a2 x (step 2 result)
```

`a1` and `a2` are the retain fractions in `CompLaB.xml`,
`<retain_fraction>`, one per stage: the fraction of a stage's own optimum
that every later stage in the chain is required to preserve. This case ships
all three at 1.0, so the chain is strict, no stage may give up any part of an
earlier one, and the run carries no fitted parameter at all while still
removing every degree of freedom that plain FBA left open. Lowering a retain
fraction below 1 lets the organism trade growth for excretion, which real
cells do and a genome-scale model will not do on its own unless told to; that
value is then a calibration against a measured yield, not something read off
the network, and belongs in any report of the result.

`<stage_reactions>` names the three reactions by their model names,
`Biomass_Ecoli_core`, `EX_ac_e`, `EX_for_e`, rather than by column number.
That matters because a column index stops meaning anything the moment the
model file is regenerated, and a wrong index would silently pin the wrong
reaction; a wrong name instead stops the run at start-up with the candidate
names it did find.

## 5. What happens each step

The solver repeats this 500 times:

1. **Stream and collide** each of the four species on its own D3Q7 lattice,
   and diffuse the biomass field the same way.
2. **Apply the boundary conditions**, re-pinning glucose and oxygen to 1.0 at
   the left face.
3. **Build the uptake bounds**, per voxel, per substrate, from the local
   concentration, capped by `<fba_maximum_uptake_flux>`.
4. **Solve the three-stage chain**, in every voxel holding biomass: maximise
   biomass, fix it at its retain fraction, maximise acetate export, fix that
   at its retain fraction, maximise formate export. Each stage after the
   first starts from the previous stage's basis, so in practice it finishes
   in a handful of pivots rather than a full solve from scratch; expect the
   whole chain to cost around 1.6 times a single solve, not three times.
5. **Read growth and every exchange flux off the final stage's solution**,
   not the first stage's, scale by the local biomass and the time step, and
   write the increments into the change lattices.
6. **Add the increments**, clamped so a step cannot draw a voxel's substrate
   below what it actually holds.

**Growth and the fluxes must come from the same solution vector.** With every
retain fraction at 1.0 the growth reported by the final stage equals the
first stage's own optimum, so this distinction is invisible in the shipped
case; the moment a retain fraction drops below 1 the two differ, and reading
growth from the wrong stage would add biomass at one rate while removing
substrate at a rate belonging to a different, inconsistent solution.

## 6. What comes out

```
output/
  glc_0000100.vti  ...  glc_0000500.vti     concentration of glucose
  o2_*.vti  ac_*.vti  for_*.vti             the other three species
  Ecoli_*.vti                                the biomass field
  rate_glc_*.vti  rate_o2_*.vti  ...        the reaction rate as a field, mol/L/s
  summary.csv                                one row every 50 steps
  run.log                                    the whole run, including the checks
```

| Quantity | What it should do |
|---|---|
| Total glc, total o2 | Fall relative to what pure diffusion alone would give, because the organism is consuming them |
| Total ac, total for | Rise from 0 and keep rising; both are closed at the left and made only by the reaction |
| Total biomass | Rises wherever both glucose and oxygen reach the biofilm |
| Growth (`rate_glc` scaled by yield near the patch) | Positive where the chain is unconstrained by local supply |
| Reproducibility | Identical, to the last printed digit, across repeated runs at the same bounds; that is the entire point of the chain |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | About 1.6 times a single-stage FBA case of the same size, not 3 times, because stages 2 and 3 warm-start |

## 7. What to check

1. **Acetate and formate appear.** Both start at zero and are closed at the
   left, so any amount present was made by the organism; if they stay at
   zero the later stages of the chain never ran.
2. **Growth is read from the final stage.** With every retain fraction at 1.0
   this happens to agree with the first stage's own optimum; it is worth
   knowing this before lowering any fraction below 1.
3. **No field goes negative.** A staged solve draws from the same voxel more
   than once, so this also checks that the retained fraction from an earlier
   stage was correctly subtracted from what a later stage was allowed.
4. **The answer is reproducible.** Run the case twice; every flux should
   agree to the last digit. If it does not, a stage is missing from the
   chain or a retain fraction was not actually applied.
5. **Cross-check against COBRApy**, the Python constraint-based modelling
   package used in case 10, by changing `<reaction_type>` to `cobrapy` and
   the two `enable_fba_*` flags to match. The two back ends walk the same
   chain through independent code, so any disagreement is a defect in one of
   the couplings rather than a property of either solver, which makes this
   the strongest check available on this path.

`postprocess.py` runs checks 1 and 3 for you: it reports every field's total
change, flags any negative value, and scans the log for infeasible or
non-converging solves.

## 8. What this case demonstrates

**Lexicographic optimisation**, the chain of linear programs that turns an
arbitrary member of an optimal set into a single, reproducible answer, with
the retain fraction as the one explicit, honestly labelled place a modeller
can choose to trade an earlier objective for a later one.

**What it deliberately leaves out.** No flow, no abiotic reaction, no change
to the pore space, and no thermodynamic gate; case 19 is the case built
around that instead. Case 22 keeps the single-stage linear program per
substrate but replaces the ambiguity problem with a different one, which
carbon source the organism is actually living on, and the two techniques
compose: a `<multi_step>` block can be added inside each of case 22's
per-source solves, at the cost of both stage counts multiplying together.

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
./scripts/setup_case.sh 21_multistep_fba run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
Palabos v2.3.0, and GLPK's development package (`libglpk-dev` on
Debian/Ubuntu, `glpk` from Homebrew on macOS). This case's `<simulation_mode>`
block sets `<enable_fba_glpk>true</enable_fba_glpk>`, which is the CMakeLists
flag `-DENABLE_GLPK=ON`:

```bash
cmake -B build -S . \
      -DCMAKE_BUILD_TYPE=Release \
      -DPALABOS_ROOT=/path/to/palabos-v2.3.0 \
      -DENABLE_GLPK=ON
cmake --build build -j
```

**When you must recompile.** `-DENABLE_GLPK=ON` only makes the GLPK solver
available inside the executable; `<enable_fba_glpk>` in `CompLaB.xml` is what
switches it on for a given run, and that line, along with every stage
reaction, retain fraction and bound in `<multi_step>`, can be edited without a
rebuild. A build without the CMake flag compiles cleanly and then refuses at
start-up, because turning the XML path on is not enough if GLPK was never
linked in. Editing `input/e_coli_core.xml` or `input/geometry.dat` also needs
no rebuild; only a change to `CMakeLists.txt` or the solver sources does.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Shared with case 09. Standard library only. |
| **run** | `CompLaB.xml` | What the solver reads, including the `<multi_step>` chain. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| | `input/e_coli_core.xml` | The metabolic model of section 4, read natively as SBML. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives and infeasible or non-converging solves. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
