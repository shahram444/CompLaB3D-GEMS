# 09 - FBA through GLPK: growth predicted by a metabolic model, not prescribed

## 1. The scenario

Examples 05 to 08 tell the organism how fast to grow: a Monod law with a
`mu_max` you chose by hand. This case does not. It gives the organism a
metabolic model, a stoichiometric matrix of its internal reactions, and asks a
linear program, in every voxel and at every step, how fast that network *can*
run given what diffusion has delivered there. Flux balance analysis, FBA for
short, is exactly that: instead of writing down a rate, you write down a
network and an objective, and the rate falls out of an optimisation.

A reactive-transport modeller can read this as the metabolic equivalent of a
local equilibrium or kinetic sub-model coupled to transport, except the
sub-model here is a small optimisation problem rather than a closed-form rate
law. A computational biologist will recognise the linear program itself: it is
the same flux balance analysis run on a genome-scale model in a well-mixed
flask, just solved separately in every voxel of a pore space instead of once
for the whole system.

Growth is the answer to that program, not an input to it. Change the substrate
supply and the growth rate changes because the optimum moves, which is the
point of flux balance analysis and the reason it costs what it costs. The
linear program is solved by GLPK, an open-source linear-program solver, called
directly from the C++ solver. No interpreter, no file round-trip. Example 10
runs the identical case through COBRApy, a Python package for constraint-based
metabolic modelling, as a cross-check: the two must agree to every printed
digit, because they are handed the same numbers.

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
domain and reads back values that nothing updates.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `S`, `O`, `P` |
| **Diffusivity** | 5e-10 m2/s for all three species, in pore and in biofilm alike |
| **Chemistry** | none written by hand. `<enable_kinetics>false</enable_kinetics>`, `<enable_abiotic_kinetics>false</enable_abiotic_kinetics>` |
| **Biology** | `Toybug`, one organism, `<reaction_type>glpk</reaction_type>`, attached biofilm seeded on its own material number, spreading by finite-difference diffusion |
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

`input/toy_model.xml` is a metabolic model deliberately small enough to check
by hand: three metabolites, four reactions. A metabolic model is a
stoichiometric matrix, the table of how many molecules of each metabolite
each reaction consumes or produces, together with an objective, the one
combination of reaction rates (fluxes) the linear program is asked to
maximise.

```
    metabolites      S_c    O_c    P_c

    reactions
      0  EX_S        S_c ->                    exchange reaction; negative
                                                flux imports the donor from
                                                outside the modelled cell
      1  EX_O        O_c ->                    negative flux imports the
                                                acceptor
      2  EX_P        P_c ->                    positive flux exports the
                                                product
      3  BIO         S_c + 0.5 O_c -> P_c      the objective: this flux IS
                                                growth
```

The linear program, in words before the algebra: maximise the objective flux
(reaction 3, growth), subject to every internal metabolite being at steady
state (nothing accumulates inside the modelled cell), and subject to each
flux staying between a lower and an upper bound. The exchange reactions are
where the outside world enters: their bounds are the uptake limits built from
what diffusion has actually delivered to that voxel, so a voxel starved of O
cannot report a flux for O beyond what is locally present.

```
    maximise   c' v          the objective flux, reaction 3
    subject to S v  = 0      every internal metabolite at steady state
               lb <= v <= ub  uptake bounded by the LOCAL concentrations
```

Because the objective reaction consumes S and 0.5 O for one unit of flux, and
because the exchange reactions are tied to it by the steady-state constraint,
the whole four-reaction system reduces to one line you can check with a
pencil: growth is capped by whichever substrate runs out first,

```
    growth  =  min( Vs , 2 * Vo )
```

where Vs and Vo are the uptake bounds built each step from the local S and O
concentrations. That number comes from the stoichiometry, not from this code,
which is what makes it a test rather than a demonstration.

## 5. What happens each step

The solver repeats this 500 times:

1. **Stream and collide** each species on its own D3Q7 lattice, advancing
   diffusion by one step, and diffuse the biomass field the same way.
2. **Apply the boundary conditions**, so S and O are re-pinned to 1.0 at the
   left face.
3. **Build the uptake bounds**, per voxel, per substrate: a Michaelis-Menten
   term on the local concentration, `V = Vmax * C / (Kc + C)`, capped by
   `<fba_maximum_uptake_flux>`.
4. **Solve the linear program**, in every voxel that holds biomass, with
   those bounds as the exchange constraints. This is the expensive step: one
   full linear-program solve per occupied voxel per iteration, which is why
   this rate path costs far more than a hand-written rate law.
5. **Read the growth rate off the objective flux** and the substrate rates
   off the exchange fluxes, scale both by the local biomass and the time
   step through `<biomass_molar_mass>` (which converts the model's mmol per
   gram dry weight per hour into the solver's mol per litre per second), and
   write the increments into the change lattices.
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
| Growth (from `rate_` fields near the patch) | Positive where the linear program is unconstrained, roughly bounded by `min(Vs, 2*Vo)` |
| S : O consumption ratio | Near 1 : 0.5, the stoichiometry of the objective reaction, wherever the program is unconstrained |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Noticeably slower than a compiled rate law of the same size, because a linear program is solved in every occupied voxel every step |

## 7. What to check

1. **Growth responds to supply.** Raise `<fba_maximum_uptake_flux>` on the
   donor and growth should rise until the acceptor becomes binding, exactly
   where `min(Vs, 2*Vo)` says it will.
2. **Substrates are consumed in the ratio the stoichiometry gives.** The
   objective reaction takes S and O in 1 : 0.5, so wherever the linear
   program is unconstrained the two draws should hold that ratio.
3. **No substrate goes negative.** The uptake bound is built from what is
   locally available, and the increment is clamped again before it is
   applied.
4. **The linear program solves in every voxel that holds biomass.** A failed
   solve is reported in the log, not silently treated as zero growth.

`postprocess.py` reads `output/summary.csv` and `run.log`, reports every
field's change, flags any negative value, and prints the lines from the log
that mention the linear program, so checks 3 and 4 are largely done for you.
Checks 1 and 2 need a comparison across runs or a look at the fields, which is
yours to do.

## 8. What this case demonstrates

**Growth as the answer to an optimisation rather than an input to one.** The
same coupling pattern (build a bound from local transport, hand it to a
solver, read a rate back) recurs in examples 10, 11, 12 and 16; what changes
between them is only who or what turns the bounds into a rate.

A linear program per voxel per step is the most expensive rate path in this
repository by some margin. Example 11 replaces it with a network fitted to
its answers, at a small fraction of the cost; both are trained from sweeps of
exactly this kind of solver.

**What it deliberately leaves out:** flow, an abiotic reaction, and any
change to the pore space. Example 10 is this same case through a different
linear-program solver, COBRApy, and the two must agree to every printed
digit; GLPK is roughly fifty times faster. Example 12 adds a second, hand
written rate path (maintenance) on top of this same organism.

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
./scripts/setup_case.sh 09_fba_glpk run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
Palabos v2.3.0, and GLPK's development package (`libglpk-dev` on
Debian/Ubuntu, `glpk` from Homebrew on macOS). GLPK is a small, fast
open-source solver for linear programs, and this case's `<reaction_type>` is
`glpk`, so the executable must be built with the flag that links it in:

```bash
cmake -B build -S . \
      -DCMAKE_BUILD_TYPE=Release \
      -DPALABOS_ROOT=/path/to/palabos-v2.3.0 \
      -DENABLE_GLPK=ON
cmake --build build -j
```

**When you must recompile.** `-DENABLE_GLPK=ON` only makes the GLPK solver
available in the executable; `<enable_fba_glpk>true</enable_fba_glpk>` in
`CompLaB.xml` is what switches it on for a given run, and that line can be
edited without a rebuild. But the two have to agree: a build without the flag
compiles cleanly and then refuses at start-up with `<enable_fba_glpk> is true
but this executable was built without GLPK`, because turning the XML path on
is not enough if the flag was never passed to CMake. Editing
`input/toy_model.xml` or `input/geometry.dat` also needs no rebuild; only a
change to `CMakeLists.txt` or the solver sources does.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Exports a metabolic model into the flat format the solver reads. This case ships that file already built, so the script does nothing by default; it shows the step you will need with a model of your own. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| | `input/toy_model.xml` | The metabolic model of section 4, in the flat `<Metabolic_Model>` format. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
