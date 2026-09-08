# 18 - Graph network: one network over the reaction graph, returning every species rate at once

## 1. The scenario

Anaerobic oxidation of methane coupled to sulphate reduction, the reaction
that consumes most of the methane produced in marine sediment before it ever
reaches the water column:

```
    CH4  +  SO4(2-)   ->   HS(-)  +  HCO3(-)  +  H2O
```

Four species, one reaction, and all four rates returned by a single
evaluation of a **graph neural network**: a network built to operate on a
graph of nodes and edges rather than on a fixed-length vector, here a
bipartite graph with one node per species and one per reaction, edges
connecting a species to the reactions it takes part in. The network is read
from `input/aom.gnn` at start-up, not compiled in, the same arrangement as
the symbolic law in example 17.

The network runs by **message passing**: in each of a fixed number of
rounds, every node collects a message from each of its graph neighbours,
combines it with its own current state, and updates that state. After a few
rounds, a reaction node's state has been shaped by every species that feeds
into it, and a species node's state carries information from every reaction
it participates in, without either having been given the other's identity
directly, only the graph structure connecting them.

**What distinguishes this path from the symbolic one is the stoichiometry.**
The stoichiometric matrix is handed to the trainer as structure, not learned
from data. The network predicts one **extent** per reaction (a single number
describing how far the reaction has proceeded) and the species rates are
formed by multiplying that extent through the known stoichiometric matrix:

```
    r  =  S xi                    S = ( -1  -1  +1  +1 )'
```

so the four returned rates stay in the ratio -1 : -1 : +1 : +1 exactly,
without that ratio ever having been a training target. In example 17 you get
the same guarantee, but only because the other species were written as
multiples of the fitted one by hand. Here it is a property of the readout
itself.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
    |....................  ####  ..........................  |
CH4 |...................  ####  ..........................  |
4e-3|[ANME][ANME].........  ########  ......................  |closed:
SO4 |[ANME][ANME].........                                     |HS, HCO3
6e-3                                                           |made inside
held|                                                         |
    +---------------------------------------------------------+
      no flow: Peclet = 0. CH4 and SO4 diffuse in from x=0
      [ANME] = the seeded biomass patch, 108 of the 3168 open voxels,
               standing for the whole methane/sulphate consortium

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition for
`CH4` and `SO4`; `HS` and `HCO3` are held closed on both ends, genuine
bounce-back, so the sulphate-methane transition zone develops inside the
domain rather than being imposed at a boundary. The biomass patch sits on
material number 3, at every `x` where `x % 4 == 0`, the same 108 voxels used
throughout this staggered-slot geometry. `ANME` is one biomass pool standing
for the whole consortium; the real system is an archaeon and a
sulphate-reducing bacterium in syntrophy, and nothing here resolves the two.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, four species: `CH4`, `SO4`, `HS`, `HCO3` |
| **Diffusivity** | 2.0e-9 / 1.0e-9 m2/s for CH4 (pore / biofilm); 1.0e-9 / 5.0e-10 for SO4; 1.5e-9 / 7.5e-10 for HS; 1.2e-9 / 6.0e-10 for HCO3 |
| **Chemistry** | none compiled in; a graph network read from `input/aom.gnn` |
| **Biology** | `ANME`, one population, seeded on material 3 at initial density 1.0e-5, cellular automaton biomass movement, rate routed through the graph-network path |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 1000 advection-diffusion steps, output every 200 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **CH4** | Dirichlet, held at 4.0e-3 | Neumann, zero gradient | 2.0e-3 |
| **SO4** | Dirichlet, held at 6.0e-3 | Neumann, zero gradient | 4.0e-3 |
| **HS** | closed | closed | 1.0e-4 |
| **HCO3** | closed | closed | 5.0e-4 |
| **ANME** (biomass) | closed | closed | 1.0e-5 in the seeded patch, 0 elsewhere |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, a reservoir. *Neumann with a zero gradient*
means nothing flows across that face, an open outlet neither substrate
reaches. *Closed* means the same thing on both ends: the two products, and
biomass, stay inside the domain entirely.

## 4. The reaction

`input/aom.gnn` is a text file. Its header is the part worth reading before
anything else:

```
    provenance fitted 2026-09-01 from training/aom_samples.csv
    provenance samples 400, readout extent, scaled mse 1.0591e-02, R 0.97196
    provenance stoichiometric residual 0 of the rate norm
    version 1
    units per_hour
    species   CH4 SO4 HS HCO3
    reactions R1
    rounds 2          message-passing rounds
    width 6           hidden width per node
    growth 1
    readout extent    one extent per reaction, species rates formed as S xi
    stoich 4 1        the matrix shape: 4 species, 1 reaction
    -1                then one coefficient per line, species-major
    -1
    1
    1
    da 1
    trainmin  2.08e-05  5.20e-05  1.17e-06  1.08e-04
    trainmax  4.99e-03  7.93e-03  2.00e-03  3.98e-03
```

The keyword `stoich` is followed by two integers, the number of species and
the number of reactions, and the matrix entries then follow one per line.
`units per_hour` is not decoration: it declares that the fitted rates are
per hour, and carries a factor of 3600 into the conversion the loader
applies. Then come the scaling blocks and the weights. The whole thing
travels as one file: the fit, its provenance, its structure and the box it
is valid in.

**The species names are the binding.** `aom.gnn` declares `CH4 SO4 HS HCO3`
and the substrates in `CompLaB.xml` carry exactly those names. A network
fitted to these cannot be applied to a simulation that calls them something
else; the run stops at start-up naming the offending substrate. As on the
symbolic path, `<half_saturation_constants>` is unused and set to `0 0 0 0`
to say so; the network carries whatever saturation behaviour it learned.

**`trainmin` and `trainmax` are enforced at every evaluation.** Every initial
and boundary value in section 3 sits inside that box on purpose. Outside it,
a neural network does not fail, it returns a confident number, and there is
nothing in the output to distinguish it from a good one. This is the single
most important thing to know about any fitted rate path.

**The Damkohler number is an edge input.** `da` enters the network alongside
the concentrations, which is what lets one fitted network span
reaction-limited and transport-limited conditions rather than being valid
only at the ratio it was swept at.

| tag | value | what it does |
|---|---|---|
| `<graphnet><enabled>` | `true` | builds and enables the graph-network path |
| `<network_file>` | `input/aom.gnn` | which network to read |
| `<reaction_type>` | `graphnet` | routes this organism through it |
| `<solver_type>` | `CA` | biomass spreads by cellular automaton |
| `<decay_coefficient>` | `0.0` | no decay, so growth is the only biomass term |

## 5. What happens each step

The solver repeats this 1000 times:

1. **Stream and collide** `CH4`, `SO4`, `HS` and `HCO3` on their D3Q7
   lattices, which advances diffusion by one step. There is no velocity
   field to advect them with.
2. **Apply the boundary conditions**, re-pinning `CH4` to 4.0e-3 and `SO4`
   to 6.0e-3 at the left face.
3. **Evaluate the graph network** in every open voxel that holds biomass:
   clamp each input to `trainmin`..`trainmax`, run the message-passing
   rounds, read one extent off the reaction node, and form the four species
   rates as `r = S xi`.
4. **Apply the biomass increment** in place.
5. **Run the cellular-automaton rule.** Wherever a voxel's biomass now
   exceeds `<maximum_biomass_density>`, spill the excess into an open
   neighbouring voxel.

The stability caveat is on step 3: the network has no conservation law
forcing non-negativity on its own, so the increment is applied as evaluated,
with no positivity clamp. A time step too large for the local concentrations
produces a `[NEG!]` warning exactly as the other reaction paths do.

## 6. What comes out

```
output/
  CH4_0000200.vti  ...  CH4_0001000.vti     concentration of methane
  SO4_*.vti  HS_*.vti  HCO3_*.vti             the other three species
  ANME_*.vti                                   biomass density of the consortium
  rate_CH4_*.vti  rate_SO4_*.vti  rate_HS_*.vti  rate_HCO3_*.vti   the reaction rate as a field, mol/L/s
  summary.csv                                  one row every 100 steps
  run.log                                      the whole run, including the [GNN] block
```

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Total biomass | Rises from the seeded patch and levels off |
| Total CH4, total SO4 | Fall near the colony, replenished from the left boundary |
| Total HS, total HCO3 | Both rise from their initial values, since both are closed and made only where biomass is |
| HS produced against HCO3 produced | Track each other one for one; their difference should stay near zero |
| Clamp count in the `[GNN]` block | Low, since the boundary and initial values sit inside `trainmin`..`trainmax` deliberately |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Seconds on one core. This is a small case on purpose |

## 7. What to check

1. **The provenance lines are echoed at start-up**, and they are the network
   you meant to use.
2. **The stoichiometric ratio holds.** CH4 and SO4 consumed one for one, HS
   and HCO3 produced one for one against them. With the extent readout this
   is exact, so any drift is a reader bug rather than a fit-quality issue,
   which is precisely what makes it a useful check.
3. **The run stays inside `trainmin`..`trainmax`.** The clamp count is
   reported; a large one means the simulation has wandered outside what the
   network knows.
4. **The sulphate-methane transition develops inside the domain.** Methane
   and sulphate both from the left, products held in; the reaction zone
   should be a band, not a layer plated onto the inlet.
5. **No `[NEG!]` warnings.**

`postprocess.py` reports the change in every field's total and flags any
negative minimum, covering check 5. It also runs the balance check between
`HS` and `HCO3`, reporting a PASS or FAIL residual relative to the largest
change, which covers check 2 directly. It greps `run.log` for the `[GNN]`
block so checks 1 and 3 are read straight off that output. Check 4 is a look
at the `.vti` output. `CH4` and `SO4` are each fed from a Dirichlet
boundary, so neither total is closed, and this case's `<conserve>`
declaration, if any, would only ever name the two closed products.

## 8. What this case demonstrates

**One network evaluation returning every species rate at once, with the
stoichiometry guaranteed by the readout rather than by hand.** Where example
17 gets a correct ratio because a person wrote the substrate lines as
multiples of the fitted rate, this case gets it because the network can only
ever produce a rate vector that lies in the span of the stoichiometric
matrix it was given as structure. That guarantee costs nothing at run time
and cannot be violated by a bad fit, only by a bad stoichiometric matrix.

Where a network like this comes from: the case ships `input/aom.gnn`, so
`offline.sh` is optional. You supply two tables, `training/aom_stoich.csv`,
the stoichiometric matrix given to the trainer as structure rather than
learned, and `training/aom_samples.csv`, 400 samples of the rate. The number
that matters is not the correlation, it is the ratio between the species
rates, which must come out as the stoichiometry without ever having been a
training target; with the default extent readout that ratio is exact by
construction, so the trainer's own "predictions off the stoichiometric
subspace" line should be at machine precision. The retrained file is written
as `input/aom_retrained.gnn` and is not installed automatically.
`training/xval_gnn.py` cross-checks the C++ reader in the solver against the
Python trainer's own forward pass over the same network and inputs, drawing
half its inputs outside the training box so the clamp is compared as well
as the arithmetic, two independent implementations of the same forward pass
being two chances to be wrong rather than one.

**What it deliberately leaves out.** No flow, no abiotic reaction, no
geometry evolution, and, importantly for this particular reaction, no
thermodynamic control. Anaerobic methane oxidation runs within a few
kilojoules per mole of equilibrium, close enough that the free energy
available genuinely limits the rate, and a network fitted without that term
will keep predicting a rate in conditions where the reaction cannot proceed
at all. Example 19 is this reaction with a thermodynamic factor added, and
it is the right companion to this case. The consortium is also one biomass
pool, so nothing here resolves the syntrophic exchange between the two
partner organisms.

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
./scripts/setup_case.sh 18_graph_network run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. No optional solver is needed. Cases 09, 12, 21 and 22
additionally need GLPK, and case 10 needs Python with COBRApy; this one
needs neither, so the plain cmake line is enough:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** This case has no chemistry compiled in: the
network lives in `input/aom.gnn` and is read at start-up, so editing it, or
`CompLaB.xml`, or the geometry, never needs a rebuild. A rebuild is only
needed if the solver source itself changes, such as the network reader in
`src/complab3d_graphnet.hh`.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and the seeded biomass patch, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Optional. Trains a graph network from the two shipped CSVs and reports its fit and its stoichiometric residual. The case ships a network, so nothing has to be trained before running it. Runs once, before the build. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/aom.gnn` | The fitted network of section 4: structure, weights, training box and provenance in one file. |
| | `input/geometry.dat` | The pore space and the patch. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and checks the stoichiometric ratio of section 7. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |
| **train** | `training/train_graphnet.py` | Trains a network from a stoichiometric matrix and a sample table, and writes the `.gnn`. |
| | `training/aom_stoich.csv` | The stoichiometric matrix, given to the trainer as structure rather than learned. |
| | `training/aom_samples.csv` | The 400 samples `input/aom.gnn` was fitted from. |
| | `training/xval_gnn.py` | Runs the C++ reader and the Python trainer over the same network and reports the largest disagreement. Needs numpy and a compiler; no Palabos. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top; this
case does not ship a `kinetics/` header of its own, since its rate law comes
from `input/aom.gnn` instead.
