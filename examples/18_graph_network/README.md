# 18 — one network over the reaction graph, returning every species rate at once

## 1. What this example does

Anaerobic oxidation of methane coupled to sulphate reduction — the reaction that
consumes most of the methane produced in marine sediment before it ever reaches
the water column:

```
    CH4  +  SO4(2-)   ->   HS(-)  +  HCO3(-)  +  H2O
```

Four species, one reaction, and **all four rates returned by a single
evaluation** of a message-passing network over the bipartite graph of species
and reactions. The network is read from `input/aom.gnn` at start-up, not
compiled in — the same arrangement as the symbolic law in example 17.

**What distinguishes this path from that one is the stoichiometry.** The
stoichiometric matrix is handed to the trainer as *structure*, not learned from
data. The network predicts one **extent** per reaction and the species rates are
formed as

```
    r  =  S xi                    S = ( -1  -1  +1  +1 )'
```

so the four returned rates stay in the ratio −1 : −1 : +1 : +1 exactly, without
that ratio ever having been a training target. In example 17 you get the same
guarantee, but only because you wrote the other species as multiples of the
fitted one by hand. Here it is a property of the readout.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 4 of them: `CH4`, `SO4`, `HS`, `HCO3` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `ANME` — a graph network over the species–reaction graph |
| **Biomass** | `ANME` — attached biofilm (seeded on its own material number), cellular automaton — biomass stays put until a voxel fills, then spills into neighbours |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 1000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

Methane and sulphate are supplied from the left face and the products are held
in — both `HS` and `HCO3` have `closed` faces, genuine bounce-back — so the
sulphate–methane transition zone develops inside the domain rather than being
imposed at a boundary.

`ANME` is one biomass pool standing for the whole consortium. The real system is
an archaeon and a sulphate-reducing bacterium in syntrophy; nothing in this case
resolves the two, and the network was fitted to the consortium's net rates.

## 3. What is inside the `.gnn` file

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

**[v1.3]** Note the shape of the `stoich` block. The keyword is followed by two
integers, the number of species and the number of reactions, and the matrix
entries then follow one per line. Writing the coefficients on the keyword line
instead makes the parser in `src/complab3d_graphnet.hh` read the first two of
them as the shape, the matrix read then fails, and the loader reports the file
as "incomplete or inconsistent" without naming the line at fault. `units
per_hour` is not decoration either: it declares that the fitted rates are per
hour, and carries a factor of 3600 into the conversion the loader applies.

Then the scaling blocks and the weights. The whole thing travels as one file:
the fit, its provenance, its structure and the box it is valid in.

**The species names are the binding.** `aom.gnn` declares `CH4 SO4 HS HCO3` and
the substrates in `CompLaB.xml` carry exactly those names. A network fitted to
these cannot be applied to a simulation that calls them something else — the run
stops at start-up naming the offending substrate. As on the symbolic path,
`<half_saturation_constants>` is unused and set to `0 0 0 0` to say so; the
network carries whatever saturation behaviour it learned.

**`trainmin` and `trainmax` are enforced at every evaluation.** Every initial and
boundary value in this case sits inside that box on purpose. Outside it a neural
network does not fail — it returns a confident number, and there is nothing in
the output to distinguish it from a good one. This is the same warning as
example 11 and it is the single most important thing about any fitted rate path.

**The Damköhler number is an edge input.** `da` enters the network alongside the
concentrations, which is what lets one fitted network span reaction-limited and
transport-limited conditions rather than being valid only at the ratio it was
swept at.

| tag | value | what it does |
|---|---|---|
| `<graphnet><enabled>` | `true` | builds and enables the graph-network path |
| `<network_file>` | `input/aom.gnn` | which network to read |
| `<reaction_type>` | `graphnet` | routes *this organism* through it |
| `<solver_type>` | `CA` | biomass spreads by cellular automaton |
| `<decay_coefficient>` | `0.0` | no decay, so growth is the only biomass term |

## 4. Where a network like this comes from

The case ships `input/aom.gnn`, so `offline.sh` is optional. It is how that file
was made.

You supply two tables:

```
    training/aom_stoich.csv     the stoichiometric matrix — STRUCTURE, given
                                to the trainer rather than learned
    training/aom_samples.csv    400 samples of the rate
```

**The number that matters is not the correlation.** It is the ratio between the
species rates, which must come out *as* the stoichiometry without ever having
been a training target. With the default `--readout extent` that ratio is exact
by construction, so the last line the trainer prints — *predictions off the
stoichiometric subspace* — should be at machine precision. Anything larger means
the writer and the reader disagree about the file format, not that the fit is
poor.

`--readout species` reproduces what every `.gnn` written before v1.2 means: each
species reads its own rate off its own node, and the ratio is only approximately
held. Use it to read an old file, not to make a new one.

The trainer writes its own fit quality into the file's `provenance` lines, so a
network can always be asked where it came from and how well it did:

```bash
grep '^provenance' input/aom.gnn
```

The retrained file is written as `input/aom_retrained.gnn` and is **not**
installed; `CompLaB.xml` still points at the shipped one.

`training/xval_gnn.py` is the other check worth knowing about. Two
implementations of the same forward pass — the C++ reader and the Python trainer
— written at different times in different languages, are two chances to be
wrong. It runs both over the same network and the same inputs and reports the
largest disagreement, drawing half its inputs outside the training box so the
clamp is compared as well as the arithmetic.

## 5. What to check

1. **the provenance lines are echoed at start-up**, and they are the network you
   meant to use;
2. **the stoichiometric ratio holds.** CH₄ and SO₄ consumed one for one, HS and
   HCO₃ produced one for one against them. With the extent readout this is exact,
   so any drift is a reader bug and not a fit quality issue — which is precisely
   what makes it a useful check;
3. **the run stays inside `trainmin`..`trainmax`.** The clamp count is reported;
   a large one means the simulation has wandered outside what the network knows;
4. **the sulphate–methane transition develops inside the domain.** Methane from
   the left, sulphate from the left, products held in — the reaction zone should
   be a band, not a layer plated onto the inlet;
5. **no `[NEG!]` warnings.** The network has no conservation law forcing
   non-negativity; the increment is clamped against what is locally present.

## 6. What this case does not do

No flow, no abiotic reaction, no geometry evolution, and — importantly for this
particular reaction — **no thermodynamic control**. Anaerobic methane oxidation
runs within a few kJ mol⁻¹ of equilibrium, close enough that the free energy
available genuinely limits the rate, and a network fitted without that term will
keep predicting a rate in conditions where the reaction cannot proceed at all.
**Example 19 is this reaction with the thermodynamic factor**, and it is the
right companion to this case.

The consortium is one biomass pool, so nothing here resolves the syntrophic
exchange between the two partners.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 18_graph_network run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and reports porosity, refusing to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Optional. Trains a graph network from the two shipped CSVs and reports its fit and its stoichiometric residual. The case ships a network, so nothing has to be trained before running it. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and checks the stoichiometric ratio of section 5. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `input/aom.gnn` | The fitted network of section 3 — structure, weights, training box and provenance in one file. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

### The offline code this case ships

| File | What it is |
|---|---|
| `training/train_graphnet.py` | Trains a network from a stoichiometric matrix and a sample table, and writes the `.gnn`. |
| `training/aom_stoich.csv` | The stoichiometric matrix, given to the trainer as structure rather than learned. |
| `training/aom_samples.csv` | The 400 samples `input/aom.gnn` was fitted from. |
| `training/xval_gnn.py` | Runs the C++ reader and the Python trainer over the same network and reports the largest disagreement. Needs numpy and a compiler; no Palabos. |

### A note on the domain

`x = 0` and `x = nx-1` carry the boundary conditions named per substrate. The
solver gives the other four faces nothing — not a wall, not a symmetry plane, not
periodicity — so `preprocess.py` draws an inert wall there, and `NY` and `NZ` are
two larger than the pore space they hold. The layers are **added**, not taken out
of the pore space, so porosity and every count are what the case declares.

> **Why shared code is copied here rather than referenced.** So that this folder
> *is* the procedure. The cost is real and worth stating: a fix to a shared tool
> has to be applied to every case that carries it, and `tests/check_repo.sh`
> fails if a copy drifts from `tools/`.
