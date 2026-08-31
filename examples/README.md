# Examples

Eighteen cases, arranged so that each one adds a single thing to the one before
it. If you are new to the model, run them in order until something breaks in a
way you do not understand — that is the piece worth reading about.

## How to run one

Every example holds its **own** pre-processing, offline and post-processing code
— not a pointer to a shared tool, the actual script for that case. Assemble it
and the directory is the whole pipeline:

```bash
./scripts/setup_case.sh 13_precipitation run/mycase
cd run/mycase
./pipeline.sh
```

`pipeline.sh` runs four steps, and each one is a file in that directory you can
read and run on its own:

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds that case's pore space, reports porosity, and refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Only in the seven cases that need it: the model export, the surrogate fit, the symbolic search, the network training. |
| **run** | `CompLaB.xml` | What the solver reads. |
| **post** | `postprocess.py` | Reads the summary CSV and the log, reports every field's change, flags negatives, and runs that case's balance check. Standard library only. |

Step by step instead:

```bash
python3 preprocess.py
./offline.sh                                    # where the case has one
cmake -B build -S . && cmake --build build -j
./build/complab CompLaB.xml 2>&1 | tee output/run.log
python3 postprocess.py
```

An example folder holds **everything that case needs**: its pore space, the
metabolic model it reads, the training code that produced whatever it was
trained on, its rate-law file, its chemistry, and its own pre- and
post-processing. Nothing has to be fetched from elsewhere in the tree.

`setup_case.sh` therefore does very little — it lays the two shared kinetics
defaults down, copies the case folder whole on top, and adds the solver
sources:

| | Pass | Where it comes from |
|---|---|---|
| 1 | the default chemistry, so the case always compiles | `config/kinetics/*.default.hh` |
| 2 | the case, whole | `examples/NN_case/` — `CompLaB.xml`, `input/`, `models/`, `training/`, `kinetics/`, and its four scripts |
| 3 | the solver sources | `src/`, `CMakeLists.txt` |

### What is inside a case folder

| | Holds |
|---|---|
| `input/` | What the solver reads: the geometry, and the `.sym` or `.gnn` or metabolic model where the case has one |
| `models/` | The genome-scale model, for the cases that sweep or export one |
| `training/` | The actual training code that produced what the case was trained on, plus the data it was trained on |
| `kinetics/` | This case's own chemistry, where it differs from the shared default |

**Duplication inside `examples/` is deliberate.** Several cases carry the same
geometry, and three carry the same model exporter. The point is that a case
folder is the whole procedure rather than a set of pointers. The cost is real: a
fix to a trainer has to be applied to every case that carries it, and
`tests/check_repo.sh` fails if a copy drifts out of step with `tools/`. Outside
`examples/`, the one-copy rule still holds.

Ten of the eighteen carry a chemistry header of their own; the rest use the
shared defaults unchanged.

## The eighteen

| | Case | What it adds | Offline work first |
|---|---|---|---|
| 01 | `flow_only` | Navier–Stokes on a pore space, nothing else | — |
| 02 | `diffusion_only` | one solute, transported, no reaction | — |
| 03 | `abiotic_kinetics` | a chemical rate law, no organisms | — |
| 04 | `equilibrium` | aqueous speciation by continued fractions | — |
| 05 | `biotic_cellular_automaton` | biomass, spread by the automaton | — |
| 06 | `biotic_finite_difference` | the same, by finite difference | — |
| 07 | `biotic_lattice_boltzmann` | the same, on a D3Q7 lattice | — |
| 08 | `two_microbes` | two populations competing for one substrate | — |
| 09 | `fba_glpk` | growth from a genome-scale linear program | **B1** |
| 10 | `fba_cobrapy` | the same through the reference implementation | **B1** |
| 11 | `surrogate` | a fitted network in place of the linear program | **B1 → B2** |
| 12 | `mixed_reaction_types` | one population on FBA beside one on kinetics | **B1** |
| 13 | `precipitation` | FeS fills the pore and seals it | — |
| 14 | `dissolution` | calcite is eaten away and the pore reopens | — |
| 15 | `precip_and_dissolution` | both at once, on different phases | — |
| 16 | `complete_pipeline` | geometry, chemistry, biology and metabolism together | **B1** |
| 17 | `symbolic_law` | the rate law read from a text file, no rebuild | — |
| 18 | `graph_network` | the whole coupled rate vector from one evaluation | — |

The offline column points into [`../pipelines/`](../pipelines/). Twelve of the
eighteen need nothing prepared at all — 17 and 18 included, because both ship
the file they read.

## Where the six rate paths appear

| Rate path | Example |
|---|---|
| Compiled kinetics | 03, 05–08 |
| Flux balance, GLPK | 09, 12, 16 |
| Flux balance, COBRApy | 10 |
| Surrogate network | 11 |
| Symbolic law | **17** |
| Graph network | **18** |

17 and 18 are driven by a file rather than a build setting, which is why neither
needs offline work to run: the `.sym` and the `.gnn` both ship. Point any of
05–08 at one and it will use that instead of its compiled kinetics.

## Which one answers my question?

**"Does my geometry work?"** → 01. If the flow does not percolate, nothing else
will either.

**"Is my chemistry right?"** → 03, with the organisms switched off. A rate law
that misbehaves is much easier to see without biomass moving underneath it.

**"Why is my pore not clogging?"** → 13. Compare its `max_precipRho` with yours;
the usual cause is a fill density derived from the wrong molar volume.

**"Why does my dissolution lose mass?"** → 14, run on one process. Then read the
MPI limitation in [`../pipelines/C_run/README.md`](../pipelines/C_run/README.md).

**"Is the surrogate worth it?"** → run 09 and 11 on the same geometry and
compare the biomass fields. That difference is what the approximation costs you.

**"Symbolic law or graph network?"** → 17 and 18 are the same question asked
twice. 17 keeps the stoichiometry right because you wrote it that way; 18 keeps
it right structurally, at 20× the cost per evaluation. With four species, use
17. With twenty, use 18.

## Which geometry each case runs on

Four geometries serve all eighteen, and they are described in
[`../config/geometry/README.md`](../config/geometry/README.md).

| Geometry | Cases |
|---|---|
| `slot_bare.dat` — wall and pore only | 01–04 |
| `slot_one_biofilm.dat` — one seeded population | 05–07, 09–12, 16–18 |
| `slot_two_biofilms.dat` — two populations facing each other | 08 |
| `slot_mineral_phase.dat` — a declared solid phase on the walls | 13–15 |

All four are 24 × 12 × 12, which runs in seconds. Scale up only once the case
does what you expect at this size.

## A note on `CompLaB.xml`

Each case has its own, and they are all genuinely different — that is the point
of having eighteen. None of them is a reference: every tag is documented once, in
[`../config/CompLaB.reference.xml`](../config/CompLaB.reference.xml), and
nowhere else.
