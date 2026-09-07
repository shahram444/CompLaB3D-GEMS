# Examples

Twenty-two cases, arranged so that each one adds a single thing to the one
before it. If you are new to the model, run them in order until something breaks
in a way you do not understand — that is the piece worth reading about.

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
| **offline** | `offline.sh` | Only in the eight cases that have one: the model export, the surrogate fit, the symbolic search, the network training, the gate curve. |
| **run** | `CompLaB.xml` | What the solver reads. |
| **post** | `postprocess.py` | Reads the summary CSV and the log, reports every field's change, flags negatives, and runs that case's balance check. Standard library only. |

Step by step instead:

```bash
python3 preprocess.py
./offline.sh                                    # where the case has one
cmake -B build -S . && cmake --build build -j    # see the flag note below
./complab CompLaB.xml 2>&1 | tee output/run.log
python3 postprocess.py
```

**[v1.3] That bare `cmake` line is not enough for every case.** Both back-end
options in `CMakeLists.txt` are OFF by default, so add the one the case needs:
`-DENABLE_GLPK=ON` for cases 09, 11, 12 and 16, and `-DENABLE_COBRAPY=ON` for
case 10. Every other case builds with no flag. Each case's own `pipeline.sh`
already passes the right one, so this only bites a hand build.

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

Eleven of the twenty-two carry a chemistry header of their own — 03, 05–08, 12,
13–15, 19 and 20 — and the rest use the shared defaults unchanged.

## The twenty-two

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
| 19 | `thermodynamic_gate` | the rate multiplied by the free energy available for it | — |
| 20 | `upscaling` | one resolved aggregate reduced to a continuum rate | — |
| 21 | `multistep_fba` | a chain of linear programs, so the by-product fluxes are determined | — |
| 22 | `cybernetic_switching` | the organism changes carbon source when the first runs out | — |

The offline column points into [`../pipelines/`](../pipelines/). Fourteen of the
twenty-two need nothing prepared at all, and 17, 18, 19 and 20 are among them because
all four ship the file they read. 21 and 22 are among them for a different
reason: both point `<model_filename>` straight at `input/e_coli_core.xml`, the
SBML as it is distributed, which the solver reads natively, so there is no flat
model to export first. 17, 18 and 19 have an offline STEP all the
same, and each is optional in a different way: 17 and 18 show how their shipped
law and network were made, while 19's reports where the gate closes **before**
the solver starts — the one thing about that case which cannot be checked
afterwards.

20 has an offline step of a different kind, and it is not part of its
`pipeline.sh` on purpose. One run measures the effectiveness factor at one
aggregate size and one bulk composition; `offline/upscale.py` sweeps both and
fits the curve a continuum model actually needs. That is one solver run per
point, so it is a command you choose to give rather than one a pipeline gives
for you.

## Where the six rate paths appear

| Rate path | Example |
|---|---|
| Compiled kinetics | 03, 05–08 |
| Flux balance, GLPK | 09, 12, 21, 22 |
| Flux balance, COBRApy | 10 |
| Surrogate network | 11, 16 |
| Symbolic law | **17** |
| Graph network | **18** |
| *(the thermodynamic factor, on top of any of them)* | **19**, **20** |

([v1.3] Case 16 was listed on the GLPK row above; it belongs on the surrogate
row. Its `CompLaB.xml` sets `<enable_fba_glpk>false</enable_fba_glpk>`,
`<enable_surrogate>true</enable_surrogate>` and
`<reaction_type>surrogate</reaction_type>`: it calls GLPK once only, at
start-up, to train the surrogate, and every transport step after that evaluates
the trained network instead. It therefore still needs `-DENABLE_GLPK=ON` at
build time, for that start-up training alone.)

17 and 18 are driven by a file rather than a build setting, which is why neither
needs offline work to run: the `.sym` and the `.gnn` both ship. Point any of
05–08 at one and it will use that instead of its compiled kinetics.

19 is not a seventh row of that table. It is a factor between zero and one,
computed per voxel from the local free energy, that multiplies whatever the
chosen path returned. Add `<thermodynamics>` and a `.thm` file to any of the
cases above and the factor applies there too.

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

**"Why do my acetate and formate fluxes change when I change solver?"** → 21.
The growth optimum is unique; the by-product fluxes that carry the same
electrons are not, so a single program returns whichever end its pivoting rule
reached. 21 pins them one at a time and the answer stops moving.

**"My organism eats both substrates at once, and real cells do not."** → 22. One
program has no way to prefer glucose over acetate, because both raise the same
objective. 22 solves one program per carbon source with the others' uptake shut,
and blends the answers by how much carbon each source is supplying.

**"What rate do I give my column model?"** → 20. It reduces one resolved
aggregate to the effectiveness factor a continuum grid block needs, and reports
it against the two dimensionless groups that set it.

## Which geometry each case runs on

Each case's `preprocess.py` writes its own `input/geometry.dat`; the shapes they
draw fall into five families, described in
[`../config/geometry/README.md`](../config/geometry/README.md).

| Geometry | Cases | Size |
|---|---|---|
| staggered slot, wall and pore only | 01–04 | 24 × 26 × 8 |
| the same slot with one seeded population | 05–07, 09–12, 16–18, 21, 22 | 24 × 26 × 8 |
| the same slot with two populations facing each other | 08 | 24 × 26 × 8 |
| the same slot with a declared solid phase on the grain faces | 13–15 | 24 × 26 × 8 |
| a single aggregate in open water | 19, 20 | 32 × 32 × 10, 32 × 34 × 10 |

The slot cases run in seconds; the two aggregate cases take minutes. Scale up
only once a case does what you expect at these sizes.

**Why `ny` and `nz` are two larger than the shape drawn.** The solver conditions
only `x = 0` and `x = nx-1`; the other four faces get nothing — not a wall, not a
symmetry plane, not periodicity. Every `preprocess.py` therefore **adds** an
inert wall layer on those four faces. The pore space itself is untouched, so
porosity and every voxel count are what the case declares.

## Every case README has the same shape

Not by convention — the "What is simulated" table in each one is generated from
that case's own `CompLaB.xml` and `input/geometry.dat`, so it cannot describe
something the case does not do.

| Section | What is in it |
|---|---|
| 1 | what the case is for, in plain language, and the reaction it runs |
| 2 | **What is simulated** — geometry, flow, solute transport, abiotic and biotic reaction, how biomass is solved, geometry evolution, run length |
| 3 … | the physics of this case: the rate law, the equations, every constant with its unit and its source |
| *n* | **What to check** — what a correct run looks like, in the order worth checking it |
| last | what the case does **not** do, and which case does that instead |

then the folder contents: the pipeline, what the solver reads, what it inherits.

## A note on `CompLaB.xml`

Each case has its own, and they are all genuinely different — that is the point
of having twenty-two. None of them is a reference: every tag is documented once, in
[`../config/CompLaB.reference.xml`](../config/CompLaB.reference.xml), and
nowhere else.
