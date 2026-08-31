# Examples

Eighteen cases, arranged so that each one adds a single thing to the one before
it. If you are new to the model, run them in order until something breaks in a
way you do not understand — that is the piece worth reading about.

## How to run one

An example directory is not runnable as it stands, and deliberately so: it holds
only what is unique to that case. Assemble a working directory instead:

```bash
./scripts/setup_case.sh 13_precipitation run/mycase
cd run/mycase
cmake -B build -S . && cmake --build build -j
./build/complab CompLaB.xml
```

`setup_case.sh` works in three passes, each overwriting the last, so a case can
always override anything it inherits:

| | Pass | Where it comes from |
|---|---|---|
| 1 | the default chemistry, so the case always compiles | `config/kinetics/*.default.hh` |
| 2 | whatever this case borrows | the paths listed in the case's own `case.files` |
| 3 | whatever this case owns | `examples/NN_case/kinetics/`, `CompLaB.xml`, `input/` |

Six of the eighteen carry a chemistry header of their own; the rest use the
shared ones unchanged. Every case borrows a geometry, the three FBA cases also
borrow the toy metabolic model, and case 18 borrows the trained graph network
from pipeline stage B4.

### `case.files`

Each case has one. It is two columns — the path from the repository root, and
where that file lands in the assembled directory — with a comment above saying
what each borrowed file is:

```
config/geometry/slot_one_biofilm.dat  input/geometry.dat
models/toy_model.xml                  input/toy_model.xml
```

Read it when you want to know what a case is actually made of. Edit it if you
want a case to start from a different geometry permanently; to do it once, just
replace `input/geometry.dat` in the assembled directory instead.

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
