# 18 — a graph network, the whole rate vector at once

> **Status.** The graph-network rate processor
> (`src/complab3d_processors_graphnet.hh`) is implemented and covered by the
> regression suite — including a cross-check showing the C++ evaluation and the
> Python trainer agree to 1.36 × 10⁻¹⁸ — but it is **not yet called from
> `src/complab.cpp`**, so the shipped solver does not read `<network_file>`.
> This case is complete and correct as configuration and is the specification
> the wiring has to satisfy. See
> [the note in the top-level README](../../README.md#known-limitations).

Anaerobic oxidation of methane coupled to sulfate reduction:

```
CH4 + SO4²⁻  →  HS⁻ + HCO3⁻ + H2O
```

Four species, one reaction, and **all four rates plus growth returned by one
evaluation** of a message-passing network over the bipartite graph of species
and reactions.

## What this case shows, that 17 does not

The stoichiometry is given to the trainer as **structure**, not learned. The
network is told that R1 consumes CH4 and SO4 one-for-one and produces HS and
HCO3 one-for-one, and it fits rates within that constraint.

The consequence is worth stating plainly: **the returned rates stay in the ratio
−1 : −1 : +1 : +1 without that ratio ever having been a training target.** In
[example 17](../17_symbolic_law/) you get the same guarantee only because you
wrote the other species as multiples of the fitted one by hand. Here it is
structural, so it holds at every point of the input space including the ones you
never checked.

That is what buys the extra cost — 2.5 µs/voxel against 0.11 for the symbolic
law. For four species and one reaction the hand-written multiples of 17 are
fine. For twenty species coupled through eight reactions they are not, and this
is the path.

## The network

`input/aom.gnn`, trained by
[`pipelines/B_offline_models/B4_graph_network/`](../../pipelines/B_offline_models/B4_graph_network/)
from 400 samples of a dual-Monod law. Its own header records what it is:

```
provenance fitted 2026-08-29 from aom_samples.csv
provenance samples 400, scaled mse 2.2148e-02, R 0.96787
provenance source law: dual-Monod, vmax 4e-3 mol/L/h, K_CH4 1e-3, K_SO4 5e-4 mol/L
units per_hour
species CH4 SO4 HS HCO3
reactions R1
rounds 2   width 6   growth 1
stoich 4 1
  -1  -1  +1  +1
```

Because it was fitted to a law with a known closed form, it can be checked
against that law rather than only against held-out samples.
`tests/check_aom.py` does exactly that on 300 samples the network never saw, and
reports R ≈ 0.95–0.96 per species.

## The training box is enforced

`trainmin` and `trainmax` travel **inside** the `.gnn` file and are applied at
every evaluation:

| Species | trainmin | trainmax | This case's initial | Inlet |
|---|---|---|---|---|
| `CH4` | 2.1e-5 | 5.0e-3 | 2.0e-3 | 4.0e-3 |
| `SO4` | 5.2e-5 | 7.9e-3 | 4.0e-3 | 6.0e-3 |
| `HS` | 1.2e-6 | 2.0e-3 | 1.0e-4 | — |
| `HCO3` | 1.1e-4 | 4.0e-3 | 5.0e-4 | — |

Every initial and boundary value sits inside that box **on purpose**. Outside
it, a neural network does not fail — it returns a confident number that is
wrong, with no indication that anything happened. The clamp is what turns that
silent failure into a counted one, and the closing report gives the count.

This is the single concrete way the graph-network and symbolic paths are safer
than the surrogate, which enforces nothing.

## Names are the binding

`aom.gnn` declares `species CH4 SO4 HS HCO3`, and the four substrates in
`CompLaB.xml` carry exactly those names, in that order. A network fitted to
these cannot be applied to a simulation that calls them something else: the run
stops at start-up naming the offending substrate rather than binding to the
wrong lattice.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 18_graph_network run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space: a staggered slot, 24 × 24 × 6, with the biomass patches seeded. Reports porosity and refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Retrains the network from the stoichiometry and samples that ship with this case, writes `input/aom_retrained.gnn`, and cross-validates it. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's balance check. Standard library only. |
| | `pipeline.sh` | Runs the four in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/aom.gnn` | **The trained network** — weights, stoichiometry and the enforced training box, in one text file read at start-up. |
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |

### The offline code this case ships

| File | What it is |
|---|---|
| `training/aom_samples.csv` | The 400 samples `input/aom.gnn` was fitted from. |
| `training/aom_stoich.csv` | The stoichiometric matrix — the STRUCTURE the trainer is given rather than asked to learn. |
| `training/train_graphnet.py` | Trains the message-passing network and writes the `.gnn`. |
| `training/xval_gnn.py` | The parity self-test: evaluates the same network through the Python trainer and the C++ loader and reports the disagreement. It takes no arguments and checks the shipped network, not one you just trained. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

> **Why the training code is copied here rather than shared.** So that this
> folder *is* the procedure. The cost is real and worth stating: a fix to a
> trainer has to be applied to every case that carries it, and
> `tests/check_repo.sh` fails if a copy drifts from `tools/`.
