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

## What to run

```bash
./scripts/setup_case.sh 18_graph_network run/gnn
cd run/gnn
cmake -B build -S . && cmake --build build -j
./build/complab CompLaB.xml
```

| | |
|---|---|
| **Inputs** | `input/geometry.dat` (shared, `slot_one_biofilm`), `input/aom.gnn` (shared, from stage B4) |
| **Offline work first** | none to run this case — the trained network ships. To train your own, stage **B4** |
| **Outputs** | VTI fields for `CH4`, `SO4`, `HS`, `HCO3`, `ANME`; the log, carrying the network's provenance and the clamp count |
| **Runtime** | seconds, at 24 × 24 × 6 |

## What to look for

**The stoichiometric ratio in the output.** Integrate the four fields over the
domain between two writes. CH4 and SO4 consumed should match HS and HCO3
produced, one for one, to within the transport error. If they do not, the graph
in `stoich` is wrong — not the fit.

**The clamp count at the end.** Rising through the run means the simulation is
walking out of the box the network was trained in. Concentrations drift as the
colony grows, so a case that starts inside the box does not necessarily stay
there.

**The `[GNN]` provenance block in the log.** It carries the training command,
the sample count and the fit quality, so the result is traceable to the data
that produced it.

## Training your own

```bash
cd pipelines/B_offline_models/B4_graph_network
STOICH=mystoich.csv SAMPLES=mysamples.csv OUT=mynetwork.gnn ./run.sh
```

You supply two CSVs — the stoichiometric matrix and the samples — and a
Damköhler number per reaction, which enters as an edge input so one network can
span reaction-limited and transport-limited regimes. `run.sh` trains, then
immediately cross-validates against the samples.

The number that matters in that output is **not** the correlation. It is the
ratio between species rates, which should come out close to the stoichiometry
without ever having been a training target. If it does not, the graph is wrong.
