# 17 — a symbolic rate law, read from a text file

> **Status.** The symbolic rate processor
> (`src/complab3d_processors_symbolic.hh`) is implemented and covered by the
> regression suite, but it is **not yet called from `src/complab.cpp`**, so the
> shipped solver does not read `<expressions_file>`. This case is complete and
> correct as configuration and is the specification the wiring has to satisfy —
> it will run as written the moment that call site exists. See
> [the note in the top-level README](../../README.md#known-limitations).

The same colony as [example 05](../05_biotic_cellular_automaton/), with one
thing changed: **the rate law is not compiled into the program.** It is four
lines of algebra in `input/growth.sym`, parsed once at start-up, printed into
the log, and walked as an expression tree per voxel.

## What this case shows

In 05, the law lives in `defineKinetics.hh`. Changing a half-saturation constant
means editing C++, rebuilding, and re-running — and the number that produced a
result is buried in a header nobody keeps beside the output.

Here the law is a file. A different law is a different file, with no rebuild,
and the law is echoed into the log so **the result carries the law that made
it**. That is the whole argument for this path.

## The law

```
version 1
units   per_hour
vars    acetate o2 Bug

rate    growth  = 0.35 * acetate / (0.05 + acetate) * o2 / (0.01 + o2)
rate    acetate = -2.5 * growth * Bug
rate    o2      = -5.0 * growth * Bug

range   acetate 1e-6 5.0e-3
range   o2      1e-6 2.0e-3
```

Dual Monod on acetate and oxygen. Yield 0.4 mol biomass per mol acetate, two
oxygen per acetate.

## Three things to copy into your own file

**1. Write the `units` line even when you mean `per_second`.** Omit it and the
file is read as per-second. A law fitted from an hourly sweep is then 3600 times
too fast, and nothing warns you — the run simply finishes with the substrate
gone. This is the single most likely way to get a wrong answer on this path.

**2. Carry the biomass explicitly in the substrate lines.** `growth` is a
*specific* rate, 1/time. A substrate rate is *absolute*, mol/L/time. That is why
the organism's name appears in `vars` and is multiplied in above. Leave it out
and the substrate disappears at the same speed whether one cell is present or a
million.

**3. Write a `range` line per variable.** They are enforced at every evaluation,
not advisory, and the closing report says how many evaluations were clamped. A
variable with no range is trusted everywhere — which is exactly the failure mode
the surrogate path has and this one does not.

## The other species are multiples, not separate fits

`acetate` and `o2` are written as multiples of the fitted `growth`, so the
stoichiometric ratios cannot drift apart during the run. Fitting each species
independently would let them.

This is the difference between this case and [18](../18_graph_network/), where
the stoichiometry is supplied to the trainer as structure and the same guarantee
comes for free.

## Names are the binding

Every name in `vars` must match a `<name_of_substrates>` or the microbe's
`<name_of_microbes>` **exactly**:

| In `growth.sym` | In `CompLaB.xml` |
|---|---|
| `acetate` | `<name_of_substrates>acetate</name_of_substrates>` |
| `o2` | `<name_of_substrates>o2</name_of_substrates>` |
| `Bug` | `<name_of_microbes>Bug</name_of_microbes>` |

A mismatch stops the run at start-up naming the offending variable, rather than
binding to the wrong lattice.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 17_symbolic_law run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space: a staggered slot, 24 × 24 × 6, with the biomass patches seeded. Reports porosity and refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Runs the symbolic search on the sample table that ships with this case, and writes `input/growth_discovered.sym` **without installing it** — compare it with the shipped law first. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's balance check. Standard library only. |
| | `pipeline.sh` | Runs the four in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `input/growth.sym` | **The rate law itself** — four lines of algebra, read at start-up. A different law is a different file, with no rebuild. |

### The offline code this case ships

| File | What it is |
|---|---|
| `training/fit_symbolic.py` | The genetic-programming search over expression trees. Returns a Pareto set — one expression per node count — not one answer. |
| `training/growth_samples.csv` | 400 samples of the dual-Monod law `input/growth.sym` describes, with 2% noise on the growth column, so the search has something realistic to work on. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

> **Why the training code is copied here rather than shared.** So that this
> folder *is* the procedure. The cost is real and worth stating: a fix to a
> trainer has to be applied to every case that carries it, and
> `tests/check_repo.sh` fails if a copy drifts from `tools/`.
