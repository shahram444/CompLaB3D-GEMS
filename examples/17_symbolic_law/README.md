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

## What to run

```bash
./scripts/setup_case.sh 17_symbolic_law run/sym
cd run/sym
cmake -B build -S . && cmake --build build -j
./build/complab CompLaB.xml
```

| | |
|---|---|
| **Inputs** | `input/geometry.dat` (shared, `slot_one_biofilm`), `input/growth.sym` (this case's own) |
| **Offline work first** | none — unlike the surrogate, nothing has to be prepared |
| **Outputs** | VTI fields for `acetate`, `o2`, `Bug`; the log, carrying the law and the clamp count |
| **Runtime** | seconds, at 24 × 24 × 6 |

## What to look for

The `[SYM]` block near the top of the log. It prints the provenance comment, the
resolved units, every expression as parsed, and the valid range of each
variable. Read it before believing the run: an expression that parsed
differently from how you read it shows up here, and nowhere else.

At the end, the clamp count. A large fraction of evaluations clamped means the
simulation is spending its time outside the box the law was fitted in, and the
answer is an extrapolation whatever the numbers look like.

## Where the file comes from

`input/growth.sym` is written out by hand here so the case is readable without
running anything. In a real study it comes from
[`pipelines/B_offline_models/B3_symbolic_law/`](../../pipelines/B_offline_models/B3_symbolic_law/):
a genetic-programming search over expression trees returns a **Pareto set** —
one expression per node count, not one answer. You take the elbow: the shortest
expression whose error is acceptable. Then you finish the file by hand, adding
the other species as multiples and the `range` lines.

`pipelines/B_offline_models/B3_symbolic_law/expected/` holds three real outputs
of that search, including `abiotic.sym` — a law with no organism at all, pointed
at by `<abiotic_file>` instead, which fires in every fluid voxel rather than
only where biomass is.
