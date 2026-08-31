# Example 14: dissolution and pore reopening

Acid enters from the left and eats into four calcite grains. As each
grain thins, its outer voxels drop below `reopen_fraction` of full
density and reopen to flow.

## What to expect

- the acid front eating into the **upstream** faces of the grains
  first, because that is where the acid arrives
- Ca2+ leaving on the right
- the porosity in the log **rising in steps**, the mirror of
  example 13
- **mass conserved**: every mole of calcite removed appears as a mole
  of Ca2+ in the water

## The switches this case exercises

`<dissolution><enabled>true`, one `<phase0>` block with a non zero
`<initial_fill>`, and `defineDissolutionRate()` in
`defineAbioticKinetics.hh`.

## Notes

`<initial_fill>27.1</initial_fill>` is what makes the ORIGINAL grains
dissolvable rather than only material that precipitated first. A
precipitate phase would use `<initial_fill>0</initial_fill>` and
`<is_precipitate>true</is_precipitate>` instead.

Nothing precipitates in this case. `defineAbioticRxnKinetics` is
empty on purpose, so the only thing moving mineral is the dissolution
hook and the mass balance is unambiguous.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 14_dissolution run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space: a staggered slot, 24 × 24 × 6, with the reactive phase lining the grain faces. Reports porosity and refuses to continue if it does not percolate. Standard library only. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's balance check. Standard library only. |
| | `pipeline.sh` | Runs the four in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |

### This case's own chemistry

| File | What it is |
|---|---|
| `kinetics/defineAbioticKinetics.hh` | This case's own chemistry, compiled in on top of the shared default. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

> **Why the training code is copied here rather than shared.** So that this
> folder *is* the procedure. The cost is real and worth stating: a fix to a
> trainer has to be applied to every case that carries it, and
> `tests/check_repo.sh` fails if a copy drifts from `tools/`.
