# Example 03: abiotic kinetics, A + B to C

A second order irreversible reaction with no biology. A enters from
the left, B from the right, and C forms where they overlap.

## What to expect

- a peak in C near the middle of the domain, with A and B drawn down
  around it
- the reaction front stationary, because there is no flow to push it
- **mass conserved**: total A + total C constant, total B + total C
  constant

## The switches this case exercises

`enable_abiotic_kinetics`, and a hand written rate law in
`defineAbioticKinetics.hh`.

## Notes

Turn `<Peclet>` up to 1 and the front moves downstream. That single
change is the cheapest way to see transport and reaction competing.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 03_abiotic_kinetics run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space: a staggered slot, 24 × 24 × 6, wall and pore only. Reports porosity and refuses to continue if it does not percolate. Standard library only. |
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
