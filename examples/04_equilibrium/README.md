# Example 04: aqueous speciation, the carbonate system

Four species over two components, re-equilibrated in every pore voxel
at every step. Acid enters from the left, so the speciation should
shift along the domain.

## What to expect

- the solver reporting the tableau it parsed: 2 components, 4 species
- H2CO3 dominant where the acid is, CO3 dominant where it is not
- **totals conserved**: total carbonate, summed over HCO3, CO3 and
  H2CO3, changes only through transport, never through speciation
- convergence of the Anderson accelerated iteration reported each step

## The switches this case exercises

`<equilibrium><enabled>true`, plus the `<components>`,
`<stoichiometry>` and `<logK>` tableau.

## Notes

This is the most expensive part of the code by a wide margin, which
is why every other example leaves it off. Even on this 3456 voxel
domain you will see it in the step time.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 04_equilibrium run/mycase
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

### The offline code this case ships

| File | What it is |
|---|---|
| `training/makeEquilibrium.py` | Builds the components / stoichiometry / log K tableau the speciation solver reads. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

> **Why the training code is copied here rather than shared.** So that this
> folder *is* the procedure. The cost is real and worth stating: a fix to a
> trainer has to be applied to every case that carries it, and
> `tests/check_repo.sh` fails if a copy drifts from `tools/`.
