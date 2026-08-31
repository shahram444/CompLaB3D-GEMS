# Example 15: precipitation and dissolution together

Both directions on the same mineral. FeS forms where the iron and
sulfide fronts meet and seals its voxels; acid arriving from the left
eats it back out again and those voxels reopen.

## What to expect

- porosity **falling, then rising** as the acid front catches up with
  the seal
- the porosity settling rather than oscillating. Chattering means
  `reopen_fraction` is too close to 1
- **mass conserved across both directions**: iron and sulfur are only
  moved between the water and the mineral, never created

## The switches this case exercises

`<precipitation>` and `<dissolution>` both enabled, on one mineral,
with `<is_precipitate>true</is_precipitate>` linking them.

## Notes

This is the case that justifies `reopen_fraction`. Set it to 0.99 and
rerun: you should see the porosity flip back and forth every update
interval, with a flow solve wasted on each flip. That is the failure
the hysteresis exists to prevent.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 15_precip_and_dissolution run/mycase
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
