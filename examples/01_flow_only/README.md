# Example 01: flow and a conservative tracer

Navier Stokes plus advection diffusion, no chemistry and no biology.
The tracer is held at 1 on the inlet face and 0 on the outlet, and
reacts with nothing.

## What to expect

- the flow solver converging, then the tracer sweeping left to right
- a monotone profile: concentration falls with x, never rises
- no `[NEG!]` warnings
- every z plane identical, because z is periodic and the geometry is
  extruded

## The switches this case exercises

`biotic_mode` false, `<Peclet>` non zero, Dirichlet solute boundaries.

## Notes

This is the smoke test. If it fails, the build or the geometry is
wrong and nothing further is worth debugging.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 01_flow_only run/mycase
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

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

> **Why the training code is copied here rather than shared.** So that this
> folder *is* the procedure. The cost is real and worth stating: a fix to a
> trainer has to be applied to every case that carries it, and
> `tests/check_repo.sh` fails if a copy drifts from `tools/`.
