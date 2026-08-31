# Example 02: diffusion only, no flow

The example 01 tracer with `<Peclet>0</Peclet>`, so the flow solver
never runs and transport is pure diffusion.

## What to expect

- a **straight line** from 1 at the inlet to 0 at the outlet once it
  converges. That is the analytic steady state for diffusion with no
  reaction, and it is the check worth making
- convergence reported before `ade_max_iT`
- no z dependence

## The switches this case exercises

`<Peclet>0</Peclet>`, which switches the flow solver off entirely.

## Notes

Curvature in the profile means a reaction is firing somewhere it
should not, or a boundary is not being held. Both are worth chasing
here, where there is nothing else going on to hide them.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 02_diffusion_only run/mycase
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
