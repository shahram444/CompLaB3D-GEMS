# Example 13: precipitation and pore clogging

FeS forms where an iron front and a sulfide front meet, accumulates
in place because it is declared immobile, and seals its voxels once
it reaches the full mineral density.

## What to expect

- FeS building up in a band across the middle
- the porosity reported in the log **falling in steps**, one per
  update interval
- the flow re-solved after each conversion, and the velocity rising
  in the remaining open path
- **mass conserved**: every mole of FeS formed removes one mole of
  Fe2+ and one of HS-

## The switches this case exercises

`<precipitation><enabled>true`, `<immobile>true</immobile>` on the
mineral substrate, and the rate law in `defineAbioticKinetics.hh`.

## Notes

`max_precipRho` is 48.9 mol/L, which is 1000 divided by mackinawite's
molar volume of 20.44 cm3/mol. Derive it from your own mineral rather
than copying this number.

If the domain seals completely the run stops and reports a percolation
limit. That is the code working, not failing.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 13_precipitation run/mycase
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
