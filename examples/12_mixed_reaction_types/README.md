# Example 12: combined reaction types, FBA plus hand written kinetics

The example 09 organism with `<reaction_type>glpk_and_kinetics`, so
the linear program and `defineKinetics.hh` both run and their rates
add. The kinetics file contributes maintenance: a donor drain that
yields no growth.

## What to expect

- the same growth rate as example 09, near 8 mmol/gDW/h, because
  maintenance does not add biomass
- the donor drawn down **faster** than in example 09, by exactly the
  maintenance term
- the metabolic summary showing 1 GLPK microbe, and the kinetics
  path also active

## The switches this case exercises

`<reaction_type>glpk_and_kinetics</reaction_type>`, both
`<enable_kinetics>` and `<enable_fba_glpk>` true, and both
`<maximum_uptake_flux>` and `<fba_maximum_uptake_flux>` present.

## Notes

Diff this case's donor profile against example 09's. The difference
IS the maintenance term, and it is the cheapest way to confirm the
two paths are genuinely adding rather than one silently winning.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 12_mixed_reaction_types run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space: a staggered slot, 24 × 24 × 6, with the biomass patches seeded. Reports porosity and refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Exports the model for the FBA organism. The population on compiled kinetics needs nothing prepared. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's balance check. Standard library only. |
| | `pipeline.sh` | Runs the four in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `input/toy_model.xml` | The metabolic model, in the directory the solver reads it from. |

### The offline code this case ships

| File | What it is |
|---|---|
| `training/extractMM.py` | Exports an SBML or BiGG model to the tabular form the GLPK path reads, and answers the first question worth asking: can this model grow at all when nothing constrains it? |
| `training/makeKinetics.py` | Generates a `defineKinetics.hh` from a reaction list, instead of writing the C++ by hand. |

### This case's own chemistry

| File | What it is |
|---|---|
| `kinetics/defineKinetics.hh` | This case's own chemistry, compiled in on top of the shared default. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

> **Why the training code is copied here rather than shared.** So that this
> folder *is* the procedure. The cost is real and worth stating: a fix to a
> trainer has to be applied to every case that carries it, and
> `tests/check_repo.sh` fails if a copy drifts from `tools/`.
