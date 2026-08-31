# Example 09: flux balance analysis through GLPK

One organism on a four reaction toy model, small enough that the
answer is a formula:

```
growth = min( Vs, 2 * Vo ),   V = Vmax * C / (Kc + C)
```

With Vmax 10 and 4 and both concentrations well above their half
saturation constants, the acceptor limits and growth should settle
near **8 mmol/gDW/h**.

## What to expect

- the model loading as 3 metabolites, 4 reactions
- the metabolic summary showing `GLPK` with 1 microbe
- the growth rate approaching 8, from the stoichiometry rather
  than from this code
- product P accumulating at the same rate biomass is made, since
  the biomass reaction makes exactly one P per unit growth
- no `[NEG!]` warnings

## The switches this case exercises

`<enable_fba_glpk>true`, `<reaction_type>glpk</reaction_type>`, and the
`<model_filename>` / `<exchange_reaction_indices>` /
`<fba_maximum_uptake_flux>` group.

## Notes

Examples 09 and 10 must agree. That comparison is the strongest
check available on the FBA coupling, because the two solvers share
nothing but the linear program.

Drop `<fba_maximum_uptake_flux>` to `4. 10. 0.` and the donor
becomes limiting instead, so growth should settle near 4. That
is a second independent check for one line of editing.

**A safer way to write the exchange mapping.** `<exchange_reaction_indices>` is
positional: the numbers are correct only for this exact model file.
Insert one reaction into it and every index after that points somewhere
else, with nothing to complain.

`<exchange_reaction_names>` names the reactions instead, and they are
resolved against the model at start-up, so a wrong one stops the run and
prints the near matches. It needs an SBML model; the matrix format that
`toy_model.xml` uses does not carry reaction names, which is why this
case still uses indices. **Example 16 shows the named form**, on a real
genome-scale model.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 09_fba_glpk run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space: a staggered slot, 24 × 24 × 6, with the biomass patches seeded. Reports porosity and refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Exports the metabolic model to the tabular form GLPK reads, then checks it can grow unconstrained. |
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

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

> **Why the training code is copied here rather than shared.** So that this
> folder *is* the procedure. The cost is real and worth stating: a fix to a
> trainer has to be applied to every case that carries it, and
> `tests/check_repo.sh` fails if a copy drifts from `tools/`.
