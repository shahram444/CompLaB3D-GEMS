# Example 16: the whole pipeline, from the XML alone

The four things that used to be manual, all in one file. The metabolic
model comes out of the bundle; the exchange reactions are named rather
than numbered; the surrogate network is fitted during start-up from
that model; and the run writes its own summary CSV.

There is no preparation step. No extractMM.py, no training scripts, no
post-processing to find out whether anything was conserved.

## What to expect

- `[MODEL] unpacking bundled models/e_coli_core.xml.gz`, then
  `revision check: 72 species, 95 reaction tags, as the manifest expects`
- `resolved 3 exchange reaction name(s) against e_coli_core`
- `[SRG] no weights file; training one now. This happens once.`
  followed by a few hundred linear programs, then a held-out R^2 above
  0.99
- **`output/ecoli.srg` written.** Run the case again: it loads instead
  of training, and start-up is instant
- `output/summary.csv`, one row per 100 steps
- at the end, `[SRG] runtime network: N evaluation(s), 0 clamped`. A
  non-zero clamp count means the simulation left the box the network
  was fitted on

## The switches this case exercises

`<model_source>`, `<exchange_reaction_names>`, `<surrogate>` with
`<train_if_missing>`, and `<diagnostics>`. Needs `-DENABLE_GLPK=ON`,
because training sweeps a linear program.

## Notes

**This run has no GLPK microbe**, so there is no persistent linear
program to sweep. A temporary one is built from this microbe's own
model for training and released as soon as the fit is done. The log
says so.

**Why the geometry is still a file here.** `<generate>` in `<domain>`
can build the pore space from the XML, and example 01 shows the tags
commented in place. It writes pore, solid and wall only -- it cannot
seed the microbe material numbers this case needs, so the seeded
geometry is shipped.

**Comparing against the FBA it replaces.** Set `<enable_fba_glpk>true`
and `<reaction_type>glpk</reaction_type>`, delete the `<surrogate>`
block, and rerun. On the reference run the two agreed to four parts in
ten thousand in final biomass, and the surrogate was about 150 times
faster.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 16_complete_pipeline run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space: a staggered slot, 24 × 24 × 6, with the biomass patches seeded. Reports porosity and refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Pre-fetches the metabolic model, so the run uses a file you hold rather than whatever a remote database serves today. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's balance check. Standard library only. |
| | `pipeline.sh` | Runs the four in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |

### The metabolic model

| File | What it is |
|---|---|
| `models/NOTICE.md` | Where these models came from and under what terms. |
| `models/e_coli_core.xml.gz` | The *E. coli* core genome-scale model — 95 reactions, the standard small test model. |
| `models/manifest.txt` | Checksums, so a model swapped underneath you is detected rather than silently used. |

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
