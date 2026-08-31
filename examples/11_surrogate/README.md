# Example 11: surrogate model, a trained network instead of FBA

The shipped Geobacter network standing in for a flux balance solve.
Two inputs, acetate and Fe(III) uptake; one output, growth rate.

## What to expect

- growth rates in the low 1e-2 1/h range: the shipped network's whole
  output range is 0 to 0.0527 1/h
- a step time far below examples 09 and 10 on the same domain, which
  is the entire reason to use a surrogate
- both substrates drawn down together

## The switches this case exercises

`<enable_surrogate>true`, `<reaction_type>surrogate</reaction_type>`,
and `surrogateModel.hh` compiled in.

## Notes

**The Vmax values matter here in a way they do not elsewhere.** The
shipped network has never seen a Fe(III) uptake above 0.5 mmol/gDW/h.
`<fba_maximum_uptake_flux>` is set to 0.4 to stay inside that. Raise
it and the run leaves the training box, where a network extrapolates
badly and returns confident nonsense.

**There are now two ways to get a network in, and this case shows the
older one.** Here the weights are compiled into `surrogateModel.hh`, so
changing them means editing that file and rebuilding.

The alternative is `<surrogate><weights_file>`, which reads the weights
at run time and can also FIT them during start-up from a metabolic
model. **Example 16 shows that.** Both routes are supported, and a run
with no `<surrogate>` block behaves exactly as this one does.

To fit your own network offline, see `surrogate_training/` in the parent
folder. `inspectSurrogate.py` prints any network's valid range,
recovered by inverting the mapminmax scaling.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 11_surrogate run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space: a staggered slot, 24 × 24 × 6, with the biomass patches seeded. Reports porosity and refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Sweeps the linear program, fits the network, verifies the exported header reproduces the trainer, installs it, and plots the response surface. The longest preparation here, and the only one ending in a rebuild. |
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
| `training/README_SURROGATE_TRAINING.md` | The full training procedure, in prose. |
| `training/example_ecoli_core_training.csv` | A sweep already done, so the fit can be run without waiting for the linear program. |
| `training/exportSurrogateHeader.m` | The MATLAB exporter. |
| `training/extractMM.py` | Exports an SBML or BiGG model to the tabular form the GLPK path reads, and answers the first question worth asking: can this model grow at all when nothing constrains it? |
| `training/generateTrainingData.py` | Sweeps the linear program over a grid of uptake bounds. This is the slow step, and the reason the surrogate exists. |
| `training/inspectSurrogate.py` | Plots and evaluates the response surface. **Run this before any production run** — this path enforces no training range, so a fitted box that returns zero growth over half its area looks exactly like one that does not. |
| `training/testExporterParity.m` | Checks the two exporters agree. |
| `training/trainSurrogate.m` | The MATLAB trainer, for parity with the Python one. |
| `training/trainSurrogate.py` | Fits the network to that sweep and writes the C++ header. |
| `training/verifyExport.py` | Checks the exported header reproduces the Python trainer to machine precision. Run it every time; a transcription error here is invisible in the weights. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

> **Why the training code is copied here rather than shared.** So that this
> folder *is* the procedure. The cost is real and worth stating: a fix to a
> trainer has to be applied to every case that carries it, and
> `tests/check_repo.sh` fails if a copy drifts from `tools/`.
