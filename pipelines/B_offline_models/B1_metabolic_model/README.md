# B1 — Preparing the metabolic model

**In:** a genome-scale metabolic model, SBML or BiGG.
**Out:** for GLPK, a tabular export; for COBRApy, the model file unchanged.
**Run:** `./run.sh`, which calls `tools/extractMM.py`.
**Needed for:** flux balance analysis on either back end, and for B2.

## Where models come from

[BiGG Models](http://bigg.ucsd.edu) is the usual source. Three are bundled under
`models/` for the examples: *E. coli* core, iJO1366, and STM v1.0. Pick a model
for the organism you are actually simulating; a model for the wrong organism
will still solve and still give you a number.

## The two back ends want different things

**COBRApy** reads the SBML or BiGG file directly. Nothing to do here: point
`<model_file>` at it and move on.

**GLPK** needs the tabular export, because the C++ side does not parse SBML:

```bash
python3 ../../../tools/extractMM.py \
       ../../../models/e_coli_core.xml \
       --objective Biomass_Ecoli_core \
       --out ecoli_core.tab
```

The export holds the stoichiometric matrix, the objective vector and the fixed
internal bounds. It does **not** hold the exchange bounds, because those are the
part the simulation rewrites at every voxel.

## Check it before you use it

```bash
grep -c '<reaction' ecoli_core_mm.xml     # what came out
```

reports the reaction and metabolite counts, names the objective, and solves the
model once at a generous set of bounds. **If that solve returns zero growth,
stop here.** A model that cannot grow unconstrained will never grow in a voxel,
and you will spend a long time blaming the transport.

## Which exchanges will the simulation control?

The exchange reactions named in `<Vmax>` and `<Kc>` in `CompLaB.xml` must exist
in the model, in the order the species are declared. Getting this wrong is the
most common failure on this path, so the loader checks the names at start-up and
refuses to run rather than silently bounding the wrong reaction.
