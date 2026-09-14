# Preparing the complete pipeline

**In:** a genome-scale model.
**Out:** a converted model the solver reads at start-up.
**Run:** `MODEL=yourmodel.xml OBJECTIVE=<biomass reaction id> ./offline.sh`
**Needed for:** this case only if you bring your own model. It ships everything.

---

## What this case is

Every other example isolates one thing. This one runs the whole chain end to
end, so it is the case to copy when you are starting a study of your own rather
than learning one feature.

```
    A   geometry            a pore space, generated or from an image
    A   initial fields      what is dissolved where, at the start
    B   this step           the metabolic model, converted
    C   the run             flow, transport, reaction, biomass, geometry update
    D   postprocess         slices, histories, the mass balance report
```

The offline stage is the shortest part of that chain and the only one that has
to happen before anything else, because the solver reads the converted model at
start-up and will not start without it.

---

## Running it

```bash
MODEL=yourmodel.xml OBJECTIVE=Biomass_Ecoli_core ./offline.sh
```

It calls `training/extractMM.py`, the same converter examples 09, 10 and 12 use.
See example 09's training README for what the tool does, what the output
contains and why to bind exchange reactions by name rather than by index.

`pipeline.sh` runs this step for you as part of the chain, so you only come here
directly when you want to convert a model without running anything else.

---

## What is different here

**Nothing about the conversion.** It is the same tool with the same flags.

What is different is that **a mistake here surfaces four stages later**. In
example 09 a badly bound exchange reaction shows up immediately, because the
case does nothing else. Here it shows up as a biomass field that grows in the
wrong place after a geometry update, and by then there are three other things it
could plausibly have been.

So the discipline for this case is to check the conversion on its own terms
before running the chain:

```bash
grep -c '<reaction'  input/converted_model.xml
grep '<objective'    input/converted_model.xml
grep '<nrxn>'        input/converted_model.xml
```

and then run stage C alone once, with a short step count, and read the start-up
log to confirm the exchange reactions bound to the substrates you meant.

---

## The mass balance report is the check that matters

Stage D writes a balance report, and this is the case where it earns its place.
A conversion that went wrong in a way the start-up log did not catch, an
exchange bound to the wrong species, say, shows up there as a species that does
not close.

```bash
python3 postprocess.py
```

It takes no arguments. It reads `output/summary.csv` and the solver log from the
case folder it is run in, and prints a verdict. `run.sh` runs it for you and tees
the output into `output/checks.log`.

For a conservation check on a sum you name, the shared tool does that instead:

```bash
python3 ../../tools/postprocess/postprocess.py . --conserve "<a sum that should be closed>"
```

A closed balance does not prove the model is right. An unclosed one proves
something is wrong, and it is nearly always either the exchange binding or a
sum named in `--conserve` that is not actually closed in this geometry.

---

## Related

- Example 09 for the converter itself.
- The `pipelines/` directory, which is the same five stages as standalone scripts
  you can run against your own case rather than against this example.
