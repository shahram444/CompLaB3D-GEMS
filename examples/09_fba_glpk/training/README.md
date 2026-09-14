# Preparing a metabolic model for the compiled linear program

**In:** a genome-scale model, SBML `.xml`, COBRA `.mat`, or `.json`.
**Out:** a flat XML file the solver reads at start-up, plus a table of exchange
reactions and their indices.
**Run:** `MODEL=yourmodel.xml OBJECTIVE=<biomass reaction id> ./offline.sh`
**Needed for:** this case only if you bring your own model. It ships
`input/toy_model.xml` and runs without ever coming here.

---

## What this step is and is not

It is a **format conversion**, not a fit. Nothing is trained, nothing is
estimated, and running it twice on the same model gives the same file. What it
does is take a model written for COBRApy or the MATLAB toolbox and flatten it
into a form a C++ solver can read at start-up without linking a parser for
three different formats.

```
    yourmodel.xml  (SBML, from BiGG / ModelSEED / KBase)
             |
             |   extractMM.py
             v
    input/converted_model.xml     the stoichiometric matrix, the bounds,
                                  the objective, the reaction names
             |
             |   <model_filename> in CompLaB.xml
             v
    GLPK solves it once per voxel per step, in process
```

---

## The shipped toy model

```bash
cat input/toy_model.xml
```

Four reactions, one column of `<S>` each, small enough to read in full. That is
the point of it: before you put a model with two thousand reactions through
this, look at one where you can check every number by eye.

---

## Converting your own

```bash
MODEL=yourmodel.xml OBJECTIVE=Biomass_Ecoli_core ./offline.sh
```

A gzipped model is expanded for you. The objective is the reaction the linear
program maximises, almost always the biomass reaction; if you leave it out the
script uses whatever objective the model already carries, which is usually but
not always what you want.

Then two edits in `CompLaB.xml`:

```xml
<model_filename>input/converted_model</model_filename>   <!-- no .xml -->
<exchange_reaction_names>EX_ac_e EX_o2_e</exchange_reaction_names>
```

**Use the names, not the indices.** `<exchange_reaction_indices>` also works and
the script prints the table for it, but an index that has silently moved because
someone revised the model is a whole afternoon of wondering why the organism is
eating the wrong thing. A name that has moved stops the run.

---

## Checking it worked

The converted file is text. Three things are worth a look before you trust it:

```bash
grep -c '<reaction' input/converted_model.xml    # as many as the model had
grep '<objective' input/converted_model.xml      # the reaction you asked for
grep '<nrxn>' input/converted_model.xml
```

Then run the case and read the start-up log. The solver prints the exchange
reactions it bound and the bounds it will hold them at. If an exchange you named
is missing, it says so and stops rather than guessing.

---

## The two that catch people

**The uptake bounds are set by the local concentration, not by the model.** The
model's own exchange bounds are overwritten every step by what the solver finds
in that voxel. Bounds in the file only matter for the reactions you did not name
as exchanges.

**An infeasible linear program is not an error in the model.** A voxel with no
donor left has nothing to eat, and the program says so. The solver treats that
as zero growth and moves on. A run where every voxel is infeasible from the
first step usually means the exchange reactions were bound the wrong way round,
so uptake is being forced rather than allowed.

---

## Related

- Example 10 does the same thing through COBRApy instead of the compiled GLPK
  path, and is the one to use if your model needs anything COBRApy can do that a
  flat matrix cannot.
- Example 11 sweeps this linear program once, offline, and fits a network to the
  answers, so the program does not have to be solved again inside the run.
- Example 12 runs this path alongside compiled kinetics in the same simulation.
