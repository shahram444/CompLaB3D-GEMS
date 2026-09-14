# Preparing a metabolic model for the COBRApy back end

**In:** a genome-scale model, SBML `.xml`, COBRA `.mat`, or `.json`.
**Out:** the same flat XML the compiled path reads, and a working Python
interpreter the solver can call into.
**Run:** `./offline.sh` to check the environment; `MODEL=yourmodel.xml
OBJECTIVE=<biomass reaction id> python3 training/extractMM.py ...` to convert.
**Needed for:** this case only if you bring your own model. It ships
`input/toy_model.xml`.

---

## How this differs from example 09

The chemistry is identical. What changes is **who solves the linear program**.

| | example 09 | this case |
|---|---|---|
| solver | GLPK, linked into the executable | COBRApy, through an embedded Python interpreter |
| per-voxel cost | lower | higher, a Python call each time |
| what you can do | a stoichiometric matrix and bounds | anything COBRApy can do: pFBA, loopless, gene knockouts, a custom objective |
| what can break | the model file | the Python environment |

Use example 09 unless you need something in the second column. Use this one when
you do, and accept that it is slower.

---

## What offline.sh checks

It does not convert anything. It answers the two questions that account for
nearly every failure on this path:

```bash
./offline.sh
```

**Is COBRApy importable, and from which interpreter?** It prints the version and
the path. The path matters: the solver embeds an interpreter, and if that is not
the same environment your shell uses, `import cobra` succeeds at the prompt and
fails inside the run.

**Is the model file structurally sound?** It parses `input/toy_model.xml` and
checks that the stoichiometric matrix really has metabolites times reactions
entries. A matrix with the wrong number of entries is read silently by anything
that does not count, and produces fluxes that look plausible and are wrong.

---

## Converting your own model

Same tool as example 09:

```bash
python3 training/extractMM.py yourmodel.xml -o input/converted_model.xml \
        --objective Biomass_Ecoli_core -f
```

Then in `CompLaB.xml`:

```xml
<model_filename>input/converted_model</model_filename>
<exchange_reaction_names>EX_ac_e EX_o2_e</exchange_reaction_names>
```

Names rather than indices, for the reason given in example 09's training README:
an index that moved because the model was revised will not stop the run, it will
just feed the organism the wrong thing.

---

## The environment problem, and how to see it coming

The solver runs an **embedded** interpreter. That interpreter does not
necessarily inherit what your shell sees.

On a cluster with EasyBuild modules this bites almost every time, because
`EBPYTHONPREFIXES` is honoured by the module system and ignored by an embedded
interpreter. The symptom is `ModuleNotFoundError` for something COBRApy depends
on, often `mpmath`, thrown from inside a run that started fine. The fix is to
set `PYTHONPATH` explicitly in the job script rather than relying on the module:

```bash
export PYTHONPATH=/path/to/site-packages:$PYTHONPATH
```

`./offline.sh` prints `sys.executable`, so compare that with the interpreter the
solver was built against before blaming the model.

---

## Checking it worked

Run the case and read the start-up log. It names the interpreter, the COBRApy
version and the exchange reactions it bound. Then look at the closing report:
the number of linear programs solved and how many came back infeasible. A large
infeasible count in the first steps usually means the exchange bounds were
applied with the wrong sign, so uptake is being forced rather than permitted.

---

## Related

- Example 09, the same chemistry through the compiled GLPK path.
- Example 11, which sweeps the program offline and fits a network to it, so the
  program does not have to be solved inside the run at all.
- Example 21, which calls the program more than once per step.
