# Checking the energetics before the solver reads them

**In:** `input/aom.thm`, a file of measured energetics.
**Out:** `output/thermo_curve.csv` and a table on the terminal.
**Run:** `./offline.sh` from the case directory.
**Needed for:** nothing. The case runs without it. Skip it once and you will
run it every time afterwards.

---

## Why this step exists

Every other rate path has an offline stage: a metabolic model to fetch, a
network to train, an expression to fit. This path has neither a fit nor a
training set. It has a file of numbers taken from measurements, and the question
that corresponds to training is a different one:

**does the threshold those numbers describe sit anywhere near the
concentrations this run will actually visit?**

That question is trivial to answer beforehand and impossible to answer
afterwards, which is the whole argument for this step.

A gate that is **shut everywhere** looks in the output exactly like a reaction
that never happened. A gate that is **open everywhere** looks exactly like no
gate at all. Neither is distinguishable from the fields after the fact. Both are
obvious from a single table printed before the run starts.

---

## Running it

```bash
./offline.sh
```

which runs `offline/thermo_curve.py`. Nothing beyond the standard library is
needed, so this works on a login node with no modules loaded.

It reads the energetics, walks the composition across the range the boundary
conditions will feed in, and prints the factor at each point together with the
concentration at which it crosses.

---

## Reading the table

Three outcomes, and only one of them means carry on.

**The factor spans the range.** High where the inlet composition is fresh,
falling as the products accumulate, reaching zero somewhere inside the domain.
This is the case worth simulating: the threshold is doing something and the
simulation will show you where.

**The factor is 1 everywhere.** The reaction is thermodynamically comfortable
across the whole range you are feeding. The run will be identical to one with no
gate at all, so either widen the composition range or accept that this
particular question does not need the gate.

**The factor is 0 everywhere.** Nothing will happen. Check the sign of the
standard free energy first: it is the usual cause, and a sign error there is
silent in every other respect.

The CSV is written so you can plot it rather than squint at the table:

```bash
python3 -c "
import csv, sys
rows = list(csv.DictReader(open('output/thermo_curve.csv')))
print(rows[0].keys())
print(len(rows), 'rows')"
```

---

## The other script in this directory

`offline/upscale.py` is a different thing and is not run by `offline.sh`. It
sweeps aggregate size and bulk composition, runs the solver at each point, and
turns the results into an upscaling law: one effectiveness factor per condition,
the ratio of the mean rate inside an aggregate to the rate the bulk composition
would give on its own.

That is a campaign, not a preparation step. It takes many runs rather than
seconds, and it belongs with example 20 rather than here. It sits in this
directory because the gate is what makes the effectiveness factor interesting:
without it the aggregate interior is merely slower, and with it the interior can
stop entirely, which is a different shape of answer.

---

## What to check after the run

The closing report says how many evaluations the gate was applied to and how the
factor was distributed. Compare that spread against the table this step printed.
They should agree. If the run reports a gate that is shut everywhere when this
step predicted a spread, the composition inside the domain went somewhere the
boundary conditions did not suggest, which is usually a product accumulating
because nothing is carrying it away.

---

## Related

- Example 20, which uses the effectiveness factor `upscale.py` produces.
- The thermodynamics section of this case's `README.md`, for what the file's
  fields mean and where the numbers come from.
