# Discovering the rate law

**In:** a table of concentrations and the rate they produced.
**Out:** a `.sym` text file the solver reads at start-up, complete and runnable.
**Run:** `python3 training/make_rate_law.py`, or `./offline.sh` from the case
directory.
**Needed for:** this case only if you want your own law. It ships
`input/growth.sym` and runs without ever coming here.

---

## What is in this directory

| | |
|---|---|
| `make_training_data.py` | writes the table, from a law you can read |
| `growth_samples.csv` | 400 rows: acetate, o2, and the growth rate they gave |
| `growth_samples_truth.json` | the law behind those rows, and the noise and seed |
| `make_rate_law.py` | asks questions about your table and runs the search |
| `fit_symbolic.py` | the search itself, if you would rather pass flags |

---

## The loop this case closes

```
    make_training_data.py          a law we chose, sampled, with noise added
              |
              v
    growth_samples.csv    +    growth_samples_truth.json    the answer key
              |
              |   the search is shown the table and not the key
              v
    a list of formulas, one per length
              |
              |   you choose one
              v
    input/growth_discovered.sym    and the key says whether it is right
```

That is the whole point. A search fitted to real data can only be judged on how
well it fits. A search fitted to data drawn from a law we wrote down can be
judged on whether it found **that law**, which is a much harder test and the
only one that means anything when the method is what is being demonstrated.

---

## The sampling range is the thing to get right

A half-saturation constant can only be recovered from data that reaches it.

Below its constant, a Monod term is a straight line in that variable. Every
value of the constant fits a straight line equally well. So a table that stops
short of the constant cannot identify it, however many rows it has and however
carefully the search is run. The search will return a product of the two
concentrations, that product will fit well, and nothing in the fit will tell you
anything is missing.

This case was originally guilty of exactly that. It sampled acetate to 0.0049
mol/L against a half-saturation constant of 0.05, a tenth of the way there, and
a plain product already fitted to within 4 per cent.

**The concentrations were not the problem.** 1e-6 to 5e-3 mol/L is what this
case runs at and what a pore actually holds. The constants were: 0.05 mol/L for
acetate and 0.01 for oxygen are one to two orders of magnitude above any
measured aerobic heterotroph. They are now

```
    growth = 0.35 * acetate / (5.0e-4 + acetate) * o2 / (1.0e-4 + o2)     per hour
```

which is defensible on its own terms and puts both constants inside the sampled
range. The same table now runs from well below each constant to ten and twenty
times above it, and the curve is covered end to end.

`make_training_data.py` prints that ratio every time it runs, and refuses
quietly with a warning if the sampling stops short of three times either
constant.

---

## Making the table

```bash
python3 make_training_data.py                       # the shipped table
python3 make_training_data.py --n 2000 --noise 0.05 # more rows, noisier
python3 make_training_data.py --seed 7 -o other.csv
```

It writes the CSV and a companion `_truth.json` holding the law, the constants,
the noise level, the seed and the sampled range, so the key cannot drift away
from the data it belongs to.

---

## Running the search

Two ways. Answer questions:

```bash
python3 make_rate_law.py
```

It reads any table, lists its columns back at you, asks which one is the rate
and which are the inputs, asks for the reaction and the yield, prints a summary
before it starts, and prints the equivalent command at the end so the run can be
repeated without answering anything.

Or pass the flags yourself, which is what `offline.sh` does:

```bash
python3 fit_symbolic.py --data growth_samples.csv --target growth \
        --inputs acetate,o2 \
        --reaction "acetate -1  o2 -2" --yield "acetate 0.4" --biomass Bug \
        --units per_hour --pop 600 --gens 60 --seed 1 --verify \
        --out ../input/growth_discovered.sym
```

---

## What the search is told, and what it is not

It is told four operators, `+ - * /`, and your column names. Nothing else. It
does not know what saturation is, it cannot reason that a saturating process
needs a division, and it never learns. It builds formulas at random, checks them
against your table, throws away what fits badly and breeds what fits well.

Division appears in the answer because formulas with division match the data,
not because anything understood why.

---

## What comes back, and how to choose

A list, one formula per length, where length is how many pieces a formula has
and a number, a name or an operator is one piece. A longer formula always fits
better, because more numbers can bend to more of the data, which is exactly why
you should not take the best-fitting one.

Read down and watch the error. Early on, each extra piece cuts it a lot: that
piece is describing something real. Later it barely moves: those pieces are
describing the noise in your table. Take the last length where the error still
fell by something worth having, with `--pick <length>`.

The automatic choice is a rule of thumb and says so. On this data it lands short.

---

## Checking it repeats

`--verify` runs the whole search a second time and reports which lengths agreed.
It compares what the formulas **compute**, not how they are spelled, because the
same law comes back written two ways more often than not, and it judges the gap
between runs against each formula's own distance from your data rather than a
fixed tolerance.

Expect the short and middle lengths to agree and the longest not to. That is not
a defect. A formula carrying more constants than 400 rows can pin down has no
single best answer, so two runs land in different places, and those are the same
formulas nobody should be quoting.

---

## What the file it writes contains

```
units    per_hour
vars     acetate o2 Bug

reaction acetate -1 o2 -2
yield    acetate 0.4
biomass  Bug

rate     growth = ...

range    acetate ...
range    o2      ...
```

plus a header carrying the command, the seed, the thread setting and the library
versions. It runs as it stands. The substrate lines are **not** in it and do not
need to be: the solver derives one per species from the reaction at start-up, so
their ratios are arithmetic rather than typing and cannot disagree with the
chemistry.

---

## Installing it

The file is written to `input/growth_discovered.sym` and is **not** installed.
`CompLaB.xml` still points at `input/growth.sym`. Compare them first, then:

```xml
<expressions_file>input/growth_discovered.sym</expressions_file>
```

---

## The one that catches people

**The `range` lines are enforced, not advisory.** Every evaluation clamps each
input to the range the table covered, and the closing report says how often that
happened. A high clamp percentage means the simulation went somewhere your table
never did, and the law is being held at the edge of what it was fitted on rather
than trusted beyond it. That is a signal to widen the table, not to ignore.

---

## Related

- Example 18 learns the same kind of rate without writing a formula down.
- `pipelines/B_offline_models/B3_symbolic_law/` is the same search as a
  standalone stage you can point at your own data.
