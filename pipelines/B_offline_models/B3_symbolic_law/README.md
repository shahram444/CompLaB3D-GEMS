# B3 — Finding a symbolic rate law

**In:** a table with one column per variable and one target column.
**Out:** a `.sym` text file the solver reads at start-up.
**Run:** `./run.sh`, which calls `tools/fit_symbolic.py`.
**Needed for:** the symbolic path only.

## Where the data comes from

Anywhere. Flux balance output, a laboratory measurement, a field experiment. The
search does not care, and this is the only rate path where that is true — the
surrogate can only be fitted to flux balance output, because it is fitted to
reproduce a linear program.

## Running the search

Two ways. Answer questions:

```bash
python ../../../tools/make_rate_law.py
```

It reads any training table, lists its columns back at you, asks which one is
the rate and which are the inputs, asks for the chemistry, and prints the
equivalent command at the end so the run can be repeated without answering
anything.

Or type that command yourself:

```bash
python ../../../tools/fit_symbolic.py \
       --data samples.csv --target growth --inputs acetate \
       --reaction "acetate -1  Fe3 -4.8  Fe2 +4.8  HCO3 +1.845" \
       --yield "acetate 0.4" --biomass Bug \
       --units per_hour \
       --pop 600 --gens 60 --depth 6 --seed 1 --verify \
       --out discovered.sym
```

`--pop` and `--gens` are the two that matter. If the list comes back with no
entry whose error drops sharply, the search was too small rather than the data
unfittable: a population of 400 over 60 rounds stalls at 13 per cent on a
dual-substrate law that a population of 600 recovers exactly by round 50.

`--verify` runs the whole search a second time and reports which lengths agreed.
It doubles the run time and it is worth it before a number goes in a paper.

## What comes back is a list, not an answer

One formula per length, where length is how many pieces a formula has and a
number, a name or an operator is one piece. **Choosing from it is your job**,
and it is the step that makes this method different from curve fitting:

```
length  typical error  formula
    1        78.50 %   0.188686
    3        37.80 %   acetate * 276.846
    5         7.59 %   acetate / (acetate + 0.000916122)
    7         0.00 %   (0.00051 / (-0.0006 - acetate)) - -0.85
    9         0.00 %   ((0.85 * acetate) - 8.0e-15) / (acetate + 0.0006)
```

Read down and watch the error. Take the last length where it still fell by
something worth having — here 7, where it collapses and nothing beyond improves
on it. That formula is Monod's law written in an unfamiliar arrangement, which
you can verify by putting both over a common denominator. Take a length other
than the automatic choice with `--pick 7`.

## What it writes

A file that runs as it stands:

```
units    per_hour
vars     acetate Bug

reaction acetate -1  Fe3 -4.8  Fe2 +4.8  HCO3 +1.845
yield    acetate 0.4
biomass  Bug

rate     growth = ((0.00051 / (-0.0006 - acetate)) - -0.85)

range    acetate 2.59305e-05 0.00389437
```

plus a header carrying the command, the seed, the thread setting and the library
versions.

The search fits **one** output column. The other species come from the reaction:
before the first step the solver computes, for each species,

```
coefficient  =  ( its number / |number of the yield species| ) / yield
```

and evaluates that coefficient times the fitted rate times the biomass. Their
ratios are then arithmetic on the reaction rather than four numbers typed by
hand that nothing checks. The log prints every derived line marked
`<- from the reaction, not the file`.

**The `range` lines are enforced at every evaluation, not advisory**, and the
search writes them from the span of your table. Without them a voxel outside the
fitted data would get an extrapolated answer; with them it is held at the edge
of what the law was fitted on, and the closing report says how often that
happened.

## For a reaction with no organism

Drop the yield and the biomass, and name the fitted rate `extent`, which is how
fast the reaction itself turns:

```bash
python ../../../tools/fit_symbolic.py \
       --data samples.csv --target rate --inputs Fe,HS \
       --reaction "Fe -1  HS -1  FeS +1" --rate-name extent \
       --units per_second --out abiotic.sym
```

Each species rate is then its number times `extent`, with no biomass factor,
which is what an abiotic law needs. `expected/abiotic_stoich.sym` shows the same
law written both ways, derived and by hand, side by side.

## The older behaviour still works

Leave `--reaction` out and the file holds only the fitted line, and you write
the other species yourself as multiples of it. Every `.sym` written before the
reaction keywords existed still loads untouched. A species may have a `rate`
line **or** a place in the `reaction` line, never both, because two sources of
truth for one number is what this removes.
