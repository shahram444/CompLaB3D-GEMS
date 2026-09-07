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

```bash
python ../../../tools/fit_symbolic.py \
       --data samples.csv --target growth \
       --vars acetate \
       --pop 600 --gens 60 --depth 6 --seed 1 \
       --out discovered.sym
```

`--pop` and `--gens` are the two that matter. If the Pareto set comes back with
no entry whose error drops sharply, the search was too small rather than the
data unfittable: a population of 400 over 60 generations stalls at 13 per cent
on a dual-substrate law that a population of 600 recovers exactly by generation
50.

## What comes back is a list, not an answer

One expression per node count. **Choosing from it is your job**, and it is the
step that makes this method different from curve fitting:

```
nodes  typical error  expression
    1        78.50 %  0.188686
    3        37.80 %  acetate * 276.846
    5         7.59 %  acetate / (acetate + 0.000916122)
    7         0.00 %  (0.00051 / (-0.0006 - acetate)) - -0.85
    9         0.00 %  ((0.85 * acetate) - 8.0e-15) / (acetate + 0.0006)
```

Take the elbow — here seven nodes, where the error collapses and nothing beyond
improves on it. That seven-node expression is Monod's law written in an
unfamiliar arrangement, which you can verify by putting both over a common
denominator.

## Finish the file by hand

The search fits **one** output column. The other species have to be written as
multiples of it, which is what holds them in exact stoichiometric ratio:

```
units  per_hour
vars   acetate

rate   growth  = ((0.00051 / (-0.0006 - acetate)) - -0.85)
rate   Fe3     = -4.8   * growth
rate   Fe2     =  4.8   * growth
rate   HCO3    =  1.845 * growth

range  acetate 2.59305e-05 0.00389437
```

**Write the `range` line.** It is enforced at every evaluation, not advisory,
and without it a voxel outside the fitted data gets an extrapolated answer.
