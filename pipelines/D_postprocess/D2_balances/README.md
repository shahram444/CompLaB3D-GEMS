# D2 — Checking the run is right

**In:** the VTI output and the log.
**Out:** a report you should read before believing anything else.
**Run:** `./run.sh`, or `tools/postprocess/postprocess.py`.

```bash
# the report, with no conservation check asked for
python ../../../tools/postprocess/postprocess.py run/mycase

# name a sum that must stay constant, and it is checked
python ../../../tools/postprocess/postprocess.py run/mycase --conserve "Fe+FeS" --tol 1e-9
```

`--conserve` takes an expression in the `<name_of_substrates>` names and may be
repeated. Without it the report still gives every field's change, but nothing is
asserted to be conserved, because which totals are closed depends on the
boundary conditions and only you know that.

## What it checks, and what to do when one fails

**Mass balance.** Every mole of mineral formed must remove the matching moles of
dissolved reactant, and every mole dissolved must appear in the water. The
report gives the residual as a fraction of the total. Anything above about
`1e-10` is worth understanding; a residual that grows with time is a real
problem.

*If dissolution is on and you ran in parallel, expect a residual and see the
limitation in [`../../C_run/README.md`](../../C_run/README.md) before looking
for a bug in your chemistry.*

**Clamped evaluations.** On the symbolic and graph-network paths, the report
gives how many rate evaluations were clamped to the fitted range. A large
fraction means the simulation is spending its time outside the data the law was
fitted from, and the answer is to refit over a wider range, not to remove the
clamp.

**Porosity history.** Should fall monotonically for precipitation, rise for
dissolution, and do neither for a run with both disabled. A porosity that moves
in a run with no geometry process enabled means something is writing into the
mask.

**Clogging.** For a precipitation run, whether the outflow reached zero, and at
what time. `Q_out = 0` with open path remaining means the flow solver did not
converge, not that the pore sealed.

## Comparing two runs

```bash
python ../../../tools/postprocess/postprocess.py run/glpk      --output cmp/glpk
python ../../../tools/postprocess/postprocess.py run/surrogate --output cmp/surrogate
diff cmp/glpk/summary.csv cmp/surrogate/summary.csv
```

There is no built-in comparison mode. Run the tool on each folder and compare
the two summary tables, which is the same thing and leaves you looking at the
numbers rather than at a verdict.

This is how the surrogate is validated against the linear program it replaces:
same case, same geometry, same seed, two rate paths, and the difference in the
biomass field is the error the approximation costs you.
