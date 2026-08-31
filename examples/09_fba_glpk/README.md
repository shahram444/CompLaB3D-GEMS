# Example 09: flux balance analysis through GLPK

One organism on a four reaction toy model, small enough that the
answer is a formula:

```
growth = min( Vs, 2 * Vo ),   V = Vmax * C / (Kc + C)
```

With Vmax 10 and 4 and both concentrations well above their half
saturation constants, the acceptor limits and growth should settle
near **8 mmol/gDW/h**.

## What to expect

- the model loading as 3 metabolites, 4 reactions
- the metabolic summary showing `GLPK` with 1 microbe
- the growth rate approaching 8, from the stoichiometry rather
  than from this code
- product P accumulating at the same rate biomass is made, since
  the biomass reaction makes exactly one P per unit growth
- no `[NEG!]` warnings

## The switches this case exercises

`<enable_fba_glpk>true`, `<reaction_type>glpk</reaction_type>`, and the
`<model_filename>` / `<exchange_reaction_indices>` /
`<fba_maximum_uptake_flux>` group.

## Notes

Examples 09 and 10 must agree. That comparison is the strongest
check available on the FBA coupling, because the two solvers share
nothing but the linear program.

Drop `<fba_maximum_uptake_flux>` to `4. 10. 0.` and the donor
becomes limiting instead, so growth should settle near 4. That
is a second independent check for one line of editing.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake -DENABLE_GLPK=ON .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
