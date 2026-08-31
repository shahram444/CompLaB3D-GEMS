# Example 05: biotic kinetics, CA biomass solver

One organism, Monod growth on one donor, biomass transported by CA.
Examples 05, 06 and 07 differ only in `<solver_type>`.

## What to expect

- biomass rising from the seeded patch, then levelling off as
  growth and decay balance
- the donor drawn down around the colony and replenished from the
  left boundary
- a product plume spreading from the colony
- no `[NEG!]` warnings: the donor must never go negative

## The switches this case exercises

`<solver_type>CA</solver_type>`, the cellular automaton: biomass that exceeds the maximum density
     spills into neighbouring voxels, which is how a biofilm spreads

## Notes

`CA_method` picks how the excess is shared. `fraction` splits it in
proportion to the space available in each neighbour.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
