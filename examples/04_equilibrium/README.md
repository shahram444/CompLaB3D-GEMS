# Example 04: aqueous speciation, the carbonate system

Four species over two components, re-equilibrated in every pore voxel
at every step. Acid enters from the left, so the speciation should
shift along the domain.

## What to expect

- the solver reporting the tableau it parsed: 2 components, 4 species
- H2CO3 dominant where the acid is, CO3 dominant where it is not
- **totals conserved**: total carbonate, summed over HCO3, CO3 and
  H2CO3, changes only through transport, never through speciation
- convergence of the Anderson accelerated iteration reported each step

## The switches this case exercises

`<equilibrium><enabled>true`, plus the `<components>`,
`<stoichiometry>` and `<logK>` tableau.

## Notes

This is the most expensive part of the code by a wide margin, which
is why every other example leaves it off. Even on this 3456 voxel
domain you will see it in the step time.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
