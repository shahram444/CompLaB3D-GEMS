# Example 14: dissolution and pore reopening

Acid enters from the left and eats into four calcite grains. As each
grain thins, its outer voxels drop below `reopen_fraction` of full
density and reopen to flow.

## What to expect

- the acid front eating into the **upstream** faces of the grains
  first, because that is where the acid arrives
- Ca2+ leaving on the right
- the porosity in the log **rising in steps**, the mirror of
  example 13
- **mass conserved**: every mole of calcite removed appears as a mole
  of Ca2+ in the water

## The switches this case exercises

`<dissolution><enabled>true`, one `<phase0>` block with a non zero
`<initial_fill>`, and `defineDissolutionRate()` in
`defineAbioticKinetics.hh`.

## Notes

`<initial_fill>27.1</initial_fill>` is what makes the ORIGINAL grains
dissolvable rather than only material that precipitated first. A
precipitate phase would use `<initial_fill>0</initial_fill>` and
`<is_precipitate>true</is_precipitate>` instead.

Nothing precipitates in this case. `defineAbioticRxnKinetics` is
empty on purpose, so the only thing moving mineral is the dissolution
hook and the mass balance is unambiguous.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
