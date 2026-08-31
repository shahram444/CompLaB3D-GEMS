# Example 13: precipitation and pore clogging

FeS forms where an iron front and a sulfide front meet, accumulates
in place because it is declared immobile, and seals its voxels once
it reaches the full mineral density.

## What to expect

- FeS building up in a band across the middle
- the porosity reported in the log **falling in steps**, one per
  update interval
- the flow re-solved after each conversion, and the velocity rising
  in the remaining open path
- **mass conserved**: every mole of FeS formed removes one mole of
  Fe2+ and one of HS-

## The switches this case exercises

`<precipitation><enabled>true`, `<immobile>true</immobile>` on the
mineral substrate, and the rate law in `defineAbioticKinetics.hh`.

## Notes

`max_precipRho` is 48.9 mol/L, which is 1000 divided by mackinawite's
molar volume of 20.44 cm3/mol. Derive it from your own mineral rather
than copying this number.

If the domain seals completely the run stops and reports a percolation
limit. That is the code working, not failing.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
