# Example 16: a rate law read from a file

Example 11's system with the rate law moved out of C++ and into
input/rates.sym, read at start-up.

## What to expect

- growth of the same order as example 11: same organism, same two
  substrates
- **Fe(III) drawn down 4.8 times faster than acetate, in molar terms.**
  That ratio is the reaction stoichiometry, and nobody tells the solver
  it at run time -- it falls out of the three lines in rates.sym. If it
  comes back at 4.8, the file parsed and every name bound correctly
- the whole rate law echoed into the log before iteration 0

## The switches this case exercises

`<symbolic>` with `<expressions_file>`, and
`<reaction_type>symbolic</reaction_type>` on the organism.

## Notes

Change 0.0527 in rates.sym and rerun without rebuilding: that is the
capability this case exists to show. Push a concentration outside a
`range` line and the end-of-run report tells you what fraction of
evaluations were clamped, and warns above 5%.

`fitting/` holds the worked material behind this: the sweeps, the fitted
files, and a checker that holds one of them against the law it came from.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
