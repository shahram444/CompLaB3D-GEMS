# Example 03: abiotic kinetics, A + B to C

A second order irreversible reaction with no biology. A enters from
the left, B from the right, and C forms where they overlap.

## What to expect

- a peak in C near the middle of the domain, with A and B drawn down
  around it
- the reaction front stationary, because there is no flow to push it
- **mass conserved**: total A + total C constant, total B + total C
  constant

## The switches this case exercises

`enable_abiotic_kinetics`, and a hand written rate law in
`defineAbioticKinetics.hh`.

## Notes

Turn `<Peclet>` up to 1 and the front moves downstream. That single
change is the cheapest way to see transport and reaction competing.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
