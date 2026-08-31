# Example 01: flow and a conservative tracer

Navier Stokes plus advection diffusion, no chemistry and no biology.
The tracer is held at 1 on the inlet face and 0 on the outlet, and
reacts with nothing.

## What to expect

- the flow solver converging, then the tracer sweeping left to right
- a monotone profile: concentration falls with x, never rises
- no `[NEG!]` warnings
- every z plane identical, because z is periodic and the geometry is
  extruded

## The switches this case exercises

`biotic_mode` false, `<Peclet>` non zero, Dirichlet solute boundaries.

## Notes

This is the smoke test. If it fails, the build or the geometry is
wrong and nothing further is worth debugging.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
