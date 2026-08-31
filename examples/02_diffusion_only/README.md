# Example 02: diffusion only, no flow

The example 01 tracer with `<Peclet>0</Peclet>`, so the flow solver
never runs and transport is pure diffusion.

## What to expect

- a **straight line** from 1 at the inlet to 0 at the outlet once it
  converges. That is the analytic steady state for diffusion with no
  reaction, and it is the check worth making
- convergence reported before `ade_max_iT`
- no z dependence

## The switches this case exercises

`<Peclet>0</Peclet>`, which switches the flow solver off entirely.

## Notes

Curvature in the profile means a reaction is firing somewhere it
should not, or a boundary is not being held. Both are worth chasing
here, where there is nothing else going on to hide them.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
