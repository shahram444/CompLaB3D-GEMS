# Example 08: two microbes competing, two different solvers

Two organisms on opposite walls sharing one donor. Bug1 uses the
cellular automaton and grows fast; Bug2 uses lattice Boltzmann and is
more efficient at low substrate.

## What to expect

- both populations growing, then their **ratio going flat**
- Bug1 ahead, because it sits nearer the inlet where the donor is
  plentiful and its higher mu_max pays off
- a donor shadow between the two colonies

## The switches this case exercises

Per microbe `<solver_type>` and `<reaction_type>`, two `<microbeN>`
material numbers, `bfilm_count` of 2.

## Notes

Swap the two `<solver_type>` lines and rerun. The ratio should barely
move: the solver decides how biomass spreads, not how much of it
there is. If it moves a lot, that is worth understanding before you
trust either solver on a real problem.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
