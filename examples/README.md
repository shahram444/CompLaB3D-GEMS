# CompLB3D example suite

Fifteen tiny cases, one per capability. Each runs in seconds on a laptop, each
is self-contained, and each says what it should produce.

They are deliberately small. The point is to prove a feature works and to show
what the switches do, not to be physically interesting. Every domain is
24 × 24 × 6 with z periodic, which makes it quasi-2D: every z-plane sees
identical conditions, so a result that depends on z is a bug rather than
physics.

```bash
python3 makeGeometry.py                       # writes input_shared/
python3 makeExamples.py                       # writes the 15 case folders
python3 validateExamples.py ../CompLB3D       # structural check, one second
./runAllExamples.sh ../CompLB3D               # build and run everything
./runAllExamples.sh ../CompLB3D 09 10         # or just these two
```

---

## The cases

| | case | what it exercises |
|---|---|---|
| 01 | `flow_only` | Navier–Stokes plus advection–diffusion, no chemistry. The smoke test. |
| 02 | `diffusion_only` | `<Peclet>0`, so no flow. **Has an analytic answer.** |
| 03 | `abiotic_kinetics` | `defineAbioticKinetics.hh`: A + B → C |
| 04 | `equilibrium` | the speciation tableau, on the carbonate system |
| 05 | `biotic_cellular_automaton` | Monod growth, CA biomass solver |
| 06 | `biotic_finite_difference` | the same case, FD solver |
| 07 | `biotic_lattice_boltzmann` | the same case, LBM solver |
| 08 | `two_microbes` | two organisms, two *different* solvers, competing |
| 09 | `fba_glpk` | flux balance analysis through GLPK. **Has an analytic answer.** |
| 10 | `fba_cobrapy` | the same model through COBRApy. Must agree with 09. |
| 11 | `surrogate` | the trained network instead of an LP |
| 12 | `mixed_reaction_types` | `glpk_and_kinetics`: both paths, rates added |
| 13 | `precipitation` | mineral growth and pore clogging |
| 14 | `dissolution` | mineral loss and pore reopening |
| 15 | `precip_and_dissolution` | both directions on one mineral |

Cases 05, 06 and 07 are **the same case with one line changed**. So are 09 and
10. Those pairs are the useful ones: they isolate a single decision.

---

## What is worth checking, not just running

`runAllExamples.sh` checks that each case runs, exits cleanly and never
reports a negative concentration. That catches crashes and instabilities. It
does **not** check the numbers, and four of these cases have a right answer
that does not come from this code:

**02** — with no reaction and one face held at 1, the other at 0, the steady
profile is a **straight line**. Curvature means a reaction is firing somewhere
it should not, or a boundary is not being held.

**03, 13, 14, 15** — **mass balance.** One A plus one B makes exactly one C;
one mole of calcite removed appears as exactly one mole of Ca²⁺. Any drift is a
bug, not chemistry.

**09 and 10** — the toy model is four reactions, so the answer is a formula:

```
growth = min( Vs, 2 · Vo ),    V = Vmax · C / (Kc + C)
```

With Vmax 10 and 4 and both concentrations well above their half-saturation
constants, the acceptor limits and growth should settle near **8 mmol/gDW/h**.
That number comes from the stoichiometry. And 09 and 10 must agree with each
other, because the two solvers receive an identical linear program and share
nothing else.

---

## Why each folder carries its own kinetics headers

CompLaB compiles its rate laws in. `defineKinetics.hh` and
`defineAbioticKinetics.hh` are `#include`d, not read at run time, so two cases
with different chemistry need different binaries. `runAllExamples.sh` copies
each case's headers in and rebuilds, restoring your own on exit — including on
Ctrl-C.

That is not this suite being awkward; it is how the code is built. The 2D suite
worked the same way: every `example/scalability` case shipped its own
`defineKinetics.hh`.

The shipped headers implement the 95-substrate uranium network. Dropped into a
2-substrate case they index `C[0]` to `C[94]` and read past the end of the
array, so every example replaces them. `validateExamples.py` checks each
header's own `C.size()` guard against the case's substrate count, which catches
exactly the mistake of forgetting to.

Three cases need a different cmake configuration: 09 and 12 need
`-DENABLE_GLPK=ON`, 10 needs `-DENABLE_COBRAPY=ON`. The runner groups them, so
the tree is reconfigured three times rather than fifteen.

---

## validateExamples.py

A one-second structural check, run before you build anything. It verifies:

- every XML is well formed, and **no comment contains a double hyphen** — XML
  forbids it and this codebase has tripped over it repeatedly
- every tag name appears in the set the **source actually reads**, extracted
  mechanically from CompLB3D rather than typed out. A tag renamed upstream
  shows up as an unknown tag instead of as silence
- the vector-length and presence rules that terminate a run: `Kc` and `Vmax`
  lengths against the substrate count, `initial_densities` against the material
  numbers, `viscosity_ratio_in_biofilm` present for exactly the seeded microbes,
  `biomass_diffusion_coefficients` present when the solver is FD, FD rejected
  for planktonic microbes
- reaction types spelled acceptably, and the matching `enable_*` switch on
- geometry files: the right voxel count, and every material number the XML
  refers to actually present
- metabolic model files: `nmet`, `nrxn`, and the lengths of S, b, c, lb, ub;
  exchange indices inside range
- precipitation and dissolution pointing at substrates that are actually
  `<immobile>true</immobile>`
- the kinetics headers **compile**, at `-Wall -Wextra`

It is not Palabos. It builds no lattices and runs nothing.

**It was checked against deliberate breakage**, which is the only way to trust
a checker. Ten mutations were introduced one at a time — a short `Kc` vector, an
exchange index past the end of the model, `gplk` for `glpk`, the enable switch
turned off, a mistyped tag, FD on a planktonic microbe, a truncated S matrix, a
double hyphen in a comment, a kinetics guard written for other chemistry, a
syntax error in a rate law. **All ten were caught.**

---

## Extending the suite

`makeExamples.py` generates the folders from shared fragments in `_build/`.
Edit and regenerate rather than hand-patching fifteen XMLs; that is how the
suite stays consistent after the validator finds something.

A case is a `case(...)` call: a number, a name, the XML, a README, a geometry,
and optionally its own kinetics headers and a metabolic model.

---

## Not covered here

- **checkpoint restart** — `<read_NS_file>` and `<read_ADE_file>`, which resume
  from a `.chk` written by `<save_CHK_interval>`. Set `save_CHK_interval` in any
  case above, run it, then rerun with the read flags on.
- **MPI** — every case runs under `mpirun` unchanged; the domain is too small
  for the scaling to mean anything.
- **the scalability benchmarks** from the 2D suite (`example/scalability`,
  four cases). Those are performance measurements rather than capability
  demonstrations, and are still to be ported.

---

## Standing caveat

None of this has been compiled against real Palabos. The XML, the geometry, the
model files and the kinetics headers are all validated — tag names, lengths,
indices, well-formedness, and the headers genuinely compile — but the first real
build is on Tahoma.
