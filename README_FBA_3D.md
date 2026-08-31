# 3D flux-balance verification case

The 3D version of the two-species verification case from CompLaB 2D
(`2D orginal/example/indimeshComparison/sim1`). It is the natural regression test
for the restored FBA solver, because the expected answer comes from published 2D
results rather than from this code.

## The biology

Two bacteria, neither of which can grow alone:

| | eats | excretes | cannot |
|---|---|---|---|
| *E. coli* K-12 (iJO1366) | lactose | acetate | make its own methionine |
| *S. enterica* LT2 (iRR1083) | acetate | methionine | eat lactose |

Each supplies what the other is missing. Their biomass settles at a fixed ratio
set by the metabolic stoichiometry of the two genome-scale models — **not** by
the geometry. That invariance is what the original sim1…sim6 series tested: six
different spatial seedings, same converged ratio.

## Files

```
CompLaB.xml                    the case, ready to run
input/geometry.dat             25 x 25 x 5, extruded from inputDomain1.dat
input/iJO1366_rGEM_c.xml       E. coli   313 metabolites, 432 reactions
input/iRR1083_rGEM_c.xml       S. enterica 224 metabolites, 346 reactions
CompLaB_reference_template.xml the general reference file, kept for lookup
```

## Running it

```bash
cd build && cmake -DENABLE_GLPK=ON .. && make
cd ..              # CompLaB.xml must be in the working directory
./complab          # or:  srun ./complab
```

## Running the same case on COBRApy instead

`CompLaB_cobrapy.xml` is the COBRApy twin of `CompLaB.xml` — same domain, same
geometry, same model files, same numbers. It differs by exactly four lines:

```xml
<enable_fba_glpk>false</enable_fba_glpk>
<enable_fba_cobrapy>true</enable_fba_cobrapy>
...
<reaction_type>cobrapy</reaction_type>     <!-- in both microbe blocks -->
```

(`<glpk_lpsolver>` is dropped, being GLPK-only.) Running both and comparing the
converged biomass ratio is the strongest check available on the FBA coupling:
the two solvers receive the identical linear programme, so a disagreement means
one of the two couplings is wrong, not that the biology is ambiguous.

```bash
python3 -m pip install cobrapy
cd build && cmake -DENABLE_COBRAPY=ON .. && make
cd .. && cp CompLaB_cobrapy.xml CompLaB.xml && ./complab
```

Two things COBRApy needs that GLPK does not:

- **`src/complab3d_cobrapy.py` must be on disk** at `<src_path>` when the
  executable runs. It is imported by name at run time, not compiled in.
- **`cobrapy` installed in the Python CMake linked against** — not merely on
  your `$PATH`.

COBRApy does **not** read SBML itself. It gets the same `input/*_rGEM_c.xml`
stoichiometry GLPK gets and rebuilds a `cobra.Model` in memory from it at
start-up, with anonymous numbered reactions. Both paths therefore depend on
`extractMM.py`. Expect one to two orders of magnitude slower per voxel, because
cobra rebuilds the problem through optlang on every `optimize()` call.

## What to check

- `[OK] XML configuration loaded and validated`
- the metabolic summary showing `FBA (GLPK) - 2 microbe(s)`
- both models loading with the metabolite and reaction counts above
- biomass of both species rising, then their **ratio going flat**
- no `[NEG!]` concentration warnings

Then compare the converged *E. coli* : *S. enterica* ratio against the 2D
result. If they agree, the 3D FBA coupling is doing what the published one did.

A second check worth running: swap the two `<microbeN>` material numbers

```xml
<microbe0>3 4</microbe0>      <!-- was 8   -->
<microbe1>8</microbe1>        <!-- was 3 4 -->
```

which reproduces `sim6`. The ratio must not move. That is the actual verification
— one run only tells you the code produced *a* number.

## Why quasi-2D

The domain is the original 25×25 mesh extruded 5 voxels deep, with z periodic, so
every z-plane sees identical conditions and the 3D answer **must** reproduce the
2D one. Any difference is a bug in the port, not new physics. Once it passes,
thicken z or use a real 3D geometry to test something genuinely three-dimensional.

There is no flow (`<Peclet>0</Peclet>`) — deliberately, so a failure points at the
FBA coupling rather than at the Navier-Stokes solver.

## Three things I had to fix to make this case work

These are corrections to the shipped code, not to the case, and they are in the
updated `src/complab3d_metabolic.hh` and `src/complab.cpp` that go with this.

**The `b` vector in both model files was the wrong length.** `extractMM.py` built
it over *reactions* (432 and 346) when the loader wants one entry per *metabolite*
(313 and 224). Harmless in 2D — `b` is all zeros for a steady-state FBA and
`nrxn > nmet`, so the wrapper read in bounds by luck — but the 3D loader now
checks the length and would have rejected both models. The copies in `input/` are
corrected; nothing physical changes.

**`fix_concentration` and `fix_lower_bounds` are written as words.** The 2D files
say `no yes no yes`; my parser read them as integers, so every 2D input file would
have thrown and silently fallen back to all-false. It accepts both spellings now.

**The finite-difference solver corrupted its own scratch lattices.** `ptr_fd_lattices`
is laid out `[B…, copies…, mask]` with the mask last, but
`updateLocalMaskNtotalLattices3D` assumes the cellular-automaton layout with the
mask second-to-last. It was therefore reading the last *copy* lattice as the mask
and writing to it. This case uses `solver_type FD`, exactly as the 2D original
does, so it would have hit this immediately. Fixed with a correctly shaped vector;
the same bug exists in the 2D code at `complab.cpp:878,882` and is worth fixing
there too, since the published verification runs went through it.

## Caveats

**Not yet compiled against real Palabos**, like everything else in this port. The
XML, the geometry and the model files are all validated — lengths, indices,
well-formedness — but the first real build is still on Tahoma.

**`<biomass_molar_mass>` is set to 24.6 g/mol**, the generic CH₁.₈O₀.₅N₀.₂ formula
weight. The 2D code had no such parameter because it kept biomass in gDW/L
directly; CompLB3D uses mol/L throughout, so the conversion has to happen
somewhere. If the 2D biomass units were already gDW/L, set this to 1.0 to
reproduce the 2D numbers exactly — worth checking before you compare ratios.

**`<fba_maximum_uptake_flux>` is my reading of the 2D intent.** The 2D file gave
only `substrate_lower_bounds` (−15, −25, −8, −8) and no explicit Vmax, so the 2D
code used the magnitude of those bounds as the Monod Vmax. I have written them
into the proper Vmax tag *and* kept them as the hard floor, which reproduces that
behaviour. If you meant the bounds purely as clamps, set the Vmax entries to 0 and
the organisms become supply-limited instead.
