# tools — the pipeline, end to end

Everything CompLB3D needs on the way in, and everything it does not produce on
the way out. No step here says "now write a script".

```
   YOUR DATA                        CompLB3D                    YOUR ANSWER
   ─────────                        ────────                    ───────────

   CT scan / .npy / raw ─┐
   or nothing at all ────┼─► geometry.py ──► geometry.dat ─┐
                         │                                 │
   chemistry, as          │                                 │
   reactions ────────────┼─► makeKinetics.py ──► *.hh ─────┼─►  complab  ─┐
                         │                                 │              │
   chemistry, as          │                                 │              │
   formulas ─────────────┼─► makeEquilibrium.py ──► XML ───┤              │
                         │                                 │              │
   genome-scale model ───┴─► extractMM.py ──► matrix ──────┘              │
                             (in the parent folder)                        │
                                                                           │
   surrogate_training/  ──► a fitted network ─────────────────────────────┘
                                                                           │
                            postprocess.py ◄──────── output/*.vti ◄────────┘
                                  │
                                  └──► summary.csv, profiles.csv,
                                       massbalance.txt, plots/*.png
```

| tool | closes the gap |
|---|---|
| `geometry.py` | pore space in, from an image stack, an array, or a recipe — and checked before you spend compute on it |
| `makeKinetics.py` | rate laws without writing C++, with a mass-balance self-test |
| `makeEquilibrium.py` | the speciation tableau derived from formulas, not assembled by hand |
| `postprocess.py` | **numbers out.** The solver writes `.vti` snapshots and nothing else — no CSV, no totals, no porosity curve |
| `vtireader.py` | a `.vti` reader that needs only numpy, so post-processing runs on the machine that produced the data |

---

## geometry.py — pore space in

```bash
# build one
geometry.py create random --nx 100 --ny 100 --nz 100 --porosity 0.35 \
    --walls y -o input/geometry.dat --vti

# or import a segmented CT stack, one image per z-slice
geometry.py import ct_slices/ --threshold 128 -o input/geometry.dat
geometry.py import volume.npy --downsample 2 -o input/geometry.dat
geometry.py import scan.raw --raw-dims 512 512 400 --raw-dtype uint16 -o input/geometry.dat

# add a biofilm patch
geometry.py seed input/geometry.dat --nx 100 --ny 100 --nz 100 \
    --code 3 --box 10 30 1 5 0 99 -o input/seeded.dat

# ALWAYS do this before a long run
geometry.py inspect input/geometry.dat --nx 100 --ny 100 --nz 100
```

Recipes: `channel`, `spheres`, `cylinders`, `random`, `fracture`, `layered`.

`inspect` is the one that saves days. It reports porosity per material, whether
the pore space **percolates** from inlet to outlet, how much of it is
**isolated**, and whether the domain is quasi-2D. A sealed domain or 30% dead
pore space is not something you want to discover from a converged run that
looks odd.

The file format has no header, so the loop order is the entire contract with
the solver — and the wrong order gives a file of exactly the right length
describing a transposed pore space. It is written in one place, here, and
round-trip tested on an array with three different dimensions so a transpose
cannot hide.

## makeKinetics.py — chemistry without C++

Write the reactions:

```
substrates: A B C
microbes:   Bug

abiotic R1:
    equation: A + B -> C
    k: 5.0e-2

microbe Bug:
    mu_max: 2.0e-4
    monod: A 0.05
    yield: A 0.4
    release: C 2.0
    decay: 1.0e-6
```

```bash
makeKinetics.py mychemistry.rxn -o .
g++ -O2 -std=c++11 -o selftest selftest_kinetics.cpp && ./selftest
```

You get `defineKinetics.hh`, `defineAbioticKinetics.hh`, and a self-test — with
the index comments, the sign conventions, the negative-flooring, the size guard
and the NaN checks all correct by construction.

**The self-test is the point.** It computes the conserved moieties of your
reaction network (the left null space of the stoichiometry matrix) and checks
them over 20,000 random states. That is the correct general statement of mass
balance, and it stays valid when several reactions share a substrate — where
a naive per-reaction ratio check is meaningless because the rates are sums.
It needs neither Palabos nor a working build of the solver.

It does not catch everything: reading the wrong concentration still conserves
mass. It catches wrong stoichiometry, dropped coefficients and flipped signs,
which is where generators actually go wrong.

Covers mass action, multi-Monod growth with inhibition, and first-order
dissolution. Anything else, write by hand — the generated file is ordinary
readable C++ and a good starting point.

## makeEquilibrium.py — the tableau from formulas

```
component HCO3-1
component H+1

species HCO3-1   logK 0
species H+1      logK 0
species CO3-2    logK -10.33
species H2CO3    logK 6.35
inert FeS
```

```bash
makeEquilibrium.py carbonate.spc
```

It balances elements and charge to derive each row, then prints the
`<equilibrium>` block ready to paste. It **refuses** rather than guessing when
the components are not independent, and it errors with the unbalanced elements
named when a species cannot be made from them.

## postprocess.py — numbers out

```bash
postprocess.py path/to/run --conserve "A+C" --conserve "B+C"
```

| file | what |
|---|---|
| `summary.csv` | one row per snapshot: total, mean, min, max of every field, plus porosity |
| `profiles.csv` | the x-profile of every field at the last step |
| `massbalance.txt` | your conservation checks, plus a negative-value check, with a verdict |
| `plots/*.png` | totals, porosity, profiles |

Column names come from `CompLaB.xml`, so they read `acetate`, not
`subsLattice0`. Totals are over open voxels only, so a run that seals pore
space is not compared against a moving denominator.

**Exit code 1 when a check fails**, so it drops into CI or
`runAllExamples.sh` unchanged.

---

## Verified

| check | result |
|---|---|
| `.vti` reader vs real VTK output, all 6 encodings (ascii / binary / appended × plain / zlib) | **exact**, every value |
| the same with a `UInt64` header | exact |
| geometry write→read round trip on a (3, 5, 7) array — all dimensions different, so a transpose cannot hide | **identical** |
| `.npy` and raw-binary import paths | byte-identical output |
| `inspect` on a sealed domain, and on one with a 5% isolated cavity | both flagged, exit 1 |
| `postprocess.py` on a run built to conserve exactly | drift **1.7e-16**, PASS |
| the same run with 3% of C deleted and one voxel driven negative | **both caught**, exit 1 |
| `makeKinetics.py` self-test on a 2-reaction network with shared substrates | 4 conserved moieties, imbalance **3.9e-16** |
| the same with a coefficient dropped, and with a sign flipped | **both caught**, exit 1 |
| generated kinetics vs the hand-written example 03, 50,000 random states | **worst difference 0.000e+00** |
| `makeEquilibrium.py` on the carbonate system | reproduces example 04's hand-built tableau **exactly** |
| the same with dependent components, and with an impossible species | both refused with a named reason |

## Requirements

```bash
python3 -m pip install numpy            # required
python3 -m pip install matplotlib       # optional: plots
python3 -m pip install pillow           # optional: image-stack import
python3 -m pip install scipy            # optional: faster connected components
```

`vtk` is **not** required.
