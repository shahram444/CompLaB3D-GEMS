# Reference pore geometries

Four small geometries to start from. **No example reads them.** Every case in
[`../../examples/`](../../examples/) generates its own `input/geometry.dat` with
its own `preprocess.py`, so that an example folder is the whole procedure rather
than a set of pointers; these four are here as clean starting points for a case
of your own, and as the plainest illustration of the material-number convention
below.

A geometry file is one integer per line, in `x`-fastest order, with as many lines
as `nx*ny*nz` in the case's `CompLaB.xml`. All four here are the staggered slot
at 24 × 24 × 6 = 3456 voxels: small enough to run in seconds on a laptop, large
enough that the flow field is not trivially symmetric.

| File | Materials present | What it is for |
|---|---|---|
| `slot_bare.dat` | 1 wall, 2 pore | Flow, transport and chemistry with no organisms. The shape examples 01–04 draw. |
| `slot_one_biofilm.dat` | 1, 2, 3 biofilm | One population seeded as a patch on each wall. The shape examples 05–07, 09–12 and 16–18 draw. |
| `slot_two_biofilms.dat` | 1, 2, 3 and 4 | Two populations facing each other across the slot. Example 08's shape. |
| `slot_mineral_phase.dat` | 0 mineral, 1, 2 | A declared solid phase lining the pore-facing grain surfaces, so it can dissolve or be added to. The shape examples 13–15 draw. |

**These four are unpadded, and the examples' own files are not.** The solver
conditions only `x = 0` and `x = nx-1`; the other four faces get nothing — not a
wall, not a symmetry plane, not periodicity. Each `preprocess.py` therefore
**adds** an inert wall layer on those four faces, which is why a case declaring
this shape runs at 24 × 26 × 8 = 4992 voxels. The pore space itself is untouched,
so porosity and every voxel count are what the case declares. Do the same to any
geometry of your own, or condition those faces yourself.

## What the numbers mean

The material number is read by the solver as follows, and the same convention
holds for a geometry you make yourself:

| Value | Meaning |
|---|---|
| `0` | A declared solid phase — one of the `<phaseN>` blocks under `<dissolution>`. It can be consumed, and the voxel reopens when it is gone. |
| `1` | Inert wall. Never changes. |
| `2` | Open pore, water. |
| `3`, `4`, … | Open pore holding an initial biomass patch: `3` is `microbe0`, `4` is `microbe1`, and so on. |

A voxel converted to solid by precipitation during a run takes the
`<solid_substrate>` value given in `<precipitation>`, not one of these.

## Making your own

```bash
python tools/geometry.py --nx 128 --ny 64 --nz 64 --type slot --out input/geometry.dat
python tools/geometry.py --inspect input/geometry.dat        # counts, porosity, percolation
```

`--type` also takes `sphere_pack`, `channel` and `image` (which reads a stack of
segmented TIFFs). Run it with `--help` for the full list. Inspect before you run:
a geometry whose pore space does not percolate will produce a flow field of
zeros and no error message.

## How an example picks one

It does not copy the file. It names it in its own `case.files`:

```
config/geometry/slot_one_biofilm.dat  input/geometry.dat
```

and `scripts/setup_case.sh` puts it in place when the case is assembled. To run
an example on your own geometry, assemble it first and then replace
`input/geometry.dat` in the assembled directory — the example itself stays
untouched.
