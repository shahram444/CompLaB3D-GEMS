# A1 — The pore space

**In:** a description of the geometry you want, or an image stack.
**Out:** `geometry.dat`, one integer per voxel.
**Run:** `./run.sh`, or `tools/geometry.py` directly.

## What the integers mean

| Value | Meaning |
|---|---|
| `0` | open pore: fluid flows, solutes are transported |
| `1` | wall: bounce-back, no flow, no transport |
| `2` and above | a named grain phase, referred to by `material_number` |

The third row is what makes dissolution possible. A grain labelled `2` can be
declared as a dissolvable phase in `CompLaB.xml`; a grain labelled `1` is inert
rock for ever. **If you intend to dissolve anything, label it now** — you cannot
add the distinction later without rebuilding the geometry.

## Ways to make one

```bash
# a slot pore, the simplest case that can clog
python ../../../tools/geometry.py --nx 128 --ny 64 --nz 64 \
       --type slot --width 20 --out geometry.dat

# a packed sphere pack with two grain phases
python ../../../tools/geometry.py --nx 200 --ny 200 --nz 200 \
       --type spheres --radius 12 --phases 2 --out geometry.dat

# from a segmented micro-CT stack
python ../../../tools/geometry.py --from-images scan/*.tif \
       --threshold 128 --out geometry.dat
```

## Check it before you use it

```bash
python ../../../tools/geometry.py --inspect geometry.dat
```

prints the dimensions, the porosity, the count of each label, and whether the
pore space percolates from inlet to outlet. **A domain that does not percolate
will run and produce nothing**, which is a slow way to discover a typo.
