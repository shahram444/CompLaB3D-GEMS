# A1 — The pore space

**In:** a description of the geometry you want, or an image stack.
**Out:** `geometry.dat`, one integer per voxel.
**Run:** `./run.sh`, or `tools/setup/geometry.py` directly.

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
# a channel, the simplest pore space that can clog
python ../../../tools/setup/geometry.py create channel --nx 128 --ny 64 --nz 64 \
       --aperture 20 -o geometry.dat

# a sphere pack
python ../../../tools/setup/geometry.py create spheres --nx 200 --ny 200 --nz 200 \
       --radius 12 -o geometry.dat

# from a segmented micro-CT stack, a directory of images or a .npy
python ../../../tools/setup/geometry.py import scan/ --threshold 128 \
       -o geometry.dat
```

## Check it before you use it

```bash
python ../../../tools/setup/geometry.py inspect geometry.dat --nx 128 --ny 64 --nz 64
```

prints the dimensions, the porosity, the count of each label, and whether the
pore space percolates from inlet to outlet. **A domain that does not percolate
will run and produce nothing**, which is a slow way to discover a typo.
