# Example 15: precipitation and dissolution together

Both directions on the same mineral. FeS forms where the iron and
sulfide fronts meet and seals its voxels; acid arriving from the left
eats it back out again and those voxels reopen.

## What to expect

- porosity **falling, then rising** as the acid front catches up with
  the seal
- the porosity settling rather than oscillating. Chattering means
  `reopen_fraction` is too close to 1
- **mass conserved across both directions**: iron and sulfur are only
  moved between the water and the mineral, never created

## The switches this case exercises

`<precipitation>` and `<dissolution>` both enabled, on one mineral,
with `<is_precipitate>true</is_precipitate>` linking them.

## Notes

This is the case that justifies `reopen_fraction`. Set it to 0.99 and
rerun: you should see the porosity flip back and forth every update
interval, with a flow solve wasted on each flip. That is the failure
the hysteresis exists to prevent.

## Running it

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && cd build && cmake .. && make
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case.
