# A2 — Starting fields

**In:** nothing, or a previous run.
**Out:** initial concentration and biomass arrays.
**Run:** only if you want something other than the uniform values in the XML.

Most cases do not need this stage. The `<initial>` value on each species in
`CompLaB.xml` fills the whole domain with one number, which is what you want for
a clean start.

You need this stage when:

- **biomass sits somewhere specific** — a biofilm on one wall rather than
  everywhere, which is the usual setup for the two-population cases;
- **you are restarting** from the end of an earlier run;
- **a plume enters from one face** rather than the whole inlet.

```bash
# put an iron-reducing population on the lower wall only, by marking those
# voxels with a material code the XML then gives an initial density to
python ../../../tools/setup/geometry.py seed geometry.dat --nx 128 --ny 64 --nz 64 \
       --code 3 --box 0 127 0 2 0 63 -o input/geometry_seeded.dat
```

Then point the XML at what you made:

```xml
<microbe0>
    <initial_file>input/biomass0.dat</initial_file>
</microbe0>
```
