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
# put an iron-reducing population on the lower wall only
python ../../../tools/geometry.py --field biomass --shape wall --face ymin \
       --value 1.0e-4 --like geometry.dat --out input/biomass0.dat

# restart from step 50000 of a previous run
python ../../../tools/postprocess.py --extract-restart previous_run/ \
       --step 50000 --out input/
```

Then point the XML at what you made:

```xml
<microbe0>
    <initial_file>input/biomass0.dat</initial_file>
</microbe0>
```
