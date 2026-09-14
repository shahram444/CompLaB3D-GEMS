# D1 — Fields and pictures

**In:** the VTI output of a run.
**Out:** slices, profiles, animations.
**Run:** `./run.sh`, or `tools/postprocess/postprocess.py`.

The tool takes the run folder as its one required argument and does the whole
job in a single pass: slices, histories, the totals and the plots.

```bash
# everything, written into run/mycase/output
python ../../../tools/postprocess/postprocess.py run/mycase

# the same, into a folder of your own
python ../../../tools/postprocess/postprocess.py run/mycase --output figs/

# numbers only, no plots, for a batch job with no display
python ../../../tools/postprocess/postprocess.py run/mycase --no-plots --quiet
```

There is no flag for one slice or one animation. If you want a particular
picture, take the fields with `vtireader.py` below and draw it yourself; that
is what the reader is for.

`vtireader.py` is the underlying reader if you would rather work in your own
script:

```python
from tools.vtireader import read_vti
f = read_vti('run/mycase/output/step_050000.vti')
print(f.keys())          # the species and fields present
print(f['FeS'].shape)    # nx, ny, nz
```
