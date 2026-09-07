# D1 — Fields and pictures

**In:** the VTI output of a run.
**Out:** slices, profiles, animations.
**Run:** `./run.sh`, or `tools/postprocess.py`.

```bash
# a mid-plane slice of every species at the last step
python ../../../tools/postprocess.py --dir run/mycase --slice z --at 0.5 --out figs/

# porosity and permeability against time
python ../../../tools/postprocess.py --dir run/mycase --history --out figs/

# an animation of the mineral front
python ../../../tools/postprocess.py --dir run/mycase --animate FeS \
       --slice z --at 0.5 --out mineral.gif
```

`vtireader.py` is the underlying reader if you would rather work in your own
script:

```python
from tools.vtireader import read_vti
f = read_vti('run/mycase/output/step_050000.vti')
print(f.keys())          # the species and fields present
print(f['FeS'].shape)    # nx, ny, nz
```
