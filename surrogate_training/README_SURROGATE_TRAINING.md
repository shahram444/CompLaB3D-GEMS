# Surrogate training toolkit

The missing piece of CompLaB's surrogate option: the code that **makes** a
surrogate network, not just the code that evaluates one.

`surrogateModel.hh` has shipped since CompLaB v1.0 with a trained network
hard-coded into it. The network is unmistakably MATLAB Deep Learning Toolbox
output — `IW1_1`, `LW2_1`, `tansig`, `mapminmax` gains and offsets are
genFunction's own naming — but **the training script was never distributed**.
The weights could be evaluated and never reproduced, retrained, checked, or
fitted to a different organism. This directory closes that gap.

---

## The three steps

```
   your genome-scale model (.xml / .mat / .json)
                 |
                 |  generateTrainingData.py      offline FBA sweep
                 v
   training_data.csv                              inputs -> growth rate
                 |
                 |  trainSurrogate.py   or   trainSurrogate.m
                 v
   surrogate_weights_<name>.hh                    generated C++
                 |
                 |  verifyExport.py               compile and compare
                 v
   two lines in surrogateModel.hh, and you are running
```

### Worked example

```bash
# 1. sweep the FBA offline
python3 generateTrainingData.py iAF987.xml \
    --exchange EX_ac_e   --range 1e-3 10   --log \
    --exchange EX_fe3_e  --range 1e-5 0.5  --log \
    --samples 4000 --seed 1 \
    -o geobacter_training.csv

# 2. fit and export
python3 trainSurrogate.py geobacter_training.csv \
    --name geobacter --microbe-id 0 \
    -o surrogate_weights_geobacter.hh

# 3. prove the C++ matches the fit
python3 verifyExport.py surrogate_weights_geobacter.hh
```

Then, in `surrogateModel.hh`:

```cpp
#include "surrogate_weights_geobacter.hh"
...
else if (microbeId == 0) {
    bioR[microbeId] = srgnet_geobacter::eval(Fin[microbeId]);
}
```

and in `CompLaB.xml`:

```xml
<enable_surrogate>true</enable_surrogate>
...
<reaction_type>surrogate</reaction_type>
```

Step 3 is not optional. An exporter that transposes a matrix or drops a digit
produces a header that compiles, runs, and is quietly wrong — there is no
independent reference to notice. So the trainer writes 512 reference
predictions next to the header, and `verifyExport.py` compiles the header and
checks them.

---

## MATLAB or Python

| | |
|---|---|
| `trainSurrogate.m` | MATLAB, Deep Learning Toolbox, `trainlm`. Reconstructs how the shipped network was made. |
| `trainSurrogate.py` | numpy + scipy, L-BFGS-B on analytic gradients. No licence needed. |

Both write the **same file format**, and `testExporterParity.m` proves it:
it reads a Python-generated header, re-emits it with the MATLAB emitter, and
diffs. Verified byte-identical on 1-, 2- and 3-input networks. A network
trained on one path drops into a build that used the other.

They are not bit-identical *fits* — no two gradient-based optimisers ever are —
but they minimise the same loss over the same architecture with the same
70/15/15 split and early stopping, and land in the same place to within the
scatter of the random restart.

---

## The files

| file | what it does |
|---|---|
| `generateTrainingData.py` | runs the FBA sweep, writes the CSV with its metadata |
| `trainSurrogate.py` | fits and exports (Python path) |
| `trainSurrogate.m` | fits and exports (MATLAB path) |
| `exportSurrogateHeader.m` | the MATLAB C++ emitter, kept separate so it can be tested without a toolbox |
| `verifyExport.py` | compiles the generated header and checks it against the fit |
| `testExporterParity.m` | proves the two emitters agree byte for byte (runs in Octave) |
| `inspectSurrogate.py` | reads an existing network and recovers what nobody wrote down |

---

## What `inspectSurrogate.py` recovered about the shipped network

The training range of a network is never written down anywhere — but
`mapminmax` is invertible, so it can be recovered exactly:

```
gain = 2 / (x_max - x_min),  offset = x_min   =>   x_max = offset + 2/gain
```

Run against the shipped `surrogateModel.hh`:

```
Architecture : 2 -> 10 -> 10 -> 10 -> 10 -> 1   (371 parameters)

RECOVERED TRAINING RANGE
  input 0 : 0.0008901598344 .. 9.997925644   mmol/gDW/h   (acetate)
  input 1 : 3.513032548e-05 .. 0.49974581    mmol/gDW/h   (Fe(III))
  output  : 0 .. 0.05274516479 1/h
```

**Read that upper bound on input 1.** The shipped network has never seen an
electron-acceptor uptake above **0.5 mmol/gDW/h**. If your
`<fba_maximum_uptake_flux>` for that substrate exceeds 0.5, the simulation
spends its time outside the training box, and a neural network outside its
training box does not fail — it returns confident nonsense. Ask it for
`(50, 5)` and it happily reports 0.0486 1/h.

This is worth stating plainly in the paper: **the shipped surrogate encodes one
specific organism on two specific substrates over one specific range.** It is a
demonstration of the mechanism, not a general-purpose solver option the way
GLPK is. Anyone applying it to their own system needs to retrain — which,
until now, they could not do.

`inTrainingRange()` in every generated header is the guard that was missing;
call it if a run may wander.

---

## Choosing the training range

Use the range the **simulation** will visit, not the range the organism can
tolerate. CompLaB's Michaelis–Menten bound is

```
v = Vmax * C / (Kc + C)
```

so the largest uptake any voxel can ever request is `Vmax`, i.e. your
`<fba_maximum_uptake_flux>` entry — that is the natural upper end of the sweep.
The lower end should be small but non-zero.

Sample **logarithmically** (`--log`) whenever the range spans more than about
one order of magnitude. Pore-network concentrations routinely span four, and a
linear grid puts nearly all of its points in the top decade, so the fit is
worst exactly where the biology is most sensitive.

Keep it to **1–3 inputs**. Past that the sample count needed for the same
accuracy grows fast and the fit stops being worth the trouble.

---

## When a surrogate is and is not the right choice

**Use one** when FBA cost dominates the run, few substrates actually limit
growth, and the simulation stays inside the trained box. Speed-ups of 100–1000×
versus in-line FBA are typical.

**Do not** when you need the individual internal fluxes — a surrogate
reproduces the objective, i.e. growth, and nothing else — or when the
simulation may wander outside the training range.

---

## Verified

| check | result |
|---|---|
| Python trainer → generated C++, compiled and compared at 512 points | agree to **2.7e-13** relative |
| same, 1-input net (`2 -> 8 -> 8 -> 1`) | **7.9e-14** |
| same, 3-input net with `--log-output` | **1.1e-13** |
| MATLAB emitter vs Python emitter, 1/2/3-input nets | **byte-identical** (205, 155, 210 lines) |
| `inspectSurrogate.py` vs the compiled shipped network, 1681-point grid | **3.5e-17** absolute |
| generated header alongside the shipped one, two nets in one build | compiles clean at `-Wall -Wextra`, both evaluate correctly |
| end-to-end on a real model (*E. coli* core, glucose × O₂, 900 FBA solves) | test R² **0.9994**, RMSE 0.0036 1/h |

All C++ here was compiled against a Palabos stub, like the rest of the port.
**First real build is still on Tahoma.**

## Requirements

```bash
python3 -m pip install cobrapy numpy scipy
```

MATLAB path additionally needs the Deep Learning Toolbox. `testExporterParity.m`
needs neither — it runs in GNU Octave.
