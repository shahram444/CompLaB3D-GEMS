# B2 — Fitting the surrogate network

**In:** the metabolic model from B1, and the range of uptake rates you expect.
**Out:** a header holding the trained weights, compiled into the solver.
**Run:** `./run.sh`, then rebuild.
**Needed for:** the surrogate path only.

This is the longest preparation in the repository, and the only one that ends in
a recompile. The tools live in [`tools/surrogate/`](../../../tools/surrogate/),
and there is a MATLAB path alongside the Python one; either produces the same
header.

## The one judgement you have to make

**Choose the range of uptake rates the simulation will visit.** Everything
downstream is valid only inside it. Too narrow and the run spends its time
extrapolating; too wide and the fit spreads its accuracy over states that never
occur.

The surrogate path does not enforce that range at run time. A voxel outside the
fitted box gets a confident answer and no warning, which is the reason for the
last step below.

## The four steps

```bash
# 1. sweep the linear program over the range you chose
python3 ../../../tools/surrogate/generateTrainingData.py \
       ../../../models/e_coli_core.xml.gz \
       --exchange EX_glc__D_e --range 0.001 10 --log \
       --exchange EX_o2_e     --range 0.00003 0.5 --log \
       --grid 141 \
       -o training_data.csv

# 2. fit the network
python3 ../../../tools/surrogate/trainSurrogate.py training_data.csv \
       --name geobacter --layers 10 10 10 10 --restarts 5 \
       -o surrogate_weights_geobacter.hh

# 3. put the header where the solver expects it
cp surrogate_weights_geobacter.hh ../../../src/surrogateModel.hh

# 4. rebuild
cd ../../.. && cmake --build build -j
```

`--grid 141` sweeps a full 141 x 141 lattice. `--samples N` draws N points at
random instead, which is the better choice above two substrates.

## Then look at what you fitted, before you trust it

```bash
# check the exported header reproduces the trainer, to machine precision
python3 ../../../tools/surrogate/verifyExport.py ../../../src/surrogateModel.hh

# evaluate it at a point you care about
python3 ../../../tools/surrogate/inspectSurrogate.py ../../../src/surrogateModel.hh \
       --eval 9.0 0.45
```

`verifyExport.py` matters more than it sounds. The trainer and the solver are
two separate implementations of the same arithmetic, and nothing about a header
file by itself stops them drifting apart.

Then **plot the response over the whole training box**, because two things are
invisible in a table of weights:

- **a dead region.** The network shipped with this repository returns exactly
  zero over 41.6 per cent of its own training box, below a donor uptake of 6.86.
  A simulation whose voxels sit mostly there is running a model that predicts no
  growth anywhere.
- **where your operating point sits.** If it is at the edge of the box, widen
  the sweep and refit.

## The MATLAB path

`trainSurrogate.m` and `exportSurrogateHeader.m` do the same job through
MATLAB's neural network toolbox, and `testExporterParity.m` checks the two
exporters agree. Use whichever you already have installed; the header is
identical either way.
