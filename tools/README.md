# tools

Everything here runs **outside** the solver, on a workstation, before or after a
simulation. Nothing in this directory is compiled, and nothing in it runs during
a time step.

The folders are named for the rate path you are setting up, because that is the
question you have when you come looking: *I am running method N, what do I have
to produce first?* Each rate path has its own guide in
[`docs/methods/`](../docs/methods/), and the guide's Figure 2 is a flowchart of
the same work these scripts do.

---

## Which tool, for which method

| You are running | Produce this first | With | Documented in |
|---|---|---|---|
| **Method 1**, flux balance analysis<br>`<reaction_type>glpk`, `cobrapy` | the flat model XML the solver reads | [`method_1_fba/extractMM.py`](method_1_fba/extractMM.py) | `Method_1_FBA_guide.docx` |
| **Method 2**, the surrogate network<br>`<reaction_type>surrogate` | a `.srg` file, or a generated header | [`method_2_surrogate/`](method_2_surrogate/) — four steps, see its own README | `Method_2_Surrogate_guide.docx` |
| **Method 3**, symbolic rate laws<br>`<reaction_type>symbolic` | a `.sym` file | [`method_3_symbolic/make_rate_law.py`](method_3_symbolic/make_rate_law.py) — **start here**, it asks questions | `Method_3_Symbolic_guide.docx` |
| **Method 4**, the graph network<br>`<reaction_type>graphnet` | a `.gnn` file | [`method_4_graphnet/train_graphnet.py`](method_4_graphnet/train_graphnet.py) | `Method_4_GraphNetwork_guide.docx` |
| **Method 5**, the thermodynamic gate<br>`<thermodynamics>` | a `.thm` file | nothing — you write it by hand from the literature | `Method_5_Thermodynamics_guide.docx` |
| **compiled kinetics**<br>`<reaction_type>kinetics` | `defineKinetics.hh` | [`setup/makeKinetics.py`](setup/makeKinetics.py) | the manual, `docs/` |

Method 5 has no tool on purpose. Its file is four numbers per reaction taken
from the literature at the temperature of the case, and a script that guessed
them would be worse than none.

---

## Everything else

| Folder | What is in it | When you need it |
|---|---|---|
| [`setup/`](setup/) | `geometry.py`, `makeKinetics.py`, `makeEquilibrium.py` | before any run: the pore space, the compiled rate laws, the speciation tableau |
| [`postprocess/`](postprocess/) | `postprocess.py`, `vtireader.py`, `upscale_sweep.py` | after a run: totals and balances, reading `.vti` without the vtk package, and turning a sweep of runs into an upscaling law |
| [`runtime/`](runtime/) | `complab3d_cobrapy.py` | **not a tool.** The solver imports it during the run. `scripts/setup_case.sh` copies it into the case's `src/`, and `comp_cobrapy.sh` will refuse to start without it. |

---

## The shortest path through each

```bash
# Method 1: a genome-scale model, into the form the solver reads
python3 tools/method_1_fba/extractMM.py iAF987.xml -o geobacter.xml

# Method 2: sweep the program, fit the network, check the export
cd tools/method_2_surrogate
python3 generateTrainingData.py geobacter.xml --exchange EX_ac_e --range 1e-3 10 --log ...
python3 trainSurrogate.py training_data.csv --srg mymodel.srg
python3 verifyExport.py                       # compiles the header and compares

# Method 3: answer questions about your own measurements
python3 tools/method_3_symbolic/make_rate_law.py mydata.csv
#   ... it prints the fit_symbolic.py command that reproduces the same file

# Method 4: fit a graph network to a reaction network
python3 tools/method_4_graphnet/train_graphnet.py --stoich S.csv --data samples.csv -o net.gnn

# before any of it
python3 tools/setup/geometry.py create channel --nx 128 --ny 64 --nz 64
python3 tools/setup/geometry.py inspect geom.dat      # porosity, percolation, dead pore space

# after
python3 tools/postprocess/postprocess.py run/mycase --conserve "A+C"
```

Every script answers `--help`, and every one of them has a header explaining
what it is for before it explains how to call it.

---

## A note on the copies inside `examples/`

Several examples ship their own copy of a tool under `training/` so the example
runs without reaching outside its own folder. Those copies are checked against
these originals by `tests/check_repo.sh`, which fails if one has drifted. If you
edit a tool here, run that script before committing.
