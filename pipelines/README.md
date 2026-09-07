# Pipelines — what to run, and when

Most of the work in a CompLaB3D-GEMS study happens **outside** the solver. The
solver reads files; something has to make them first. This directory is the map
of what those somethings are.

Four stages, in the order they happen:

| | Stage | When you need it | Where |
|---|---|---|---|
| **A** | **Pre-processing** | Always. Something has to say what the pore space looks like. | [`A_preprocess/`](A_preprocess/) |
| **B** | **Offline models** | Only for the rate path you chose. Four of the six need nothing here. | [`B_offline_models/`](B_offline_models/) |
| **C** | **The run** | Always. | [`C_run/`](C_run/) |
| **D** | **Post-processing** | Always, unless you enjoy reading raw VTI. | [`D_postprocess/`](D_postprocess/) |

---

## Which of these do I actually need?

Start from the rate path. That single choice decides whether stage B exists for
you at all.

| Your rate path | A: pre | B: offline | C: run | D: post |
|---|---|---|---|---|
| Compiled kinetics — you write the rate law in C++ | A1 | **nothing** | yes | D1, D2 |
| Flux balance analysis, GLPK | A1 | **B1** | yes | D1, D2 |
| Flux balance analysis, COBRApy | A1 | **B1** | yes | D1, D2 |
| Surrogate network | A1 | **B1 → B2** | yes | D1, D2 |
| Symbolic law | A1 | **B3** | yes | D1, D2 |
| Graph network | A1 | **B4** | yes | D1, D2 |

B3 and B4 produce a file the solver reads, so neither needs a rebuild — and
neither is needed at all to *try* those paths: examples 17 and 18 ship the
`.sym` and the `.gnn` they read. Stage B is where you make your own.

Two things worth reading off that table before you start.

**The surrogate is the only path with a chain.** You cannot fit a surrogate
without first having a metabolic model to sweep, so B1 comes before B2. Budget
for that: B2 is the most expensive preparation in the repository, because it
solves a linear program at every point of a grid.

**Precipitation and dissolution add nothing to this table.** They are geometry
processes, not rate paths, and they run alongside whichever rate path you chose.
If you are doing precipitation with compiled kinetics, your offline work is
still nothing at all: edit two files and start the solver.

---

## The shape of every stage

Each numbered directory holds the same three things, so you always know where to
look:

```
B3_symbolic_law/
├── README.md      what it is for, what goes in, what comes out, what to run
├── run.sh         the command, with the arguments filled in
└── expected/      what a correct result looks like, so you can tell
```

`run.sh` is meant to be read as much as executed. Every one is short, and every
one names its inputs and outputs at the top.

---

## The order, drawn out

```
   A1 geometry.dat ─────────────────────────────────────┐
                                                        │
   your rate path ──┬─ compiled kinetics ───────────────┤
                    │                                   │
                    ├─ FBA ── B1 model export ──────────┤
                    │            │                      │
                    ├─ surrogate ─┴─ B2 fit + rebuild ──┤
                    │                                   │
                    ├─ symbolic ── B3 search → .sym ────┤
                    │                                   │
                    └─ graph net ── B4 train → .gnn ────┤
                                                        │
                                                        ▼
                                              C: run the solver
                                                        │
                                                        ▼
                                        D1 fields    D2 balances
```

---

## A worked order, start to finish

Somebody doing the precipitation case with a symbolic rate law would run, in
this order:

```bash
# A1 — make the pore space
python tools/geometry.py --nx 128 --ny 64 --nz 64 --type slot --out input/geometry.dat

# B3 — find a rate law from data, offline, once
python tools/fit_symbolic.py --data mydata.csv --target growth \
       --pop 600 --gens 60 --out mylaw.sym

# C — assemble the case and run it
./scripts/setup_case.sh 13_precipitation run/mycase
cd run/mycase
cmake -B build -S . && cmake --build build -j
./complab CompLaB.xml

# D — read the result
python ../../tools/postprocess.py --dir . --report
```

Five commands. The two long ones are the fit and the run.
