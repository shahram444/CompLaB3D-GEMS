# B — Offline models

Everything here runs **before** the solver, produces a file, and then never runs
again. Whether you need any of it depends entirely on the rate path.

| | Produces | Needed for | Typical cost |
|---|---|---|---|
| [B1](B1_metabolic_model/) | a metabolic model the solver can read | flux balance analysis, and B2 | minutes |
| [B2](B2_surrogate_network/) | weights compiled into a header | the surrogate path | hours |
| [B3](B3_symbolic_law/) | a `.sym` text file | the symbolic path | minutes to hours |
| [B4](B4_graph_network/) | a `.gnn` text file | the graph-network path | minutes |

**If you are using compiled kinetics, skip this whole directory.** You write the
rate law in C++ and there is nothing to prepare.

## The one dependency

B2 needs B1 first. The surrogate is fitted to flux balance output, so the
metabolic model has to exist and be solvable before there is anything to fit.
The other three are independent of each other and of B1.

## What separates B2 from B3 and B4

B2 ends in a **rebuild**. The surrogate weights are pasted into
`src/surrogateModel.hh` and compiled in, so changing the fitted model means
recompiling the solver. B3 and B4 end in a **text file** that the solver reads
at start-up, so a new law can be tried by editing a file.

That is a real difference in how it feels to work: swapping a `.sym` file is a
few seconds, swapping a surrogate is a build.
