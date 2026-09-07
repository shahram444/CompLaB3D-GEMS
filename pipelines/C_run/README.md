# C — Running the solver

**In:** everything A and B produced.
**Out:** VTI fields, a log, and the reports.
**Run:** `./run.sh`, or the three commands below.

## Assemble, build, run

```bash
# 1. put a working directory together from an example
./scripts/setup_case.sh 13_precipitation run/mycase

# 2. build
cd run/mycase
cmake -B build -S . && cmake --build build -j

# 3. run
./complab CompLaB.xml                 # serial
mpirun -np 8 ./complab CompLaB.xml    # parallel
```

`setup_case.sh` exists because the examples do not each carry a copy of the two
kinetics headers. It lays down the shared defaults from `config/kinetics/`, then
puts the example's own chemistry on top where it has any. See
[`../../examples/README.md`](../../examples/README.md).

## Read the start-up lines before you walk away

The solver prints what it decided in the first few lines, and every one of them
is a chance to catch a mistake early:

```
geometry     128 x 64 x 64, porosity 0.312, percolates: yes
precipitation enabled, solid_substrate 7, max_precipRho 5.0e-4
dissolution  disabled
rate path    symbolic, file mylaw.sym, 4 rates, 1 variable
range        acetate 2.59e-05 .. 3.89e-03, enforced
```

If any of those lines disagrees with what you meant, stop now.

## One limitation to know before a production run

**Dissolution products are lost at MPI block boundaries.** A dissolving voxel
deposits its products into its open face neighbours; when a neighbour lies in
another MPI block the deposit is written into the local envelope and is not
communicated back. Dissolution results therefore depend on the processor count.

A run that needs a quantitative dissolution mass balance should be done on a
single process, or checked against a single-process run at reduced resolution.
Precipitation is unaffected.
