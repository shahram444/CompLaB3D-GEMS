# Step by step: from nothing to a trained PRT-DeepONet-3D

Every command, in order, from an empty directory to a number that says whether
the whole chain worked. Read `README_1_SCENARIO.md` first if you want to know
what is being simulated and why; this file assumes you do not, and tells you what
to type.

Commands are written from **inside this package directory**. `$V132` is the
CompLaB3D v1.3.2 tree, which is the folder above this one:

```bash
cd  "<...>/13 - v1.3.2, positivity clamp and batch sweep/ABC_sweep_package"
V132=".."
```

---

## What is in this package

```
README_1_SCENARIO.md         the physics, the numbers, and why each was chosen.
                             Section 8 lists all 270 runs and which split each is
                             in; sections 11 and 12 annotate every value in
                             CompLaB.xml and in the rate law.
README_2_STEP_BY_STEP.md     this file

kinetics/
  defineAbioticKinetics.hh   A + B -> C, R = k[A][B], k from PRT_KABIO
  defineKinetics.hh          every biotic rate zero, and still mandatory

xml/
  CompLaB.xml.template       one case, with @TOKEN@ placeholders

patch/
  apply_damkohler_feed_reference.py
                             makes the Damkohler banner work when nothing
                             starts in the domain. Idempotent.

tools/
  setup_build.sh             assemble a buildable tree carrying this chemistry
  make_geometries.py         the 10 sphere packs and the straight duct
  make_campaign.py           build / calibrate / status / retry
  check_campaign.py          check (the gate) / clamp (triage)
  collect_to_h5.py           .vti -> dataset.h5
  verify_dataset.py          the verdict on the CompLaB3D half
  train_prt3d.py             train, validate, test, and the final verdict
  analytic_channel.py        the duct against its closed form
  selftest.py                the whole chain on synthetic data, no solver needed
```

---

## Step 0. What you need

| | |
|---|---|
| CompLaB3D v1.3.2 | the tree this package sits in |
| Palabos 2.3.0 | unpacked somewhere; the build needs its path |
| a C++11 compiler and CMake | plus MPI if you want more than one rank |
| Python 3 | numpy and h5py. SciPy and matplotlib are optional: SciPy only speeds up one distance transform, which has an exact fallback, and matplotlib only draws the geometry previews. |
| PyTorch | only for step 11, and only if you use the standalone trainer |

```bash
python3 -c "import numpy, h5py; print('ok')"
```

---

## Step 1. Prove the pipeline before building anything

This runs every stage after the solver on synthetic data. Ten minutes, no
compiler, no cluster. If it fails, the problem is on this machine and not in your
case setup.

```bash
python3 tools/selftest.py --work /tmp/abc_selftest
```

It ends with `SELFTEST PASSED` and a list of the stages it exercised. What it
cannot test is the physics; step 10 does that.

---

## Step 2. Make the pore spaces

```bash
python3 tools/make_geometries.py --out geometries
```

Eleven directories, each holding `geometry.dat` for CompLaB and a `.npz` holding
the material map and the two distance fields for PRT. It prints the grain count
and achieved porosity of each, then round-trips every file and checks the shape,
the sealed y and z faces, that no no-dynamics voxel touches pore, and the x-slice
ordering. **All eleven should say they round-tripped.**

Look at `geometries/geom_0000/geom_0000.png`. Three orthogonal slices: white pore,
grey bounce-back, black grain interior. Black should never touch white.

`geometries/index.csv` is the manifest: gid, porosity, tortuosity, split.

---

## Step 3. Assemble and build the solver

```bash
./tools/setup_build.sh "$V132" ~/abc_build
cd ~/abc_build
cmake -B build -S . -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
cd -
```

`setup_build.sh` copies the solver, applies `patch/apply_damkohler_feed_reference.py`
**to the copy** and never to your original tree, and puts this package's two
kinetics headers at the root where `complab3d_processors_part1.hh` includes them
from. It writes `CHEMISTRY.txt` saying what went in.

The binary lands at `~/abc_build/complab`.

> **The one mistake that costs weeks.** If you build against the CompLaB defaults
> in `config/kinetics/` instead, the rate constant is compiled in and `PRT_KABIO`
> is ignored. Every case in the campaign then runs identical chemistry while its
> `params.json` records a different Damköhler. Every run succeeds. The dataset
> builds. It teaches the network that Damköhler does nothing, and nothing anywhere
> says so. Step 5 tests for exactly this, and so does step 9.

---

## Step 4. Build the campaign

```bash
python3 tools/make_campaign.py build \
    --geometries geometries --out campaign02 --complab ~/abc_build/complab \
    --throttle 12 --mem 16gb --time 02:00:00 --email you@example.edu
```

`--complab` is what makes the timestep **measured** rather than estimated: the
tool runs the solver once per geometry and Péclet for a single step and reads the
`[ADE] dt=` line it prints. That takes a couple of minutes and it is worth it,
because the timestep depends on the pore space (`README_1_SCENARIO.md`, section 3)
and both the run length and the rate law's per-step cap are sized from it.

Without `--complab` the campaign still builds, from an estimate, and says so on
the terminal, in `campaign.json` and in every `params.json`.

You get **270 cases**: 30 geometries x 3 Péclet (0.02, 0.2, 2.0) x 3 Damköhler
(0.1, 1.0, 10). Each run keeps 20 of its 21 snapshots, so the dataset is **5400
training pairs**. Each case is a
self-contained directory holding `CompLaB.xml`, `env.sh`, `params.json`,
`input/geometry.dat` and an empty `output/`.

They are already assigned to the three splits, by porosity, and the assignment
travels with the data from here to the trained model:

```
train        porosity 0.32, 0.38, 0.44, 0.47, 0.53, 0.59
             gids 0-2, 6-8, 12-17, 21-23, 27-29
             18 geometries, 162 runs, 3240 training pairs

validation   porosity 0.35, 0.50
             gids 3-5, 18-20
             6 geometries, 54 runs, 1080 training pairs

test         porosity 0.41, 0.56
             gids 9-11, 24-26
             6 geometries, 54 runs, 1080 training pairs
```

Both held-out levels in each split sit strictly inside the training range
[0.32, 0.59], so the question being asked is interpolation in porosity.

`campaign02/runs.csv` has a row per case with its gid, porosity, split, Pe, Da, rate
constant, step count and snapshot interval. `README_1_SCENARIO.md` section 8 has
the full matrix and the run-id arithmetic.

---

## Step 5. The gate

```bash
python3 tools/check_campaign.py check \
    --campaign campaign02 --complab ~/abc_build/complab
```

Fatal findings stop the submission. Warnings do not. It reads every case and then
runs one of them for a single step under `PRT_KABIO=1234.5678`, looking for that
number in the `[KIN]` line of the solver's own log.

**Do not submit until this says `no fatal findings`.** What it checks:

- the binary really does read `PRT_KABIO`
- the Damköhler banner reports the Damköhler the campaign asked for, which is also
  a check that the patch in step 3 went in
- the flow solver hit its Péclet target
- `geometry.dat` length against `<nx><ny><nz>`, the failure nothing else catches,
  because `readGeometry()` stops when it has what it wants and ignores the rest,
  giving a truncated pore space that still runs
- percolation along x, sealed y and z faces, no bare solid against pore
- `<abiotic_rate_scale>` pinned at 1.0 while `PRT_KABIO` varies, so the rate is
  not set twice
- `PRT_DT` against the timestep `params.json` recorded
- `ade_max_iT = vtk_interval * 20 + 1`, so every run holds the same 21 frames
- `<delta_P>` above zero, or the solver silently sets Pe = 0 and every case
  becomes pure diffusion
- cell Péclet, the per-step reactant fraction against the cap, the disk footprint,
  and that all three splits have geometries in them

---

## Step 6. Run it

**On one machine:**

```bash
cd campaign02
for c in runs/run_*; do ./run_one.sh "$c" ~/abc_build/complab; done
cd -
```

**On a cluster (the UGA GACRC layout):**

```bash
cd campaign02
sbatch submit.sbatch
cd -
```

A Slurm **array**, one task per case, not one big job. Array tasks are separate
processes, so one crash cannot take the rest down, and `run_one.sh` never
propagates a failure: it writes `status.json` either way. Edit the partition and
the memory in `submit.sbatch` to match your allocation.

Expect **23 to 31 core-hours** in total and about **21 GB** of output. With the
array throttled to 12 concurrent tasks that is **2 to 3 hours of wall clock**. No
single case runs longer than about 15 minutes, which is why each array task asks
for a 2 hour limit.

Each case leaves `output/*.vti`, `output/run.log`, `output/summary.csv` and
`status.json`.

---

## Step 7. See what came back

```bash
python3 tools/make_campaign.py status --campaign campaign02
python3 tools/check_campaign.py clamp --campaign campaign02
```

`status` counts done, failed and not started, and prints the wall time. `retry`
prints the `--array` list for whatever did not finish:

```bash
python3 tools/make_campaign.py retry --campaign campaign02
```

`clamp` reads the `[CLAMP]` line out of every `run.log` and ranks the runs by how
much reaction had to be refused. **Above about one per cent, that run is not the
reaction it was configured for**; rebuild it with a shorter step and run it again.

This step matters because the collector silently drops any run more than 5 per
cent negative into `failures.csv`. A Damköhler sweep sweeps toward larger k, so the
high-Da half is exactly the half at risk, and without this the first sign of
trouble is a dataset quietly smaller than the campaign that produced it.

---

## Step 8. Collect into one .h5

```bash
python3 tools/collect_to_h5.py \
    --campaign campaign02 --geometries geometries --out dataset
```

Writes `dataset/dataset.h5`, plus `report.md` and, if anything was dropped,
`failures.csv` with a reason per run.

**If you have the PRT-DeepONet-GUI, use its collector instead.** The campaign is
built to satisfy it exactly: every `params.json` carries all thirteen keys it
requires, the geometry directories are in the layout it reads, and the species
names match the `.vti` base names CompLaB writes.

```bash
python3 <PRT-DeepONet-GUI>/3D/tools/collect_complab_output.py \
    --campaign campaign02 --geometries geometries --out dataset
```

Either produces the same layout and everything downstream accepts both.

---

## Step 9. The verdict on the CompLaB3D half

```bash
python3 tools/verify_dataset.py --dataset dataset/dataset.h5
```

This is the one to read carefully. In rough order of how badly a No would hurt:

1. **Does Damköhler do anything?** For every geometry and Péclet, the two Damköhler
   cases must give different fields. If they do not, you built against the wrong
   header (step 3) and the whole campaign is one chemistry wearing sixty labels.
2. **Does Péclet do anything?** If not, `<delta_P>` was zero or the flow never
   solved.
3. **Does the pore space do anything?** If not, the geometries are not reaching
   the solver.
4. **Are the boundary conditions the ones that were asked for?** A held at 1.0 on
   x = 0 and B held at 1.0 on x = nx-1, and each falling across the domain the way
   it should. Catches a swapped boundary and a reversed x axis. It also reports the
   largest value anywhere against the feed: nothing here can concentrate, so
   anything above the feed was manufactured at a boundary. Over 1.05 x feed it
   warns, over 1.25 x feed it fails the dataset.
5. Nothing negative, C actually produced and produced inside rather than at a
   face, float16 headroom, and a split that shares no geometry.

It ends with `VERDICT:` and exits non-zero if the dataset is not usable.

---

## Step 10. The physics check

Steps 1 to 9 test that the pipeline carries numbers correctly. This tests that the
numbers are right. It uses the straight duct, because a duct has answers that can
be written down.

```bash
python3 tools/make_campaign.py build \
    --geometries geometries --out channel --channel --boundaries dirichlet \
    --pe 0 0.5 --da 0.1 10 --run-factor 2 --complab ~/abc_build/complab

cd channel && for c in runs/run_*; do ./run_one.sh "$c" ~/abc_build/complab; done; cd -

python3 tools/analytic_channel.py --run channel/runs/run_0000 \
    --geometries geometries --plot channel_check.png
```

`--run-factor 2` because the closed-form statements are steady-state statements.
At Pe = 0 the slowest diffusive mode relaxes in `L^2/(pi^2 D)` = 10.4 s, so
205 s is twenty time constants and comfortably steady. `analytic_channel.py`
reports the drift over the last snapshot interval and says so if it is not;
raise the factor if it complains. `--pe 0` turns the flow solver off entirely
(`complab.cpp:1156`), which is the case with the cleanest answer.

Four cases, one duct. They are longer than a sweep case and there are only four
of them.

`--boundaries dirichlet` is the default, and it is not optional here: an open face
has no boundary value to compare against. Holding every face makes the duct a
counter-diffusion cell, which has an answer that can be written down.
`analytic_channel.py` refuses a case built any other way rather than reporting a
meaningless number. The sweep is built the same way, so this result speaks directly
to the campaign rather than to a separate configuration.

It compares the three **reaction-free combinations** (`A - B`, `A + C`, `B + C`)
against the one-dimensional steady profile. Each of them is exact for any rate
constant at all, because the reaction is one-for-one-for-one and the three species
share a diffusivity (`README_1_SCENARIO.md`, section 9). A disagreement is
therefore a transport fault: unequal diffusivities, broken stoichiometry, or a
boundary that is not what the XML says.

At Pe = 0 they must hold to the tolerance. At Pe > 0 the comparison is approximate
and the residual grows like Pe², which is Taylor dispersion in a duct and is real
physics rather than a solver fault; the tool reports it rather than failing on it.

---

## Step 11. Train, validate and test

### The standalone way, which needs nothing but PyTorch

```bash
python3 tools/train_prt3d.py --dataset dataset/dataset.h5 --species C \
    --epochs 60 --out model
```

It honours the train / validation / test split recorded in the file, which is held
out **by porosity**: it trains on 0.32, 0.38, 0.44, 0.47, 0.53 and 0.59 (18
geometries, 162 runs, 3240 training pairs), validates on 0.35 and 0.50 (6
geometries, 54 runs, 1080 pairs) and tests on 0.41 and 0.56 (6 geometries, 54
runs, 1080 pairs). It prints those three lists
on startup, so check them against the table in `README_1_SCENARIO.md` section 8. The architecture is the published 3D one at this resolution: a five-block
geometry CNN, a parameter branch over (Pe, Da), and an eight-layer trunk over
(x, y, z, t, geodesic distance).

Train one model per species. `--species C` is the interesting one; `A` and `B` are
worth running too, and are easier.

### The GUI way, which is the one to use once you are past the first dataset

The same `.h5` goes straight into the PRT-DeepONet-GUI, or into its trainer:

```bash
python3 <PRT-DeepONet-GUI>/3D/model/train.py \
    --data dataset/dataset.h5 --species C --distance gdf --epochs 200
```

That is the code the published results came from and it carries the three feature
switches. It draws its own test split at random by geometry rather than reading
the one in the file, so pass `--test-frac 0.2` and check which geometries it took.

---

## Step 12. What the final number means

`train_prt3d.py` ends with the held-out test error against two baselines:

```
model                      the number you got
baseline, dataset mean     predict one number everywhere
baseline, per-sample mean  predict each sample's own mean everywhere
```

The **per-sample mean** is the one that matters. Beating the dataset mean only
means the model knows roughly how much product there is; beating the per-sample
mean means it has learned the **shape** of the field, which is the thing the pore
space and the two dimensionless numbers determine and the thing this whole sweep
exists to teach.

- **PASS**, comfortably below both: the chain works. CompLaB3D produced fields
  whose inputs explain their outputs, the collector carried them intact, and the
  network can read them.
- **MARGINAL**: train longer first. If it stays marginal, the sweep may be too
  narrow to be worth generalising over, rather than broken.
- **FAIL**: run `verify_dataset.py` again and read every line of it. A dataset
  that passes step 9 and fails here is a training problem; a dataset that fails
  step 9 was never going to work.

---

## If something goes wrong

| what you see | what it means |
|---|---|
| no `[KIN] abiotic A+B->C` in a run log | the binary was built against the CompLaB defaults, not this package's headers. Redo step 3. Nothing downstream is worth keeping. |
| `Damkohler: not reported` in the startup banner | the patch in step 3 did not go in. Every `<initial_concentration>` is zero, so stock v1.3.2 has no reference composition to build k_r from. |
| `[NS] Pe achieved` far from the target | the flow solve did not converge. Raise `<ns_max_iT1>`. |
| every case looks like pure diffusion | `<delta_P>` is zero or below threshold, so `complab_functions.hh` silently set Pe = 0. The gate in step 5 checks this. |
| `[CLAMP]` reports more than 1 per cent refused | the step is too long for the rate. Rebuild that case with a shorter one. |
| the collector drops runs into `failures.csv` | read the reason column. "more than 5% negative" means the clamp is not in the binary you built. |
| `verify_dataset.py` says Damköhler changes nothing | the same failure as the first row, caught after the fact. |
| `verify_dataset.py` says the pore space changes nothing | each case's `input/geometry.dat` is not its own. Rebuild the campaign. |
| `verify_dataset.py` reports something exceeding the feed | the campaign was built with `--boundaries mixed`, and an open plane is manufacturing mass. Rebuild with `--boundaries dirichlet`, which is the default. `README_1_SCENARIO.md` section 4 has the measurements. |
| `analytic_channel.py` says STOP, boundaries = mixed | build the duct case with `--boundaries dirichlet`, the default; the exact statements need every face held. |
| a geometry does not percolate | `make_geometries.py` rejects those and moves to the next seed; if you changed the porosity targets, the tightest ones may need more attempts. |
| the build cannot find Palabos | pass `-DPALABOS_ROOT=/path/to/palabos-v2.3.0`, or set it in the environment. |
| `selftest.py` fails | the problem is on this machine, not in the case setup. The step it stopped at is named. |
