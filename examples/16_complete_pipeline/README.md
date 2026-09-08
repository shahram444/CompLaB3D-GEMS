# 16 - Complete pipeline: model fetched, network trained, all from one file

## 1. The scenario

Every earlier metabolic case in this repository needs something prepared
before the solver starts: a model exported by a Python script, a surrogate
network fitted by an offline sweep, a header recompiled into the binary.
This case needs none of it. The metabolic model is unpacked from the bundle
in this folder, the exchange reactions (the reactions a metabolic model uses
to represent uptake and excretion across the cell boundary) are matched by
name against that model's own reaction list, the surrogate network, a small
feed-forward network fitted to reproduce flux balance analysis, is trained
during start-up from a sweep of the linear program, and the fit is checked
before the first transport step. Copy the folder, build, run.

The organism is *E. coli* on the `e_coli_core` model, a standard,
small-scale genome model of *Escherichia coli* central metabolism widely
used as a teaching and testing case in the metabolic-modelling community.
Glucose and oxygen diffuse in, biomass and acetate come out, growing as
attached biofilm in a diffusive pore space:

```
    glucose  +  O2   ->   biomass  +  acetate
```

Acetate is the interesting output. It is not a waste product the
configuration declares; it appears because *E. coli* overflows to
fermentation when the carbon supply outruns the oxygen supply, a well known
behaviour called overflow metabolism, and the linear program finds that on
its own. In a pore space with both delivered by diffusion from the same
face, that happens deep in the domain where oxygen has been drawn down but
glucose has not.

A reactive-transport modeller can read this case as the point where the
metabolic coupling stops needing anything hand-prepared: the model, the
exchange mapping and the reduced-order rate law are all resolved at run
time. A computational biologist will recognise `e_coli_core` immediately,
and can watch its known overflow behaviour emerge from a spatial, diffusive
setting rather than a well-mixed flask.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um
   (the pore space itself is 24 x 24 x 6; preprocess.py adds a one voxel
   wall on the y and z faces the solver conditions with nothing, so NY
   and NZ come out two larger than the pore space it declares)

   x=0                                                          x=23
    |                                                             |
glc |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> |
o2  |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> | Dirichlet,
1.0 |>>>  ====                                          >>>>>>>> | held at 0
    |               ####    ####    ####    ####    ####         |
    |               ====                                         |
    +-------------------------------------------------------------+
      grain blocks alternate between the y=0 and y=23 walls every
      four voxels along x, so diffusion has to weave to get through

      #### = solid grain (wall)     ==== = initial Ecoli patch
      porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry boundary conditions. The
other four faces are an inert wall, drawn by `preprocess.py`, because the
solver gives a face with no declared condition nothing at all, and a lattice
Boltzmann field with nothing conditioning it streams off the edge of the
domain and reads back values that nothing updates. This is the same
geometry examples 09 through 12 use.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `glucose`, `o2`, `acetate` |
| **Diffusivity** | 5e-10 m2/s for all three species, in pore and in biofilm alike |
| **Chemistry** | none written by hand. `<enable_kinetics>false</enable_kinetics>`, `<enable_abiotic_kinetics>false</enable_abiotic_kinetics>` |
| **Biology** | `Ecoli`, one organism, `<reaction_type>surrogate</reaction_type>`, trained at start-up, attached biofilm seeded on its own material number |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 400 advection-diffusion steps, output every 200 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **glucose** | Dirichlet, held at 1.0 | Dirichlet, held at 0 | 1.0 |
| **o2** | Dirichlet, held at 1.0 | Dirichlet, held at 0 | 1.0 |
| **acetate** | closed | closed | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does; here both ends are Dirichlet, so glucose and
oxygen are held at a reservoir value on the left and drawn down to zero on
the right, driving a sustained gradient across the whole domain rather than
the reflecting outlet examples 09 to 12 use. *Closed* means no flux crosses
that face at all: acetate is made inside the domain by overflow metabolism
and stays there.

Two settings in `<microbiology>` are not optional for a biofilm organism,
and leaving them out is not harmless: `<thrd_biofilm_fraction>` defaults to
0, which makes the biomass-at-or-above-threshold test true in every pore
voxel, empty ones included, and an empty voxel has no species to label, so
the mask update reaches its should-never-happen branch and prints an error,
once per empty voxel per sweep. `<maximum_biomass_density>` defaults to a
value no biofilm reaches, so the cap never engages and biomass grows without
bound in the best-supplied voxels. This case sets both:
`<maximum_biomass_density>50</maximum_biomass_density>` and
`<thrd_biofilm_fraction>0.01</thrd_biofilm_fraction>`.

## 4. The reaction

There is no model file shipped in `input/` for this case, and no
`<model_filename>` tag either. Instead, `CompLaB.xml` names where the model
comes from and lets the solver resolve it at start-up:

```xml
    <model_source>bigg:e_coli_core</model_source>
    <model_bundle>models</model_bundle>
    <model_cache>input</model_cache>
    <allow_download>false</allow_download>
```

`e_coli_core` ships gzipped in `models/`. At start-up it is unpacked into
`<model_cache>` and checked against `models/manifest.txt`, so a model that
has been revised upstream shows up in the log rather than silently
differing. `<allow_download>false</allow_download>` is the default, and it
is why this case runs on a compute node with no route to the internet.

The exchange reactions are matched by name rather than by position:

```xml
    <exchange_reaction_names>EX_glc__D_e EX_o2_e EX_ac_e</exchange_reaction_names>
```

One reaction name per substrate, in `<name_of_substrates>` order, resolved
against the model's own reaction list at start-up. This is the safe form:
the positional `<exchange_reaction_indices>` that examples 09, 10 and 12 use
cannot be checked, a wrong index produces a well-posed linear program with
the wrong biology, which runs to completion and reports plausible numbers.
A wrong name stops the run and prints the near matches.

The linear program itself, in words: maximise the biomass objective flux,
subject to every internal metabolite of `e_coli_core` being at steady state,
subject to each flux staying within a bound, and with the exchange bounds
for glucose and oxygen built each step from the local concentration through
a Michaelis-Menten term. What is different from example 09 is that no
transport step actually calls this linear program: it is only ever solved
during the one-time training sweep below, and every step of the run itself
calls a network fitted to its answers instead.

```xml
    <surrogate>
        <enabled>true</enabled>
        <weights_file>output/ecoli.srg</weights_file>
        <train_if_missing>true</train_if_missing>
        <train>
            <samples>300</samples>          points swept from the linear program
            <restarts>2</restarts>          independent fits, best kept
            <epochs>400</epochs>
            <verify_points>64</verify_points>   LP re-solves the fit never saw
            <layers>8</layers>
            <input0>glucose 0.1 10.0</input0>
            <input1>o2      0.1 20.0</input1>
        </train>
    </surrogate>
```

No weights file exists the first time, so one is trained: the linear program
is swept over the two ranges above, a network is fitted to what it returned,
and the fit is checked by re-solving the linear program at 64 points it
never saw. The result is written to `<weights_file>`, and every later run of
this folder loads it instead of retraining. Delete `output/ecoli.srg` to
force a retrain.

The training ranges and `<fba_maximum_uptake_flux>` (`10. 20. 0.`) are the
same numbers, and that is not a coincidence: the Michaelis-Menten uptake
bound means no voxel can ever request more than Vmax, so the top of each
training range is Vmax. Raise one without retraining and the run leaves the
region the network was fitted on; it will say so, and say how often, but it
will not stop, and outside that region a network does not fail, it
extrapolates. Example 11 is the longer discussion of that failure mode.

## 5. What happens each step

Once, before the first step:

1. **Resolve the model.** Unpack `models/e_coli_core.xml.gz` into
   `<model_cache>`, check it against the manifest, and match
   `<exchange_reaction_names>` against the model's own reaction list.
2. **Load or train the surrogate.** If `output/ecoli.srg` exists, load it in
   milliseconds. Otherwise sweep the linear program over the two training
   ranges (300 samples, 2 restarts), fit the network, and verify it against
   64 linear-program solves the fit never saw. This training sweep is the
   only place in this case a linear program runs; it happens once and its
   cost does not repeat on later runs.

Then, 400 times:

3. **Stream and collide** each species on its own D3Q7 lattice, advancing
   diffusion by one step, and diffuse the biomass field the same way.
4. **Apply the boundary conditions**, so glucose and o2 are re-pinned at
   both ends.
5. **Build the uptake bounds**, per voxel, per substrate, from the local
   concentration.
6. **Evaluate the network**, in every voxel that holds biomass: no linear
   program is solved during this loop, only the trained surrogate.
7. **Read the growth rate off the network's output**, scale by the local
   biomass and the time step through `<biomass_molar_mass>`, and write the
   increments into the change lattices.
8. **Add the increments**, clamped so a step cannot draw a voxel's
   substrate below what it actually holds.

## 6. What comes out

```
output/
  glucose_0000200.vti  glucose_0000400.vti     concentration of glucose
  o2_*.vti  acetate_*.vti                       the other two species
  Ecoli_*.vti                                   the biomass field
  rate_glucose_*.vti  rate_o2_*.vti  rate_acetate_*.vti   reaction rate, mol/L/s
  ecoli.srg                                     the trained surrogate network, read on later runs
  summary.csv                                   one row every 100 steps
  run.log                                       the whole run, including model resolution and the fit check
```

| Quantity | What it should do |
|---|---|
| Total glucose, total o2 | Both maintained near their Dirichlet ends, drawn down toward the interior, and further reduced near the biomass patch |
| Total acetate | Rises from 0 and keeps rising, since it is closed at both ends |
| Acetate field, spatially | Concentrated where oxygen has been drawn low but glucose has not, the overflow-metabolism signature of section 1 |
| Total biomass | Rises where both substrates are present |
| Fit-verification error (in `run.log`) | Small, reported before the first transport step |
| Out-of-training-box clamp count (in `run.log`) | Zero or small over the run |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Second run's start-up time | Markedly shorter than the first, since it loads `output/ecoli.srg` instead of retraining |

## 7. What to check

1. **The model resolves.** The start-up lines name the model, its reaction
   count, and each exchange reaction matched by name. If any of that is not
   what you meant, stop there.
2. **The fit is verified before the run starts.** The 64 held-out points are
   re-solved with the linear program and compared; a poor agreement reported
   at start-up means the network is not a stand-in for the program and
   nothing after it is worth reading.
3. **The run stays inside the training region**, and the out-of-region
   count stays small.
4. **Acetate appears where oxygen is drawn down.** That is the overflow
   metabolism of section 1, and it is the result that would be hard to get
   from a hand-written rate law.
5. **The second run is fast.** It should load `output/ecoli.srg` rather than
   retrain; if it retrains, the weights file was not written where
   `<weights_file>` says.

`postprocess.py` reads `output/summary.csv` and `run.log`, reports every
field's change, flags any negative value, and prints the log lines about how
the model was resolved and which rate paths ran, covering most of check 1.
It also states plainly that this case has enough moving parts that no single
field tells you whether it worked, and that a check reported SKIPPED is not
a pass. Checks 2 through 5 need a look at the log and the fields, which is
yours to do.

## 8. What this case demonstrates

**A complete, self-contained metabolic pipeline: model resolution, named
exchange matching, and surrogate training, all resolved at start-up from one
configuration file, with nothing prepared offline.** Copy the folder, build,
run.

**What it deliberately leaves out:** flow (`<Peclet>0</Peclet>`), an abiotic
reaction, and any change to the pore space. The surrogate is fitted over two
inputs only, glucose and o2; everything else the model could respond to is
held at whatever the training sweep held it at. A network is a stand-in for
the linear program on the axes it was swept over, and on no others. Example
11 is the same idea with a pre-fitted, compiled-in network instead of one
trained at start-up; example 09 is the plain linear program this surrogate
approximates.

## 9. How to build and run

Two files do this, and they are deliberately separate because compiling and
running belong in different places on a cluster.

| File | What it is | Where you run it |
|---|---|---|
| `COMPILE.txt` | The interactive session that builds `./complab`, step by step | Once, by hand, on an interactive node |
| | `run.sh` | The SLURM batch job: modules, pre-processing, the solver, the checks. Submit with `sbatch run.sh`. |
| | `COMPILE.txt` | The interactive session that builds `./complab`, one step at a time. |

**First time in this case folder,** open `COMPILE.txt` and follow it. It asks
for an interactive node, loads the modules, runs cmake with the flags this case
actually needs, and compiles. It also says exactly when you have to come back
and recompile, which for most cases is almost never.

**Every run after that** is one command:

```bash
sbatch run.sh
squeue -u $USER                 # watch it
tail -f output/run.log          # read it while it runs
```

`run.sh` carries a comment on every line: the SLURM header, the module loads
that must match what you compiled with, the pre-processing, the solver call and
the post-processing checks. It runs on one rank on purpose, and the note at the
bottom of the file says why.

To work in a scratch copy instead of dirtying this folder:

```bash
./scripts/setup_case.sh 16_complete_pipeline run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. `CompLaB.xml` declares
`<enable_fba_glpk>false</enable_fba_glpk>` and
`<enable_fba_cobrapy>false</enable_fba_cobrapy>`; only
`<enable_surrogate>true</enable_surrogate>` is set, and the surrogate solver
needs neither GLPK nor Python and is always compiled in, so no optional flag
is needed:

```bash
cmake -B build -S . \
      -DCMAKE_BUILD_TYPE=Release \
      -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** Nothing in this case's own configuration
requires a rebuild to change: the model source, the exchange-reaction names
and the whole surrogate training block live in `CompLaB.xml` and are read at
start-up, so editing any of them, or deleting `output/ecoli.srg` to force a
retrain, takes effect on the next run with no compile step. A rebuild is
needed only for a change to `CMakeLists.txt` or the solver sources, and it
would also be needed if this case's `<reaction_type>` were switched to a
GLPK or COBRApy path, which would then need the matching `-DENABLE_GLPK=ON`
or `-DENABLE_COBRAPY=ON`, exactly as in examples 09 and 10.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Optional. The model is unpacked and the network trained by the solver itself at start-up; this script exists for inspecting the bundled model by hand and for pinning a reproducible copy of it. |
| **run** | `CompLaB.xml` | What the solver reads, including the model source, the named exchange reactions and the surrogate training block. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| | `models/e_coli_core.xml.gz` | The metabolic model, unpacked into `input/` at start-up and checked against `models/manifest.txt`. |
| | `output/ecoli.srg` | The trained surrogate network. Written on the first run, loaded on every later one. Delete it to force a retrain. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
