# 02 - Diffusion only: a profile you can check by hand

## 1. The scenario

This case takes example 01's tracer and switches the flow solver off entirely.
Nothing is pushed through the pore space. The tracer is held at 1 on the inlet
face and 0 on the outlet face, and spreads by diffusion alone, reacting with
nothing.

This is the diffusion-cell limit of example 01: no pressure difference applied,
just a fixed concentration difference across a static medium, the setup behind
a Fick's-law diffusion coefficient measurement. A reactive-transport modeller
will recognise it as the textbook check every transport solver should pass
before it is trusted with anything harder.

That makes this the one case in the whole set whose correct answer can be
written down without running it. Diffusion with no reaction and no flow,
between two fixed values, has an analytic steady state: a straight line.

```
    d2C
    ---  =  0        ->      C(x)  =  1 - x / L
    dx2
```

Curvature in the converged profile means something is firing that should not
be: a reaction, a boundary that is not actually being held, or a transport
step that is not conservative. There is nothing else happening here to hide
it.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
    |....................  ####  ..........................  |
1.0 |....................  ####  ..........................  |0.0
held|....................  ########  ......................  |held
    |                                                         |
    +---------------------------------------------------------+
      no flow: Peclet = 0, the Navier-Stokes solve never runs
      tracer diffuses from x=0 toward x=23, dots = diffusion, not advection

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition. The
other four faces are an inert wall, drawn by `preprocess.py`, because the
solver gives those faces nothing at all and a face with no condition reads
back values that nothing updates.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so the Navier-Stokes solve never runs |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one species: `tracer` |
| **Diffusivity** | 5e-10 m2/s, in pore and in biofilm alike |
| **Chemistry** | none. Both `<enable_kinetics>` and `<enable_abiotic_kinetics>` are false |
| **Biology** | none. `<biotic_mode>false</biotic_mode>`, no biomass field allocated |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 2000 advection-diffusion steps, output every 200 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **tracer** | Dirichlet, held at 1.0 | Dirichlet, held at 0.0 | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no matter
what the interior does. With both ends held, this is a fixed head-difference
diffusion cell: the tracer is fed at the inlet and drained at the outlet, so
its total is set by the boundaries, not by the run.

## 4. The physics

No flow, no reaction: transport is the diffusion equation alone, solved on the
D3Q7 advection-diffusion lattice (seven discrete directions per voxel: the six
face neighbours plus staying put). With a fixed value at each end and nothing
else acting, the steady state is the straight line given in section 1. Nothing
in this case's chemistry is compiled in: no `kinetics/` header ships with it,
because there is no reaction for one to define.

## 5. What happens each step

The solver repeats this for each of the 2000 steps:

1. **Stream and collide** the tracer on its D3Q7 lattice, which advances
   diffusion by one step. There is no velocity field to advect it with.
2. **Apply the boundary conditions**, re-pinning the tracer to 1.0 at the left
   face and 0.0 at the right.

There is a stability caveat worth knowing, though it is not about
instability in the usual sense: `<ade_converge_iT>` is read from
`CompLaB.xml` but never checked against, so the advection-diffusion loop
always runs the full `<ade_max_iT>` = 2000 steps and never stops early even
once the profile has settled. Steadiness is something you confirm by looking
at the record, not something the solver reports; a total still moving at the
last row of `summary.csv` means the run needs to be longer, not that anything
is wrong.

## 6. What comes out

```
output/
  tracer_0000200.vti  ...  tracer_0002000.vti      concentration of the tracer
  summary.csv                                       one row every 200 steps
  run.log                                            the whole run
```

There is no `rate_*.vti` field here: nothing reacts, so no rate exists to
write. There is no `nsLattice_*.vti` either: the flow solver never runs.

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Total tracer | Rises from 0 toward a steady value as the profile fills in, then flattens |
| Profile in x, once settled | A straight line from 1 at the inlet to 0 at the outlet |
| Profile in y, near a grain block | Varies: the staggered grain blocks force diffusion to go around them |
| Profile in z | Uniform. The geometry does not vary in z inside the walls |
| Minimum tracer | Zero or above; a conservative species going negative is a real failure |
| Wall clock | Seconds on one core, and noticeably faster than example 01 with the flow solver off |

## 7. What to check

1. **The profile is straight** in x between the two boundary values, once the
   run has settled. This is the whole point of the case.
2. **The tracer total has stopped changing** between the last rows of
   `output/summary.csv`. Because the loop never stops early (see section 5),
   steadiness has to be read off the record rather than trusted from a
   convergence flag.
3. **No `[NEG!]` warnings.** A conservative tracer between 0 and 1 cannot
   leave that range.
4. **No dependence on z.** The geometry is uniform in z inside its walls, so
   the profile must be too; a gradient across z would mean the boundary
   treatment differs between the two z faces. There is dependence on y near a
   grain block, and that is expected: `preprocess.py` builds a staggered slot
   whose grain blocks alternate between the two y walls every four voxels
   along x, so diffusion has to go around them, and the concentration
   necessarily varies with y close to a block. That is a property of the
   geometry, not a sign of a boundary problem.

`postprocess.py` runs check 3 for you (it flags any negative field) and prints
the flow-solver log lines so you can confirm they show `Peclet = 0` was
honoured and the solve was skipped, as check 4's counterpart. It also reports
that no conservation balance is possible here: the tracer is held at both
ends, so no total in this case is closed. Checks 1, 2 and 4 are read from the
`.vti` output and `summary.csv`.

## 8. What this case demonstrates

**A transport-only case with a known analytic answer.** Because the steady
state is a straight line derivable from the diffusion equation alone, this is
the case to trust the solver against before adding anything that has no
independent answer to check it with.

**What it deliberately leaves out:** flow, reaction, biology, and any change
to the pore space. Example 01 is the same tracer with the flow switched back
on; example 03 adds the first reaction, a chemical rate law with no organism
in it.

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
./scripts/setup_case.sh 02_diffusion_only run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. No optional solver is needed. Cases 09, 12, 21 and 22
additionally need GLPK, and cases 10 and 16 need Python with COBRApy; this one
needs neither, so the plain cmake line is enough:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** This case has no chemistry of its own to edit:
both kinetics paths are switched off in `CompLaB.xml`, and no `kinetics/`
header ships in this folder (`setup_case.sh` still lays down the shared
defaults so the executable builds, but neither is ever called here). A
rebuild is only needed if you change the solver source itself. Editing
`CompLaB.xml` or `input/geometry.dat` never needs one.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
