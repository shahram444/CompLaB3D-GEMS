# 01 - Flow only: steady flow through a pore space, and a conservative tracer carried by it

## 1. The scenario

Water is pushed through a pore space by a pressure difference, and once the flow
has settled down a tracer is released at the inlet and swept through by that
flow. Nothing reacts. Nothing grows. The pore space never changes shape.

If you have run a column experiment to measure hydraulic conductivity and then
injected a conservative tracer to check the flow field with a breakthrough
curve, this is that experiment, done at the pore scale instead of the core
scale. A reactive-transport modeller will recognise the two-stage structure:
solve the flow first, then advect a solute on top of the converged velocity
field.

This is the smoke test of the whole example set. If this case does not build,
run and produce a sane tracer profile, nothing built on top of it is worth
debugging. Read it first.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
flow|>>>>>>>>>>>>>>>>>>  ####  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> |flow
1.0 |>>>>>>>>>>>>>>>>>>  ####  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> |0.0
held|>>>>>>>>>>>>>>>>>>  ########  >>>>>>>>>>>>>>>>>>>>>>>>>> |held
    |                                                         |
    +---------------------------------------------------------+
      pressure drop pushes flow left -> right
      tracer released at x=0, held at 1, carried downstream and held at 0 at x=23

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition, for both
the flow and the tracer. The other four faces are an inert wall, drawn by
`preprocess.py`, because the solver gives those faces nothing at all and a face
with no condition reads back values that nothing updates.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet number of 1 |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one species: `tracer` |
| **Diffusivity** | 5e-10 m2/s, in pore and in biofilm alike |
| **Chemistry** | none. Both `<enable_kinetics>` and `<enable_abiotic_kinetics>` are false |
| **Biology** | none. `<biotic_mode>false</biotic_mode>`, no biomass field allocated |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 400 advection-diffusion steps, output every 100 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **tracer** | Dirichlet, held at 1.0 | Dirichlet, held at 0.0 | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no matter
what the interior does. Held at both ends, the tracer here is fed from a
reservoir at the inlet and drained into one at the outlet, so its total is set
by the boundaries, not by anything the interior does.

## 4. The physics

Flow is a steady Navier-Stokes solve on the D3Q19 lattice Boltzmann scheme
(nineteen discrete velocity directions per voxel), with BGK collision. `<tau>`
sets the kinematic viscosity:

```
    nu     =  (tau - 1/2) / 3          tau = 0.8   ->   nu = 0.1
    omega  =  1 / tau                                   omega = 1.25
```

You do not hand the solver a pressure drop directly. You hand it a target
*Peclet number*, the ratio of how fast the tracer would be carried by advection
to how fast it would spread by diffusion alone over the domain's characteristic
length, and the solver works backwards to the pressure drop that achieves it on
this particular geometry:

```
    1. solve the flow at a seed pressure drop        <delta_P> = 1e-6, a guess
    2. measure the permeability of this pore space    from that solution
    3. compute the pressure drop that hits the        using the measured
       target Peclet on this geometry                 permeability
    4. solve the flow again at that pressure drop
```

The seed value is discarded once it has done its job. Permeability depends on
the pore space you happened to build, so a pressure drop that gives Pe = 1 in
one geometry gives something else in another; asking for the Peclet number
directly is what makes runs on two different geometries comparable. Setting
`<Peclet>0</Peclet>`, as in example 02, skips this whole procedure and the flow
solver with it.

Once the flow has converged, the tracer is advected on it and diffuses with its
own coefficient. It reacts with nothing, so its only sources are the two
boundaries.

## 5. What happens each step

The run has two phases.

1. **Flow, solved to convergence.** The Navier-Stokes lattice is iterated,
   re-checking convergence every `<ns_update_interval>` steps, up to a cap of
   `<ns_max_iT1>` = 20000 iterations against a tolerance of
   `<ns_converge_iT1>` = 1e-6.
2. **Permeability measured, pressure drop corrected, flow re-solved** to the
   tighter cap `<ns_max_iT2>` = 2000 and tolerance `<ns_converge_iT2>` = 1e-4,
   as described above.
3. **The tracer is released** and, for each of the 400 advection-diffusion
   steps: stream and collide on the D3Q7 lattice, advected by the converged
   velocity field; re-apply the two Dirichlet faces.

The stability caveat is on the flow solve, not the tracer: if the first solve
hits the `<ns_max_iT1>` cap without meeting its tolerance, it has not converged,
and the permeability measured from it, the corrected pressure drop, and every
tracer step built on that velocity field all inherit the error. A run that hits
the cap is not a run to trust.

## 6. What comes out

```
output/
  nsLattice_*.vti           the converged velocity field
  tracer_0000100.vti  ...  tracer_0000400.vti      concentration of the tracer
  summary.csv                one row every 40 steps
  run.log                    the whole run, including the [NS] convergence lines
```

There is no `rate_*.vti` field here: nothing reacts, so no rate exists to write.

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Flow convergence | Reached well inside the 20000-iteration cap on the first solve |
| Achieved Peclet | Close to the target of 1, within a fraction of a percent |
| Total tracer | Rises from 0 as the pore space fills, then levels off |
| Tracer profile in x | Falls monotonically from the inlet face to the outlet face |
| Minimum tracer | Zero or above; a conservative species going negative is a real failure |
| Wall clock | Seconds on one core. This is a small case on purpose |

## 7. What to check

1. **The flow converged**, not hit the iteration cap. `run.log` carries the
   `[NS]` lines that say which happened.
2. **The achieved Peclet number is close to the target.** Look for the
   achieved-versus-target line in the log.
3. **No `[NEG!]` warnings.** A tracer that only advects, diffuses and is pinned
   at two Dirichlet faces has no mechanism that can drive it negative; a warning
   here would point at the transport step itself.
4. **The tracer profile is monotone in x**, falling from the inlet toward the
   outlet, in the `tracer_*.vti` files.

`postprocess.py` runs check 3 for you (it flags any negative field) and prints
the `[NS]` lines from the log so you can read checks 1 and 2 off directly. It
also reports that no conservation balance is possible here: the tracer is held
at both ends, so no total in this case is closed, and any conservation check
would fail for reasons that mean nothing. Check 4 is a look at the `.vti`
output.

## 8. What this case demonstrates

**Driving flow to a target Peclet number instead of a fixed pressure drop.**
One number, `<Peclet>`, and the solver finds the pressure drop that achieves it
on whatever pore space it is handed, which is what makes runs on different
geometries comparable at all.

**What it deliberately leaves out:** any solute that reacts, any organism, and
any change to the pore space. Example 02 removes the flow entirely and keeps
only the transport, to isolate what diffusion alone does. Example 03 adds the
first reaction, a chemical rate law with no organism in it.

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
./scripts/setup_case.sh 01_flow_only run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer, and
Palabos v2.3.0. No optional solver is needed. Cases 09, 12, 21 and 22 additionally
need GLPK, and cases 10 and 16 need Python with COBRApy; this one needs neither,
so the plain cmake line is enough:

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
