# 08 - Two microbes: two populations competing for one substrate

## 1. The scenario

Two populations share one pore space and one donor, and they are not the
same organism. `Bug1_fast` grows quickly but needs a lot of donor to do it.
`Bug2_efficient` grows more slowly but keeps going at donor concentrations
where the first has stalled. Each is seeded in its own patches, and from
there they compete for the same shared substrate field.

A reactive-transport modeller will recognise this as the classic r-strategist
versus K-strategist trade-off: a fast, wasteful grower that wins where
resources are plentiful, against a slow, efficient one that wins where they
are scarce. Which one dominates where depends on where each sits relative to
the donor source, not just on the parameters in its own rate law.

This is also the case that proves the per-microbe settings in this solver are
genuinely independent. The two organisms differ in three ways at once:
maximum growth rate, half-saturation constant, and which biomass solver
moves them, `CA` for one and `LBM` for the other, running side by side in the
same simulation.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
    |[1]..[1]..[1]..[1]..[1]..[1]#### ......................  |
donor|..[2]..[2]..[2]..[2]..[2]..[2]#....................... |zero
0.5 |....................  ########  ......................  |gradient
held|                                                         |
    +---------------------------------------------------------+
      flow: Peclet = 0.6, pushes donor and product left -> right
      [1] = Bug1_fast, against the y=0 wall, at every x % 4 == 0
      [2] = Bug2_efficient, against the y=ny-1 wall, at every x % 4 == 3
      both interleaved along the whole length, 108 voxels each

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Both populations are seeded in patches spread along the entire length of
the domain, one against each of the two staggered walls, interleaved with
each other every four voxels along x. The single patch closest to the donor
source at `x = 0` belongs to `Bug1_fast`, which is why it has the initial
advantage the XML comments call out. Only the two end faces `x = 0` and
`x = 23` carry a boundary condition for the solutes; the other four faces
are an inert wall, drawn by `preprocess.py`.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet number of 0.6 |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, two species: `donor`, `product` |
| **Diffusivity** | 5e-10 m2/s for both species, in pore and in biofilm alike |
| **Chemistry** | `kinetics/defineKinetics.hh`, compiled into the executable, byte-identical to examples 05 to 07 |
| **Biology** | `Bug1_fast` on material 3, cellular automaton; `Bug2_efficient` on material 4, lattice Boltzmann at 1e-10 m2/s |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 1000 advection-diffusion steps, output every 200 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **donor** | Dirichlet, held at 0.5 | Neumann, zero gradient | 0.1 |
| **product** | closed | closed | 0 |
| **Bug1_fast** / **Bug2_efficient** (biomass) | closed | closed | 1.0 in each population's own patch, 0 elsewhere |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, a reservoir that keeps both colonies fed.
*Neumann with a zero gradient* means nothing flows across that face, an open
outlet the donor never reaches. *Closed* means the same thing on both ends:
the product is made inside the domain and stays there, and neither
population's biomass can cross a domain face.

Unlike the earlier biotic cases, this one runs with flow on, at a target
Peclet number of 0.6 rather than 1. Raising it to 1 pushes the grid Peclet
number, the Peclet number evaluated on a single lattice cell rather than
over the whole domain, to 2.6, above the value of 2 at which advection
overshoots what diffusion can smooth out within one cell; the case's own
`CompLaB.xml` records that such a run produced product concentrations of
-4.3e-4 against a peak of 4.5e-3, flagged `[NEG!]` at every diagnostic
interval. At 0.6 the grid Peclet number is 1.6 and the Mach number, the
ratio of the flow speed to the lattice's own speed of sound, is 0.20, both
inside the limits the solver checks at start-up.

## 4. The reaction

The same `kinetics/defineKinetics.hh` as examples 05 to 07, written for any
number of populations: it loops over every organism and takes each one's
parameters from a table indexed in `CompLaB.xml` order.

```
                                 [S]
    growth[m]  =  mu_max[m] * ───────────── * B[m]           gDW L-1 s-1
                               Ks[m] + [S]

    donor consumed   =   SUM over m of  growth[m] / Y[m]
    product made     =   SUM over m of  growth[m] * fP[m]
```

| symbol | `Bug1_fast` | `Bug2_efficient` | unit | meaning |
|---|---|---|---|---|
| mu_max | 2.0e-4 | 1.0e-4 | 1/s | maximum specific growth rate |
| Ks | 5.0e-2 | 5.0e-3 | mol/L | half-saturation constant: the donor concentration at which growth runs at half its maximum rate |
| Y | 0.4 | 0.4 | gDW/mol | yield: biomass made per mole of donor consumed |
| fP | 2.0 | 2.0 | mol/gDW | product released per unit of biomass made |
| kd | 1.0e-6 | 1.0e-6 | 1/s | first-order decay rate |

`Bug1_fast` grows twice as fast but needs ten times more donor to reach half
that rate; `Bug2_efficient` is the opposite trade. Both organisms draw from
the same `donor` field, and their two contributions are summed into one
donor increment before it is applied, so competition is resolved on the
shared budget rather than by whichever organism the solver happens to
evaluate first. As with examples 05 to 07, this compiled kinetics path
applies the summed increment as computed, with no clamp against what the
voxel actually holds; if the two organisms' combined draw exceeds the local
donor, the voxel goes negative and the solver reports `[NEG!]`.

Only `<solver_type>` and `<reaction_type>` are declared per microbe in
`CompLaB.xml`; there is no `<mu_max>` tag anywhere in the repository, and
both `mu_max` and this rate law's `Ks` come from the per-microbe table in
`kinetics/defineKinetics.hh`. `<half_saturation_constants>` is present per
microbe in the XML and carries the same two numbers, but it feeds the
flux-balance and surrogate rate paths, which build their own uptake bounds
from it, and a hand-written rate law like this one does not read it.

## 5. What happens each step

The solver repeats this 1000 times:

1. **Solve the flow** to the target Peclet number, as in example 01, once
   before the reaction loop begins.
2. **Stream and collide** `donor` and `product` on their D3Q7 lattices,
   advected by the converged flow and diffusing with their own coefficient.
3. **Apply the boundary conditions**, re-pinning `donor` to 0.5 at the left
   face.
4. **Evaluate the rate law** in every open voxel that holds either
   population's biomass: read the local donor concentration, compute each
   organism's `growth`, and sum both contributions into one donor increment
   and one product increment.
5. **Apply each population's own biomass increment**, then move biomass by
   whichever solver that population uses: the cellular-automaton rule for
   `Bug1_fast`, streaming on its own D3Q7 lattice for `Bug2_efficient`.

The stability caveat is the Peclet number described in section 3: too high
for this geometry and advection overshoots what a single lattice cell's
diffusion can smooth, driving `product` negative at every interval. The
donor side carries the usual compiled-kinetics caveat as well, no positivity
clamp on the summed reaction increment.

## 6. What comes out

```
output/
  donor_0000200.vti  ...  donor_0001000.vti        concentration of the donor
  product_*.vti                                     the product plume
  Bug1_fast_*.vti  Bug2_efficient_*.vti              biomass density, one field per population
  rate_donor_*.vti  rate_product_*.vti               the reaction rate as a field, mol/L/s
  nsLattice_*.vti                                    the converged flow field
  summary.csv                                        one row every 100 steps, one column per population
  run.log                                             the whole run, including the checks
```

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Both biomass totals | Rise from their seeded patches as each population grows |
| `Bug1_fast` relative to `Bug2_efficient` | Ahead nearer the inlet, where donor is plentiful; the gap should narrow or reverse further from the inlet |
| Total donor | Falls where either colony is active, replenished from the left boundary |
| Total product | Rises from 0, since it is closed and made only where biomass is |
| Minimum donor | Zero or above, unless the Peclet number has been raised past the limit in section 3 |
| Achieved Peclet | Close to the target of 0.6 |
| Wall clock | Seconds on one core, more than examples 05 to 07 since the flow is now solved too |

## 7. What to check

1. **Which population wins, and where.** The fast one should dominate near
   the inlet where donor is plentiful; the efficient one should hold on
   further in.
2. **The donor is drawn down by both populations**, and never goes negative
   at this Peclet number.
3. **The two biomass totals separately.** `output/summary.csv` carries a
   column per population, so the competition can be read as two curves
   rather than one sum.
4. **No `[NEG!]` warnings.** One at every interval means the Peclet number
   is too high for this geometry, as described in section 3; check
   `run.log` for a `[STABILITY]` block before believing any output at a
   Peclet number higher than the one shipped here.

`postprocess.py` runs check 2 for you (it flags any negative field) and
reports the change in `product`. Checks 1 and 3 are read from the two
biomass columns in `summary.csv`, and check 4 from `run.log`.

## 8. What this case demonstrates

**Per-microbe settings that are genuinely independent, not global.**
`<solver_type>` and `<reaction_type>` are declared once per organism, and a
run may mix a cellular automaton with a lattice-Boltzmann biomass field, or
flux-balance analysis with a hand-written rate law, in the same step; this
case exercises the first of those combinations. It also demonstrates two
populations competing for one shared, clamped resource, where the outcome
depends on both physiology and position, not on physiology alone.

**What it deliberately leaves out:** any change to the pore space, and any
rate path besides compiled kinetics for either organism. Example 12 mixes a
kinetics-driven population with a flux-balance one in a single run; examples
05 to 07 are the one-organism version of this case's chemistry, each with a
different biomass solver.

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
./scripts/setup_case.sh 08_two_microbes run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. No optional solver is needed. Cases 09, 12, 21 and 22
additionally need GLPK, and cases 10 and 16 need Python with COBRApy; this
one needs neither, so the plain cmake line is enough:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** The chemistry in this case is a C++ header
compiled into the executable, not data read at start-up. So editing
`defineKinetics.hh`, including any organism's `mu_max`, `Ks`, `Y`, `fP` or
`kd`, means rebuilding. Editing `CompLaB.xml`, including `<Peclet>`,
`<solver_type>` or `<biomass_diffusion_coefficients>`, or
`input/geometry.dat`, does not.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and both seeded biomass patches, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **build** | `kinetics/defineKinetics.hh` | This case's chemistry, compiled in. Byte-identical to examples 05, 06 and 07. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space and both patches. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
