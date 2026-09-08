# 06 - Biotic finite difference: the same organism, biomass that spreads by diffusion

## 1. The scenario

This is example 05 with one line changed: `<solver_type>FD</solver_type>`
instead of `CA`. The chemistry is byte-identical, the same rate law, the same
parameters, the same geometry, the same seeded patch, so everything that
differs between the two runs is the biomass solver and nothing else.

Instead of holding biomass in place until a voxel fills and then spilling the
excess, `FD` (finite difference) spreads it continuously by diffusion on the
same grid, with its own diffusion coefficient. A reactive-transport modeller
will recognise this as treating biomass the way a mobile, loosely attached
population is often modelled: as a diffusing field rather than a solid that
only advances by displacement.

The two solvers are meant to be compared directly. Both grow the colony by
the same amount, because the growth law is the one thing that has not
changed; what differs is only where that biomass ends up, tightly held in the
seeded patch under `CA`, or smoothed out over a much larger footprint under
`FD`. That is also a warning about how to read the output: a spreading colony
has a falling peak density even while its total keeps rising, so a check that
only watches the peak would call this population shrinking when it is not.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
    |....................  ####  ..........................  |
donor|...................  ####  ..........................  |zero
0.5 |@@@@@@@@@@~~~~~~~~~..  ########  ......................  |gradient
held|@@@@@@@@@@~~~~~~~~~..                                    |
    |                                                         |
    +---------------------------------------------------------+
      no flow: Peclet = 0. donor diffuses in from x=0
      @@@@@ = the seeded patch, 108 of the 3168 open voxels
      ~~~ = biomass spreading beyond the patch by its own diffusion

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition for
the solutes. The other four faces are an inert wall, drawn by
`preprocess.py`. The seeded biomass patch is the same as example 05: material
number 3, at every `x` where `x % 4 == 0`, in the three voxel rows next to the
wall block there, through the full depth in z, 108 voxels in total, carved
out of open pore rather than added to it.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, two species: `donor`, `product` |
| **Diffusivity** | 5e-10 m2/s for both species, in pore and in biofilm alike |
| **Chemistry** | `kinetics/defineKinetics.hh`, compiled into the executable, byte-identical to example 05 |
| **Biology** | `Bug`, one population, seeded on material 3 at initial density 1.0, biomass spread by finite difference at 3e-13 m2/s |
| **Geometry change** | none. The pore space is fixed for the whole run; biomass diffuses into voxels but does not remove or add solid |
| **Run length** | 1000 advection-diffusion steps, output every 200 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **donor** | Dirichlet, held at 0.5 | Neumann, zero gradient | 0.1 |
| **product** | closed | closed | 0 |
| **Bug** (biomass) | closed | closed | 1.0 in the seeded patch, 0 elsewhere |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, a reservoir that keeps the colony fed.
*Neumann with a zero gradient* means nothing flows across that face, an open
outlet the donor never reaches. *Closed* means the same thing on both ends:
the product is made inside the domain and stays there, and biomass, which
this solver diffuses only within the domain, never crosses a face at all.

## 4. The reaction

Identical to example 05, one organism, one donor, Monod-limited growth,
first-order decay:

```
                              [S]
    growth  =  mu_max  *  ───────────  *  B                gDW L-1 s-1
                           Ks + [S]

    donor consumed   =   growth / Y            product made  =  growth * fP
    decay             =   kd * B                (first order, always on)
```

| symbol | value | unit | meaning |
|---|---|---|---|
| mu_max | 2.0e-4 | 1/s | maximum specific growth rate |
| Ks | 5.0e-2 | mol/L | half-saturation constant on the donor: the concentration at which growth runs at half its maximum rate |
| Y | 0.4 | gDW/mol | yield: biomass made per mole of donor consumed |
| fP | 2.0 | mol/gDW | product released per unit of biomass made |
| kd | 1.0e-6 | 1/s | first-order decay rate |

The rate law is byte-identical to `kinetics/defineKinetics.hh` in example 05,
so any difference between the two runs is the biomass solver, not the
chemistry. See example 05's README for the code excerpt.

`FD` is the only one of the three biomass solvers that requires
`<biomass_diffusion_coefficients>`; here it is set to 3e-13 m2/s, in pore and
in biofilm alike, roughly a thousand times slower than the solute
diffusivities above, since biomass diffuses far more slowly than a small
dissolved molecule. `FD` also has no relaxation time to satisfy, unlike the
`LBM` solver in example 07, so this diffusivity carries without the numerical
constraint that path has.

## 5. What happens each step

The solver repeats this 1000 times:

1. **Stream and collide** `donor` and `product` on their D3Q7 lattices, which
   advances diffusion by one step. There is no velocity field to advect them
   with.
2. **Apply the boundary conditions**, re-pinning `donor` to 0.5 at the left
   face.
3. **Evaluate the rate law** in every open voxel that holds biomass: read the
   local donor concentration and biomass density, compute `growth`, and write
   the donor, product and biomass increments.
4. **Apply the biomass increment** in place, `B += (growth - kd*B) * dt`.
5. **Diffuse the biomass field** on the same grid, one finite-difference step
   at `<biomass_diffusion_coefficients>`, spreading it into neighbouring
   voxels continuously rather than only once a voxel is full.

The stability caveat is on step 3, as in example 05: the compiled kinetics
path applies each reaction increment as computed, with no clamp against what
a voxel actually holds, so a time step or a rate constant too large for the
local donor will drive it negative and produce a `[NEG!]` warning. This FD
biomass path is also the one path in this repository not exercised by any
other shipped example, so it is worth checking a small run of your own
parameters before relying on it for anything larger.

## 6. What comes out

```
output/
  donor_0000200.vti  ...  donor_0001000.vti     concentration of the donor
  product_*.vti                                  the product plume
  Bug_*.vti                                       biomass density, spreading beyond the seeded patch
  rate_donor_*.vti  rate_product_*.vti            the reaction rate as a field, mol/L/s
  summary.csv                                     one row every 100 steps
  run.log                                         the whole run, including the checks
```

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Total biomass | Rises by close to the same amount as example 05, since the growth law is unchanged |
| Occupied voxel count | Grows well beyond the 108 seeded voxels, as biomass diffuses outward |
| Peak biomass density | Falls relative to example 05, because the same total is spread over more voxels |
| Total donor | Falls near the colony, replenished from the left boundary |
| Total product | Rises from 0, since it is closed and made only where biomass is |
| Minimum donor | Zero or above. A negative value is a real failure |
| Wall clock | Seconds on one core, comparable to example 05 |

## 7. What to check

1. **The total biomass rises**, this is the growth. The peak falls, and that
   is the biomass spreading, not dying.
2. **It spreads from the seeded patch outward** smoothly, with no jumps,
   since diffusion has no threshold the way the automaton does.
3. **The donor is drawn down where the biomass is**, and replenished from the
   left boundary.
4. **No `[NEG!]` warnings.**
5. **Run 05, 06 and 07 and compare the totals**, not just this one run in
   isolation. 05 and 06 should agree closely, since the chemistry is the same
   file and neither solver streams biomass with the flow; 07 differs more,
   because it does.

`postprocess.py` runs checks 1, 3 and 4 for you: it reports the change in
`product` and whether it is monotone, and it flags any negative field. Checks
2 and 5 are read from the `.vti` output and the biomass total in
`summary.csv` across the three runs.

## 8. What this case demonstrates

**The finite-difference biomass solver**, one of three interchangeable ways
this solver can move biomass. Where the cellular automaton in example 05
holds biomass rigid until a voxel is full, `FD` treats it as a continuously
diffusing field from the first step, which is the right choice when the
population is motile or loosely enough attached that growth redistributes it
smoothly rather than in jumps.

**What it deliberately leaves out:** any change to the pore space, any flow,
and any coupling of biomass movement to the flow field, which is what the
`LBM` solver in example 07 adds. Example 05 is this same case with the
cellular automaton instead; example 08 adds a second population competing
for the same donor.

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
./scripts/setup_case.sh 06_biotic_finite_difference run/mycase
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
`defineKinetics.hh`, including changing `mu_max`, `Ks`, `Y`, `fP` or `kd`,
means rebuilding. Editing `CompLaB.xml`, including `<solver_type>` or
`<biomass_diffusion_coefficients>`, or `input/geometry.dat`, does not.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and the seeded biomass patch, the same as example 05, so the two solvers are comparable. Reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **build** | `kinetics/defineKinetics.hh` | This case's chemistry, compiled in. Byte-identical to examples 05, 07 and 08. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space and the patch. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
