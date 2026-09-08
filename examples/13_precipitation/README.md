# 13 - Precipitation: a mineral that grows, fills its voxels, and blocks the flow

## 1. The scenario

Every case before this one runs on a pore space that never changes. This one
does not. Two dissolved species enter from opposite ends, meet in the middle,
and react to form a mineral. The mineral accumulates in the voxel where it
formed, and once a voxel holds enough of it, that voxel is converted to
**solid**. The flow field is then re-solved around the new geometry, so the
precipitate changes the transport that produced it.

The reaction is iron monosulphide precipitation, the mineral that forms
wherever dissolved iron meets sulphide in anoxic sediment:

```
    Fe(2+)  +  HS(-)   ->   FeS(s)  +  H(+)
```

A reactive-transport modeller will recognise this as pore clogging by mineral
scale: the textbook reason an injection well loses permeability over time.
Iron enters from the left face, sulphide from the right. Neither is present
at the start. They diffuse and advect toward each other and the mineral forms
in a band where they overlap, not at either boundary and not uniformly, but
wherever transport happens to bring them together. Nothing in the
configuration says where that band should sit.

The feedback is the point of the case. Mineral fills a voxel, the voxel
seals, the flow finds another path, the reactants meet somewhere else.
Porosity, the fraction of the domain that is open pore rather than solid,
falls as the run proceeds, and `output/summary.csv` records it every
interval.

The proton on the right of the reaction above is not carried as a substrate
in this case: the three transported fields are `Fe2`, `HS` and `FeS`, and the
rate law in `kinetics/defineAbioticKinetics.hh` writes only those three. Do
not go looking for `H` in the output.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
Fe2+|>>>>>>>>>>>>>>>>>>>>  ####  <<<<<<<<<<<<<<<<<<<<<<<<<<<<< |HS-
5.0 |>>>>>>>>>>>>>>>>>>>>  ####::.....................<<<<<<<< |5.0
held|>>>>>>>>>>>>>>>>>>>>  ########  <<<<<<<<<<<<<<<<<<<<<<<<< |held
    |                                                         |
    +---------------------------------------------------------+
      flow: Peclet = 1, pushes both species left -> right
      #### = solid grain block   :: = one-voxel reactive mineral face
      FeS forms and seals a voxel where the two fronts overlap

      grain blocks alternate between the y=0 and y=25 walls every 4 voxels
      along x; the reactive face is the one-voxel layer of each block that
      touches open pore, 72 voxels in total

      porosity 0.6346, 3168 open voxels, 72 declared reactive-phase voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition. The
other four faces are an inert wall, drawn by `preprocess.py`, because the
solver gives those faces nothing at all and a face with no condition reads
back values that nothing updates. The reactive layer itself is one voxel
thick, on the pore-facing side of each grain block only: mineral buried
inside a grain can never be reached by the water, so tracking it there would
cost solve time and change nothing.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346, plus 72 voxels of a declared reactive phase |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet number of 1 |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `Fe2`, `HS`, `FeS`. `FeS` is immobile: it never streams, only accumulates its reaction term |
| **Diffusivity** | 5e-10 m2/s for `Fe2` and `HS`, in pore and in biofilm alike; 0 for `FeS` |
| **Chemistry** | `kinetics/defineAbioticKinetics.hh`, compiled into the executable |
| **Biology** | none. `<biotic_mode>false</biotic_mode>`, no biomass field is allocated |
| **Geometry change** | precipitation: a mineral fills open voxels and converts them to solid |
| **Run length** | 2000 advection-diffusion steps, output every 200 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **Fe2** | Dirichlet, held at 5. | Neumann, zero gradient | 0 |
| **HS** | Neumann, zero gradient | Dirichlet, held at 5. | 0 |
| **FeS** | closed | closed | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, which is a reservoir. *Neumann with a zero
gradient* means nothing flows across that face, an open outlet the species
never reaches. *Closed* means the same thing on both ends: FeS is made inside
the domain and stays there, which is also what lets its inventory be tracked
against `max_precipRho` below.

## 4. The reaction

Second order, one term, in `kinetics/defineAbioticKinetics.hh`:

```
    R  =  k_FeS * [Fe2+] * [HS-]                  mol L-1 s-1

    k_FeS = 1.0    L mol-1 s-1

    d[Fe2+]/dt = -R          d[HS-]/dt = -R          d[FeS]/dt = +R
```

This is a **kinetic** precipitation law, not an equilibrium one: the rate is
proportional to the product of the two concentrations at all times, with no
solubility product and no supersaturation term. A real mineral's rate would
carry a **saturation index** Omega, the ratio of the ion activity product to
the solubility product, with Omega greater than 1 meaning the solution is
**supersaturated** (more dissolved mineral than equilibrium allows, which is
what drives precipitation) and Omega less than 1 meaning it is undersaturated.
`R = k A (Omega - 1)^n`, with `A` the reactive surface area, is the usual
form; it would go in the same function, in the same place. This case skips it
because the point here is the geometry feedback, not the thermodynamics.

The rate law itself is fifteen lines of C++:

```cpp
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    if (C.size() < 3 || subsR.size() < 3) return;   // a short substrate list
    if (mask < 2) return;                           // solid or wall voxel

    const double Fe2 = std::max(C[0], 0.0);         // floor tiny negatives
    const double HS  = std::max(C[1], 0.0);
    const double R = ExamplePrecip::k_FeS * Fe2 * HS;   // mol/L/s

    subsR[0] = -R;
    subsR[1] = -R;
    subsR[2] = +R;                                  // the mineral, immobile
```

The guards at the top are the part to copy when writing a new rate law: the
first stops it reading past the end of a shorter-than-expected concentration
vector, the second keeps chemistry out of solid grains and walls. The floor
at zero matters because the transport solver can leave a value fractionally
below zero, and a rate law that multiplies two such values turns a rounding
artefact into a reaction.

`FeS` is declared `<immobile>true</immobile>` with a diffusion coefficient of
zero, so its lattice never streams: it accumulates exactly in the voxel where
it formed and stays there. A mineral that diffused would smear the band away
before it ever reached a sealing density.

## 5. What happens each step

The solver repeats this 2000 times:

1. **Stream and collide** `Fe2` and `HS` on their D3Q7 lattices, advected by
   the current velocity field and diffusing at the same time.
2. **Apply the boundary conditions**, re-pinning `Fe2` to 5. at the left face
   and `HS` to 5. at the right.
3. **Evaluate the rate law** in every open voxel: read `Fe2` and `HS` there,
   compute `R = k * Fe2 * HS`, and write the three increments into the change
   lattices.
4. **Add the increments** to the concentrations, `C += R * dt`.
5. **Every `<update_interval>` (100) steps, check the geometry.** Any voxel
   whose `FeS` inventory has reached `<max_precipRho>` (48.9 mol/L, see
   below) is converted to solid; `<surface_only>true</surface_only>` means
   only a voxel already touching an existing solid surface is eligible, so
   the mineral cannot **nucleate**, form a first solid seed, out in the open
   water. If any voxel converted, the flow field is re-solved on the new
   geometry.

Step 4 is an explicit update with no clamp against what is actually available
in the voxel. That is the stability limit of this case: if `k * [Fe2] * dt`
approaches 1, a voxel can be driven negative. The solver reports that as a
`[NEG!]` warning, and the fix is a smaller time step or a smaller `k`.

The **48.9 mol/L** threshold itself comes from the mineral, not from the
model. A voxel is completely full of solid mineral when its FeS density
equals the mineral's own molar density, which is 1000 divided by the
**molar volume**, the volume one mole of the solid occupies:

```
    max_precipRho  =  1000 / Vmolar          mol L-1, Vmolar in cm3 mol-1

    mackinawite    87.91 g mol-1
                    4.30 g cm-3     ->   20.44 cm3 mol-1   ->   48.9 mol L-1
```

If you change minerals, recompute `<max_precipRho>` the same way.

| tag | value | what it controls |
|---|---|---|
| `<solid_substrate>` | `2` | which substrate index is the mineral (FeS) |
| `<max_precipRho>` | `48.9` | mol/L in a completely full voxel |
| `<surface_only>` | `1` | heterogeneous nucleation: growth only on wall-adjacent voxels |
| `<perm_ratio>` | `0` | a filled voxel becomes an impermeable wall, not a low-permeability one |
| `<update_interval>` | `100` | steps between geometry re-checks and flow re-solves |

`<update_interval>` is a cost trade: the geometry check is cheap, the flow
re-solve that follows a change is not. A value of 100 means the flow field
can lag the geometry by up to 100 steps, which is acceptable while sealing is
slow and stops being acceptable once it accelerates.

## 6. What comes out

```
output/
  Fe2_0000200.vti  ...  Fe2_0002000.vti     concentration of Fe2+
  HS_*.vti                                    concentration of HS-
  FeS_*.vti                                    the mineral inventory, mol/L per voxel
  rate_Fe2_*.vti  rate_HS_*.vti  rate_FeS_*.vti   the reaction rate as a field, mol/L/s
  nsLattice_*.vti                              the velocity field, re-solved whenever the geometry changes
  summary.csv                                  one row every 100 steps, including porosity
  run.log                                      the whole run, including the checks
```

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Porosity | Falls, in steps of at most `<update_interval>` iterations, from 0.6346 |
| Total Fe2, total HS | Both rise from 0, because both are fed from a reservoir |
| Total FeS | Rises from 0 and keeps rising while any voxel remains open near the reaction band |
| `rate_FeS` peak | Positive, largest in the band where the two fronts overlap |
| Flow re-solves | More than the initial one, appearing in the log as soon as a voxel first seals |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Seconds to low minutes on one core. This is a small case on purpose |

## 7. What to check

1. **Porosity falls monotonically.** It is the headline result; read it from
   the `porosity` column of `output/summary.csv` rather than by eye in the
   log. This case only seals, it cannot reopen anything, so any rise is a
   bug.
2. **The mineral forms in a band, not at a boundary.** Iron enters left and
   sulphide right; a mineral layer plated onto either inlet face means one
   supply is overwhelming the other and the case is not showing what it
   intends to.
3. **Fe2 and HS are consumed one for one, and FeS gains what they lose.** The
   stoichiometry is 1 : 1 : 1.
4. **The flow re-solves after a change.** Look for the flow solver reporting
   again mid-run; if it never does, no voxel ever reached `max_precipRho`.
5. **No `[NEG!]` warnings.** A second-order rate can overshoot in a step
   where one reactant is nearly exhausted, which is what a smaller `k` or
   `<ade_dt>` is for.

`postprocess.py` reports the porosity trend, the change in every field's
total, and flags any negative minimum for you, which covers checks 1, 3 and
5. It also runs the closed-species report on `FeS`, confirming its total only
accumulates. Checks 2 and 4 are read from the `.vti` output and the log.
There is no `<conserve>` tag declared for this case: `Fe2` and `HS` are each
fed from a Dirichlet boundary, so neither total is conserved and naming
either in a mass-balance check would guarantee a failure that means nothing.

## 8. What this case demonstrates

**Geometry that changes because of the chemistry running on it.** A mineral
substrate is not just a number carried by the transport solver, it is also a
trigger: crossing `<max_precipRho>` converts a voxel from pore to solid, and
the flow is re-solved on the result. Porosity itself becomes an output that
moves during the run, not a fixed property of the input geometry.

**What it deliberately leaves out.** The conversion is one way: nothing here
reopens a sealed voxel, the mineral can only grow. Example 14 runs the
reverse, a declared mineral phase dissolving and its voxels returning to
pore; example 15 runs both at once, which is where the two paths have to
agree about who owns a voxel. There is also no aqueous speciation here:
`Fe2` and `HS` are total dissolved concentrations, not the free ions an
activity model would give; example 04 adds the speciation solver.

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
./scripts/setup_case.sh 13_precipitation run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. No optional solver is needed. Cases 09, 12, 21 and 22
additionally need GLPK, and case 10 needs Python with COBRApy; this one needs
neither, so the plain cmake line is enough:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** The chemistry in this case is a C++ header
compiled into the executable, not data read at start-up. So editing
`defineAbioticKinetics.hh`, including changing `k_FeS`, means rebuilding.
Editing `CompLaB.xml`, including `<max_precipRho>` or `<update_interval>`, or
`input/geometry.dat`, does not.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the staggered slot, lines the pore-facing grain surfaces with the reactive phase, reports porosity, refuses to continue if the slot does not percolate. Standard library only. |
| **build** | `kinetics/defineAbioticKinetics.hh` | This case's chemistry, compiled in. Its `defineDissolutionRate()` is present but returns zero: this case does not dissolve anything, and an unused entry point still has to link. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top; the
kinetics header above then replaces the shared default it laid down, so the
law described in section 4 is what gets compiled in.
