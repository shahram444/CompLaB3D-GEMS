# 03 - Abiotic kinetics: a chemical reaction with no organisms in it

## 1. The scenario

Two dissolved species enter a pore space from opposite sides. They diffuse
toward each other, meet somewhere in the middle, and where they overlap they
react to make a third species that cannot leave.

If you have ever set up a counter-diffusion experiment in a gel plug or a
diffusion cell, this is that experiment. A reactive-transport modeller will
recognise it as the simplest possible reactive mixing problem: two solutes
supplied from fixed-concentration reservoirs at either end, a second-order
reaction where they overlap, and a product that accumulates in a band.

Nothing is alive here. No biomass field is allocated, no organism appears in any
equation, and the pore space never changes shape. This is the case to read first
if your chemistry is mineral reactions, redox couples, or sorption rather than
microbiology.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
 A  |>>>>>>>>>>>>>>>>>>  ####  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< | B
 1.0|>>>>>>>>>>>>>>>>>>  ####  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< |1.0
 held|>>>>>>>>>>>>>>>  ########  <<<<<<<<<<<<<<<<<<<<<<<<<<<< |held
    |                                                         |
    +---------------------------------------------------------+
      A diffuses right  ->        <-  B diffuses left
                        C is made where they overlap

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry boundary conditions. The other
four faces are an inert wall, drawn by `preprocess.py`, because the solver gives
those faces nothing at all and a face with no condition reads back values that
nothing updates.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species |
| **Diffusivity** | 5e-10 m2/s for all three species, in pore and in biofilm alike |
| **Chemistry** | `kinetics/defineAbioticKinetics.hh`, compiled into the executable |
| **Biology** | none. `<biotic_mode>false</biotic_mode>`, no biomass field allocated |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 1500 advection-diffusion steps, output every 200 |

Boundary conditions, which are the whole experiment:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **A** | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **B** | Neumann, zero gradient | Dirichlet, held at 1.0 | 0 |
| **C** | closed | closed | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no matter
what the interior does, which is a reservoir. *Neumann with a zero gradient*
means nothing flows across that face, which is an open outlet the species never
reaches. *Closed* means the same thing on both ends: C is made inside the domain
and stays there.

## 4. The reaction

One second-order reaction, one order in each reactant:

**A + B -> C**

```
    R  =  k · [A] · [B]                  k = 5.0e-2  L mol-1 s-1

    dA/dt = -R          dB/dt = -R          dC/dt = +R
```

That is the entire chemistry, and it is fifteen lines of C++ in
`kinetics/defineAbioticKinetics.hh`:

```cpp
    if (C.size() < 3 || subsR.size() < 3) return;   // a short substrate list
    if (mask < 2) return;                           // solid or wall voxel

    const double A = std::max(C[0], 0.0);           // floor tiny negatives
    const double B = std::max(C[1], 0.0);
    const double R = ExampleAbiotic::k_AB * A * B;  // mol/L/s

    subsR[0] = -R;   subsR[1] = -R;   subsR[2] = +R;
```

The first two lines are guards and are the part to copy when you write your own
rate law. The first stops the law reading past the end of a concentration vector
shorter than it expects. The second keeps chemistry out of solid grains and
walls. The floor at zero matters because the transport solver can leave a value
fractionally below zero, and a rate law that multiplies two such values turns a
rounding artefact into a reaction.

## 5. What happens each step

The solver repeats this 1500 times:

1. **Stream and collide** each species on its own D3Q7 lattice, which advances
   diffusion by one step.
2. **Apply the boundary conditions**, so A is re-pinned to 1.0 at the left face
   and B at the right.
3. **Evaluate the rate law** in every open voxel: read A and B there, compute
   `R = k·A·B`, and write the three increments into the change lattices.
4. **Add the increments** to the concentrations, as `C += R·dt`.

Step 4 is an explicit update with no clamp against what is actually available in
the voxel. That is the stability limit of this case, and it is worth knowing:
if `k·[A]·dt` approaches 1, a voxel can be driven negative. The solver reports
that as a `[NEG!]` warning, and the fix is a smaller time step or a smaller `k`.

## 6. What comes out

```
output/
  A_0000200.vti  ...  A_0001500.vti        concentration of A
  B_*.vti  C_*.vti                          the other two species
  rate_A_*.vti  rate_B_*.vti  rate_C_*.vti  the reaction rate as a field, mol/L/s
  summary.csv                               one row every 100 steps
  run.log                                   the whole run, including the checks
```

Numbers to expect, as orders of magnitude rather than exact values:

| Quantity | What it should do |
|---|---|
| Total A, total B | Both **rise** from 0, because both are fed from a reservoir |
| Total C | Rises from 0 and keeps rising; C is the only closed species |
| `rate_C` peak | Positive, largest in the mixing band near the middle |
| `rate_A` / `rate_C` | **Exactly -1.000000** in every voxel, from the 1:1:1 stoichiometry |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Seconds on one core. This is a small case on purpose |

The stoichiometry ratio is the sharpest check in the whole example set: it is
recovered from the output fields alone, with no knowledge of the rate law, and
it comes out exact to six decimals.

## 7. What to check

1. **All three totals rise, and C rises.** A and B are each fed from a
   reservoir, so the 1:1:1 stoichiometry constrains the per-step reaction
   increments, not the reported totals. C is the only closed species, so C is
   the total whose rise is entirely the reaction.
2. **No `[NEG!]` warnings** in `run.log`. One means the time step or the rate
   constant is too large for this chemistry.
3. **`rate_A` divided by `rate_C` is -1** wherever the reaction runs. This is
   the stoichiometry read straight out of the output.
4. **The reaction stops where the reactants run out**, not at a fixed time. Look
   at where the `rate_C` band sits and confirm it is where A and B overlap.

`postprocess.py` runs checks 1, 2 and 4 for you and prints a verdict. Note that
this case declares no `<conserve>` tag, so the solver itself checks no closed
sum. If you add one, do not name A or B in it: a species fed from a Dirichlet
boundary is not conserved and never will be, so such a sum would guarantee a
mass-balance failure that means nothing.

## 8. What this case demonstrates

**The abiotic reaction path.** One switch, `<enable_abiotic_kinetics>`, and one
header file you write yourself. It fires in every open voxel whether or not
anything is alive there, which is what distinguishes it from the biotic path
where a rate is proportional to the biomass present.

It also demonstrates the two things every reaction case in this set relies on:
that the rate is a field you can look at rather than a number in a log, and that
a case can be checked against its own stoichiometry without trusting the solver.

**What it deliberately leaves out:** flow, organisms, and any change to the pore
space. Example 13 uses this same abiotic path to precipitate a mineral and seal
the pore space; example 14 to dissolve one. Example 06 is the equivalent case
with an organism in it.

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
./scripts/setup_case.sh 03_abiotic_kinetics run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer, and
Palabos v2.3.0. No optional solver is needed. Cases 09 to 12 and 16 additionally
need GLPK or Python with COBRApy; this one does not, so the plain cmake line is
enough:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** The chemistry in this case is a C++ header
compiled into the executable, not data read at start-up. So editing
`defineAbioticKinetics.hh`, including changing `k_AB`, means rebuilding. Editing
`CompLaB.xml` or `input/geometry.dat` does not.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **build** | `kinetics/defineAbioticKinetics.hh` | This case's chemistry, compiled in. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
