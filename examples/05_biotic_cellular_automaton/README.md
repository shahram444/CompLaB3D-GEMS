# 05 - Biotic cellular automaton: the first organism, and biomass that spills when a voxel fills

## 1. The scenario

This is the first case with something alive in it. One population, seeded in
a patch, consumes a dissolved donor, grows, and releases a product:

**donor -> biomass + product**

The rate is Monod in the donor (growth rises with donor concentration and
saturates at high concentration, the same functional form as
Michaelis-Menten enzyme kinetics), so growth is fast where the donor is
plentiful and slows as it runs out. Decay is first order and always on, so
the population approaches a steady value rather than growing without bound.
A reactive-transport modeller will recognise this as the standard substrate
limited growth model used for a biofilm colonising a pore surface.

What is specific to this case is how the biomass moves once it grows. With
`<solver_type>CA</solver_type>` it does not diffuse at all. `CA` stands for
cellular automaton: a model where each voxel's state updates from a local
rule rather than a continuous transport equation. Here the rule is the
classical biofilm one: a voxel accumulates biomass until it reaches
`<maximum_biomass_density>`, and only then does the excess spill into its
open neighbours. Nothing moves before a voxel fills.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
    |....................  ####  ..........................  |
donor|...................  ####  ..........................  |zero
0.5 |@@@@@@@@@@..........  ########  ......................  |gradient
held|@@@@@@@@@@..........                                     |
    |                                                         |
    +---------------------------------------------------------+
      no flow: Peclet = 0. donor diffuses in from x=0
      @@@@@ = the seeded biomass patch, 108 of the 3168 open voxels,
              against the wall blocks at x % 4 == 0

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition. The
other four faces are an inert wall, drawn by `preprocess.py`, because the
solver gives those faces nothing at all and a face with no condition reads
back values that nothing updates. The biomass patch itself sits on material
number 3, at every `x` where `x % 4 == 0`, in the three voxel rows next to
the wall block there, through the full depth in z: 108 voxels in total,
carved out of open pore rather than added to it, so the porosity is
unchanged from the abiotic cases.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, two species: `donor`, `product` |
| **Diffusivity** | 5e-10 m2/s for both species, in pore and in biofilm alike |
| **Chemistry** | `kinetics/defineKinetics.hh`, compiled into the executable |
| **Biology** | `Bug`, one population, seeded on material 3 at initial density 1.0, biomass moved by the cellular automaton |
| **Geometry change** | none. The pore space is fixed for the whole run; biomass fills voxels but does not remove or add solid |
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
this solver moves only voxel to voxel, never crosses a domain face at all.

## 4. The reaction

One organism, one donor, Monod-limited growth, first-order decay:

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
| Ks | 5.0e-2 | mol/L | half-saturation constant on the donor: the donor concentration at which growth runs at half its maximum rate |
| Y | 0.4 | gDW/mol | yield: biomass made per mole of donor consumed |
| fP | 2.0 | mol/gDW | product released per unit of biomass made |
| kd | 1.0e-6 | 1/s | first-order decay rate |

These live in `kinetics/defineKinetics.hh`, next to the equation that uses
them, rather than in `CompLaB.xml`. The `<half_saturation_constants>` and
`<maximum_uptake_flux>` tags in the XML feed the flux-balance and surrogate
rate paths, which build their own uptake bounds from them; a hand-written
kinetics file is not obliged to read them, and this one does not.

The core of the rate law, evaluated for every organism in every open voxel:

```cpp
    if (C.size() < 2 || subsR.size() < 2) return;   // a short substrate list
    if (mask < 2) return;                            // no biology in solid or wall voxels

    const double S = std::max(C[0], 0.0);

    for (std::size_t m = 0; m < B.size() && m < bioR.size(); ++m) {
        const double Bm = std::max(B[m], 0.0);
        if (Bm <= 0.0) continue;

        const double growth = mu_max[p] * monod(S, Ks[p]) * Bm;   // gDW/L/s

        bioR[m]  += growth - kd[p] * Bm;
        subsR[0] -= growth / Y[p];                                // donor consumed
        subsR[1] += growth * fP[p];                               // product released
    }
```

The guards at the top are the part to copy when writing a new rate law: the
first stops it reading past the end of a shorter-than-expected concentration
vector, the second keeps chemistry out of solid grains and walls. The floor
at zero on both `S` and `Bm` matters because the transport and biomass
solvers can each leave a value fractionally below zero, and a rate law that
multiplies such a value into a growth term turns a rounding artefact into a
reaction.

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
5. **Run the cellular-automaton rule.** Wherever a voxel's biomass now exceeds
   `<maximum_biomass_density>`, spill the excess (`<CA_method>fraction</CA_method>`
   moves only the amount above the maximum, not the whole voxel) into an open
   neighbouring voxel.

The stability caveat is on step 3, not step 5: the compiled kinetics path
applies each reaction increment as computed, with no clamp against what a
voxel actually holds. A time step or a rate constant large enough to consume
more donor than a voxel holds will drive that voxel negative. The solver
reports that as a `[NEG!]` warning, and the fix is a smaller `<ade_dt>` or a
smaller rate constant, not a positivity clamp, because this path has none.

## 6. What comes out

```
output/
  donor_0000200.vti  ...  donor_0001000.vti     concentration of the donor
  product_*.vti                                  the product plume
  Bug_*.vti                                       biomass density of the colony
  rate_donor_*.vti  rate_product_*.vti            the reaction rate as a field, mol/L/s
  summary.csv                                     one row every 100 steps
  run.log                                         the whole run, including the checks
```

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Total biomass | Rises from the seeded patch and levels off as growth and decay balance |
| Occupied voxel count | Stays at 108, the seeded patch, unless a voxel hits `<maximum_biomass_density>` and the automaton fires |
| CA trigger count | Zero, or small, at these parameters and this run length; a run long enough or seeded dense enough will eventually fill a voxel and spill |
| Total donor | Falls near the colony, replenished from the left boundary |
| Total product | Rises from 0, since it is closed and made only where biomass is |
| `rate_donor` | Negative wherever the colony is active; zero elsewhere |
| Minimum donor | Zero or above. A negative value is a real failure |
| Wall clock | Seconds on one core. This is a small case on purpose |

## 7. What to check

1. **Biomass rises and then levels off** as growth and decay balance.
2. **The donor is drawn down around the colony** and replenished from the
   left boundary; check the donor minimum, not its mean, since a flat total
   can hide a patch that is not being reached at all.
3. **A product plume spreads from the colony.** It is made only where biomass
   is, so its shape is a map of where the reaction actually ran.
4. **No `[NEG!]` warnings.** The compiled kinetics path applies each
   increment as computed, with no positivity clamp, so a time step or rate
   constant too large for the local donor will drive it negative.
5. **Report the total biomass, not the peak.** For this solver the two move
   together, since biomass never leaves the seeded voxels until one fills;
   for the `FD` and `LBM` solvers in examples 06 and 07 they do not, because
   biomass there is already spreading out.

`postprocess.py` runs checks 1, 3 and 4 for you: it reports the change in
`product` and whether it is monotone, and it flags any negative field. Checks
2 and 5 are read from the `.vti` output and the biomass total in
`summary.csv`.

## 8. What this case demonstrates

**The compiled kinetics reaction path on an organism, and the cellular
automaton for biomass.** The rate fires in proportion to the biomass
present, which is what distinguishes it from the abiotic path in example 03,
where the rate fires in every open voxel whether or not anything is alive
there. The automaton itself demonstrates the simplest of the three ways this
solver can move biomass: nothing happens until a voxel is full, and then the
excess spills.

**What it deliberately leaves out:** any change to the pore space, any flow,
and any spreading of biomass short of a full voxel. Examples 06 and 07 run
this exact same chemistry, the same geometry, the same seeded patch, and
change only `<solver_type>`, to `FD` and `LBM`, so the biomass spreads
continuously instead. Example 08 adds a second population competing for the
same donor.

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
./scripts/setup_case.sh 05_biotic_cellular_automaton run/mycase
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
`<maximum_biomass_density>`, or `input/geometry.dat`, does not.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and the seeded biomass patch, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **build** | `kinetics/defineKinetics.hh` | This case's chemistry, compiled in. Shared byte-for-byte with examples 06, 07 and 08. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space and the patch. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
