# 07 - Biotic lattice Boltzmann: biomass on its own lattice, and the relaxation time that constrains it

## 1. The scenario

This is example 05 with `<solver_type>LBM</solver_type>`. Biomass is given
its own D3Q7 advection-diffusion lattice and transported exactly like a
dissolved species, with the same streaming and collision machinery and the
same kind of boundary conditions a solute gets, and it would advect with the
flow if there were any.

The population here is still the attached biofilm seeded in example 05, on
material number 3, not a planktonic population riding a current: `<Peclet>0
</Peclet>` means no flow field is built at all, and both of `Bug`'s
boundaries are closed, so biomass cannot leave the domain by any of three
independent routes. What distinguishes this case from 05 and 06 is only the
biomass solver, lattice Boltzmann rather than cellular automaton or finite
difference, and what that solver adds is a genuine numerical constraint the
other two do not have: a relaxation time.

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
      @@@@@ = the seeded patch, 108 of the 3168 open voxels,
              transported on its own D3Q7 lattice, unable to leave: both
              biomass boundaries are closed and there is no flow to ride

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition for
the solutes. The other four faces are an inert wall, drawn by
`preprocess.py`. The seeded biomass patch is the same as examples 05 and 06:
material number 3, at every `x` where `x % 4 == 0`, in the three voxel rows
next to the wall block there, through the full depth in z, 108 voxels in
total.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion and there is no velocity field for biomass to ride |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, two species: `donor`, `product` |
| **Diffusivity** | 5e-10 m2/s for both species, in pore and in biofilm alike |
| **Chemistry** | `kinetics/defineKinetics.hh`, compiled into the executable, byte-identical to examples 05 and 06 |
| **Biology** | `Bug`, one population, seeded on material 3 at initial density 1.0, biomass on its own D3Q7 lattice at 1e-10 m2/s |
| **Geometry change** | none. The pore space is fixed for the whole run |
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
the product is made inside the domain and stays there, and here, unlike CA
or FD, `closed` on the biomass boundary is what installs bounce-back on the
biomass lattice itself, one more reason biomass cannot leave.

## 4. The reaction

Identical to examples 05 and 06, one organism, one donor, Monod-limited
growth, first-order decay:

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

The rate law is byte-identical to `kinetics/defineKinetics.hh` in example 05;
see that case's README for the code excerpt. What is different in this case
is how the resulting biomass moves, and that movement carries its own
constraint, described next.

**The relaxation time, which is this solver's real constraint.** `LBM`
carries a diffusivity as a lattice relaxation time, not directly:

```
                   D
    tau   =   ───────────  *  ( tau_ref - 1/2 )   +   1/2
               D_ref
```

with `D_ref` the reference solute diffusivity of 5e-10 m2/s used elsewhere in
this case, and `tau_ref` the same 0.8 used for the flow lattice. At this
case's biomass diffusivity of 1e-10 m2/s, that gives tau = 0.56.

Example 06's biomass diffusivity, 3e-13 m2/s, would put tau at just above
0.5 by that same formula, and a BGK relaxation time that close to 0.5 does
not diffuse slowly, it rings: spurious oscillation dominates the signal and
produces negative biomass in a large fraction of voxels. **That is why this
case cannot reuse example 06's diffusivity.** It uses a value four orders of
magnitude larger instead, one that is no longer physically the same
population's biomass diffusivity, purely so the relaxation time is stable;
that trade-off is exactly the cost of putting biomass on this solver. The
solver itself audits every relaxation time that will actually be used,
skipping species marked `<immobile>` and biomass on `CA` or `FD` whose tau is
never touched, and refuses to start below 0.51, so this cannot be walked
into by accident.

That constraint is the reason slow biomass usually belongs on `CA` or `FD`
instead: neither has a relaxation time to satisfy, `CA` needs no diffusivity
at all, and `FD` carries a diffusivity like example 06's without difficulty,
as shown in that case's README.

## 5. What happens each step

The solver repeats this 1000 times:

1. **Stream and collide** `donor`, `product` and `Bug` on their own D3Q7
   lattices, which advances diffusion by one step. There is no velocity
   field, so nothing here advects.
2. **Apply the boundary conditions**, re-pinning `donor` to 0.5 at the left
   face and bouncing `Bug` back at both ends.
3. **Evaluate the rate law** in every open voxel that holds biomass: read the
   local donor concentration and biomass density, compute `growth`, and write
   the donor, product and biomass increments.
4. **Apply the biomass increment**, which then streams and collides on the
   biomass lattice on the next step like any other conserved field.

The stability caveat is the relaxation time described in section 4: too
close to 0.5 and it produces negative biomass that has nothing to do with
the chemistry. The donor side of the reaction has the same caveat as
examples 05 and 06: the compiled kinetics path applies each increment as
computed with no positivity clamp, so a time step or rate constant too
large for the local donor drives it negative and reports `[NEG!]`.

## 6. What comes out

```
output/
  donor_0000200.vti  ...  donor_0001000.vti     concentration of the donor
  product_*.vti                                  the product plume
  Bug_*.vti                                       biomass density, on its own lattice
  rate_donor_*.vti  rate_product_*.vti            the reaction rate as a field, mol/L/s
  summary.csv                                     one row every 100 steps
  run.log                                         the whole run, including the checks
```

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Conserved biomass total, `Bug_total + Bug_held` | Rises, by roughly the same amount as examples 05 and 06, since the chemistry is unchanged |
| Reported `Bug_total` alone | Can fall even while the conserved total rises, since `computeDensity()` on a bounce-back wall ignores populations in flight there; `Bug_held` is exactly that difference |
| `Bug_min` | Zero or above. A negative value here points at the relaxation time, not the chemistry |
| Total donor | Falls near the colony, replenished from the left boundary |
| Minimum donor | Zero or above; a negative value here is the usual kinetics stability issue |
| Wall clock | Seconds on one core, somewhat more than examples 05 and 06 for the extra lattice |

## 7. What to check

1. **The conserved total, not the visible one.** `Bug_total + Bug_held` in
   `output/summary.csv` is the quantity that should hold steady under pure
   transport and rise only from the reaction. The visible `Bug_total` alone
   can fall while the real total rises, because of what `computeDensity()`
   misses at a bounce-back wall, described above.
2. **No negative biomass anywhere.** `Bug_min` in the record. A negative
   value means the relaxation time is too close to 0.5, not that the
   chemistry is wrong.
3. **No `[NEG!]` warnings on the donor.**
4. **Compare the total against 05 and 06.** Same chemistry, so growth should
   be of the same order of magnitude. It need not match to the last digit,
   because the biomass lattice streams and the other two solvers do not; this
   is a sanity check on scale, not a digit-for-digit match.

`postprocess.py` runs check 3 for you and reports the change in `product`.
Checks 1, 2 and 4 are read from `summary.csv`, using the `_held` column
described above.

## 8. What this case demonstrates

**Biomass transported on its own lattice, and the relaxation-time constraint
that comes with it.** It is the right choice when biomass genuinely moves
with a flow field, a planktonic population in a flowing pore space, and
overkill when it does not, as here, where the case exists to show what the
solver costs rather than to need it. It also demonstrates the `_held`
diagnostic column: the difference between what a bounce-back wall reports
and what it actually holds, which matters most exactly on the field that
spends the most time near a wall.

**What it deliberately leaves out:** any flow for biomass to actually ride,
and any change to the pore space. Examples 05 and 06 are this same case with
the cellular automaton and finite-difference solvers instead, neither of
which carries this relaxation-time constraint; example 08 adds a second
population competing for the same donor.

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
./scripts/setup_case.sh 07_biotic_lattice_boltzmann run/mycase
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
means rebuilding. Editing `CompLaB.xml`, including
`<biomass_diffusion_coefficients>`, or `input/geometry.dat`, does not; the
solver itself checks the resulting relaxation time at start-up and refuses to
run if it falls below 0.51.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and the seeded biomass patch, the same as examples 05 and 06. Reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **build** | `kinetics/defineKinetics.hh` | This case's chemistry, compiled in. Byte-identical to examples 05, 06 and 08. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space and the patch. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change including the `_held` correction, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
