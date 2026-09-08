# 12 - Mixed reaction types: flux balance plus a compiled law, one organism

## 1. The scenario

A metabolic model, a stoichiometric matrix of an organism's internal
reactions solved by flux balance analysis (FBA), describes what an
organism's network *can* do. It does not describe everything the organism
actually does. The classic omission is maintenance: the donor a cell
consumes simply to stay alive, to hold its membrane potential, turn over
proteins, run its pumps, which yields no biomass and appears in no
flux-balance objective. Leave it out and the model over-predicts growth
wherever substrate is scarce, which in a pore space is most of the domain.

A reactive-transport modeller can read maintenance as the microbial
equivalent of a background sink term, a rate that runs regardless of growth.
A computational biologist will recognise it as a standard correction to a
plain FBA objective, usually added as a fixed lower bound on an ATP-drain
reaction inside the model; here it is added outside the model instead, as a
second, independently written rate path on the same organism.

This case runs one organism through two rate paths at once. Flux balance
analysis supplies the growth, exactly as in example 09; a compiled rate law
in `kinetics/defineKinetics.hh` supplies the maintenance draw the linear
program knows nothing about. The two rates are added. That is what the
combined reaction types in this code are for:

```
    glpk_and_kinetics        cobrapy_and_kinetics        surrogate_and_kinetics
```

Each names a metabolic back end plus the compiled kinetics header, per
organism. This case uses `glpk_and_kinetics`.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um
   (the pore space itself is 24 x 24 x 6; preprocess.py adds a one voxel
   wall on the y and z faces the solver conditions with nothing, so NY
   and NZ come out two larger than the pore space it declares)

   x=0                                                          x=23
    |                                                             |
S,O |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> |
1.0 |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> | Neumann,
    |>>>  ====                                          >>>>>>>> | open outlet
    |               ####    ####    ####    ####    ####         |
    |               ====                                         |
    +-------------------------------------------------------------+
      grain blocks alternate between the y=0 and y=23 walls every
      four voxels along x, so diffusion has to weave to get through

      #### = solid grain (wall)     ==== = initial Toybug patch
      porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry boundary conditions. The
other four faces are an inert wall, drawn by `preprocess.py`, because the
solver gives a face with no declared condition nothing at all, and a lattice
Boltzmann field with nothing conditioning it streams off the edge of the
domain and reads back values that nothing updates. This is the same
geometry and the same one-organism patch examples 09 and 10 use.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `S`, `O`, `P` |
| **Diffusivity** | 5e-10 m2/s for all three species, in pore and in biofilm alike |
| **Chemistry** | `kinetics/defineKinetics.hh`, compiled into the executable, alongside FBA |
| **Biology** | `Toybug`, **one** organism, `<reaction_type>glpk_and_kinetics</reaction_type>`, attached biofilm seeded on its own material number |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 500 advection-diffusion steps, output every 100 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **S** (donor) | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **O** (acceptor) | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **P** (product) | closed | closed | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, which is a reservoir. *Neumann with a zero
gradient* means nothing flows across that face, an open outlet the species
never actually reaches in this run. *Closed* means the same thing on both
ends: P is made inside the domain by the reaction and stays there.

## 4. The reaction

Two rates on one organism, added together.

From the linear program, the same four-reaction toy model examples 09 and 10
use, with the uptake bounds rebuilt every step from the local concentrations.
In words: maximise the objective flux (growth) subject to every internal
metabolite being at steady state and every flux staying within its bound.

```
    maximise   c' v          the objective flux (growth)
    subject to S v  = 0
               lb <= v <= ub  with  V = Vmax * C / ( Kc + C )
```

From the compiled kinetics, a maintenance drain, first order in biomass and
independent of the reaction network, running as long as any donor is
present:

```
    r_maint  =  m_S * B                     mol of donor per litre per second

    m_S = 2.0e-6    mol donor per gDW per second, yielding no biomass
```

capped at what the voxel actually holds, so a step cannot draw a substrate
negative.

The sign discipline is the whole difficulty of a combined type. The kinetics
header in this folder writes `subsR`, the substrate rates, and deliberately
leaves `bioR`, the biomass rate, at zero:

```cpp
    subsR[0] = -std::min(drain, S);

    // bioR deliberately left at zero: growth is the FBA path's job.
```

Growth belongs to the flux-balance path. Writing it in both places
double-counts it, and nothing complains: the run completes, the numbers look
plausible, and the organism grows at twice the rate the model says. That is
the classic mistake with combined reaction types, and the reason the comment
sits in the header itself rather than only here.

Two Vmax tags, because the two paths read different units, is worth stating
on its own:

| tag | read by | units |
|---|---|---|
| `<maximum_uptake_flux>` | the compiled kinetics path | mol per mol cells per second |
| `<fba_maximum_uptake_flux>` | the flux-balance path | mmol gDW-1 h-1 |

One number cannot mean both things, so a combined microbe carries both tags.
Here `<maximum_uptake_flux>` is `0. 0. 0.`, the maintenance law does not use
it, it uses its own `m_S`, while `<fba_maximum_uptake_flux>` is `10. 4. 0.`,
the same bounds as example 09.

Maintenance and decay are not the same thing, and both are present here.
Maintenance is a substrate draw that produces no biomass; `<decay_coefficient>`,
set to `0.01`, is a first-order biomass loss that returns nothing to the
solutes. A cell pays the first while alive and suffers the second
regardless.

## 5. What happens each step

The solver repeats this 500 times:

1. **Stream and collide** each species on its own D3Q7 lattice, advancing
   diffusion by one step, and diffuse the biomass field the same way.
2. **Apply the boundary conditions**, so S and O are re-pinned to 1.0 at the
   left face.
3. **Build the uptake bounds** for the flux-balance path, per voxel, per
   substrate, from the local concentration.
4. **Solve the linear program**, in every voxel that holds biomass. This is
   the expensive step, one full linear-program solve per occupied voxel per
   iteration, the same cost as example 09 pays for the same organism.
5. **Evaluate the compiled kinetics law** in the same voxel: the maintenance
   drain, computed from the local biomass alone, independent of the linear
   program's result.
6. **Add the two rates.** The flux-balance path's growth increment and the
   kinetics path's substrate drain both land in the same change lattices for
   S, O and P, so they accumulate rather than overwrite each other.
7. **Add the increments**, clamped so a step cannot draw a voxel's substrate
   below what it actually holds, regardless of which path (or both) drew it
   down.

## 6. What comes out

```
output/
  S_0000100.vti  ...  S_0000500.vti        concentration of the donor
  O_*.vti  P_*.vti                          the other two species
  Toybug_*.vti                              the biomass field
  rate_S_*.vti  rate_O_*.vti  rate_P_*.vti  the COMBINED reaction rate as a field, mol/L/s
  summary.csv                               one row every 50 steps
  run.log                                   the whole run, including the checks
```

| Quantity | What it should do |
|---|---|
| Growth, compared to example 09 | **Lower**, on the same geometry and the same bounds, because the donor is now also spent on maintenance |
| The gap versus example 09 | Widens deep in the domain, where diffusion delivers less substrate and maintenance takes a larger share of a smaller supply |
| S consumption vs O consumption | S falls faster than the 1 : 0.5 objective-reaction ratio would predict alone, because only S carries a maintenance term; O still follows the objective's stoichiometry |
| Total biomass | Rises more slowly than example 09's, for the same reason |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Close to example 09's; the added kinetics evaluation is cheap next to the linear-program solve |

## 7. What to check

1. **Growth is lower than example 09's**, on the same geometry and the same
   bounds, because the donor is now being spent on maintenance as well.
2. **The gap widens where substrate is scarce.** Deep in the domain, where
   diffusion delivers little, maintenance takes a larger share of a smaller
   supply, which is the whole reason to include it.
3. **The donor drops faster than the acceptor.** Only S carries a
   maintenance term, so its balance against biomass no longer follows the
   objective reaction's stoichiometry, while O still does.
4. **Set `m_S` to zero and example 09 comes back.** That is the cleanest
   test that the two paths are being added rather than one silently
   replacing the other.
5. **No `[NEG!]` warnings.** The maintenance draw is capped at what the
   voxel holds and the solver clamps again downstream.

`postprocess.py` reads `output/summary.csv` and `run.log`, reports every
field's change, flags any negative value, and states plainly that this
case's point is the shared substrate budget: both paths draw from the same
voxel under the same clamp, so a negative minimum is a solver bug worth
reporting rather than a tuning problem. That covers checks 3 and 5. Checks
1, 2 and 4 need a comparison against example 09, which is yours to run.

## 8. What this case demonstrates

**Two rate paths on one organism, added rather than chosen between, and the
sign discipline that keeps them from double-counting growth.** The combined
reaction types (`glpk_and_kinetics`, `cobrapy_and_kinetics`,
`surrogate_and_kinetics`) let a metabolic back end supply what its network
describes while a compiled header supplies what it does not, on the same
biomass field, in the same voxel, every step.

**What it deliberately leaves out:** flow, an abiotic reaction, and any
change to the pore space. The maintenance law here is the simplest one that
makes the point: constant per biomass, on one substrate. A real maintenance
term usually carries its own half-saturation constant and often a
temperature dependence; both would go in the same header, in the same
place. Example 09 is this same organism and geometry with the maintenance
term removed, the baseline this case is checked against.

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
./scripts/setup_case.sh 12_mixed_reaction_types run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
Palabos v2.3.0, and GLPK's development package (`libglpk-dev` on
Debian/Ubuntu, `glpk` from Homebrew on macOS). The organism's
`<reaction_type>` is `glpk_and_kinetics`, so, exactly as in example 09, the
executable must be built with the flag that links GLPK in:

```bash
cmake -B build -S . \
      -DCMAKE_BUILD_TYPE=Release \
      -DPALABOS_ROOT=/path/to/palabos-v2.3.0 \
      -DENABLE_GLPK=ON
cmake --build build -j
```

**When you must recompile.** `-DENABLE_GLPK=ON` only makes GLPK available in
the executable; `<enable_fba_glpk>true</enable_fba_glpk>` in `CompLaB.xml`
switches it on for a given run, and that line can be edited without a
rebuild. But the two have to agree: a build without the flag compiles
cleanly and then refuses at start-up, because turning the XML path on is not
enough if the flag was never passed to CMake. Separately, this case's
`kinetics/defineKinetics.hh` is a C++ header compiled into the executable,
so any edit to it, including changing `m_S`, means rebuilding.
`setup_case.sh` copies it over the shared default at the case root, which is
the copy that actually compiles. Editing `CompLaB.xml` or
`input/geometry.dat` needs no rebuild.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Exports the metabolic model into the flat format the solver reads. This case ships that file already built, so the script does nothing by default. |
| **build** | `kinetics/defineKinetics.hh` | This case's maintenance law, compiled in over the shared default. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| | `input/toy_model.xml` | The metabolic model of section 4, in the flat `<Metabolic_Model>` format. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top;
`kinetics/defineKinetics.hh` then replaces the default it laid down.
