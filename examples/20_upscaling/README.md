# 20 - Upscaling: one aggregate reduced to two numbers a continuum model needs

## 1. The scenario

A pore-scale run like case 19 resolves everything happening inside a single
microbial aggregate, voxel by voxel. A column-, basin- or reservoir-scale
reactive-transport model cannot afford to do that: it carries one
concentration per grid block, possibly a metre across, and needs a single
volumetric rate to close its mass balance for that block. The obvious way to
get that rate is to evaluate the aggregate's own rate law at the block's
concentration. That is wrong, because the volume average of a rate that
saturates and shuts off is not the same as the rate evaluated at the
volume-averaged concentration, and the gap between the two is exactly the
information a coarse model is missing.

The standard fix is the effectiveness factor, written eta: the aggregate's
actually-measured mean rate, averaged over every voxel inside it, divided by
the rate the same law would give a single, well-supplied cell sitting right at
the bulk composition outside. Eta says how much of the aggregate's rated
capacity the aggregate is actually delivering; a continuum model multiplies
its block-scale rate estimate by eta and gets back the answer the resolved
simulation would have given.

Eta on its own is not enough to apply elsewhere, because it depends on how
big the aggregate is and how starved its interior is, both of which change
from one grid block to the next. The classical tool for the size dependence
is the Thiele modulus, a dimensionless number that compares how fast the
reaction consumes a substrate against how fast diffusion can resupply it
across the aggregate's own radius. A small Thiele modulus means diffusion
easily keeps up and eta stays close to 1; a large one means the interior runs
out of substrate faster than diffusion refills it and eta falls well below 1.

This case is the same aggregate as case 19, posed so that eta is a genuine
steady-state property rather than a snapshot of a still-changing transient,
and it reports eta against both the Thiele modulus and the thermodynamic gate
at the bulk composition, so the measurement can be read back later at sizes
and compositions this one run never visited.

## 2. The picture

```
   32 x 34 x 10 voxels at 10 um per voxel  =  320 x 340 x 100 um

   x=0                                                        x=31
    |                                                            |
CH4 |>>>>>>>>>>>>>>>>>       ,-''--,        <<<<<<<<<<<<<<<<<<<< | CH4
SO4 |>>>>>>>>>>>>>>>>>      /  ANME \       <<<<<<<<<<<<<<<<<<<< | SO4
HS  |>>>>>>>>>>>>>>>>>     |  sphere |      <<<<<<<<<<<<<<<<<<<< | HS
HCO3|>>>>>>>>>>>>>>>>>      \ r=60um/       <<<<<<<<<<<<<<<<<<<< | HCO3
    |                        `--,,-'                             |
    +------------------------------------------------------------+
      all four species held at the SAME value on both faces
      "the aggregate sits in a bulk of composition C" is the whole setup

      the sphere = biofilm, material 3, radius 6 voxels (60 um),
      896 voxels, biomass held fixed for the whole run (frozen catalyst)
      open water = material 2 (the bulk), porosity 0.7529, 8192 open voxels
```

The other four faces are inert wall, added by `preprocess.py` the same way as
in case 19: the solver conditions nothing there, so a face left open would
stream a lattice Boltzmann field off the edge of the block and read back an
envelope nothing updates.

## 3. What goes in

| | |
|---|---|
| **Domain** | 32 x 34 x 10 voxels, dx = 10 um, 8192 open voxels, porosity 0.7529 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `CH4`, `SO4`, `HS`, `HCO3` |
| **Chemistry** | `kinetics/defineKinetics.hh`, the same dual Monod law as case 19, compiled in |
| **Thermodynamic gate** | on, from `input/aom.thm`, identical file to case 19 |
| **Biology** | one organism, `ANME`, biomass held fixed. `<freeze_biomass>true</freeze_biomass>` |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Upscaling** | on. Reports eta, the Thiele modulus and the classical sphere result every interval |
| **Run length** | 24000 advection-diffusion steps, VTK output every 4000 |

Boundary conditions, changed from case 19 in exactly the way this case needs:

| Species | Left face (x = 0) | Right face (x = 31) | Initial |
|---|---|---|---|
| **CH4** | Dirichlet, 1.0e-2 mol/L | Dirichlet, 1.0e-2 mol/L | 1.0e-2 |
| **SO4** | Dirichlet, 2.8e-2 mol/L | Dirichlet, 2.8e-2 mol/L | 2.8e-2 |
| **HS** | Dirichlet, 1.0e-5 mol/L | Dirichlet, 1.0e-5 mol/L | 1.0e-5 |
| **HCO3** | Dirichlet, 2.3e-3 mol/L | Dirichlet, 2.3e-3 mol/L | 2.3e-3 |

*Dirichlet on both ends, at the same value,* is what turns the open water into
a genuine bulk: a single, fixed exterior composition, exactly the one number a
continuum grid block carries. Holding the products at their reservoir value
too does not force the answer; the gate still closes because sulfide and
bicarbonate build up inside the aggregate relative to that fixed exterior, and
that internal gradient is free to develop however the chemistry dictates.

**Two more settings needed for eta to mean anything, on top of the Dirichlet
change.** `<freeze_biomass>true</freeze_biomass>` still lets the reaction
consume and produce every step, but discards the biomass increment so the
catalyst distribution itself never moves. Without it, biomass keeps growing
faster at the rim, where it is better fed, than in the interior, and the
measured eta keeps drifting for a reason that has nothing to do with
transport; every run without this setting reported not steady at five
diffusion times. And the run has to be long enough for the interior
concentration profile to relax, which happens on the diffusive timescale, the
aggregate radius squared over the diffusivity: at 60 um and 1.5e-10 m2/s that
is about 24 seconds, so the run covers five of them, 24000 steps at 5 ms each.

## 4. The reaction

Identical to case 19, section 4: the same AOM stoichiometry, the same dual
Monod intrinsic rate law with the same `mu_max`, `K_CH4`, `K_SO4` and `Y`, and
the same thermodynamic gate equation reading the same `input/aom.thm`. Nothing
about the chemistry changes between the two cases; only the boundary
conditions, the freeze on biomass and the run length do, because those three
are what a steady-state measurement of eta requires and case 19's transient
setup does not provide.

## 5. What happens each step

The solver repeats this 24000 times:

1. **Stream and collide** each of the four species on its own D3Q7 lattice.
2. **Apply the boundary conditions**, re-pinning all four species to the same
   Dirichlet value at both faces.
3. **Evaluate the intrinsic rate and the thermodynamic gate**, exactly as in
   case 19, in every voxel holding biomass.
4. **Add the increments to the four concentration lattices only.** The
   biomass increment is computed as usual and then discarded, which is what
   `<freeze_biomass>` means: the catalyst distribution from the last biomass
   update stays fixed for the rest of the run.
5. **Sample the effectiveness factor.** The upscaling diagnostic reads the
   rate increments the solver is about to apply, before they are applied,
   averages them over the aggregate, and divides by the rate the same law
   gives at the current bulk concentration.

**The steadiness check is the solver's own, not a postprocessing guess.**
It compares the last three diagnostic intervals of eta and reports the run as
not steady if they disagree by more than 1 percent. Reading eta from a run
still reporting not steady is reading a transient value, not the steady-state
quantity the effectiveness factor is defined to be.

## 6. What comes out

```
output/
  CH4_0004000.vti  ...  CH4_0024000.vti     concentration of methane
  SO4_*.vti  HS_*.vti  HCO3_*.vti           the other three species
  ANME_*.vti                                 the biomass field, unchanging by design
  rate_CH4_*.vti  rate_SO4_*.vti  ...       the reaction rate as a field, mol/L/s
  dG_ANME_*.vti                              free energy of the reaction, kJ/mol,
                                              negative where it releases energy
  FT_ANME_*.vti                              the gate, 0 to 1, the factor the
                                              rate was multiplied by
  summary.csv                                one row every 1000 steps
  upscaling.csv                              one row every 1000 steps: bulk
                                              concentration, rate at bulk,
                                              measured mean rate, eta, the
                                              Thiele modulus, the classical eta,
                                              the gate at bulk, steadiness
  run.log                                    the whole run, including the
                                              [UPSCALE] report
```

`dG_ANME` and `FT_ANME` are written as zero in solid, wall or outer-column
voxels and nowhere else, the same convention as case 19: a zero there always
means no chemistry, never a reaction sitting at exactly zero free energy.

| Quantity | What it should do |
|---|---|
| eta at the first interval | Exactly 1, or extremely close to it. Every voxel starts at the bulk composition, so the measured mean rate must equal the rate at the bulk |
| eta over the run | Falls from 1 and settles below it, as the internal profile develops |
| Steadiness | Reported not steady early on, reported steady once the last three intervals agree to within about 1 percent |
| Thiele modulus | A single number, order 1 for this aggregate size and diffusivity |
| Classical eta at that Thiele modulus | Printed beside the measured eta, from the exact sphere solution for first-order kinetics with no energy limit |
| Measured eta versus classical | Below the classical curve, because the thermodynamic gate holds the interior back in a way no classical, energy-blind treatment captures |
| Gate at the bulk composition | Close to but below 1; even the exterior sits somewhat inside the thermodynamically favourable region |

## 7. What to check

1. **Eta is 1 at iteration zero.** This is the measurement checking itself:
   every voxel is still at the bulk composition, so the numerator and
   denominator of eta must be the same quantity. If this is not 1, nothing
   later in the record is meaningful, and the fault is in how the rate is
   read, not in the physics.
2. **The run reports steady before you quote eta.** An effectiveness factor
   read off a transient is not the quantity the case is defined to measure.
3. **Measured eta against the classical curve.** Read the direction of the
   gap, not just its size: a measured value below the classical one says the
   thermodynamic gate, not diffusion alone, is what is limiting the interior.
4. **One run is one point.** A single solver run gives eta at one aggregate
   radius and one bulk composition. Building a usable eta(Thiele modulus,
   gate) relationship needs a sweep over both, which `offline/upscale.py`
   runs as a series of separate solver invocations; it is not part of
   `run.sh` because each point is its own multi-minute run.

`postprocess.py` runs checks 1 and 2 for you, in that order, and stops with a
clear failure message if either does not pass; it prints the numbers for
check 3 and tells you which side of the classical curve they fall on.

## 8. What this case demonstrates

**The effectiveness factor and the Thiele modulus**, the two numbers a
continuum model needs to stand in for a resolved aggregate, measured under
conditions, a fixed bulk, a frozen catalyst and a long enough run, that make
eta an actual steady-state property rather than a number still in motion. The
gap between the measured eta and the classical sphere result is exactly the
thing a continuum model gets wrong if it applies textbook diffusion-and-
reaction theory to a reaction with an energy threshold.

**What it deliberately leaves out.** The same list as case 19: no flow, no
abiotic reaction, no change to the pore space, and no aqueous speciation. It
also does not, on its own, give you a usable relationship for a different
aggregate size or a different bulk composition; that is what the sweep in
`offline/upscale.py` is for, and one run of this case is one point on that
curve, not the curve itself.

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
./scripts/setup_case.sh 20_upscaling run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. No optional solver is needed; upscaling and the
thermodynamic gate are both part of the core solver, not build flags:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** `CompLaB.xml`, `input/geometry.dat` and
`input/aom.thm` are all read at start-up, so editing any of them, including
the boundary conditions, `<freeze_biomass>` or the energetics, needs no
rebuild. `kinetics/defineKinetics.hh` is compiled in, so editing the intrinsic
rate law itself does.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the aggregate and the walled box around it, and reports the equivalent-sphere radius the solver uses for the Thiele modulus. Standard library only. |
| **build** | `kinetics/defineKinetics.hh` | The ungated intrinsic rate law, identical to case 19, compiled in. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The aggregate and the bulk around it. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| | `input/aom.thm` | The energetics, identical to case 19. Read at start-up and echoed into the log. |
| **post** | `postprocess.py` | Checks that eta starts at 1, checks steadiness, and reports eta against the classical curve. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |
| | `offline/upscale.py` | Sweeps aggregate radius and bulk sulfide, one solver run per point, and fits the correction the classical curve needs once the energy limit is present. Not part of `run.sh`. |
| | `offline/expected_sweep.csv` | Reference points from a prior sweep, to check your own against. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
