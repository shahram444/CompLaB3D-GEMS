# 19 - Thermodynamic gate: a reaction that stops where the energy runs out

## 1. The scenario

Every rate law in this repository so far answers one question: how fast can
this organism go, given what substrate is here. A Monod term, a linear
program, a fitted network all do that and stop there. None of them asks a
second question that matters just as much in a diffusion-limited pore: is
there enough energy in the reaction, at the concentrations actually sitting in
this voxel, for the organism to run it at all.

The energy in a reaction is its Gibbs free energy change, written dG here. A
reaction with dG below zero releases energy and can, in principle, be run by a
cell; a reaction with dG at or above zero cannot, no matter how much substrate
is present. dG is not fixed: it depends on the local mix of reactants and
products through the reaction quotient Q, the ratio of product activities to
reactant activities. Push products up and reactants down and dG rises, even
though the reaction written on paper has not changed. A cell also cannot spend
every last joule of that energy on itself: some has to go toward making ATP,
the molecule a cell uses to store and move chemical energy internally, and a
reaction that does not clear that threshold yields nothing usable no matter
how favourable it looks on paper.

This case sets up exactly that situation for a single microbial aggregate, a
compact clump of cells rather than a thin film, 120 um across, oxidising
methane with sulfate. Reactant supply comes from open water outside the
aggregate; the products it makes have nowhere to go but out, by diffusion,
through the same crowded interior the reactants came in by. A dual Monod rate
law, the same shape used in examples 05 through 08, looks only at the two
reactants and keeps the reaction running at full tilt at the centre of the
aggregate long after the accumulated products have made that reaction
energetically pointless there. The thermodynamic gate is the correction: a
factor between 0 and 1 that the rate is multiplied by, computed fresh in every
voxel from the local dG, so the aggregate can develop a reacting shell around
a core that has genuinely stopped.

If your background is reactive transport, this is a local-equilibrium type
correction applied to a kinetic rate rather than a full speciation solve. If
your background is microbial physiology or bioenergetics, this is the same
energy bookkeeping you already know, run separately in every voxel of a
resolved pore space instead of once for a well-mixed culture.

## 2. The picture

```
   32 x 32 x 10 voxels at 10 um per voxel  =  320 x 320 x 100 um

   x=0                                                        x=31
    |                                                            |
CH4 |>>>>>>>>>>>>>>>>>       ,-''--,        <<<<<<<<<<<<<<<<<<<< |
SO4 |>>>>>>>>>>>>>>>>>      /  ANME \       <<<<<<<<<<<<<<<<<<<< |
1.0e|>>>>>>>>>>>>>>>>>     |  sphere |      <<<<<<<<<<<<<<<<<<<< |
held|>>>>>>>>>>>>>>>>>      \ r=60um/       <<<<<<<<<<<<<<<<<<<< |
    |                        `--,,-'                             |
    +------------------------------------------------------------+
      CH4, SO4 diffuse in ->     <- HS, HCO3 diffuse out
      the gate closes deepest in the sphere, where products pile up

      the sphere = biofilm, material 3, radius 6 voxels (60 um),
      120 um across, centred in x and y, spanning the full depth in z
      open water = material 2, porosity 0.7500, 7680 open voxels
```

Only `x = 0` and `x = 31` carry boundary conditions. The other four faces (`y
= 0`, `y = 31`, `z = 0`, `z = 9`) are an inert wall, drawn by `preprocess.py`,
because the solver gives a face with no declared condition nothing at all: a
lattice Boltzmann field left unconditioned streams off the edge of the block
there and reads back an envelope nothing updates. Those wall layers are
*added* around the water and the aggregate, not carved out of them, so the
120 um sphere and the porosity above are exactly what the case declares.

## 3. What goes in

| | |
|---|---|
| **Domain** | 32 x 32 x 10 voxels, dx = 10 um, 7680 open voxels, porosity 0.7500 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `CH4`, `SO4`, `HS`, `HCO3` |
| **Chemistry** | `kinetics/defineKinetics.hh`, an ordinary dual Monod law, compiled in |
| **Thermodynamic gate** | `<thermodynamics><enabled>true</enabled>`, energetics read from `input/aom.thm` |
| **Biology** | one organism, `ANME`, cellular automaton biofilm, seeded as the sphere above |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Upscaling** | on. Reports the effectiveness factor and Thiele modulus each interval; see case 20 |
| **Run length** | 12000 advection-diffusion steps, VTK output every 500 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 31) | Initial |
|---|---|---|---|
| **CH4** | Dirichlet, held at 1.0e-2 mol/L | Neumann, zero gradient | 1.0e-2 |
| **SO4** | Dirichlet, held at 2.8e-2 mol/L | Neumann, zero gradient | 2.8e-2 |
| **HS** | Dirichlet, held at 1.0e-5 mol/L | Neumann, zero gradient | 1.0e-5 |
| **HCO3** | Dirichlet, held at 2.3e-3 mol/L | Neumann, zero gradient | 2.3e-3 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, which is a reservoir: here it is the open-water
composition around the aggregate, held the same for every species including
the two products. *Neumann with a zero gradient* means nothing is forced to
flow across that face, an open outlet. All four species can reach it, but the
right-hand boundary sits far enough from the aggregate that little of the
methane or sulfate supplied at the left actually gets there before the
reaction consumes it. Every species starts at its own left-boundary value
everywhere, so the run begins with the whole domain already at the open-water
composition and only the aggregate's own reaction moves it away from there.

## 4. The reaction

Anaerobic oxidation of methane with sulfate as the electron acceptor,
catalysed by ANME archaea in syntrophic partnership with sulfate-reducing
bacteria. This one organism in the model, `ANME`, stands in for the whole
consortium: nothing here resolves the electron transfer between the two
partners, and that omission reappears below as a single tuned constant.

```
    CH4  +  SO4(2-)   ->   HS(-)  +  HCO3(-)  +  H2O
```

The intrinsic rate law, `kinetics/defineKinetics.hh`, is an ordinary dual
Monod term in the two reactants and knows nothing about the two products or
about energy:

```
                              [CH4]                [SO4]
    growth  =   mu_max  ·  ───────────────  ·  ───────────────  ·  B     gDW/L/s
                            K_CH4 + [CH4]       K_SO4 + [SO4]

                growth
    r       =   ──────                                                  mol/L/s
                  Y
```

with `mu_max = 1.0e-1 /s`, `K_CH4 = 1.0e-3 mol/L`, `K_SO4 = 5.0e-4 mol/L`,
`Y = 0.25 gDW/mol`, decay off. `mu_max` here is chosen to be far faster than a
real ANME-SRB consortium, which doubles on a timescale of months. What is
actually being reproduced in a short run is not that rate constant but the
Thiele modulus, the ratio of reaction speed to diffusion speed defined
properly in case 20, which is the dimensionless group that decides whether an
interior gradient forms at all.

The thermodynamic gate multiplies that intrinsic rate by a factor `F_T`
between 0 and 1, computed at every voxel from `input/aom.thm`:

```
    dG   =  dG0  +  R T ln Q

    F_T  =  max( 0 ,  1 - exp( (dG + m dG_ATP) / (chi R T) ) )

    realised rate  =  r x F_T
```

Read left to right: `dG0` is the reaction's standard free energy at the
file's own temperature, the energy it would release under a fixed reference
composition. `R T ln Q` corrects that for the composition actually sitting in
this voxel, where `Q` is the reaction quotient, products over reactants. `m
dG_ATP` is the energy the cell has to set aside to make its ATP, the wage the
reaction has to clear before any of it can be spent; `m` is the moles of ATP
made per mole of reaction and `dG_ATP` is the free energy of making one mole
of it. `chi R T` is the energy yardstick, built from the gas constant `R`, the
temperature `T` and the average stoichiometric number `chi`, that turns the
exponent into a plain dimensionless number. `F_T` itself is dimensionless and
always between 0 and 1, so the realised rate can only be equal to or smaller
than the intrinsic one, never larger.

The energetics that go into this, read from `input/aom.thm`:

| symbol | value | unit | meaning |
|---|---|---|---|
| temperature | 277.15 | K | 4 C, the seep temperature this dG0 is quoted at |
| dG0 | -33.12 | kJ/mol | standard free energy of the reaction at that temperature |
| atp | 0.25 | mol ATP / mol reaction | ATP made per mole of reaction, m above |
| dGatp | 50.0 | kJ/mol | free energy of making one mole of ATP |
| chi | 1.0 | - | average stoichiometric number |
| cmin | 1e-9 | mol/L | activity floor, so a species driven to zero cannot send ln Q to minus infinity |
| dGloss | 25.7 | kJ/mol | interspecies electron-transfer loss between ANME and SRB, added to dG |

`m x dGatp = 0.25 x 50.0 = 12.5 kJ/mol` is the energy threshold: the reaction
must clear that before the organism can conserve any of it as ATP.

**Where dGloss comes from, and why it matters that it is honest about not
being measured.** The consortium in this model is one lumped organism, so
nothing here computes the actual energy lost moving electrons from the ANME
cell to its sulfate-reducing partner. `dGloss` stands in for that loss as a
constant, and it is the only tuned number in the case. It is fixed against the
one measurement available: the growth efficiency of these consortia has been
measured at 1.0 to 2.75 percent of the methane carbon. With no transfer loss
at all, the yield equation below returns 5.18 percent at the aggregate rim,
which is too high; a loss of 25.7 kJ/mol brings that down to 2.75 percent, the
top of the measured band. That is a single degree of freedom fixed against a
measurement that is not the one being predicted, which is the honest form a
fitted parameter can take. Set it to 0 in a copy of the file to recover the
assumption of direct, lossless contact between the two partners, and both the
yield and the gate's position move together.

Optionally, the yield itself can be computed from the energetics rather than
fixed at a constant, which this file turns on:

```
    dG_dis  =  200 + 18 (6 - C)^1.8

    Y       =  dG / -( dG_ana + dG_dis )
```

with `dGana = -30.0 kJ/mol`, `carbon = 1`, `reduction = 8.0`, `yieldmax =
0.05` capping what is reported.

## 5. What happens each step

The solver repeats this 12000 times:

1. **Stream and collide** each of the four species on its own D3Q7 lattice,
   advancing diffusion by one step, and update the biomass field by the same
   cellular-automaton rule used in example 05: biomass stays put until a voxel
   fills, then spills into a neighbour.
2. **Apply the boundary conditions**, re-pinning all four species to their
   Dirichlet values at the left face.
3. **Evaluate the intrinsic rate**, the dual Monod law of section 4, in every
   voxel holding biomass.
4. **Evaluate the thermodynamic gate**, reading the local concentrations of
   all four species, computing `dG` and `F_T` from `input/aom.thm`, and
   multiplying the rate from step 3 by `F_T`.
5. **Add the increments** from the gated rate to the four concentration
   lattices and to the biomass field.

**Why the deepest voxel settles on the threshold instead of crossing it.**
The gate sits inside its own feedback loop. As `F_T` falls the reaction slows,
so it stops making the very products that were closing it, and the deepest
point in the aggregate settles asymptotically on the energy threshold rather
than passing through it. Expect the run's minimum `F_T` to be small but not
exactly zero for that reason; a report of exactly zero everywhere in the core
means transport, not the local chemistry, has pushed the composition past the
threshold. `F_T` also falls sharply rather than gradually with depth: with
`chi = 1` the whole transition from open to shut happens over only a few
kJ/mol of `dG`, so expect a front rather than a fade.

## 6. What comes out

```
output/
  CH4_0000500.vti  ...  CH4_0012000.vti     concentration of methane
  SO4_*.vti  HS_*.vti  HCO3_*.vti           the other three species
  ANME_*.vti                                 the biomass field
  rate_CH4_*.vti  rate_SO4_*.vti  ...       the reaction rate as a field, mol/L/s
  dG_ANME_*.vti                              free energy of the reaction, kJ/mol,
                                              negative where it releases energy
  FT_ANME_*.vti                              the gate itself, 0 to 1, the factor
                                              the rate was multiplied by
  summary.csv                                one row every 500 steps
  upscaling.csv                              effectiveness factor and Thiele
                                              modulus each interval, see case 20
  run.log                                    the whole run, including the [THM]
                                              and [UPSCALE] reports
```

`dG_ANME` and `FT_ANME` are written as zero in solid, wall or outer-column
voxels, and nowhere else: a live reaction sitting at exactly zero free energy
is a coincidence of measure zero, so a zero in these files reads unambiguously
as "no chemistry here" rather than as a real value.

| Quantity | What it should do |
|---|---|
| CH4, SO4 totals | Fall relative to what pure diffusion alone would give, because the aggregate is consuming them |
| HS, HCO3 totals | Rise from their small initial value; the rise is entirely the reaction |
| `FT_ANME` | Near 1 at the rim, falling sharply with depth, small but not exactly zero at the centre |
| `dG_ANME` | Negative everywhere reacting, approaching the 12.5 kJ/mol threshold from below at depth |
| `rate_CH4` | Largest near the rim, small to negligible at the centre; the shell-and-core pattern this case exists to show |
| Species ratio | CH4, SO4, HS, HCO3 increments in 1:1:1:1, from the reaction's own stoichiometry |
| Minimum of any species | Zero or above. A negative value is a real failure |

## 7. What to check

1. **The `.thm` file is echoed at start-up**, and the reaction, the threshold
   and the temperature the solver parsed are the ones you meant. Read the
   `[THM]` lines at the top of the log before anything else.
2. **The gate acted, and closed somewhere.** A mean `F_T` near 1 everywhere
   means the gate never engaged, which is numerically the same as not having
   one; a mean near 0 everywhere means the gate shut before the reaction ever
   ran and something in the energetics is very likely wrong, most often `dG0`
   quoted at the wrong temperature.
3. **The gated and ungated runs differ.** Run the case again with
   `<thermodynamics><enabled>false</enabled>` and nothing else changed; every
   difference in the fields is the gate and only the gate.
4. **The minimum `F_T` is small but not exactly zero**, for the reason given
   in section 5. Zero everywhere in the core is a sign to look at transport,
   not the reaction.
5. **No `[NEG!]` warnings**, and the four species stay in 1:1:1:1 balance; the
   gate scales the whole rate vector at once, so it cannot break the
   stoichiometry on its own.

`postprocess.py` runs checks 1 through 4 for you: it pulls the `[THM]` block
out of the log, prints a verdict on whether the gate shut everywhere, opened
everywhere, or acted as intended, and, when both a gated and an ungated log
are present, prints the field-by-field difference between them.

## 8. What this case demonstrates

**The thermodynamic gate.** A dimensionless factor between 0 and 1 that
multiplies whichever rate path an organism already uses, computed fresh every
voxel from the local composition. It composes with a compiled rate law here,
and equally with the linear-program paths of examples 09, 10, 21 and 22, none
of which had to change to accept it.

**What it deliberately leaves out.** No flow (`<Peclet>0</Peclet>`), no
abiotic reaction, and no change to the pore space. No aqueous speciation
either: the activities in the reaction quotient are plain concentrations, with
no ion pairing, no ionic-strength correction and no carbonate equilibrium
behind `HCO3`; example 04 adds that speciation solver. The consortium is one
lumped biomass pool, so the electron-transfer loss between ANME and its
sulfate-reducing partner enters as the constant `dGloss` rather than being
computed from any solved mechanism. And the effectiveness factor this run
reports is measured while the aggregate is still growing and the surrounding
water is still being drawn down, so it is not a steady-state quantity; case 20
poses the identical aggregate so that it is.

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
./scripts/setup_case.sh 19_thermodynamic_gate run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. No optional solver is needed; the thermodynamic gate is
part of the core solver, not a build flag. The plain cmake line is enough:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** `CompLaB.xml`, `input/geometry.dat` and
`input/aom.thm` are all read at start-up, so editing any of them, including
every number in section 4 above, needs no rebuild: change the file and run
again. `kinetics/defineKinetics.hh` is a C++ header compiled into the
executable, so editing the intrinsic rate law itself, including `mu_max`,
means rebuilding.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the aggregate: a sphere of biomass in open water, in a box walled on four sides. Standard library only. |
| **offline** | `offline.sh` | Runs `offline/thermo_curve.py`, which sweeps the composition path the run will follow and reports where the gate would close, before the solver runs. |
| **build** | `kinetics/defineKinetics.hh` | The ungated intrinsic rate law, compiled in. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The aggregate and the water around it. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| | `input/aom.thm` | The energetics of section 4. Read at start-up and echoed into the log. |
| **post** | `postprocess.py` | Reads the `[THM]` block from the log and `output/summary.csv`, judges whether the gate acted, and compares the gated and ungated runs. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |
| | `offline/upscale.py` | Sweeps aggregate radius and bulk sulfide, one solver run per point; not part of `run.sh`. See case 20, which is posed for this measurement. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
