# 15 - Precipitation and dissolution: a mineral that seals a voxel and can be dissolved back out of it

## 1. The scenario

Examples 13 and 14 each run the pore space one way: 13 only seals voxels, 14
only reopens them. This one runs both directions on the same mineral, so a
voxel can precipitate shut, sit sealed while the flow reroutes around it, and
later be dissolved back open by acid arriving down a different path.

```
    precipitation    Fe(2+)  +  HS(-)      ->  FeS(s)
    dissolution      FeS(s)  +  H(+)       ->  Fe(2+)  +  HS(-)
```

The precipitation line produces no proton: the shipped
`kinetics/defineAbioticKinetics.hh` writes only the iron, the sulphide and
the mineral in its precipitation hook, and never touches the proton entry;
the reaction as configured in `CompLaB.xml` is written the same way. This is
deliberate, and it is why the two rate laws are not inverses of one another,
and why acid, once delivered, only ever falls.

Nothing is mineral at the start. Iron comes in from the left, sulphide from
the right, acid from the left; FeS forms where iron and sulphide overlap, and
is attacked wherever acid reaches it. Porosity is no longer a one-way trend,
which is the whole reason to run this case rather than the two before it.

**The two rate laws are not inverses of one another**, and this is the
design point worth understanding before reading the output. Precipitation
here is second order in the dissolved ions; dissolution is first order in
acid and in the mineral remaining. They are two separate laws that happen to
move the same mineral. That is how this code is meant to be used: there is
no solubility product anywhere in it and therefore no **saturation index**,
the ratio of the ion activity product to the solubility product that says
how far a solution sits from equilibrium, so there is nothing to drive a
single reversible rate. If detailed balance, the forward and back rates
equal at equilibrium, matters for your problem, you have to impose it
yourself, in the two laws you write.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
Fe2+|>>>>>>>>>>>>>>>>>>>>  ####  <<<<<<<<<<<<<<<<<<<<<<<<<<<<< |HS-
H+  |>>>>>>>>>>>>>>>>>>>>  ####::.....................<<<<<<<< |5.0
held|>>>>>>>>>>>>>>>>>>>>  ########  <<<<<<<<<<<<<<<<<<<<<<<<< |held
    |                                                         |
    +---------------------------------------------------------+
      flow: Peclet = 1, pushes Fe2+, H+ and HS- as marked
      #### = solid grain block   :: = the reactive face: seals under FeS,
             later reopens if acid reaches it faster than the seal reforms

      grain blocks alternate between the y=0 and y=25 walls every 4 voxels
      along x; the reactive face is 72 voxels of the 3168 open voxels

      porosity starts at 0.6346 and moves BOTH ways during the run
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition. The
other four faces are an inert wall, drawn by `preprocess.py`. The geometry is
the same staggered slot examples 13 and 14 use; nothing here is mineral at
the start, so the reactive face begins as open pore rather than solid, and
becomes solid only once a voxel there seals.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346, plus 72 voxels of a declared reactive phase |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet number of 1 |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `Fe2`, `HS`, `H`, `FeS`. `FeS` is immobile: it never streams, only accumulates its reaction term |
| **Diffusivity** | 5e-10 m2/s for `Fe2`, `HS` and `H`, in pore and in biofilm alike; 0 for `FeS` |
| **Chemistry** | `kinetics/defineAbioticKinetics.hh`, compiled into the executable |
| **Biology** | none. `<biotic_mode>false</biotic_mode>`, no biomass field is allocated |
| **Geometry change** | both directions: precipitation fills and seals a voxel, dissolution consumes it back open |
| **Run length** | 3000 advection-diffusion steps, output every 200 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **Fe2** | Dirichlet, held at 5. | Neumann, zero gradient | 0 |
| **HS** | Neumann, zero gradient | Dirichlet, held at 5. | 0 |
| **H** | Dirichlet, held at 1e-2 | Neumann, zero gradient | 0 |
| **FeS** | closed | closed | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, which is a reservoir. *Neumann with a zero
gradient* means nothing flows across that face, an open outlet the species
never reaches. *Closed* means the same thing on both ends: `FeS` is made and
consumed only inside the domain.

## 4. The reaction

Both live in `kinetics/defineAbioticKinetics.hh`, in the two separate hooks
the solver calls.

**Precipitation**, `defineAbioticRxnKinetics()`, a homogeneous reaction
between two dissolved species:

```
    Rp  =  kp_FeS * [Fe2+] * [HS-]            mol L-1 s-1        kp_FeS = 1.0 L mol-1 s-1

    d[Fe2+]/dt = -Rp      d[HS-]/dt = -Rp      d[FeS]/dt = +Rp
```

**Dissolution**, `defineDissolutionRate()`, a surface reaction handed the
concentration in the water touching the mineral, because a solid voxel has
no water of its own:

```
    Rd  =  kd_FeS * [H+] * m                  mol L-1 s-1        kd_FeS = 5.0e-2 L mol-1 s-1

    dm/dt = -Rd     d[Fe2+]/dt = +Rd     d[HS-]/dt = +Rd     d[H+]/dt = -Rd
```

`m` is the mineral inventory left in that voxel, so the rate falls as the
voxel empties rather than running flat out and stopping in one step. Because
a proton is consumed by dissolution and none is supplied by precipitation in
this header, acid is a genuine reactant here, delivered from the left face
and consumed as it works. Where it runs out, dissolution stops and whatever
has sealed stays sealed.

**How the two paths share a voxel.** The two configuration blocks,
`<precipitation>` and `<dissolution>`, are independent, which is why both can
be on at once. What connects them is `<phase0>`:

```xml
    <phase0>
        <name>FeS_precipitate</name>
        <material_number>0</material_number>
        <substrate>3</substrate>
        <full_density>48.9</full_density>
        <initial_fill>0</initial_fill>          nothing there at the start
        <is_precipitate>true</is_precipitate>   voxels sealed by precipitation
    </phase0>
```

`<is_precipitate>true</is_precipitate>` is the link. A voxel sealed by the
precipitation path is stamped as belonging to this phase, which is what
makes it eligible for the dissolution path later. Without that stamp the two
halves would not know they were operating on the same mineral, and a sealed
voxel would be permanent, as it is in example 13.

`48.9 mol L-1` is mackinawite's molar density, and it appears in both
blocks for the same reason: `<max_precipRho>` is the density a voxel must
reach to seal, and `<full_density>` is the density it is measured against
when reopening. They are the same physical number and must be written the
same in both places:

```
    max_precipRho = full_density = 1000 / Vmolar

    mackinawite   87.91 g mol-1 / 4.30 g cm-3  =  20.44 cm3 mol-1  ->  48.9 mol L-1
```

| tag | value | what it controls |
|---|---|---|
| `<solid_substrate>` | `3` | which substrate index holds the mineral |
| `<max_precipRho>` / `<full_density>` | `48.9` | the same mineral density, in both blocks |
| `<surface_only>` | `1` (both) | growth and attack only on wetted, wall-adjacent voxels |
| `<perm_ratio>` | `0` | a sealed voxel is an impermeable wall |
| `<reopen_fraction>` | `0.9` | the hysteresis gap described below |
| `<update_interval>` | `100` (both) | steps between geometry re-checks and flow re-solves |

**`<reopen_fraction>` earns its keep here and nowhere else in this set of
cases.** With both directions active a voxel can sit right at the sealing
threshold. On a single threshold it would flip open and shut every update
interval, and every flip re-solves the flow field. `0.9` puts a gap between
sealing at full density and reopening below 0.9 times full, so a marginal
voxel settles instead of chattering. Porosity should move in steps and then
hold, not oscillate at the update interval.

## 5. What happens each step

The solver repeats this 3000 times:

1. **Stream and collide** `Fe2`, `HS` and `H` on their D3Q7 lattices, advected
   by the current velocity field and diffusing at the same time.
2. **Apply the boundary conditions**, re-pinning `Fe2` and `H` at the left
   face and `HS` at the right.
3. **Evaluate the precipitation rate** in every open voxel: `Rp = kp_FeS *
   Fe2 * HS`, and write the increments for `Fe2`, `HS` and `FeS`.
4. **Evaluate the dissolution rate** on every open voxel touching a stamped
   `FeS_precipitate` voxel: `Rd = kd_FeS * H * m`, using the acid
   concentration in the water touching the mineral, and write the increments
   for `H`, `Fe2`, `HS` and the mineral inventory.
5. **Add the increments** to the concentrations and to the mineral
   inventory.
6. **Every `<update_interval>` (100) steps, check the geometry both ways.**
   A voxel whose `FeS` inventory has reached `<max_precipRho>` is sealed and
   stamped `is_precipitate`; a stamped voxel whose inventory has fallen below
   `<reopen_fraction>` of `<full_density>` is reopened. If any voxel changed
   state, the flow field is re-solved on the new geometry.

The stability caveat is sharper here than in examples 13 or 14: two rate
laws are drawing on the same mineral inventory in the same interval, which
is the situation most likely to overdraw it. The increments in step 5 are
applied as computed, with no clamp against what is actually present, so a
time step or rate constant too large produces a `[NEG!]` warning exactly as
in the earlier cases.

## 6. What comes out

```
output/
  Fe2_0000200.vti  ...  Fe2_0003000.vti     concentration of Fe2+
  HS_*.vti  H_*.vti                          the other two dissolved species
  FeS_*.vti                                   the mineral inventory, mol/L per voxel
  rate_Fe2_*.vti  rate_HS_*.vti  rate_H_*.vti  rate_FeS_*.vti   the reaction rate as a field, mol/L/s
  nsLattice_*.vti                             the velocity field, re-solved whenever the geometry changes
  summary.csv                                 one row every 100 steps, including porosity
  run.log                                     the whole run, including the checks
```

Numbers to expect, as orders of magnitude and directions rather than exact
values, **with the shipped rate constants** (see the sensitivity note in
section 7):

| Quantity | What it should do |
|---|---|
| Porosity | Moves both ways from 0.6346, in steps of at most `<update_interval>` iterations, and settles rather than oscillating |
| Total FeS | Rises from 0, may peak once dissolution catches up with precipitation, and can fall afterward |
| Total Fe2, total HS | Both fed and consumed; net direction depends on which process dominates locally |
| Total H (acid) | Falls whenever dissolution is running, since nothing in this header replaces it |
| Flow re-solves | More than one, from both a seal and a reopening |
| Minimum of any species | Zero or above at a well-chosen rate constant; see section 7 for what a poorly chosen one looks like |
| Wall clock | Seconds to low minutes on one core. This is a small case on purpose |

## 7. What to check

1. **Porosity moves both ways, and settles.** It should fall as the mineral
   forms and rise where acid reaches it, in steps at the update interval,
   and then hold. Oscillation with a period of exactly `<update_interval>`
   is the chattering `<reopen_fraction>` exists to prevent, and is the
   failure this case is built to expose.
2. **Iron and sulphide balance against the mineral.** Both are 1 : 1 with
   FeS in both directions, so the sum of what is dissolved and what is
   mineral changes only through the boundaries.
3. **Acid only ever falls.** It is consumed by dissolution and produced by
   nothing in this header.
4. **Mineral appears before any of it dissolves.** `<initial_fill>` is 0, so
   at step 0 there is nothing to attack; a dissolution rate reported before
   any precipitation has occurred means a phase is being stamped that
   should not be.
5. **The flow re-solves more than once.** Both a seal and a reopening
   trigger one.
6. **No `[NEG!]` warnings.**

`postprocess.py` reports the porosity trend, the change in every field's
total, and flags any negative minimum, covering checks 2, 3 and 6. It also
runs the closed-species report on `FeS`. It reports only the **net** change
in porosity between the first and last row of `summary.csv`: a run that
oscillates and returns to its starting value prints "unchanged" there, which
is why check 1 has to be read from the porosity column itself, plotted
against iteration, rather than taken from the printed verdict. Checks 4 and
5 are read from the log.

**A note on the shipped rate constant.** `kinetics/defineAbioticKinetics.hh`
ships `kp_FeS = 1.0`, the same value example 13 uses for precipitation
alone. With both directions active in the same voxel, that value is large
enough to drive `[NEG!]` warnings here and fail check 6: precipitation can
draw a voxel's `Fe2` or `HS` down faster than the explicit update in step 5
tolerates once dissolution is also feeding the same reservoir back in on the
next interval. Turning `kp_FeS` down to 0.1 removes the warnings but never
lets a voxel reach `<max_precipRho>` at all, so nothing ever seals and
dissolution never has anything to act on; the run is stable and shows
nothing. A value around 0.3 is the middle ground: it seals voxels, produces
no `[NEG!]` warnings, and lets `FeS` peak above the sealing threshold before
declining monotonically as dissolution catches up. If you change `kp_FeS`,
this is the tradeoff you are on: too high overdraws the reservoir in a
single interval, too low never triggers the geometry change this case
exists to show. Do not edit the shipped header to "fix" this without
knowing which behaviour you want; a copy in a scratch case is the place to
experiment.

## 8. What this case demonstrates

**Two independent geometry-changing rate laws sharing one mineral, and the
hysteresis needed to keep them from fighting over a single voxel every
update interval.** Porosity is not a monotone quantity here the way it is in
examples 13 and 14; it is a genuine output of the run, and which direction
it moves in any given region depends on which front, iron and sulphide or
acid, gets there first.

**What it deliberately leaves out.** No biology and no aqueous speciation:
`Fe2`, `HS` and `H` are total dissolved concentrations, not free ions, so
there is no FeHS+ complex, no sulphide protonation equilibrium, and no
ionic-strength correction, all of which control where real mackinawite
forms. Example 04 adds the speciation solver. And, as section 1 says, the
forward and back laws are independent: they do not share a solubility
product, so the mineral in this case does not have a well-defined
equilibrium, only whichever steady state the two rate constants and the
transport happen to produce.

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
./scripts/setup_case.sh 15_precip_and_dissolution run/mycase
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
`defineAbioticKinetics.hh`, including `kp_FeS` or `kd_FeS`, means
rebuilding. Editing `CompLaB.xml`, including `<reopen_fraction>` or
`<update_interval>`, or `input/geometry.dat`, does not.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the staggered slot, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **build** | `kinetics/defineAbioticKinetics.hh` | This case's chemistry, both hooks implemented, compiled in. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and reports the net change in porosity. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top; the
kinetics header above then replaces the shared default it laid down, so both
laws described in section 4 are what get compiled in.
