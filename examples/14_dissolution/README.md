# 14 - Dissolution: acid eating the grains away, and the pore space reopening

## 1. The scenario

Example 13 seals pore space. This one opens it. Acid is delivered from the
left face, eats into a mineral coating on the grain surfaces, and the voxels
it empties are converted back to open pore, after which the flow is
re-solved through the wider channel.

```
    CaCO3(s)  +  H(+)   ->   Ca(2+)  +  HCO3(-)
```

Calcite dissolution in an acidified pore space is the textbook case for
this: it is fast, it is strongly transport-limited, and its feedback runs
the other way from precipitation. A dissolving channel widens, which lets
more acid through, which widens it faster. That instability is why
dissolution fronts finger rather than advance flat, and why a pore-scale
model is worth running for it at all.

The mineral does not have to have precipitated first in this case:
`<initial_fill>` starts the declared grain surfaces completely solid, at
full mineral density, and they are eaten away from there.

The bicarbonate on the right of the reaction above is not carried as a
substrate in this case: the three fields are `H`, `Ca` and `calcite`, and the
rate law in `kinetics/defineAbioticKinetics.hh` writes only `Ca` and `H`. Do
not go looking for `HCO3` in the output.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
H+  |>>>>>>>>>>>>>>>>>>>>  ####  ...........................  |zero
1e-2|>>>>>>>>>>>>>>>>>>>>  ####::..........................  |gradient
held|>>>>>>>>>>>>>>>>>>>>  ########  ......................  |
    |                                                         |
    +---------------------------------------------------------+
      flow: Peclet = 1, pushes acid left -> right
      #### = solid grain block, full of calcite at t = 0
      :: = the one-voxel reactive face; acid eats it and the voxel reopens

      grain blocks alternate between the y=0 and y=25 walls every 4 voxels
      along x; the reactive face is the pore-facing layer of each block,
      72 voxels in total

      porosity 0.6346, 3168 open voxels, 72 declared reactive-phase voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition. The
other four faces are an inert wall, drawn by `preprocess.py`. The geometry is
the same staggered slot example 13 uses, and the 72 voxels of declared phase
are the same one-voxel layer on each grain's pore-facing side; the rest of
each block stays inert wall. Mineral buried inside a grain can never be
reached by water, so declaring it there would cost solve time and change
nothing, and would make the reactive surface area of the case a fiction.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346, plus 72 voxels of a declared reactive phase |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet number of 1 |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `H`, `Ca`, `calcite`. `calcite` is immobile: it never streams, only accumulates its reaction term |
| **Diffusivity** | 5e-10 m2/s for `H` and `Ca`, in pore and in biofilm alike; 0 for `calcite` |
| **Chemistry** | `kinetics/defineAbioticKinetics.hh`, compiled into the executable |
| **Biology** | none. `<biotic_mode>false</biotic_mode>`, no biomass field is allocated |
| **Geometry change** | dissolution: a mineral is consumed and its voxels reopen to pore |
| **Run length** | 2000 advection-diffusion steps, output every 200 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **H** | Dirichlet, held at 1e-2 | Neumann, zero gradient | 1e-4 |
| **Ca** | closed | closed | 0 |
| **calcite** | closed | closed | 0 in open pore; 27.1 in the declared grain surface |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, which is a reservoir. *Neumann with a zero
gradient* means nothing flows across that face, an open outlet acid never
reaches. *Closed* means the same thing on both ends: `Ca` is made inside the
domain and stays there, and `calcite` never leaves the voxel it occupies at
all.

## 4. The reaction

**Only material inside a `<phaseN>` block can dissolve.** That restriction is
deliberate and it is not a convenience.

Consider the alternative rule, that a solid voxel with no mineral left
becomes pore. An inert quartz grain has never held any mineral. It is
already below any threshold you could write. On the first step the entire
grain pack would dissolve. Declaring phases makes inert material safe by
construction rather than by a guard that someone has to remember to write.

```xml
    <phase0>
        <name>calcite_grain</name>
        <material_number>0</material_number>   which voxels in geometry.dat
        <substrate>2</substrate>               which substrate holds its inventory
        <full_density>27.1</full_density>      mol/L when completely solid
        <initial_fill>27.1</initial_fill>      starts solid
    </phase0>
```

`<full_density>` comes from the mineral's **molar volume**, the volume one
mole of the solid occupies, the same way `<max_precipRho>` does in example
13:

```
    full_density  =  1000 / Vmolar          mol L-1, Vmolar in cm3 mol-1

    calcite       100.09 g mol-1
                    2.71 g cm-3   ->   36.9 cm3 mol-1   ->   27.1 mol L-1
```

**The rate law itself.** `defineDissolutionRate()` is a separate hook from the precipitation one,
in the same header. First order in acid and first order in the mineral
remaining:

```
    R  =  k_calcite * [H+] * m                mol L-1 s-1

    k_calcite = 1.0e-2   L mol-1 s-1 per mol L-1 of remaining mineral
    m                    the mineral inventory still in that voxel

    dm/dt = -R           d[Ca2+]/dt = +R           d[H+]/dt = -R
```

The `m` factor is what makes the rate fall as a voxel empties, instead of
running at full speed until the inventory hits zero and then stopping in one
step. This case ships no thermodynamic term either: a full kinetic law would
carry a **saturation index**, the ratio of the ion activity product to the
solubility product, and run at zero once the solution reaches equilibrium
with the mineral. This one runs at whatever `[H+] * m` gives until the acid
or the mineral itself runs out.

**Which water the rate sees.** A solid voxel has no water of its own that
moves. So the rate law is handed the concentration in the water touching the
mineral, averaged over that voxel's open neighbours, not a concentration
stored in the solid voxel itself. A declared voxel with no open neighbour is
unreachable and does not dissolve at all, which is the same physical
statement as the geometry note in section 4.

| tag | value | what it controls |
|---|---|---|
| `<reopen_fraction>` | `0.9` | a voxel reopens once it falls below 0.9 times full density |
| `<surface_only>` | `1` | only wetted voxels dissolve |
| `<update_interval>` | `100` | steps between geometry re-checks and flow re-solves |

**`reopen_fraction` is a hysteresis knob.** Reopening at exactly full density
would flicker a voxel between solid and pore on the arithmetic noise of a
single step, and each flip costs a flow re-solve. 0.9 means a voxel has to
have genuinely lost a tenth of its inventory before the geometry changes.

## 5. What happens each step

The solver repeats this 2000 times:

1. **Stream and collide** `H` and `Ca` on their D3Q7 lattices, advected by
   the current velocity field and diffusing at the same time.
2. **Apply the boundary conditions**, re-pinning `H` to 1e-2 at the left
   face.
3. **Evaluate the dissolution rate** on every open voxel touching a declared
   phase: average the local water composition over the open neighbours,
   compute `R = k_calcite * [H+] * m`, and write the increments for the
   mineral and for `Ca` and `H`.
4. **Add the increments** to the concentrations and to the mineral
   inventory.
5. **Every `<update_interval>` (100) steps, check the geometry.** Any
   declared-phase voxel whose mineral inventory has fallen below
   `<reopen_fraction>` (0.9) times `<full_density>` is converted back to
   open pore. If any voxel reopened, the flow field is re-solved on the new
   geometry.

The stability caveat is the same as in example 13: the increments in step 4
are applied as computed, with no clamp against what is actually present. The
`m` factor in the rate makes an overdraw of the mineral itself unlikely, but
`H` can still be driven negative by a large `<ade_dt>` or a large
`k_calcite`; the fix is the same, a smaller one or the other.

## 6. What comes out

```
output/
  H_0000200.vti  ...  H_0002000.vti          concentration of acid
  Ca_*.vti                                    concentration of Ca2+
  calcite_*.vti                               the remaining mineral inventory, mol/L per voxel
  rate_H_*.vti  rate_Ca_*.vti  rate_calcite_*.vti   the reaction rate as a field, mol/L/s
  nsLattice_*.vti                             the velocity field, re-solved whenever the geometry changes
  summary.csv                                 one row every 100 steps, including porosity
  run.log                                     the whole run, including the checks
```

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Porosity | Rises, in steps of at most `<update_interval>` iterations, from 0.6346 |
| Total calcite | Falls from its starting inventory as the grains are eaten back |
| Total Ca | Rises from 0, one mole released for every mole of calcite lost |
| The dissolution front | Advances further along the high-flow paths than along the slow ones, not as a flat line |
| Flow re-solves | More than the initial one, appearing in the log as soon as a voxel first reopens |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Seconds to low minutes on one core. This is a small case on purpose |

## 7. What to check

1. **Porosity rises monotonically.** Only dissolution is enabled, so nothing
   can close a voxel; any fall is a bug.
2. **The calcium balance closes.** Calcium is created only by dissolution and
   destroyed by nothing, and it enters through no boundary: `total Ca == the
   calcite removed from the declared phase`. Both of calcium's faces are
   closed, genuine bounce-back, exactly no-flux, rather than zero-gradient,
   because a zero-gradient face is an open outflow that copies whatever is
   inside it outward and would quietly manufacture mass in exactly this kind
   of balance. This is the check to trust before any of the others.
3. **The front is not flat.** Dissolution should reach further along the
   high-flow paths, because those deliver more acid. A flat front means the
   flow field is not being re-solved, or transport is diffusion-dominated
   and the case is not showing what it is for.
4. **Acid is consumed one for one with calcium released.** The stoichiometry
   is 1 : 1.
5. **The mineral inventory only falls**, and no voxel goes below zero. The
   rate's `m` factor should make that impossible before the clamp is needed.
6. **No `[NEG!]` warnings.**

**Run this case on one process.** A dissolving voxel deposits its products
into its open face neighbours, and where a neighbour lies in a different MPI
block, the deposit is written into that block's local envelope and never
communicated back. The result then depends on the processor count, with the
error concentrated at block interfaces. If the balance in check 2 fails on
several ranks and passes on one, that is this limitation and not the
chemistry.

`postprocess.py` runs the calcium balance of check 2 for you and reports it
as PASS or FAIL, and it flags any negative minimum, covering check 6. The
porosity trend, the front shape and the stoichiometry are read from
`summary.csv` and the `.vti` output. This case declares no `<conserve>` tag
in `CompLaB.xml`; `H` is fed from a Dirichlet boundary, so its total is not
conserved and naming it in a mass-balance check would guarantee a failure
that means nothing.

## 8. What this case demonstrates

**Geometry that reopens because the chemistry consumed what was blocking
it.** A declared solid phase is not permanent: falling below
`<reopen_fraction>` of `<full_density>` converts a voxel back to pore, and
the flow is re-solved through the wider channel. Porosity itself becomes an
output that moves during the run, in the opposite direction from example 13.

**What it deliberately leaves out.** Nothing precipitates here: a voxel that
reopens stays open. Example 15 enables both paths at once, which is where
the two have to agree about who owns a voxel and what happens when a voxel
is being filled and emptied in the same interval. The acid is a bare `H`
field with no carbonate system behind it, no bicarbonate, no CO2, no pH
buffering, so the front advances faster than a real carbonate system would
allow; example 04 adds the aqueous speciation solver that would supply the
buffering.

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
./scripts/setup_case.sh 14_dissolution run/mycase
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
`defineAbioticKinetics.hh`, including changing `k_calcite`, means rebuilding.
Editing `CompLaB.xml`, including `<reopen_fraction>` or `<update_interval>`,
or `input/geometry.dat`, does not.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the staggered slot and coats the pore-facing grain surfaces with the declared phase, refuses to continue if the slot does not percolate. Standard library only. |
| **build** | `kinetics/defineAbioticKinetics.hh` | This case's chemistry, compiled in. Its `defineAbioticRxnKinetics()` is present but returns zero: this case runs no homogeneous reaction, and an unused entry point still has to link. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and checks the calcium balance of section 7. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top; the
kinetics header above then replaces the shared default it laid down, so the
law described in section 4 is what gets compiled in.
