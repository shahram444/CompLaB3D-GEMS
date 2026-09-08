# 04 - Equilibrium: aqueous speciation solved in every voxel, every step

## 1. The scenario

The three cases before this one move species around and react them on a
timescale. This one adds a reaction that is already finished everywhere: acid
base speciation, which equilibrates far faster than transport can carry
anything, and so is not a rate at all but a constraint the local composition
must satisfy at every instant.

The chemistry is the carbonate system, in a slab of water with acid entering
from the inlet face: bicarbonate, protons, carbonate and carbonic acid, related
by two mass-action equilibria. A reactive-transport modeller will recognise
this as the standard closure used whenever a set of fast reactions is replaced
by an equilibrium constraint instead of being integrated in time, the same
approach behind PHREEQC-style speciation, applied here inside a pore-scale
transport solver rather than as a separate post-hoc step.

Two components, total carbonate and protons, are carried as transported
fields. Two more species, carbonate and carbonic acid, are not transported at
all: they are reconstructed from the two components in every open voxel, every
step, by solving the tableau below.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
    |....................  ####  ..........................  |
HCO3|....................  ####  ..........................  |zero
2e-3|....................  ####  ..........................  |gradient
 H  |....................  ########  ......................  |zero
1e-6|                                                         |gradient
held|                                                         |
    +---------------------------------------------------------+
      no flow: Peclet = 0. HCO3 and H diffuse in from x=0
      CO3 and H2CO3 are made in place wherever HCO3 and H overlap

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition. The
other four faces are an inert wall, drawn by `preprocess.py`, because the
solver gives those faces nothing at all and a face with no condition reads
back values that nothing updates.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, four species: `HCO3`, `H`, `CO3`, `H2CO3` |
| **Diffusivity** | 5e-10 m2/s for all four species, in pore and in biofilm alike |
| **Chemistry** | aqueous speciation, `<equilibrium><enabled>true</enabled>`, solved every voxel every step |
| **Biology** | none. `<biotic_mode>false</biotic_mode>`, no biomass field allocated |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 600 advection-diffusion steps, output every 100 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **HCO3** (component) | Dirichlet, held at 2e-3 | Neumann, zero gradient | 1e-3 |
| **H** (component) | Dirichlet, held at 1e-6 | Neumann, zero gradient | 1e-7 |
| **CO3** (complex) | closed | closed | 0 |
| **H2CO3** (complex) | closed | closed | 0 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, a reservoir. *Neumann with a zero gradient*
means nothing flows across that face, an open outlet the species never
reaches. *Closed* means the same thing on both ends: a complex is made inside
the domain by speciation and stays there, since it has no boundary source or
sink of its own.

## 4. The equilibrium

Two mass-action equilibria over two components, total carbonate (`HCO3`) and
protons (`H`):

```
    HCO3-              <->   CO3(2-)  +  H+          log K = -10.33
    HCO3-  +  H+       <->   H2CO3                   log K =  +6.35
```

Written in `CompLaB.xml` as a stoichiometric matrix, one row per species, one
column per component, with a formation constant per row:

```
    components          HCO3    H       log K
    HCO3  (component)     1     0        0
    H     (component)     0     1        0
    CO3                   1    -1      -10.33
    H2CO3                 1     1      +6.35
```

A component's own row is the identity with log K 0. Every other species is
written as a combination of components, which is what lets the solver carry
only two transported fields and reconstruct the rest. In every open voxel,
at every step, the solver finds the composition that satisfies this system
exactly, given whatever the transport step just left behind: not a Newton
iteration on the full system, but continued-fraction iteration accelerated to
converge without needing a Jacobian or a good starting guess, since it is
being restarted from scratch in every voxel every step. It is iterated to a
residual tolerance or a maximum number of iterations, whichever comes first,
and the fraction of voxel-steps that failed to converge is reported once, at
the end of the run.

Aqueous speciation here is ideal and isothermal. There is no Debye-Huckel or
Davies activity correction, no temperature dependence, no gas phase, no redox
couple and no mineral saturation index. The log K values supplied are
conditional constants at whatever ionic strength and temperature they were
measured at, and it is the user's responsibility that they apply here.

## 5. What happens each step

The solver repeats this for each of the 600 steps:

1. **Stream and collide** the two component fields, `HCO3` and `H`, on their
   D3Q7 lattices, which advances diffusion by one step. There is no velocity
   field to advect them with.
2. **Apply the boundary conditions**, re-pinning `HCO3` to 2e-3 and `H` to
   1e-6 at the left face.
3. **Solve the equilibrium tableau** in every open voxel: given the local
   totals of `HCO3` and `H` after transport, find the composition of all four
   species that satisfies both mass-action equations simultaneously, and
   write `CO3` and `H2CO3` from that solution.

The stability caveat is on step 3: the solve starts fresh in every voxel every
step from whatever transport left behind, with no memory of the previous
step's solution. Before the first transport step even runs, the initial
composition itself is checked for feasibility; an infeasible start is a
tableau error, not a numerical one, and needs the tableau fixed rather than a
tighter tolerance.

## 6. What comes out

```
output/
  HCO3_0000100.vti  ...  HCO3_0000600.vti      the transported components
  H_*.vti
  CO3_0000100.vti  ...  CO3_0000600.vti         the reconstructed complexes
  H2CO3_*.vti
  summary.csv                                    one row every 100 steps
  run.log                                         the whole run, including the convergence report
```

There is no `rate_*.vti` field here. Speciation redistributes a fixed total
between species; it does not create or destroy mass, so there is no rate to
write a field for.

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Equilibrium convergence | Very close to 100 % of voxel-steps, reported once at the end of the run |
| Initial feasibility | Satisfied immediately, in a handful of iterations, before step 1 |
| Total `H2CO3` | Rises from 0 as acid drives the second equilibrium toward it |
| Total `CO3` | Not necessarily monotone: it can rise, then partly reverse as the acid front arrives and pushes carbonate toward `H2CO3` |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Noticeably slower than a transport-only case of the same size, since the tableau is solved in every open voxel every step |

## 7. What to check

1. **The convergence report is (very close to) 100 %.** Any voxel-step that
   fails leaves a composition that satisfies no equilibrium, and the error
   propagates through every later step.
2. **The initial composition is feasible**, reported separately at start-up
   before the first transport step.
3. **The complexes are formed in place.** `CO3` and `H2CO3` have closed
   boundaries, neither fed nor drained, so their totals change only by
   speciation; `postprocess.py` reports them separately for exactly that
   reason.
4. **`CO3` need not be monotone.** It can rise, then partly reverse as the
   acid front arrives and pushes carbonate toward `H2CO3`. That is the
   chemistry, not an instability.

`postprocess.py` runs check 3 for you, reporting the change in `CO3` and
`H2CO3` and whether each is monotone, and it flags any negative value.
Checks 1 and 2 are read from the convergence lines in `run.log`.

## 8. What this case demonstrates

**Aqueous speciation as a constraint solved everywhere, every step, rather
than a rate integrated in time.** It is the most expensive path in the
solver by a wide margin, precisely because it runs a full equilibrium solve in
every open voxel at every step rather than evaluating a closed-form rate.

**What it deliberately leaves out:** activity corrections, temperature
dependence, a gas phase, redox chemistry and mineral saturation, as stated in
section 4. It also leaves out flow, biology, and any change to the pore
space. Example 03 is the equivalent kinetic case, a rate law with no
equilibrium constraint; examples 13 to 15 add a mineral phase that can
precipitate or dissolve; example 19 gates a rate by the free energy actually
available for it.

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
./scripts/setup_case.sh 04_equilibrium run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. No optional solver is needed. Cases 09, 12, 21 and 22
additionally need GLPK, and cases 10 and 16 need Python with COBRApy; this one
needs neither, so the plain cmake line is enough:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** This case has no chemistry compiled in: no
`kinetics/` header ships with it, since the tableau lives entirely in
`CompLaB.xml` as data, not as code. Editing the tableau, the boundary
conditions, or `input/geometry.dat` never needs a rebuild. A rebuild is only
needed if you change the solver source itself.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| | `training/makeEquilibrium.py` | An optional helper: builds an `<equilibrium>` tableau from chemical formulas, checking that the components are independent and every row balances. Not run automatically; the tableau it would produce is already written into `CompLaB.xml`. |
| **run** | `CompLaB.xml` | What the solver reads, including the equilibrium tableau. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
