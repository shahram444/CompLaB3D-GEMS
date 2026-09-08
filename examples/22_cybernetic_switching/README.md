# 22 - Cybernetic switching: an organism that changes carbon source when the first runs out

## 1. The scenario

Case 21 fixes which member of an optimal set a linear program returns. This
case fixes a different problem: which substrate the organism is actually
living on, when more than one carbon source will support growth.

Offered glucose and acetate together, E. coli does not eat both at once. It
eats the glucose, excreting acetate as a by-product while doing so, and only
turns around and starts consuming that same acetate once the glucose is
gone. A microbiologist recognises this immediately as carbon catabolite
repression, the general phenomenon by which a cell preferentially consumes
one available carbon source and keeps the machinery for a second one
switched off until the first runs out. Flux balance analysis on its own
cannot reproduce it: handed both carbon sources at once, the linear program
of case 09 simply takes whichever combination of both maximises growth, from
the very first step, with no preference for one over the other. A Monod term
cannot reproduce it either, and for a sharper reason: acetate is a product in
the glucose-eating regime and a substrate in the acetate-eating regime, and a
concentration on its own does not say which regime a voxel is in.

`<cybernetic>` treats the organism as a competition between separate growth
options, one full linear-program solution per available carbon source, and
blends them by how much carbon each option would actually bring in given the
concentrations sitting in that voxel right now. Glucose dominates the blend
while it is present; as it depletes locally, the blend continuously shifts
toward the acetate solution, no threshold and no extra parameter required
beyond the same kinetic rates a Monod law would already need.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um
   (the pore space itself is 24 x 24 x 6; preprocess.py adds a one voxel
   wall on the y and z faces the solver conditions with nothing, so NY
   and NZ come out two larger than the pore space it declares)

   x=0                                                          x=23
    |                                                             |
glc |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> |
o2  |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>> |
1.0 |>>>  ====                                          >>>>>>>> | Neumann,
    |               ####    ####    ####    ####    ####         | open outlet
    |               ====             ac starts at 0.4 EVERYWHERE  |
    +-------------------------------------------------------------+
      near the inlet: glucose plentiful, acetate excreted (u_glc near 1)
      deep in the biofilm: glucose scarce, acetate consumed (u_ac near 1)
      that spatial split is the pore-scale version of the switch

      #### = solid grain (wall)     ==== = initial Ecoli patch
      porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry boundary conditions. The
other four faces are an inert wall, drawn by `preprocess.py`, for the same
reason as every other case in this set: a face the solver conditions with
nothing streams a lattice Boltzmann field off the edge of the block and reads
back values nothing updates.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `glc`, `o2`, `ac` |
| **Chemistry** | none written by hand. Growth and every exchange flux come from the metabolic model |
| **Biology** | `Ecoli`, one organism, `<reaction_type>glpk</reaction_type>`, `e_coli_core.xml` read natively as SBML, attached biofilm spreading by finite-difference diffusion, two competing growth options per voxel |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 500 advection-diffusion steps, output every 100 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **glc** (donor) | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **o2** (acceptor) | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **ac** (product and substrate) | closed | Neumann, zero gradient | 0.4 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, which is a reservoir. *Neumann with a zero
gradient* means nothing is forced to flow across that face, an open outlet.
*Closed* on the left means acetate has no external source; the whole 0.4
starting concentration is already in the domain when the run begins, so both
growth options have something to work with from the first step, rather than
the switch waiting on acetate to first be produced from nothing.

## 4. The reaction

`input/e_coli_core.xml` is the same E. coli core metabolic model as case 21:
72 metabolites, 95 reactions, with `EX_glc__D_e`, `EX_o2_e` and `EX_ac_e` as
the three exchange reactions this case watches. Rather than one linear
program handed both carbon sources at once, the organism is solved once per
source, each time with the other source's uptake shut:

```
    glucose only:   max biomass,  ac uptake forbidden   ->  v_glc, growth mu_glc
    acetate only:   max biomass,  glc uptake forbidden   ->  v_ac,  growth mu_ac
```

Each option is scored by how much carbon it would actually bring in at the
local concentrations, using an ordinary Michaelis-Menten kinetic rate for
each source:

```
    r_k  =  k_k C_k / (K_k + C_k)               the unregulated kinetic rate
    u_k  =  n_k r_k / SUM_j n_j r_j              the cybernetic weight
```

`n_k` is the number of carbon atoms the source carries, 6 for glucose and 2
for acetate, which is what makes the competition genuinely about how much
carbon is coming in rather than about which molecule is more concentrated.
The weights `u_k` are non-negative and always sum to one. The two flux
vectors are then blended directly, one full-voxel solution:

```
    v  =  SUM_k  u_k v_k                          growth, and every exchange flux
```

**Shutting a source forbids only its uptake, never its release.** An organism
solved as running on glucose alone must still stay free to export the acetate
that the acetate-only option will later consume; closing that exchange
outright would remove the very thing being switched to, and the second regime
would never have anything to work with once glucose ran out. In the
configuration this is the lower bound on the shut source's uptake being
raised to zero, with its upper bound, the export side, left untouched.

**The blend conserves mass exactly, not approximately.** Every flux vector
that solves the linear program satisfies `S v = 0`, the internal metabolites
balance, and a weighted average of vectors that each satisfy that equation
still satisfies it, because the equation is linear. Averaging the output of
two independently solved models is not generally a safe thing to do; here it
is, and that is the reason.

The competition itself is set in `CompLaB.xml`:

```xml
<cybernetic>
    <sources>                 glucose  acetate </sources>
    <substrate_ids>              0        2    </substrate_ids>
    <carbon_number>              6        2    </carbon_number>
    <uptake_kmax>              10.0      4.4   </uptake_kmax>
    <uptake_half_saturation>    0.05     0.05  </uptake_half_saturation>
    <weight_floor>            1e-3            </weight_floor>
</cybernetic>
```

`<weight_floor>` skips solving any option whose weight has fallen below it
that step, which is most of the saving described in section 5; setting it to
0 solves every source every time and is the reference case.

## 5. What happens each step

The solver repeats this 500 times:

1. **Stream and collide** each of the three species on its own D3Q7 lattice,
   and diffuse the biomass field the same way.
2. **Apply the boundary conditions**, re-pinning glucose and oxygen to 1.0 at
   the left face.
3. **Build the uptake bounds** for each substrate from the local
   concentration, as in case 09.
4. **Score each carbon source**, computing `r_k` and `u_k` from section 4 at
   the local glucose and acetate concentrations, skipping any source whose
   weight is below `<weight_floor>`.
5. **Solve one linear program per surviving source**, each with every other
   source's uptake bound forced to zero. Because the two programs differ only
   in which bounds are open, each solve warm-starts from the previous
   solve's basis and finishes quickly.
6. **Blend the flux vectors** by their weights, write the blended growth and
   exchange fluxes into the change lattices, scaled by the local biomass and
   the time step.
7. **Add the increments**, clamped so a step cannot draw a voxel's substrate
   below what it actually holds.

**Where the weights sum to zero rather than one.** In a voxel with no carbon
present at all, both `r_k` values are zero, both weights are zero, and the
organism is correctly skipped entirely for that step. That is the right
answer, not a failed solve, and should not be counted among infeasible
programs when checking the log.

## 6. What comes out

```
output/
  glc_0000100.vti  ...  glc_0000500.vti     concentration of glucose
  o2_*.vti  ac_*.vti                         the other two species
  Ecoli_*.vti                                the biomass field
  rate_glc_*.vti  rate_o2_*.vti  rate_ac_*.vti   the reaction rate as a field, mol/L/s
  summary.csv                                one row every 50 steps
  run.log                                    the whole run, including the checks
```

| Quantity | What it should do |
|---|---|
| Total glc, total o2 | Fall relative to what pure diffusion alone would give, because the organism is consuming them |
| `rate_ac` near the inlet | Positive: glucose is plentiful there, so the glucose-only option dominates and acetate is excreted |
| `rate_ac` deep in the biofilm | Negative: glucose is scarce there, so the acetate-only option dominates and acetate is consumed |
| Total biomass | Rises wherever either carbon source is available in useful amounts |
| Cybernetic weights, `u_glc + u_ac` | Exactly 1 wherever any carbon is present, 0 and 0 together where none is |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Close to twice a single-stage FBA case of the same size, minus whatever `<weight_floor>` skips |

## 7. What to check

1. **Look for the sign change in `rate_ac` across the biofilm.** Positive
   near the inlet, negative deep inside, is the switch actually happening
   spatially. If it is positive everywhere, glucose is reaching everywhere
   and the domain is too small, or the biomass too sparse, to develop the
   gradient the case is built to show.
2. **The two weights sum to one wherever carbon is present.** Both zero
   together, with the organism skipped, is the correct answer when no carbon
   is present at all and should not be mistaken for a failed solve.
3. **Mass should close**, the same `summary.csv` check as case 09; if it does
   not, the blend is being applied to fluxes and growth taken from
   inconsistent solves.
4. **Combining with case 21.** `<multi_step>` and `<cybernetic>` compose:
   each of the two options here can itself be solved as a lexicographic
   chain rather than a single program, at the cost of the two stage counts
   multiplying together, two sources times three stages is six linear
   programs per voxel per step, so try any such combination on a small
   domain first.

`postprocess.py` reports every field's total change, flags any negative
value, and scans the log for infeasible or non-converging solves, which
covers checks 2 and 3 for you; check 1 needs a look at the `rate_ac` field
itself.

## 8. What this case demonstrates

**Cybernetic switching**, a carbon-weighted blend of separately solved
per-substrate optima, which is how an organism's preference for one carbon
source over another, and the delay before it switches to a second one, is
represented without adding any threshold or extra state.

**What it deliberately leaves out.** No flow, no abiotic reaction, no change
to the pore space, and no thermodynamic gate; case 19 is the case built
around that. It also leaves out any memory of past exposure: the weights here
depend only on the concentrations sitting in a voxel right now, with no
enzyme level carried from one step to the next, so this reproduces the
spatial pattern of switching but not the temporal lag and flat growth plateau
a real diauxic culture shows in a well-mixed flask while it re-tools its
enzymes. That simplification is adequate exactly when the switch is fast
compared with transport, which is the pore-scale regime this case is built
for, and would not be adequate for reproducing a batch-culture growth curve.

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
./scripts/setup_case.sh 22_cybernetic_switching run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
Palabos v2.3.0, and GLPK's development package (`libglpk-dev` on
Debian/Ubuntu, `glpk` from Homebrew on macOS). This case's `<simulation_mode>`
block sets `<enable_fba_glpk>true</enable_fba_glpk>`, which is the CMakeLists
flag `-DENABLE_GLPK=ON`:

```bash
cmake -B build -S . \
      -DCMAKE_BUILD_TYPE=Release \
      -DPALABOS_ROOT=/path/to/palabos-v2.3.0 \
      -DENABLE_GLPK=ON
cmake --build build -j
```

**When you must recompile.** `-DENABLE_GLPK=ON` only makes the GLPK solver
available inside the executable; `<enable_fba_glpk>` in `CompLaB.xml` is what
switches it on for a given run, and that line, along with every source's
carbon number, kinetic parameters and `<weight_floor>` in `<cybernetic>`, can
be edited without a rebuild. A build without the CMake flag compiles cleanly
and then refuses at start-up, because turning the XML path on is not enough
if GLPK was never linked in. Editing `input/e_coli_core.xml` or
`input/geometry.dat` also needs no rebuild; only a change to `CMakeLists.txt`
or the solver sources does.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Shared with case 09. Standard library only. |
| **run** | `CompLaB.xml` | What the solver reads, including the `<cybernetic>` block. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| | `input/e_coli_core.xml` | The metabolic model of section 4, read natively as SBML. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives and infeasible or non-converging solves. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
