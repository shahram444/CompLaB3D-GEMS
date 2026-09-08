# 11 - Surrogate model: a fitted network in place of the linear program

## 1. The scenario

Flux balance analysis, FBA for short, answers the right question but pays a
linear program for it: an optimisation problem solved in every voxel, at
every step. That is the most expensive rate path in this repository. This
case replaces the linear program with a surrogate, a small feed-forward
neural network fitted offline to reproduce the linear program's answers over
the range of inputs it will actually see, so the run pays a few hundred
multiplications per voxel per step instead of an optimisation solve.

The organism is *Geobacter metallireducens* respiring acetate with Fe(III)
(ferric iron) as the terminal electron acceptor, the reaction that dominates
microbial iron reduction in anoxic sediments and is a standard example in
subsurface biogeochemistry:

```
    CH3COO(-)  +  8 Fe(III)  +  4 H2O   ->   2 HCO3(-)  +  8 Fe(II)  +  9 H(+)
```

Neither the growth rate nor the uptake fluxes are written down anywhere in
this case. They are what the network returns when handed the two local
uptake bounds, and the network returns them because it was fitted to a sweep
of the linear program over exactly those two bounds.

The single thing to understand before running this: a fitted network is
valid only inside the region it was trained on, and outside that region it
does not fail, it extrapolates, confidently, and returns a number that looks
like every other number. The shipped network was trained over

```
    acetate    0.00089   to  9.998    mmol gDW-1 h-1
    Fe(III)    0.000035  to  0.4997   mmol gDW-1 h-1
```

so this case sets `<fba_maximum_uptake_flux>` to `8.` and `0.4`, safely
inside. Raise the second above 0.5 and the run leaves the training region; it
will still print results, and they will not mean anything.
`training/inspectSurrogate.py` recovers those bounds from any network by
inverting its scaling transform, and the evaluator clamps to the region and
counts how often it had to.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um
   (the pore space itself is 24 x 24 x 6; preprocess.py adds a one voxel
   wall on the y and z faces the solver conditions with nothing, so NY
   and NZ come out two larger than the pore space it declares)

   x=0                                                              x=23
    |                                                                 |
acet |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>>>> |
Fe3  |>>>  ####    ####    ####    ####    ####    ####  >>>>>>>>>> | Neumann,
     |>>>  ====                                          >>>>>>>>>> | open outlet
     |               ####    ####    ####    ####    ####           |
     |               ====                                           |
     +-----------------------------------------------------------------+
      grain blocks alternate between the y=0 and y=23 walls every
      four voxels along x, so diffusion has to weave to get through

      #### = solid grain (wall)     ==== = initial Geobacter patch
      porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry boundary conditions. The
other four faces are an inert wall, drawn by `preprocess.py`, because the
solver gives a face with no declared condition nothing at all, and a lattice
Boltzmann field with nothing conditioning it streams off the edge of the
domain and reads back values that nothing updates. This is the same
geometry examples 09 and 10 use, so the surrogate can be held against the
linear program it stands in for.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, one lattice per species: `acetate`, `Fe3` |
| **Diffusivity** | 5e-10 m2/s for both species, in pore and in biofilm alike |
| **Chemistry** | none written by hand. `<enable_kinetics>false</enable_kinetics>`, `<enable_abiotic_kinetics>false</enable_abiotic_kinetics>` |
| **Biology** | `Geobacter`, one organism, `<reaction_type>surrogate</reaction_type>`, attached biofilm seeded on its own material number, spreading by finite-difference diffusion |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 500 advection-diffusion steps, output every 100 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **acetate** | Dirichlet, held at 1.0 | Neumann, zero gradient | 0 |
| **Fe3** | Dirichlet, held at 0.5 | Neumann, zero gradient | 0.5 |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, which is a reservoir. *Neumann with a zero
gradient* means nothing flows across that face, an open outlet the species
never actually reaches in this run.

## 4. The reaction

There is no stoichiometric matrix in this case at run time, only a network
that was fitted to one offline. `src/surrogateModel.hh` is a small
feed-forward network written out as plain C++, taking the two local uptake
bounds and returning a growth rate:

```
    2 inputs                 the two uptake bounds, mmol gDW-1 h-1
      -> mapminmax           scale each input to -1 .. +1 using the
                             training minimum and maximum
      -> 4 hidden layers     10 tansig (hyperbolic-tangent) neurons each
      -> 1 linear output     no squashing, so the output is unbounded
      -> mapminmax reverse   back to a growth rate
```

Every layer is one matrix multiply, one bias add and one `tanh`. The header
writes them out one line per layer, weights and all, so the forward pass is
readable rather than a call into a library:

```
    a1 = tansig( b1 + IW1_1 * x )
    a2 = tansig( b2 + LW2_1 * a1 )
    a3 = tansig( b3 + LW3_2 * a2 )
    a4 = tansig( b4 + LW4_3 * a3 )
    y  =          b5 + LW5_4 * a4        the last layer is LINEAR
```

The uptake bounds handed in are built exactly as they are on the
flux-balance path, a Michaelis-Menten term, the same half-saturation form
used to describe enzyme-limited uptake, on the local concentration, capped
by `<fba_maximum_uptake_flux>`:

```
    V  =  Vmax * C / ( Kc + C )
```

so the only difference between this case and example 09 is who turns those
two bounds into a growth rate: a linear program there, four matrix
multiplies and a `tanh` here.

## 5. What happens each step

The solver repeats this 500 times:

1. **Stream and collide** each species on its own D3Q7 lattice, advancing
   diffusion by one step, and diffuse the biomass field the same way.
2. **Apply the boundary conditions**, so acetate and Fe3 are re-pinned at
   the left face.
3. **Build the uptake bounds**, per voxel, per substrate: the same
   Michaelis-Menten term the flux-balance path uses, capped by
   `<fba_maximum_uptake_flux>`.
4. **Evaluate the network**, in every voxel that holds biomass: four matrix
   multiplies, four `tanh` calls, one linear output. No optimisation solve
   happens on this path at all, which is the entire point; where example 09
   would call GLPK here, this case does a few hundred floating-point
   multiplications instead.
5. **Read the growth rate off the network's output**, scale by the local
   biomass and the time step through `<biomass_molar_mass>`, and write the
   increment into the change lattice. The uptake fluxes for acetate and Fe3
   are then apportioned from that growth rate by the reaction
   stoichiometry, since this network was fitted on growth alone.
6. **Add the increments**, clamped so a step cannot draw a voxel's substrate
   below what it actually holds; a surrogate has no stoichiometric matrix
   enforcing non-negativity on its own, so this clamp is the only thing
   standing between it and a negative concentration.

## 6. What comes out

```
output/
  acetate_0000100.vti  ...  acetate_0000500.vti     concentration of acetate
  Fe3_*.vti                                          the electron acceptor
  Geobacter_*.vti                                    the biomass field
  rate_acetate_*.vti  rate_Fe3_*.vti                 the reaction rate as a field, mol/L/s
  summary.csv                                        one row every 50 steps
  run.log                                            the whole run, including the checks
```

| Quantity | What it should do |
|---|---|
| Total acetate, total Fe3 | Both **fall** relative to what a Dirichlet supply alone would give, because the organism is consuming them |
| Total biomass | Rises where both substrates are present, flattens or falls where either runs out |
| Acetate : Fe3 consumption ratio | Somewhere near the 1 : 8 ratio of the reaction in section 1, wherever neither bound is limiting |
| Out-of-training-box clamp count (in `run.log`) | Zero or small. A large count means the run has wandered outside the region the network was fitted on |
| Minimum of any species | Zero or above. A negative value is a real failure |
| Wall clock | Markedly faster than example 09's GLPK run on the same geometry, since no optimisation solve happens at all |

## 7. What to check

1. **The run stays inside the training region.** The evaluator clamps to it
   and counts how often; a large count means the numbers should not be
   read.
2. **Growth responds to both substrates.** Lower the acetate boundary and
   growth should fall; lower Fe3 and it should fall too. A network that
   responds to only one input was fitted on a sweep where the other never
   bound.
3. **Acetate and Fe3 are drawn together**, near the 1 : 8 ratio of the
   reaction above, wherever neither bound is limiting.
4. **No `[NEG!]` warnings.** A surrogate enforces no stoichiometry on its
   own, so the increment is clamped against what is locally present before
   it is applied; that clamp firing often is itself a signal.
5. **Against example 09.** The surrogate is only worth its speed if it
   agrees with the linear program it replaced, over the region it was
   fitted on.

`postprocess.py` reads `output/summary.csv` and `run.log`, reports every
field's change, flags any negative value, and prints the log lines
mentioning the surrogate path, so checks 1 and 4 are largely done for you.
Checks 2, 3 and 5 need a look at the fields or a comparison across runs,
which is yours to do.

## 8. What this case demonstrates

**A trained network standing in for an optimisation solve, at a fraction of
the cost, and the one discipline that makes that safe: staying inside the
region it was trained on.** The same uptake-bound construction as example 09
feeds this network instead of a linear program, which is what makes the two
directly comparable.

**What it deliberately leaves out:** flow, an abiotic reaction, any change
to the pore space, and any thermodynamic control on the rate. The network
predicts growth from two numbers and nothing else, no pH, no temperature, no
inhibition that was not in the sweep. Example 09 is the linear program this
network was fitted to reproduce; example 16 trains a surrogate of this kind
at start-up, from a bundled genome-scale model, rather than shipping a
pre-fitted header.

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
./scripts/setup_case.sh 11_surrogate run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. The surrogate solver needs neither GLPK nor Python and
is always compiled in, so no optional flag is needed:

```bash
cmake -B build -S . \
      -DCMAKE_BUILD_TYPE=Release \
      -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** The fitted network is a C++ header compiled
into the executable, not data read at start-up, so installing a newly fitted
network means copying it to `surrogateModel.hh` at the case root and
rebuilding. Note the location: `setup_case.sh` also lays a copy in `src/`,
but that copy is never compiled, `src/complab3d_processors_surrogate.hh`
includes the header as `"../surrogateModel.hh"`, so a new fit dropped into
`src/` produces bit-identical results from the old network with no error to
say so. Editing `CompLaB.xml` or `input/geometry.dat` needs no rebuild.
Refitting the network is a separate, optional, offline procedure, described
below; it is not needed to run this case.

**The offline procedure this case ships (optional, not needed to run the
case):**

```
    1.  sweep      solve the linear program on a grid of the two
                   uptake bounds                        -> training_data.csv
    2.  fit        train the network on that sweep       -> surrogate_weights.hh
                                                            surrogate_weights.srg
    3.  verify     check the exported header reproduces the trainer,
                   output by output
    4.  inspect    look at the response surface before trusting it
```

The sweep is the cost: `GRID^2` linear-program solves. `offline.sh`
demonstrates the procedure on the bundled *E. coli* core model, since the
*Geobacter* model the shipped header was fitted from (iAF987) is not
bundled, so running `offline.sh` as it stands fits a new network rather than
reproducing the one in the repository. It needs COBRApy, so `run.sh` does
not run it unless you ask it to.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Sweeps the linear program, fits the network, verifies the exported header against the trainer. Optional; a fitted header ships with the repository, and this step needs COBRApy. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| | `surrogateModel.hh` | The fitted network of section 4, compiled in, at the case root beside `CMakeLists.txt`. Ships fitted on iAF987. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.
