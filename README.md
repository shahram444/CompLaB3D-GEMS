# CompLaB3D-GEMS

**G**eometry evolution · **E**quilibrium and kinetics · **M**etabolism · **S**urrogates

A three-dimensional pore-scale reactive transport model. Lattice Boltzmann flow
and solute transport are coupled to microbial growth, abiotic chemistry and
aqueous speciation, in a pore space that **changes shape while the simulation
runs** as minerals precipitate and dissolve. The reaction rate itself is
pluggable: it can come from six different places, a single run may mix them, and
any of them can be gated by the free energy actually available at each voxel.

Built on [Palabos](https://palabos.unige.ch) 2.3.0 and on the two-dimensional
[CompLaB](https://bitbucket.org/MeileLab/complab) of Jung, Song and Meile.

---

## Contents

- [Everything it does](#everything-it-does) — the full feature inventory
- [What each rate path costs](#what-each-rate-path-costs)
- [Which rate path do I want?](#which-rate-path-do-i-want)
- [Start here](#start-here)
- [What to run, and when](#what-to-run-and-when)
- [The twenty-two examples](#the-twenty-two-examples)
- [Layout](#layout)
- [Configuration](#configuration) — including one XML with every ability in it
- [Verification](#verification)
- [Known limitations](#known-limitations)

---

# Everything it does

Eleven groups. Everything below is configured from one XML file. Two files document
that XML, and they answer different questions:
[`config/CompLaB.everything.xml`](config/CompLaB.everything.xml) has **one line
for every ability the solver has**, each with its unit and default — copy the
lines you need and delete the rest. [`config/CompLaB.reference.xml`](config/CompLaB.reference.xml)
is the same tags with the reasoning: why each exists and what goes wrong if you
get it wrong.

## 1. Flow

Flow is solved with a D3Q19 lattice Boltzmann scheme and BGK collision, run to a
convergence tolerance rather than a fixed number of steps. `<tau>` sets the
viscosity through ν = (τ − ½)/3.

You do not have to work out a pressure drop. Give a target Péclet number and the
solver measures the permeability of *your* geometry, computes the pressure drop
that hits the target on your characteristic length, and re-solves. `<delta_P>`
is only a starting guess. Setting `<Peclet>0</Peclet>` skips the flow solver
entirely, which is what a diffusion-dominated case wants.

A biofilm voxel resists flow: `<viscosity_ratio_in_biofilm>` makes it *n* times
harder to push fluid through, and zero makes it a solid wall. Because the
geometry changes as minerals form and biomass grows, the flow is re-solved every
`<ns_update_interval>` reaction steps, so a narrowing pore really does slow down.
The first solve and the later re-solves get separate iteration caps and
tolerances, because the first one is the expensive one.

## 2. Solute transport

Each dissolved species gets its own D3Q7 advection–diffusion lattice, and there
is no limit on how many you declare. Each carries two diffusivities, `<in_pore>`
and `<in_biofilm>`, applied per voxel according to what is actually in that voxel.

Boundary conditions are set per species and per side, and there are three kinds.
**Dirichlet** holds a value — that is how you inject something. **Neumann** sets
the plane to its neighbour every step, which is an outflow: mass leaves freely.
**`closed`** is bounce-back, and no flux crosses at all.

The third one matters more than it sounds. Neumann is not a synonym for closed,
and using it as one is a mass source: where the interior concentration is rising,
the boundary plane is topped up from nothing each step and streams that back in.
On example 14 that manufactured 32 % of the calcium. A species that is neither
fed nor drained — a dissolution product in a closed column, biomass that must
stay in the domain — wants `closed`.

Only `x = 0` and `x = nx-1` carry a boundary condition at all. The other four
faces get none, so a face left as open pore in the geometry has nothing defined
on it; every run counts those voxels and says so, and the shipped geometries all
draw their own confining wall.

A species marked `<immobile>true</immobile>` never advects, diffuses, streams or
collides. It sits still and accumulates its reaction term. That is how a mineral
is represented, and it is stable at any diffusivity including zero.

Nothing can go negative. Every reaction increment is clamped as
ΔC = max(RΔt, −C) before it is applied, at any step length, on every path.

## 3. The pore space

You can read a material map — an ASCII file with one integer per voxel — or have
the solver build one. `<generate>` makes a `channel`, `spheres`, `cylinders`,
`random`, `fracture` or `layered` domain from a porosity, a grain radius, an
aperture, a roughness, a layer count and a seed, and writes the result out as a
`.dat` so the run is reproducible from its own output. Segmented CT data comes in
through `<import_raw>` with a threshold and an optional inversion; TIFF and PNG
stacks go through `tools/geometry.py` first.

Whichever way it arrives, the geometry is inspected before the first step. The
solver reports the porosity, whether the pore space percolates from inlet to
outlet, and how much of it is isolated. **A sealed domain stops the run**,
because a pressure drop across one is not a well-posed problem.

Material numbers are declared, not guessed: pore, inert wall, grain interior, and
one number per attached biofilm. A microbe listed in `<material_numbers>` lives
in a biofilm; one that is not is a free-floating pool advected with the flow.

## 4. Microbiology

Any number of populations, each with its own name, initial density, decay
coefficient, kinetics and boundary conditions. Several may share one substrate,
one biofilm material number, or both — that is how competition is set up.

Biomass spreads by one of three solvers. `CA` is the cellular automaton, the
classic biofilm rule, and it has two spilling variants: `fraction` moves only the
excess above the maximum density, `half` moves half the voxel's total. `FD`
spreads by finite-difference diffusion. `LBM` treats the population as
planktonic, advected and diffused on its own D3Q7 lattice.

A pore voxel starts counting as biofilm — and starts resisting flow — above
`<thrd_biofilm_fraction>` of the maximum density. Decay is first order, per
population, in 1/s; it is also the fallback when a metabolic solve returns no
growth, so an organism that cannot grow decays rather than freezing.

Different populations in one run may use different rate paths: one on flux
balance analysis beside another on compiled kinetics is a supported
configuration, not a workaround.

## 5. Reaction — six interchangeable rate paths

All six present the same interface to the solver, so an organism picks one with
`<reaction_type>`, and any of them can be combined with compiled kinetics through
the `_and_kinetics` variants.

**Compiled kinetics** (`kinetics`) is the rate law you write yourself in
`defineKinetics.hh`. It is the cheapest by four orders of magnitude and the
easiest to read, and it needs a rebuild to change. **Abiotic kinetics** is the
same idea with no organism involved — mineral reactions, redox, sorption — in
`defineAbioticKinetics.hh`, switched on with `<enable_abiotic_kinetics>`.

**Flux balance analysis** predicts growth from a genome-scale model instead of
prescribing it. `glpk` solves the linear program in process, `cobrapy` sends it
through an embedded Python interpreter to the reference implementation. Both read
the same configuration; GLPK is the one for production and COBRApy is useful as a
cross-check.

**The surrogate network** (`surrogate`) replaces the linear program with a fitted
network — a few hundred weights that reproduce its answers at roughly 1/400 of
the cost. It can be compiled in or read at start-up from a `.srg` file.

**The symbolic law** (`symbolic`) is a short algebraic expression in a `.sym`
text file, discovered by symbolic regression from the same samples the surrogate
is fitted to. It fits slightly worse than a network and can be read, checked
against known kinetics, and quoted in a paper.

**The graph network** (`graphnet`) does message passing over the bipartite
species–reaction graph and returns the whole coupled rate vector in one
evaluation. The stoichiometric matrix is supplied as structure rather than
learned.

Nothing here needs a rebuild except compiled kinetics and a compiled surrogate. A
different `.sym`, `.gnn` or `.srg` is a different file.

### 5a. What the flux-balance paths actually do

At every voxel and every step they solve **max cᵀv subject to Sv = 0 and
ℓ ≤ v ≤ u**. The uptake bounds are Michaelis–Menten, ℓ = −V_max·C/(K_c + C),
built from `<maximum_uptake_flux>` and `<half_saturation_constants>`; a zero
V_max means the organism is limited by supply rather than by its enzymes. The
previous voxel's simplex basis warm-starts the next solve, which is most of why
the GLPK path is usable at all.

Exchange reactions are addressed by name through `<exchange_reaction_names>`,
which is resolved against the model and stops the run if a name is missing.
Positional indices still work but point silently elsewhere if the model is
revised. You can add arbitrary reaction bounds, and `<equate_bounds>` forces
reversible pairs together.

A model comes from a local SBML file, a BiGG identifier, or one of the three
shipped in `models/`. Downloads are off unless you allow them.

Units convert at exactly one boundary: mmol/gDW/h inside the linear program,
mol/L everywhere else, through each organism's `<biomass_molar_mass>`. With
speciation switched on, `<fba_concentration_basis>total</fba_concentration_basis>`
rebuilds each substrate from the equilibrium tableau — every complex that carries
it — and draws consumption back down across those complexes in proportion.

Three LP algorithms are available: simplex by default, interior point, and exact
rational arithmetic. Each microbe's problem can be dumped as `.lp` or `.mps` at
start-up when something is wrong.

### 5b. The three learned paths, side by side

They differ in what they return and in whether you can read them.

The **surrogate** is 2 inputs → four hidden layers of 10 tanh units → a linear
output, 465 parameters, fitted to the linear program's growth over a swept grid.
It returns growth, and optionally every exchange flux as well. It is not
readable, and a new one needs either a rebuild or a new `.srg`.

The **symbolic law** is an expression tree fitted to any table of data. It
returns one rate, and you extend it to the other species by writing them as
multiples — which keeps the stoichiometry exact because you wrote it that way.
**It is algebra you can quote and disagree with**, which is the whole argument
for this path.

The **graph network** is message passing on the species–reaction graph, with the
stoichiometry supplied as structure. It returns the entire coupled rate vector at
once. It is not readable.

Ranges are handled differently and it matters. The symbolic and graph-network
paths clamp every evaluation to the box they were fitted over and count how often
that happened; the closing report says so. The surrogate does neither — see
[Known limitations](#known-limitations).

#### A surrogate can return growth alone, or growth and the fluxes

The difference shows up in a run. If the network returns growth alone, the solver
computes substrate consumption from a Monod term, which is exact only where the
swept uptake bound was the binding constraint, and cannot release a product at
all. A multi-output network returns the flux the linear program actually ran, and
is correct anywhere inside its training box.

The cost per voxel is identical — one extra row in the last matrix. The
multi-output network is harder to fit, because growth is smooth while the flux
columns are piecewise linear with kinks.

Both come from the same two commands, and both load through the same
`<weights_file>` tag, so you can compare them on one case by editing one line:

```bash
python3 tools/surrogate/generateTrainingData.py MODEL \
        --exchange EX_ac_e --range 1e-3 10 --log \
        --exchange EX_o2_e --range 1e-5 0.5 --log \
        --also EX_co2_e --grid 141 -o sweep.csv     # --growth-only for the old behaviour
python3 tools/surrogate/trainSurrogate.py sweep.csv --layers 10 10 10 10 \
        -o net.hh --srg net.srg                     # both formats, one fit
python3 tools/surrogate/verifyExport.py net.hh      # checks every output, not just growth
```

#### Two formats, and why

One fit can be written twice. `net.hh` is compiled in and marginally faster;
`net.srg` is read at start-up through `<weights_file>`, so changing the network is
one line of XML instead of a rebuild. `--srg` writes both from the same fit, so
they cannot disagree by construction — and `tests/test_surrogate_parity.cpp`
checks that the two evaluators agree anyway, to 8 × 10⁻¹¹. That check earns its
place: the forward pass genuinely exists twice in this codebase, and nothing else
would catch a drift between them.

The symbolic path has a full expression language: precedence, right-associative
exponentiation, unary minus, `exp` `log` `sqrt` `min` `max` `pow` `abs`, the
constants `e` and `pi`, division by zero returning zero rather than infinity, and
malformed input **rejected with a character position** rather than guessed at.
Every expression and its declared range is printed into the log, so a result
carries the law that produced it.

### 5c. Two refinements of the linear program

Neither is a seventh rate path. Both sit inside the flux-balance paths and
change how the program is posed, not what the solver interface looks like, so
everything downstream is unchanged: transport, the thermodynamic factor and
geometry evolution all see the same thing they saw before.

**Multi-step flux balance analysis** (`<multi_step>`) fixes the fluxes a single
program leaves free. Growth has a unique optimum; the by-products that carry the
same electrons do not. On the E. coli core model at 10 mmol gDW⁻¹ h⁻¹ of glucose
and oxygen, acetate export may be anything from 9.9057 to 11.5033 at exactly the
optimal growth rate, and the simplex returns whichever end its pivoting rule
reached. The chain replaces the one program with several: maximise growth, pin
that optimum as a constraint, maximise the first by-product subject to it, pin
that too, and continue. After the last stage nothing is free, so the same model
at the same bounds returns the same numbers on any solver, every time.

```xml
<multi_step>
    <stage_reactions>  Biomass_Ecoli_core  EX_ac_e  EX_for_e  </stage_reactions>
    <retain_fraction>         1.0            1.0      1.0     </retain_fraction>
    <stage_direction>         max            max      max     </stage_direction>
</multi_step>
```

Stages are named, not indexed, and every name is checked against the model at
start-up. `<retain_fraction>` defaults to 1.0 throughout, which gives up nothing
and adds no fitted parameter; below 1.0 the organism may trade growth for
excretion, which is what measured cells do, and that fraction then belongs in
whatever the run is reported in. This is the method of Song et al. (2025); the
same device is lexicographic optimisation in DFBAlab (Gomez, Höffner & Barton
2014). Example 21 runs it.

**Cybernetic switching** (`<cybernetic>`) gives the organism a preference
between carbon sources. One program has none, because both sources raise the
same objective, so a cell offered glucose and acetate consumes both at once.
The method solves one program per source with the *other* sources' uptake shut.
Uptake only, never release, so the cell stays free to excrete what it will later
eat. The flux vectors are then blended by how much carbon each source is
supplying, from a Monod score per source.

```xml
<cybernetic>
    <sources>                glucose  acetate </sources>
    <substrate_ids>             0        2    </substrate_ids>
    <carbon_number>             6        2    </carbon_number>
    <uptake_kmax>             10.0      4.4   </uptake_kmax>
    <uptake_half_saturation>   0.05     0.05  </uptake_half_saturation>
    <weight_floor>           1e-3            </weight_floor>
</cybernetic>
```

The blend is mass-consistent rather than approximately so: each vector satisfies
Sv = 0, and a convex combination of vectors in the null space of S is itself in
that null space. `<weight_floor>` skips any source whose weight is negligible,
which is what keeps the cost near that of a single program rather than one
program per source. What it does not carry is enzyme state between steps, so
there is no lag and no diauxic plateau. That is adequate exactly when the switch
is fast compared with transport, which is the pore-scale case. Example 22 runs
it.

## 6. Thermodynamic control

All six rate paths answer one question: how fast *can* this organism run its
reaction, given what is here. None of them asks the other one: does the reaction,
at these particular concentrations, release enough energy to be run at all.

That second answer differs voxel by voxel in a pore-scale domain. A reaction can
be strongly exergonic where its substrate arrives and past its thermodynamic
threshold sixty micrometres further in, once its own products have built up
against the diffusive resistance of the aggregate they were made in. A Monod law
keeps it running there, because it never reads a product.

`<thermodynamics>` supplies what is missing. At every voxel the solver computes
the free energy at the local composition, ΔG = ΔG° + RT ln Q, subtracts what the
organism must conserve as ATP, and returns a number between zero and one:

**F_T = max(0, 1 − exp((ΔG + m·ΔG_ATP) / (χRT)))**

after Jin & Bethke (2003), and Craig (2024) equations 3.3 to 3.6. The rate the
chosen path produced is multiplied by it. Far from the threshold F_T is one and
nothing changes; at the threshold it is zero and the reaction stops.

**This is a factor, not a seventh rate path.** It multiplies the answer the
chosen path already gave, so all six inherit it and none of them had to be
rewritten to accept it. On the two flux-balance paths it scales the *solved*
fluxes rather than the uptake bounds, because a linear program handed smaller
bounds re-optimises and can return a different byproduct pattern — that would be
a modelling change rather than an energy constraint.

The same free energy also gives the growth yield, which every other path takes as
a constant. Heijnen & van Dijken (1992) estimate the energy dissipated per unit
of biomass built from the carbon source alone, and the yield follows as
Y = ΔG / −(ΔG_ana + ΔG_dis).

**The energetics live in a file, not in the XML.** A `.thm` file carries the
reaction, its standard free energy, its ATP threshold and optionally the yield
terms. Those are cited claims that get revised and quoted in a paper, so they
belong in something a reader can open, diff and cite — the same reasoning that
puts a fitted rate law in a `.sym` file. The whole file is echoed into the log, so
the result carries the energetics that made it.

Two limits worth knowing before you use it. The gate applies to biotic paths
only: an abiotic reaction is either a kinetic law you wrote, in which case it is
yours to gate, or an equilibrium, in which case the equilibrium solver has already
answered the same question exactly. And on the compiled-kinetics path the rate law
returns one combined rate vector for every organism in a voxel at once, so it can
carry only one gate; a file with more than one reaction block is refused at
start-up when that path is in use.

Worked case: [example 19](examples/19_thermodynamic_gate/).

## 7. Aqueous speciation

Speciation is set up as a components-and-species tableau: master species, a
stoichiometry row per substrate, and a formation constant per substrate. It is
solved by continued fractions in every pore voxel at every step, and it is by a
wide margin the most expensive part of the code.

It couples to the metabolic layer through the free-versus-total concentration
basis described above. A substrate whose stoichiometry row is all zeros is
ignored by the solver, which is how a mineral or a biomass pool that reacts only
kinetically is declared.

The limits are stated rather than implied: ideal activities with no Debye–Hückel
or Davies correction, no temperature dependence, no gas phase, no redox couple,
no saturation index. The log K values you supply are conditional constants at your
own ionic strength and 25 °C. A complete 95-species, 17-component uranium tableau
ships as a worked case.

## 8. Geometry evolution

These are the two capabilities that make the pore space a variable rather than a
constant. Both run alongside whichever rate path you chose.

**Precipitation closes the pore.** An immobile mineral accumulates in a pore
voxel until it reaches `<max_precipRho>`, derived from the mineral's molar volume.
`<surface_only>` restricts growth to voxels touching a surface, which is
heterogeneous nucleation and the physical case. A full voxel becomes solid, and
`<perm_ratio>0` makes it an impermeable wall. The flow is re-solved on the reduced
pore space at `<update_interval>`, so the permeability really falls; the run
reports the porosity dropping in steps and stops when the pore space seals.

**Dissolution opens it again.** Any number of `<phaseN>` blocks declare a solid,
each with a name, a material number, the substrate it releases, its full density
and its initial fill. It is consumed by the water touching it, averaged over the
voxel's open faces, and the products are deposited into the neighbouring pore with
D3Q7 weights — ¼ at rest and ⅛ in each of six directions, summing to exactly the
increment requested. A voxel reopens below `<reopen_fraction>` × full, so one
sitting on the threshold cannot flicker between states and force a flow re-solve
every interval. `<is_precipitate>` distinguishes a phase that formed during the
run from one that was there at the start.

## 9. Pore scale out — one aggregate, one continuum rate

A pore-scale run resolves an aggregate; a column- or reservoir-scale model carries
one concentration per grid block and needs one rate. `<upscaling>` measures the
bridge between them, the effectiveness factor

*η* = ⟨*r*⟩<sub>aggregate</sub> / *r*(*C*<sub>bulk</sub>)

together with the Thiele modulus *φ* = *R*√(*k*/*D*) and the classical sphere
result for comparison. The numerator is not recomputed: every rate path has
already written its result into the increment lattices, and the diagnostic samples
those after the rate processors run and before the increments are applied, so it
averages the number the solver is about to use. At iteration zero, when every
voxel is still at the bulk composition, it reports *η* = 1.0000000002 — the
measurement checking itself.

Why this needs a thermodynamic gate to be more than a textbook exercise: with
Monod kinetics alone the core of an aggregate runs slowly, but with the gate on it
**stops**, at the depth where its own products have raised Δ*G* past what the
organism can use. That is a moving internal boundary, and no first-order Thiele
analysis contains it. So *F*<sub>T</sub>(*C*<sub>bulk</sub>) is reported beside
*φ*, and that is also the answer to whether you could just evaluate the gate at
the bulk and be done: *r*(*C*<sub>bulk</sub>) already includes it, so whatever *η*
departs from 1 is precisely the error that shortcut makes.

`<freeze_biomass>` holds the catalyst while the concentration profile relaxes.
Without it *η* never settles, because biomass grows faster at the rim than in the
core and keeps moving the answer for a reason that has nothing to do with
transport. The solver checks steadiness itself and refuses to call a transient a
result. Example 20 is the case; `offline/upscale.py` sweeps radius and bulk
composition and fits the curve.

## 10. Diagnostics and output

Fields are written as VTI — velocity, every substrate, every microbe, and the
material map — named from your own `<name_of_substrates>` rather than by index.

The summary CSV gets one row per interval: porosity, and the total, mean, minimum
and maximum of every substrate and microbe, **over open voxels only**, so a run
that seals pore space is not divided by a moving denominator. Totals run over
x = 1…nx−2, excluding the ghost columns, so a total does not jump when you change
an inlet concentration.

Every field also gets a `_held` column, and it is the one to know about.
`computeDensity()` asks each cell's dynamics, and Palabos's `BounceBack` and
`NoDynamics` both answer from a stored number and ignore the populations they are
holding. So mass in flight at a wall, mass resting inside a grain, and mass in a
closed boundary plane are all absent from `_total`. `_held` is that difference,
measured from the populations themselves, and **`_total + _held` is the conserved
quantity** — it is what `<conserve>` is checked against. The distinction is not
academic: in a fully closed box with every reaction off, `_total` alone drifts by
7 % while the sum is constant to 4 × 10⁻¹³.

Each `<conserve>` line names a sum the reaction network cannot create or destroy —
one mole of Fe²⁺ removed must appear as one mole of FeS. A drift beyond tolerance
is reported *the moment it happens*. A check naming a substrate that does not
exist is reported as SKIPPED, never as PASS.

There is a verbose per-iteration mass-balance mode for when a run is going wrong,
and `<track_performance>` times the solver kernels and suppresses all output
writing.

The most useful diagnostic is the cheapest one: the geometry, the enabled
features and each organism's rate path are echoed before the first step, so you
can stop a wrong run in the first second rather than in the second week.

## 11. Running it

Parallelism is Palabos block decomposition over MPI; a serial run needs no MPI at
all. The optional back ends are opt-in at build time with `-DENABLE_GLPK=ON` and
`-DENABLE_COBRAPY=ON`, and if you switch one on in the XML that was not compiled
in, the solver **stops at start-up and names the build option** rather than
quietly ignoring you. `FBA_BULK_ONLY` restricts the FBA processors to the bulk
domain for 20–40 % more speed. The surrogate, symbolic, graph-network and
thermodynamic paths need no build option at all.

Binary checkpoints are written at an interval and reloaded on restart, so a long
run survives a queue limit.

`scripts/setup_case.sh` assembles everything a case needs into one directory,
including the solver sources, so it builds on its own.

---

## What each rate path costs

Measured on one core, for one rate evaluation at one voxel.

| Rate path | µs/voxel | Relative | What it returns |
|---|---:|---:|---|
| Compiled kinetics | 0.003 | 1× | The rate law you wrote |
| Symbolic law (`.sym`) | 0.11 | 37× | A fitted law, still readable as algebra |
| Graph network (`.gnn`) | 2.5 | 830× | The whole coupled rate vector |
| Flux balance, GLPK | 40 | 13,000× | A genome-scale optimum, with internal fluxes |
| Flux balance, COBRApy | 800 | 270,000× | The same, through the reference implementation |

Relative cost, not a portable performance claim. The surrogate reproduces the
GLPK column at roughly 1/400 of its cost.

---

## Which rate path do I want?

Read down until the first *yes*.

**Is a rate law already known?** Use **compiled kinetics**. It is cheaper than
everything else by four orders of magnitude and it is the easiest thing in the
repository to read.

**Do you need the internal fluxes, not just growth?** Use **flux balance
analysis**. It is the only path that says what the organism is doing inside.

**Is growth the only output, with one or two substrates limiting?** Use the
**surrogate network**. It reproduces the linear program's growth at about 1/400
of its cost.

**Should the result be readable and arguable?** Use the **symbolic law**. The
answer is algebra you can quote, and a new one needs no rebuild.

**Otherwise — many species coupled through many reactions?** Use the **graph
network**. It returns the whole rate vector in one evaluation, in stoichiometric
ratio.

Then, separately: **is the reaction low-energy, and are its products
transported?** If so add `<thermodynamics>` on top of whatever you chose. That
covers most anaerobic respiration in sediments — methanogenesis, sulfate
reduction, anaerobic methane oxidation, syntrophic fermentation — where a local
product build-up is the difference between a reaction that runs and one that does
not. It is one factor and it composes with every path above.

And if you chose flux balance analysis, two further questions. **Do you report a
by-product flux, not just growth?** Add `<multi_step>`, because the by-products
are not uniquely determined by the growth optimum and the number you get is
otherwise the solver's choice. **Does your organism have more than one carbon
source available at once?** Add `<cybernetic>`, because a single program will eat
both rather than preferring the better one. Neither applies to the other four
paths, and neither is a rate path itself: both change how the program is posed.

Precipitation and dissolution are not on this list. They are geometry processes
and run alongside whichever rate path you choose.

---

## Start here

```bash
git clone https://github.com/shahram444/CompLaB3D-GEMS.git
cd CompLaB3D-GEMS

# assemble a runnable case, build it, run it
./scripts/setup_case.sh 13_precipitation run/mycase
cd run/mycase
cmake -B build -S . && cmake --build build -j
./complab CompLaB.xml
```

That case is set up to form FeS where an iron front meets a sulfide front and to
seal the voxels as they fill. **It does not currently do so** — see
[Known limitations](#known-limitations) — so read it as the shortest complete
configuration to run, not as a result.

Build requirements and the optional dependencies are in
[`INSTALL.md`](INSTALL.md).

---

## What to run, and when

Most of the work happens outside the solver. **[`pipelines/`](pipelines/) is the
map of what has to be prepared, in what order.**

```
pipelines/
├── A_preprocess/       the pore space, and starting fields        always
├── B_offline_models/   the model your rate path needs         path-dependent
├── C_run/              assemble, build, run                       always
└── D_postprocess/      fields, and whether the run is right       always
```

The one thing worth knowing before you start: **four of the six rate paths need
nothing at all from stage B.** Compiled kinetics needs no preparation;
precipitation and dissolution need none either. The surrogate needs the most,
and is the only path whose preparation ends in a recompile.

Concretely: compiled kinetics needs nothing. Flux balance analysis needs **B1**,
which exports the metabolic model. The surrogate needs **B1 then B2** — sweep the
linear program, fit the network, paste it in, rebuild. The symbolic law needs
**B3**, which runs the search and hands you a Pareto set to choose from. The graph
network needs **B4**, which trains the network and writes the `.gnn`. The
thermodynamic gate needs no training at all — you write the energetics by hand
from measured values — but it has an offline *check* worth running, which reports
where the gate closes before the solver starts.

---

## The twenty-two examples

Each adds one thing to the one before it. Run them in order until something
breaks in a way you do not understand — that is the piece worth reading about.

| | Case | What it adds | Offline first |
|---|---|---|---|
| 01 | `flow_only` | Navier–Stokes on a pore space, nothing else | — |
| 02 | `diffusion_only` | one solute, transported, no reaction | — |
| 03 | `abiotic_kinetics` | a chemical rate law, no organisms | — |
| 04 | `equilibrium` | aqueous speciation by continued fractions | — |
| 05 | `biotic_cellular_automaton` | biomass, spread by the automaton | — |
| 06 | `biotic_finite_difference` | the same, by finite difference | — |
| 07 | `biotic_lattice_boltzmann` | the same, on a D3Q7 lattice | — |
| 08 | `two_microbes` | two populations competing for one substrate | — |
| 09 | `fba_glpk` | growth from a genome-scale linear program | **B1** |
| 10 | `fba_cobrapy` | the same through the reference implementation | **B1** |
| 11 | `surrogate` | a fitted network in place of the linear program | **B1 → B2** |
| 12 | `mixed_reaction_types` | one population on FBA beside one on kinetics | **B1** |
| 13 | `precipitation` | FeS fills the pore and seals it | — |
| 14 | `dissolution` | calcite is eaten away and the pore reopens | — |
| 15 | `precip_and_dissolution` | both at once, on different phases | — |
| 16 | `complete_pipeline` | geometry, chemistry, biology and metabolism together | **B1** |
| 17 | `symbolic_law` | the rate law read from a text file, no rebuild | — |
| 18 | `graph_network` | the whole coupled rate vector from one evaluation | — |
| 19 | `thermodynamic_gate` | the rate multiplied by the free energy available for it | — |
| 20 | `upscaling` | one resolved aggregate reduced to a continuum rate | — |
| 21 | `multistep_fba` | a chain of programs, so the by-product fluxes are determined | — |
| 22 | `cybernetic_switching` | the organism changes carbon source when the first runs out | — |

**Each case carries its own pipeline.** Not a pointer to a shared tool — the
actual `preprocess.py` that builds that geometry, the `offline.sh` that prepares
that model, and the `postprocess.py` that knows which of *that case's* totals
are closed and can therefore be balanced. Assemble a case and run `./pipeline.sh`,
or run the four steps yourself.

Full details, and which geometry each one runs on, in
[`examples/README.md`](examples/README.md).

---

## Layout

```
CompLaB3D-GEMS/
├── src/                the solver: Palabos data processors and headers
├── config/
│   ├── CompLaB.everything.xml    one line for every ability, with units
│   ├── CompLaB.reference.xml     the same tags, with the reasoning behind each
│   ├── kinetics/                 the shared chemistry headers
│   └── geometry/                 the four shared pore geometries
├── examples/           twenty-two cases, from flow-only to the full pipeline
├── pipelines/          what to run before and after the solver
├── tools/              the offline tools: model export, fitting, training
│   └── surrogate/      the surrogate training path, Python and MATLAB
├── tests/              the regression suite, and the structural checks
├── models/             three genome-scale models, and the FBA toy model
├── scripts/            setup_case.sh, which assembles a runnable case
└── docs/               the user guide, and how to publish this repository
```

### The offline tools

`tools/geometry.py` builds a pore space or inspects one, reporting porosity,
percolation and isolated pore. `tools/extractMM.py` converts an SBML, MATLAB or
JSON genome-scale model into the flat XML the GLPK path reads, and prints the
exchange-reaction table you fill the configuration from.
`tools/makeKinetics.py` generates a `defineKinetics.hh` from a reaction list, and
`tools/makeEquilibrium.py` builds the components/stoichiometry/log K tableau.

For the learned paths: `tools/surrogate/` sweeps the linear program, fits the
network, verifies that the export reproduces the trainer, and plots the response
surface — in Python and in MATLAB. `tools/fit_symbolic.py` is the
genetic-programming search that returns a Pareto set of rate laws.
`tools/train_graphnet.py` trains a graph network and writes the `.gnn`.

Afterwards, `tools/postprocess.py` produces slices, histories and the mass-balance
report, and `tools/vtireader.py` reads VTI output into NumPy.
`tools/complab3d_cobrapy.py` is the Python side of the COBRApy back end.

### One thing about the examples

**An example folder holds everything that case needs.** Not a pointer to a
shared tool — the metabolic model it reads, the training code that produced what
it was trained on, the data that code was run on, its rate-law file, its pore
space, and its own pre- and post-processing:

```
examples/11_surrogate/
├── CompLaB.xml
├── preprocess.py      builds this case's pore space, checks it percolates
├── offline.sh         sweep -> fit -> verify -> install -> inspect
├── postprocess.py     reads the summary, runs this case's balance check
├── pipeline.sh        the four in order
├── input/geometry.dat
├── models/            e_coli_core.xml.gz, and where it came from
└── training/          generateTrainingData.py, trainSurrogate.py,
                       verifyExport.py, inspectSurrogate.py, the MATLAB
                       equivalents, and a sweep already done
```

`scripts/setup_case.sh` lays down the two shared kinetics defaults, copies the
case folder whole on top, and adds the solver sources. Then `./pipeline.sh`.

**Duplication inside `examples/` is deliberate, and it has a cost.** Several
cases carry the same geometry; three carry the same model exporter. The point is
that a case folder is the whole procedure. The price is that a fix to a trainer
has to be applied to every case that carries it — so `tests/check_repo.sh` fails
if any copy drifts out of step with the original in `tools/`. Outside
`examples/`, nothing exists twice.

---

## Configuration

Everything the solver does is configured from one XML file, and two files in
`config/` document it. They answer different questions and neither repeats the
other.

**[`config/CompLaB.everything.xml`](config/CompLaB.everything.xml) — what can this
code do?** One line for every ability the solver has, each with its unit and its
default, grouped into the twelve blocks the XML actually has. It is not a case and
will not run as it stands, because several things in it are alternatives to one
another — you cannot both read a geometry file and generate one, and an organism
has one rate path, not seven. Read it to see the whole surface at once, then copy
the lines you need into your own `CompLaB.xml` and delete the rest.

**[`config/CompLaB.reference.xml`](config/CompLaB.reference.xml) — why is this tag
here, and what goes wrong if I get it wrong?** The same tags with the reasoning
attached. That file is the manual; the other one is the index.

Each example's own `CompLaB.xml` is a working subset, and there is no third copy
of the documentation anywhere in the tree.

Units are one convention throughout: **micrometres, mol/L, seconds, m²/s**.
Biomass is in mol/L, the same unit as the chemicals, on purpose. Flux balance
analysis is the single exception — it works in mmol/gDW/h internally because
that is what published metabolic models use, and converts at the boundary of the
linear program through each organism's `<biomass_molar_mass>`.

The blocks this repository adds to the base model, in short — the full set is in
`CompLaB.everything.xml`:

```xml
<precipitation>                          <!-- geometry closes -->
    <enabled>true</enabled>
    <solid_substrate>7</solid_substrate>
    <max_precipRho>48.9</max_precipRho>   <!-- 1000 / molar volume, cm3/mol -->
    <surface_only>1</surface_only>        <!-- heterogeneous nucleation -->
    <perm_ratio>0.0</perm_ratio>          <!-- a filled voxel is a wall -->
    <update_interval>500</update_interval>
</precipitation>

<dissolution>                            <!-- geometry opens -->
    <enabled>true</enabled>
    <reopen_fraction>0.9</reopen_fraction>
    <update_interval>200</update_interval>
    <phase0>
        <name>calcite</name>
        <material_number>0</material_number>
        <substrate>13</substrate>
        <full_density>27.1</full_density>
        <initial_fill>27.1</initial_fill>
    </phase0>
</dissolution>

<symbolic>                               <!-- a readable rate law -->
    <expressions_file>mylaw.sym</expressions_file>
</symbolic>

<graphnet>                               <!-- the whole rate vector at once -->
    <network_file>geobacter_network.gnn</network_file>
</graphnet>

<thermodynamics>                         <!-- may the reaction run here at all? -->
    <enabled>true</enabled>
    <energetics_file>input/aom.thm</energetics_file>
</thermodynamics>

<microbe0>                               <!-- where this organism's rate comes from -->
    <solver_type>CA</solver_type>
    <reaction_type>graphnet_and_kinetics</reaction_type>
</microbe0>

<diagnostics>                            <!-- prove the run conserved mass -->
    <enabled>true</enabled>
    <summary_csv>summary.csv</summary_csv>
    <interval>500</interval>
    <tolerance>1e-6</tolerance>
    <conserve>Fe2+FeS</conserve>
</diagnostics>
```

---

## Verification

Two scripts, and they check different things.

```bash
./tests/check_repo.sh     # the tree: no duplicate file, every example assembles
./tests/run_tests.sh      # the solver: seven suites, 154 checks
```

`check_repo.sh` takes a second and needs nothing installed. It is the one to run
after editing the tree — it catches a broken `case.files`, a duplicated asset, a
dead link in a README, an example that no longer assembles.

`run_tests.sh` is the real suite. The rate processors are compiled against a
minimal Palabos stub and run on lattices whose contents the test sets directly,
so the assertions are exact arithmetic identities rather than tolerance bands on
a simulation output. What it establishes:

- the expression language: precedence, right-associative exponentiation, unary
  minus, every function, division by zero returning zero, forward references
  rejected rather than silently mis-parsed;
- the `.sym` and `.gnn` loaders, including that a truncated file is refused
  rather than read past its end;
- the mass budget: no species can be drawn below zero, at any step length;
- the D3Q7 deposit: the seven populations sum to exactly the increment
  requested;
- the chemical path runs with no biomass present at all;
- the shipped `aom.gnn` reproduces the dual-Monod law it was fitted from, on
  300 samples it never saw;
- **the Python trainer and the C++ solver agree to machine precision** —
  1.36 × 10⁻¹⁸ across 50 inputs and 6 outputs, several deliberately outside the
  training box so the clamping is compared as well as the arithmetic.

---

## Working on it

[`CONTRIBUTING.md`](CONTRIBUTING.md) states the rule the tree is built on and what
to update when you change the chemistry interface.
[`docs/PUBLISHING.md`](docs/PUBLISHING.md) covers putting this on GitHub step by
step and the edit loop afterwards. [`docs/README_MANUAL.md`](docs/README_MANUAL.md)
explains how to build the user guide from its LaTeX sources, and
[`CHANGELOG.md`](CHANGELOG.md) records what each release added.

---

## Authors

**Shahram Asgari** and **Christof Meile**, Department of Marine Sciences,
University of Georgia, Athens, GA, USA — <shahram.asgari@uga.edu>

Meile Lab, University of Georgia. This work extends the two-dimensional CompLaB
v1.0 of Heewon Jung, Hyun-Seob Song and Christof Meile, whose decision to keep
the chemistry in a user-editable header rather than inside the solver is what
made every addition here possible.

Supported by the U.S. Department of Energy, Office of Science, Office of
Biological and Environmental Research, Genomic Science Program, under Award
Number DE-SC0022991.

## Licence

GNU Affero General Public License v3.0 or later. See [`LICENSE`](LICENSE).

## Citing

See [`CITATION.cff`](CITATION.cff).
