# CompLaB3D-GEMS

GEMS:
Geometry evolution, precipitation, dissolution
Equilibrium and kinetics, the base chemistry these attach to
Metabolism, flux balance analysis, GLPK and COBRApy
Surrogates, the neural, symbolic and graph-network rate laws

A three-dimensional pore-scale reactive transport model. Lattice Boltzmann flow
and solute transport are coupled to microbial growth, abiotic chemistry and
aqueous speciation, in a pore space that **changes shape while the simulation
runs** as minerals precipitate and dissolve. The reaction rate itself is
pluggable: it can come from six different places, and a single run may mix them.

Built on [Palabos](https://palabos.unige.ch) 2.3.0 and on the two-dimensional
[CompLaB](https://bitbucket.org/MeileLab/complab) of Jung, Song and Meile.

---

## Contents

- [Everything it does](#everything-it-does) — the full feature inventory
- [What each rate path costs](#what-each-rate-path-costs)
- [Which rate path do I want?](#which-rate-path-do-i-want)
- [Start here](#start-here)
- [What to run, and when](#what-to-run-and-when)
- [The sixteen examples](#the-sixteen-examples)
- [Layout](#layout)
- [Configuration](#configuration)
- [Verification](#verification)
- [Known limitations](#known-limitations)

---

# Everything it does

Nine groups. Everything in them is configured from one XML file, and every tag
named below is documented in
[`config/CompLaB.reference.xml`](config/CompLaB.reference.xml).

## 1. Flow

| Feature | What it means |
|---|---|
| **D3Q19 lattice Boltzmann, BGK collision** | Steady Stokes/Navier–Stokes flow through the pore space, solved to a convergence tolerance rather than a fixed step count. |
| **Péclet auto-calibration** | You give a target Péclet number, not a pressure drop. The solver measures the permeability of *your* geometry, works out the pressure drop that hits the target on the characteristic length, and re-solves. `<delta_P>` is only a seed. |
| **Diffusion-only mode** | `<Peclet>0</Peclet>` skips the flow solver entirely. |
| **Biofilm as a flow resistance** | `<viscosity_ratio_in_biofilm>` makes a biofilm voxel *n* times harder to push fluid through; `0` makes it a solid wall. |
| **Flow re-solved as the geometry changes** | Every `<ns_update_interval>` reaction steps, so a narrowing pore really does slow down. |
| **Two convergence budgets** | The first solve and later re-solves get their own iteration caps and tolerances — the first solve is the expensive one. |
| **Relaxation time control** | `<tau>` sets viscosity through ν = (τ − ½)/3. |

## 2. Solute transport

| Feature | What it means |
|---|---|
| **D3Q7 advection–diffusion, one lattice per species** | Any number of dissolved species, each with its own transport. |
| **Two diffusivities per species** | `<in_pore>` and `<in_biofilm>`, applied per voxel according to what is actually there. |
| **Per-species, per-side boundary conditions** | Dirichlet (hold a value — this is how you inject) or Neumann (free outflow), set independently on the inlet and the outlet for every species. |
| **Immobile species** | `<immobile>true</immobile>` marks a solid: it never advects, diffuses, streams or collides. It only accumulates its reaction term in place. This is how a mineral is represented. |
| **Mass-conserving sink clamp** | ΔC_k = max(R_k Δt, −C_k). No species can ever be drawn below zero, at any step length. |
| **Diffusivity refresh interval** | `<ade_update_interval>`, so a run with a slowly evolving biofilm does not rebuild the coefficient field every step. |

## 3. The pore space

| Feature | What it means |
|---|---|
| **Read a material map** | ASCII `geometry.dat`, one integer per voxel. |
| **Generate one** | `<generate>` builds `channel`, `spheres`, `cylinders`, `random`, `fracture` or `layered` from porosity, grain radius, aperture, roughness, layer count and a seed. The result is written out as `.dat`, so the run is reproducible from its own output. |
| **Import segmented CT data** | `<import_raw>` reads a raw binary volume with a threshold and an optional inversion. TIFF and PNG stacks go through `tools/geometry.py`. |
| **Inspected before the first step** | Porosity, whether the pore space percolates inlet to outlet, and how much of it is isolated. **A sealed domain stops the run**, because a pressure drop across one is not a well-posed problem. |
| **Declared material numbers** | Pore, inert wall (bounce-back), grain interior, and one number per attached biofilm. |
| **Attached versus planktonic** | A microbe listed in `<material_numbers>` occupies a biofilm; one that is not is a free-floating pool advected with the flow. |

## 4. Microbiology

| Feature | What it means |
|---|---|
| **Any number of populations** | Each with its own name, initial density, decay coefficient, kinetics and boundary conditions. |
| **Three biomass spreading solvers** | `CA` — cellular automaton, the classic biofilm rule; `FD` — finite-difference diffusion; `LBM` — planktonic, advected and diffused on its own D3Q7 lattice. |
| **Two CA spilling rules** | `fraction` moves only the excess above the maximum density; `half` moves half the voxel's total. |
| **Biofilm threshold** | The fraction of maximum density at which a pore voxel starts counting as biofilm — and starts resisting flow. |
| **First-order decay** | Per population, in 1/s. Also the fallback when a metabolic solve returns no growth: the organism then decays rather than freezing. |
| **Populations that compete** | Several pools may share one substrate, one biofilm material number, or both. |
| **Mixed rate paths in one run** | One population on flux balance analysis beside another on compiled kinetics. |

## 5. Reaction — six interchangeable rate paths

All six present the same interface to the solver, so they are chosen per
organism with `<reaction_type>`, and any of them can be combined with compiled
kinetics through the `_and_kinetics` variants.

| Path | `reaction_type` | Where the rate comes from | Needs a rebuild? |
|---|---|---|---|
| **Compiled kinetics** | `kinetics` | `defineKinetics.hh` — the rate law you write in C++ | yes |
| **Abiotic kinetics** | (`enable_abiotic_kinetics`) | `defineAbioticKinetics.hh` — mineral reactions, redox, sorption, with no organism involved | yes |
| **Flux balance, GLPK** | `glpk` | A genome-scale linear program solved in process, in every voxel, every step | no |
| **Flux balance, COBRApy** | `cobrapy` | The same, through an embedded Python interpreter and the reference implementation | no |
| **Surrogate network** | `surrogate` | `surrogateModel.hh`, a small network fitted offline to what the linear program returns | yes |
| **Symbolic law** | `symbolic` | A `.sym` text file of algebra, parsed at start-up | no |
| **Graph network** | `graphnet` | A `.gnn` file — message passing over the species–reaction graph | no |

Combined variants: `glpk_and_kinetics`, `cobrapy_and_kinetics`,
`surrogate_and_kinetics`, `symbolic_and_kinetics`, `graphnet_and_kinetics`.

### 5a. Flux balance analysis, in detail

| Feature | What it means |
|---|---|
| **max cᵀv s.t. Sv = 0, ℓ ≤ v ≤ u**, per voxel, per step | Growth is *predicted* from the genome, not prescribed. |
| **Two back ends** | GLPK in process (fast, recommended for production) and COBRApy in an embedded interpreter (slower, accepts the cobra ecosystem directly, useful as a cross-check). |
| **Warm-started simplex** | The previous voxel's basis seeds the next solve, which is most of why the GLPK path is usable at all. |
| **Michaelis–Menten uptake bounds** | ℓ = −V_max·C/(K_c + C), built from `<maximum_uptake_flux>` and `<half_saturation_constants>`. A zero V_max means supply-limited rather than enzyme-limited. |
| **Exchange reactions by name** | `<exchange_reaction_names>` resolves against the model and stops the run if a name is missing. Positional indices are still accepted but silently point elsewhere if the model is revised. |
| **Extra constraints** | Arbitrary reaction bounds, and `<equate_bounds>` for forcing reversible pairs. |
| **Model sources** | A local SBML file, a BiGG identifier, or one of the three shipped in `models/`. Downloads are off unless explicitly allowed. |
| **Unit conversion at one boundary only** | mmol/gDW/h inside the linear program, mol/L everywhere else, converted through `<biomass_molar_mass>`. |
| **Free versus total substrate basis** | With speciation on, `total` rebuilds each substrate from the equilibrium tableau, including every complex that carries it, and draws consumption back down across the complexes in proportion. |
| **Three LP algorithms** | Simplex (default), interior point, exact rational arithmetic. |
| **LP dump for debugging** | Write each microbe's problem as `.lp` or `.mps` at start-up. |

### 5b. The three learned paths

| | Surrogate network | Symbolic law | Graph network |
|---|---|---|---|
| **Form** | 2 inputs → 4×10 tanh → 1 linear, 465 parameters | An algebraic expression tree | Message passing on the bipartite species–reaction graph |
| **Fitted to** | The linear program's growth over a swept grid | Any table of data | Any table, with the stoichiometry supplied as *structure*, not learned |
| **Returns** | Growth rate | One rate, which you extend to the others by stoichiometry | The whole coupled rate vector at once |
| **Lives in** | `surrogateModel.hh`, compiled in | A `.sym` text file | A `.gnn` text file |
| **Readable?** | No | **Yes — it is algebra you can quote and disagree with** | No |
| **New one needs** | a rebuild | nothing | nothing |
| **Range handling** | none — see limitations | clamps and counts | clamps and counts |

The symbolic path has a full expression language: precedence, right-associative
exponentiation, unary minus, `exp` `log` `sqrt` `min` `max` `pow` `abs`, the
constants `e` and `pi`, division by zero returning zero rather than infinity,
and malformed input **rejected with a character position** rather than guessed
at. Every expression and its declared range is printed into the log, so a result
carries the law that produced it.

## 6. Aqueous speciation

| Feature | What it means |
|---|---|
| **Components-and-species tableau** | Master species, a stoichiometry row per substrate, and a formation constant per substrate. |
| **Solved in every pore voxel, every step** | By continued fractions. The most expensive part of the code by a wide margin. |
| **Coupled to the metabolic layer** | Through the free/total concentration basis described above. |
| **Species the solver ignores** | A row of zeros marks a mineral or biomass that reacts only kinetically. |
| **Stated limits** | Ideal activities (no Debye–Hückel or Davies), no temperature dependence, no gas phase, no redox couple, no saturation index. The log K values you supply are conditional constants at your own ionic strength and 25 °C. |

A complete 95-species, 17-component uranium tableau ships as a worked case.

## 7. Geometry evolution

The two capabilities that make the pore space a variable rather than a
constant. Both run alongside whichever rate path you chose.

### Mineral precipitation

| Feature | What it means |
|---|---|
| **Volume-of-Pixel filling** | An immobile mineral accumulates in a pore voxel until it reaches `<max_precipRho>`, derived from the mineral's molar volume. |
| **Surface-only growth** | `<surface_only>1` restricts precipitation to voxels touching a surface — heterogeneous nucleation, which is the physical case. |
| **Node conversion** | A full voxel becomes solid. `<perm_ratio>0` makes it an impermeable wall. |
| **Flow re-solved on the reduced pore space** | At `<update_interval>`, so permeability really falls. |
| **Clogging detection** | The run reports porosity falling in steps and stops when the pore space seals. |

### Mineral dissolution

| Feature | What it means |
|---|---|
| **Declared solid phases** | Any number of `<phaseN>` blocks, each with a name, material number, the substrate it releases, its full density and its initial fill. |
| **Consumed by the water touching it** | Averaged over the voxel's open faces. |
| **Mass-conserving release** | Products are deposited into the neighbouring pore with D3Q7 weights ¼ at rest and ⅛ in each of six directions — the seven populations sum to exactly the increment requested. |
| **Voxel reopening with hysteresis** | A voxel reopens below `<reopen_fraction>` × full, so a voxel on the threshold cannot flicker between states. |
| **Precipitate or original grain** | `<is_precipitate>` distinguishes a phase that formed during the run from one that was there at the start. |

## 8. Diagnostics and output

| Feature | What it means |
|---|---|
| **VTI fields** | Velocity, every substrate, every microbe, and the material map — named from your own `<name_of_substrates>`. |
| **Summary CSV** | One row per interval: porosity, and the total, mean, minimum and maximum of every substrate and microbe — **over open voxels only**, so a run that seals pore space is not divided by a moving denominator. |
| **Conservation checks** | Each `<conserve>` line names a sum the network cannot create or destroy (one mole of Fe²⁺ removed must appear as one mole of FeS). A drift beyond tolerance is reported *the moment it happens*. A check naming a substrate that does not exist is reported as SKIPPED, never as PASS. |
| **Ghost columns excluded** | Totals run over x = 1…nx−2, so a total does not jump when the inlet concentration changes. |
| **Per-iteration mass-balance printing** | Verbose, for when a run is going wrong. |
| **Kernel timing mode** | `<track_performance>` times the solver kernels and suppresses all output writing. |
| **Start-up echo** | The geometry, the enabled features and each organism's rate path are printed before the first step — so you can stop a wrong run in the first second rather than the second week. |

## 9. Running it

| Feature | What it means |
|---|---|
| **MPI parallel** | Palabos block decomposition. Serial runs need no MPI at all. |
| **Optional back ends are opt-in at build time** | `-DENABLE_GLPK=ON`, `-DENABLE_COBRAPY=ON`. If you switch one on in the XML and it was not compiled in, the solver **stops at start-up and names the build option** rather than quietly ignoring you. |
| **Binary checkpoints** | Written at an interval and reloaded on restart, so a long run survives a queue limit. |
| **Build options** | `ENABLE_GLPK`, `ENABLE_COBRAPY`, `ENABLE_MPI`, and `FBA_BULK_ONLY` — which restricts the FBA processors to the bulk domain for 20–40% more speed. The surrogate, symbolic and graph-network paths need no build option at all. |
| **One case, one directory** | `scripts/setup_case.sh` assembles everything a case needs, including the sources, so it builds on its own. |

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

| | Question | If yes | Why |
|---|---|---|---|
| 1 | Is a rate law already known? | **compiled kinetics** | Cheapest by four orders of magnitude, and the easiest to read |
| 2 | Do you need the internal fluxes, not just growth? | **flux balance analysis** | The only path that says what the organism is doing inside |
| 3 | Is growth the only output, with one or two substrates limiting? | **surrogate network** | Reproduces the linear program's growth at about 1/400 of its cost |
| 4 | Should the result be readable and arguable? | **symbolic law** | The answer is algebra, and a new one needs no rebuild |
| 5 | Otherwise: many species coupled through many reactions | **graph network** | The whole rate vector in one evaluation, in stoichiometric ratio |

Precipitation and dissolution are not on this list. They are geometry
processes and run alongside whichever rate path you choose.

---

## Start here

```bash
git clone https://github.com/shahram444/CompLaB3D-GEMS.git
cd CompLaB3D-GEMS

# assemble a runnable case, build it, run it
./scripts/setup_case.sh 13_precipitation run/mycase
cd run/mycase
cmake -B build -S . && cmake --build build -j
./build/complab CompLaB.xml
```

That case forms FeS where an iron front meets a sulfide front, seals the voxels
as they fill, and reports the porosity falling in steps until the pore clogs.

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

| Your rate path | Offline work needed |
|---|---|
| Compiled kinetics | none |
| Flux balance analysis | B1 — export the metabolic model |
| Surrogate network | B1, then B2 — sweep, fit, paste, rebuild |
| Symbolic law | B3 — run the search, choose from the Pareto set |
| Graph network | B4 — train, write the `.gnn` |

---

## The sixteen examples

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

Full details, and which geometry each one runs on, in
[`examples/README.md`](examples/README.md).

---

## Layout

```
CompLaB3D-GEMS/
├── src/                the solver: Palabos data processors and headers
├── config/
│   ├── CompLaB.reference.xml     every tag, annotated, in one place
│   ├── kinetics/                 the shared chemistry headers
│   └── geometry/                 the four shared pore geometries
├── examples/           sixteen cases, from flow-only to the full pipeline
├── pipelines/          what to run before and after the solver
├── tools/              the offline tools: model export, fitting, training
│   └── surrogate/      the surrogate training path, Python and MATLAB
├── tests/              the regression suite, and the structural checks
├── models/             three genome-scale models, and the FBA toy model
├── scripts/            setup_case.sh, which assembles a runnable case
└── docs/               the user guide, and how to publish this repository
```

### The offline tools

| Tool | What it does |
|---|---|
| `tools/geometry.py` | Build a pore space, or inspect one: porosity, percolation, isolated pore |
| `tools/extractMM.py` | Export an SBML or BiGG model to the tabular form the GLPK path reads, and check it can grow |
| `tools/makeKinetics.py` | Generate a `defineKinetics.hh` from a reaction list |
| `tools/makeEquilibrium.py` | Build the components/stoichiometry/log K tableau |
| `tools/surrogate/` | Sweep the linear program, fit the network, verify the export reproduces the trainer, and plot the response surface — in Python and in MATLAB |
| `tools/fit_symbolic.py` | Genetic-programming search returning a Pareto set of rate laws |
| `tools/train_graphnet.py` | Train a graph network and write a `.gnn` |
| `tools/postprocess.py` | Slices, histories, and the mass-balance report |
| `tools/vtireader.py` | Read VTI output into NumPy |
| `tools/complab3d_cobrapy.py` | The Python side of the COBRApy back end |

### One thing about the examples

**Nothing in this repository exists twice.** An example directory holds only
what is unique to it — its `CompLaB.xml`, its README, and a chemistry header
where its chemistry genuinely differs. Everything an example shares with another
example lives once, under `config/` or `models/`, and the example names it in a
short `case.files`:

```
config/geometry/slot_one_biofilm.dat  input/geometry.dat
models/toy_model.xml                  input/toy_model.xml
```

That is what removed thirty-nine duplicate files: twenty-three identical
kinetics headers, twelve identical geometries and two identical metabolic
models. Changing the chemistry interface now means editing two shared defaults
and six real overrides, not thirty-two files.

`scripts/setup_case.sh` reads `case.files`, lays the shared halves down, puts
the example's own files on top, and leaves you a directory that builds. This is
why you assemble a case rather than `cd` into it.

---

## Configuration

There is **one** annotated reference file,
[`config/CompLaB.reference.xml`](config/CompLaB.reference.xml), documenting every
tag the solver reads — with its units, its default, and whether it is required.
Each example's own `CompLaB.xml` is a working subset of it, and there is no
second copy of the documentation anywhere in the tree.

Units are one convention throughout: **micrometres, mol/L, seconds, m²/s**.
Biomass is in mol/L, the same unit as the chemicals, on purpose. Flux balance
analysis is the single exception — it works in mmol/gDW/h internally because
that is what published metabolic models use, and converts at the boundary of the
linear program through each organism's `<biomass_molar_mass>`.

The blocks this repository adds to the base model:

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

## Known limitations

Stated plainly, because finding these out yourself is expensive.

**Dissolution products are lost at MPI block boundaries.** A dissolving voxel
deposits its products into its open face neighbours; when a neighbour lies in a
different MPI block the deposit is written into the local envelope and is not
communicated back. Dissolution results therefore depend on the processor count,
with the discrepancy concentrated at block interfaces. A run needing a
quantitative dissolution mass balance should be done on a single process, or
checked against one at reduced resolution. Precipitation is unaffected.

**The COBRApy path clamps positive exchange draws to zero** and does not run the
repair loop the GLPK path uses. Expect the two back ends to agree on uptake and
growth, and to differ where a model excretes strongly. GLPK is the recommended
path for production. GLPK and COBRApy organisms cannot be mixed in one run.

**The surrogate path does not enforce its training range.** A voxel outside the
fitted box gets a confident answer and no warning. The symbolic and graph-network
paths clamp and count; the surrogate does neither. Plot the response surface with
`tools/surrogate/inspectSurrogate.py` before a production run — a large region of
a fitted box returning zero growth is common and invisible in the weights.

**Speciation is ideal and isothermal.** No Debye–Hückel or Davies correction, no
temperature dependence, no gas phase, no redox couple, no mineral saturation
index. The log K values you supply are conditional constants at your own ionic
strength and 25 °C.

**The `FD` biomass solver is not exercised by any shipped case.** Verify it on a
small run before relying on it. `CA` and `LBM` both have examples.

**The shipped `defineKinetics.hh` expects the 95-substrate uranium network** and
indexes `C[]` by hard-coded position. Edit it to match your own substrate list
before enabling compiled kinetics, or it will read past the end of the array.

---

## Working on it

| | |
|---|---|
| [`docs/PUBLISHING.md`](docs/PUBLISHING.md) | Putting this on GitHub, step by step, and the edit loop afterwards |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | The rule the tree is built on, and what to update when you change the chemistry interface |
| [`docs/README_MANUAL.md`](docs/README_MANUAL.md) | Building the user guide from its LaTeX sources |
| [`CHANGELOG.md`](CHANGELOG.md) | What this release added to the base model |

---

## Authors

| | |
|---|---|
| **Shahram Asgari** | Department of Marine Sciences, University of Georgia, Athens, GA, USA — <shahram.asgari@uga.edu> |
| **Christof Meile** | Department of Marine Sciences, University of Georgia, Athens, GA, USA |

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
