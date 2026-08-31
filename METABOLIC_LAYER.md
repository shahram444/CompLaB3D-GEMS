# The optional metabolic layer in CompLB3D

CompLB3D shipped with hand-written Monod kinetics only. This change adds back the
three metabolic solvers CompLaB 2D v1.0.0 had — GLPK flux balance analysis,
COBRApy flux balance analysis, and a pre-trained surrogate model — plus the
`extractMM.py` model converter, as **optional** modules in **separate files**.

A run that asks for none of them behaves exactly as CompLB3D did before.

---

## Quick start

**Build.** Both back ends default to OFF, so the code still compiles on a machine
with neither GLPK nor Python installed:

```bash
cmake ..                                  # kinetics + surrogate only
cmake -DENABLE_GLPK=ON ..                 # + GLPK flux balance analysis
cmake -DENABLE_COBRAPY=ON ..              # + COBRApy flux balance analysis
cmake -DENABLE_GLPK=ON -DENABLE_COBRAPY=ON ..
```

CMake locates GLPK and the Python headers itself and stops with a specific
message if it can't find them.

**Convert a metabolic model.**

```bash
python3 extractMM.py iAF987.sbml          # writes iAF987_mm.xml, prints a report
cp iAF987_mm.xml input/
```

The report tells you the objective reaction index and the exchange reaction
indices you need for the XML.

**Switch it on** in `CompLaB.xml`:

```xml
<simulation_mode>
    <enable_fba_glpk>true</enable_fba_glpk>
</simulation_mode>
...
<microbe0>
    <reaction_type>glpk</reaction_type>
    <model_filename>iAF987_mm</model_filename>
    <exchange_reaction_indices>62 121 -1</exchange_reaction_indices>
    <biomass_molar_mass>24.6</biomass_molar_mass>
    <half_saturation_constants>1e-4 2e-5 0</half_saturation_constants>
    <maximum_uptake_flux>10.0 0.5 0</maximum_uptake_flux>
</microbe0>
```

If you ask for a solver the executable wasn't built with, CompLB3D stops at
start-up and names the CMake option to turn on. It never ignores you silently.

---

## Two switches, and why

**Build switch** (`cmake -DENABLE_GLPK=ON`) — decides whether the code is put
into the program at all. This is what lets CompLB3D build on a machine that has
no GLPK.

**Run switch** (`<enable_fba_glpk>` in the XML) — decides whether a program that
*was* built with it uses it on this particular run. This is the one you change
day to day.

---

## New files

| File | Lines | What it is |
|---|---|---|
| `src/complab3d_metabolic.hh` | 918 | Config struct, XML parsing, model loading, the unit conversion, the free/total concentration machinery |
| `src/complab3d_processors_fba.hh` | 641 | `runFBA_glpk3D` and `runFBA_cobrapy3D` |
| `src/complab3d_processors_surrogate.hh` | 263 | `run_surrogate3D` |
| `src/complab3d_glpkcpp.hh` | 958 | GLPK wrapper (guarded by `COMPLAB_ENABLE_GLPK`) |
| `src/complab3d_pythonAPI.hh` | 716 | Embedded CPython bridge (guarded by `COMPLAB_ENABLE_COBRAPY`) |
| `src/complab3d_cobrapy.py` | 303 | The Python side |
| `surrogateModel.hh` | 503 | User-editable surrogate template (root, next to `defineKinetics.hh`) |
| `extractMM.py` | 396 | SBML/MATLAB/JSON → XML converter |

Changed: `src/complab.cpp`, `src/complab_functions.hh` (reaction-type vocabulary
only), `src/complab3d_processors.hh` (three includes), `CMakeLists.txt`,
`CompLaB.xml`.

The old uranium example is preserved unchanged as
`CompLaB_example_uranium_zhao.xml`.

---

## Units

CompLB3D keeps every concentration and every biomass pool in **mol/L**. Flux
balance analysis is posed in **mmol/gDW/h** against a biomass in gDW/L.

Rather than change CompLB3D's convention, each metabolic microbe carries
`<biomass_molar_mass>` (gDW per mole of cells, default 24.6 for
CH₁.₈O₀.₅N₀.₂) and the conversion happens at exactly one place:

```
B[gDW/L]  = B[mol/L] × biomass_molar_mass
dC[mol/L] = v[mmol/gDW/h] × B[gDW/L] × dt[s] / 3.6e6
dB[mol/L] = mu[1/h] × B[mol/L] × dt[s] / 3600
```

These factors appear only in `MetabolicConfig::fluxToDeltaC()` and
`growthToDeltaB()`. If a unit ever looks wrong, look there and nowhere else.
`<maximum_uptake_flux>` is the only XML tag in mmol/gDW/h.

---

## Free vs total concentration

With the equilibrium solver on, the acetate lattice holds only *free* acetate —
the rest is in CaAc⁺, FeAc⁺, MgAc and so on. Feeding the free value to a
Michaelis–Menten bound understates what the organism can reach.

```xml
<fba_concentration_basis>free</fba_concentration_basis>    <!-- default -->
<fba_concentration_basis>total</fba_concentration_basis>
```

`total` rebuilds each component's total **generically, from the equilibrium
stoichiometry matrix**, and draws consumption back down across the complexes in
proportion — element-exact. This is deliberately better than what
`defineKinetics.hh` does, which hard-codes `C[56]`, `C[81]..C[86]` and silently
corrupts itself if anyone reorders the substrate list.

Default is `free` because it reproduces CompLaB 2D exactly. **With speciation on,
`total` is almost certainly the physically correct choice** — this is the one
design decision I'd want you to make deliberately rather than inherit.

---

## Bugs fixed during the port

Each carries a `// [FIX-3D]` comment at the site naming the 2D line it came from.
65 markers in total. The ones that change results:

**Flux balance**
- `complab_processors.hh:130` — `if (fixLB && c2f < conc)` is true in essentially
  every cell (`c2f ≤ 0 ≤ conc`), so `fix_lower_bounds` silently overrode the
  Michaelis–Menten bound everywhere. The COBRApy twin at line 366 had the correct
  comparison; both are correct now.
- The COBRApy **flux index mismatch** — C++ read `vec2_flux[iB][iS]` over all
  substrates while Python returned a vector packed over that microbe's *used*
  exchanges. Correct only when every microbe used every substrate. The packing
  contract is now pinned in a comment block in both files, length-checked at run
  time, and unpacked through `expandFluxToSubstrates()`.
- `complab_processors.hh:204` — the feasibility repair only fired when
  `conc > thrd`, so a substrate already floored to zero accepted an unlimited
  negative draw. 
- **No positivity clamp existed anywhere** in the increment application path. The
  commented-out in-place kinetics code had one; the live delta code did not.
  Added, in both FBA processors and the surrogate.
- `complab_processors.hh:357` — the per-microbe bound vector was declared outside
  the microbe loop and never cleared, so microbe *i* inherited every earlier
  microbe's bounds.
- `complab_processors.hh:377` — `vec2_maxRelease[bmLattice[bloc]]` indexed the
  compacted list with an already-compacted index.
- `run_glpk`'s success set omitted 0, which is what `glp_interior`/`glp_exact`
  return on success — `lpsolver` 2 and 3 returned uninitialised flux vectors.
- `initialize_glpk`'s version check was broken twice over (`sizeof` on a pointer;
  the warning fired when the version *matched*).
- `memset(c, 0, ncol*sizeof(int))` on a `double[]` zeroed half of each buffer.
- The triplet arrays were `calloc`'d at the dense worst case `nrow*ncol+1` and
  never freed — ~66 MB per genome-scale model, now ~80 kB.
- Every `exit(EXIT_FAILURE)` inside a data processor replaced with a return code.
  One rank calling `exit()` under MPI hangs the whole job.

**COBRApy**
- `Py_DECREF(pArgs)` was commented out because `PyTuple_SetItem` steals the model
  references — leaking one tuple plus every boxed float on **every solve at every
  voxel at every timestep**. Fixed by increfing the models before the steal.
  Verified: model refcount 1 before and after 300 solves, flat RSS at 5000.
- The module and `run_FBA` were re-imported on every call; now cached.
- Model construction was O(nmets × nrxns) including zero coefficients; now sparse
  and batched.
- Infeasibility (`objective_value is None`) reached `PyFloat_AsDouble`
  unchecked. Now returns a per-microbe status the C++ side honours.

**Surrogate**
- The **index-0 hardwiring**: 2D built compacted per-voxel arrays and the shipped
  network read `Fin[0][*]` / wrote `bioR[0]`, so "microbe 0" meant "whichever
  organism was listed first in *this* voxel". With more than one surrogate
  organism the network was applied to an arbitrary one and the rest silently
  became pure decay. The signature now takes an explicit `microbeId` and the
  arrays are always full length, indexed by CompLaB.xml microbe number.
- 2D passed `flag = 0`, disabling the depletion limiter, so the estimate could
  demand more substrate than the voxel held. The limiter is unconditional now.

**extractMM.py**
- `objloc` was undefined if no reaction had objective coefficient exactly 1.0,
  causing a `NameError` 30 lines later.
- **An SBML input named `model.xml` overwrote itself** with its own output.
- `"%g"` truncated everything to 6 significant digits; now `%.17g`.
- `b` was built with `nrxn` entries but read with `nmet` — out of bounds whenever
  `nmet > nrxn`.
- Added a round-trip check that reparses the written `S` and compares.

**Pre-existing 3D bug, fixed because the new code touches the same lattices**
- `complab.cpp:993` reset `dC[]` once, above the biotic block. The abiotic block
  then accumulated into the **same** lattices with no intervening reset, so
  `update_abiotic_rxnLattices` applied the biotic contribution a **second time** —
  silently doubling every biotic reaction in any run with both
  `enable_kinetics` and `enable_abiotic_kinetics` on, which is the shipped
  uranium configuration. A reset now precedes the abiotic block.

---

## Post-delivery review: what was found and fixed

The first delivery was reviewed adversarially, module by module, with reviewers
told to prove bugs by compiling and running rather than by reading. They found
real defects — including two that meant neither FBA back end could execute a
single voxel. All are fixed and each fix was re-verified by test.

### Crashes (the code could not run at all)

- **`run_glpk` wrote into zero-length vectors.** It writes `xmin[0..ncol-1]`,
  `lambda[0..nrow-1]` and `redcosts[0..ncol-1]` but never sized them; the 2D
  caller happened to size them at its call site and the 3D port did not. First
  voxel, guaranteed segfault (proven under ASAN). `run_glpk` now sizes them
  itself — the contract belongs with the code that decides how much it writes.
- **`cfg.vec_model` was never sized.** `prep_cobrapy` writes `vec_model[j]` into
  it. Every other per-microbe vector was allocated in `initialize_metabolic`;
  this one was missed. Heap write past a zero-capacity allocation on the first
  model. Now sized before the call.

### Wrong results

- **`<maximum_uptake_flux>` was never read.** The documented Vmax tag was dead;
  the code used `<substrate_lower_bounds>` — a hard floor, not a rate — in its
  place. A user who set Vmax to 4.05e-6 and no lower bound got the
  supply-limited fallback instead. There is now a proper `vmax` field fed by
  that tag, `maxUptake` is a clamp only, and all three solvers gate on the same
  predicate.
- **The uptake branch was gated on the release ceiling.** Supplying
  `<substrate_lower_bounds>` without `<substrate_upper_bounds>` silently
  bypassed the Monod attenuation entirely, and the GLPK and COBRApy paths
  disagreed on identical input.
- **Two organisms in one voxel could drive a substrate negative.** The COBRApy
  and surrogate paths clamped per organism against an `avail` that was never
  decremented, so N organisms each took the whole supply. Both now sum first and
  clamp once, like the GLPK path. Also fixed for the `total` concentration
  basis, where a complex carrying two components was debited twice.
- **The surrogate depletion limiter was mathematically a no-op.** `v` was
  already `-maxDraw·mm` with `mm ∈ [0,1]`, so the test `v > -maxDraw` could
  never be true. Swept over thousands of states it fired zero times; it now
  fires in the 135 of 2625 states where it should, and never clamps low.
- **`MetabolicTotals` broke on any tableau with negative rows.** A proton
  condition (TOTH) is negative in alkaline water; `total()` discarded it and
  `drawDown` then divided by the free concentration instead. Measured: a request
  to remove 1 nmol/L of the H component removed 10 µmol/L of ammonia, moving the
  balance 10,000× the wrong way. The guard is also relative now, not an absolute
  1e-30 — a near-cancelling total used to amplify a 1e-10 request into a 1e-4
  change in two other substrates.
- **`GLP_EFAIL == GLP_OPT == 5`**, so a failed interior-point or exact solve was
  accepted as optimal and the caller read a growth rate the solver never wrote.
  Those failure codes are offset now and `status` is set to a value that cannot
  be mistaken for success.
- **`<objective_direction>` reached GLPK but not COBRApy** — a microbe asking to
  minimise was silently maximised. Since the two back ends cannot be mixed,
  nothing would have caught it by cross-checking. It now travels on the wire.
- **`extractMM.py` wrote `inf`**, which C++ stringstream extraction cannot
  parse: the list came back short and the run aborted with a length error that
  said nothing about the cause. It writes a ±1e30 sentinel now, and the GLPK
  wrapper treats that magnitude as a free bound.
- **Decay ran 3600× too slow.** `<decay_coefficient>` is per second everywhere
  in CompLB3D, but the fallback was pushed through the per-hour growth
  converter. There is a separate `decayToDeltaB` now.

### Robustness

- `vec_Kc` / `vec_mu` are validated before the metabolic layer indexes them. The
  existing parser prints a message on a wrong-length
  `<half_saturation_constants>` but does not stop, and pushes a one-element
  `{-99}` row when the tag is missing — the metabolic layer is the first code to
  actually read that vector, and would have read past the end of the short row
  for every substrate above the first.
- Bound-vector lengths are checked rather than silently truncated.
- Unrecognised on/off text (`ture`, `y`, `enabled`) is now an error, not a silent
  "off".
- `<reaction_type>1</reaction_type>` is accepted again — the numeric spellings
  the comment promised still worked were missing from the table.
- `exhaustionFlux` is capped, so a voxel barely above the biomass threshold no
  longer hands GLPK a column spanning twenty orders of magnitude.
- `MetabolicScopeGuard` releases the solvers on **every** exit path from main().
  Previously one of fourteen early returns did so.
- `MetabolicConfig` copying is forbidden — it owns raw `glp_prob*`, `char*` and
  `PyObject*` handles with no destructor.
- The interior-point solver no longer prints its full iteration log per voxel.
- `GLP_FX` is used when lower and upper bounds coincide (`<equate_bounds>` makes
  that reachable), instead of `GLP_DB`, which returns `GLP_EBOUND`.
- `extractMM.py` now prints the exchange-reaction index table you need for
  `<exchange_reaction_indices>` — the one manual step where an off-by-one is
  silent.
- A warning fires when only some microbes use kinetics, because
  `defineRxnKinetics` then receives a compacted biomass vector while the shipped
  `defineKinetics.hh` indexes it by global microbe number.
- The dead `infeasibleCount` (whose comment claimed the run summary reported it)
  now drives a rate-limited message, so an infeasible LP is visible.
- CMake falls back to `FindPythonLibs` below CMake 3.12 instead of failing with a
  message that reads as a missing dependency.

### Known and deliberately left

- **`<maximum_uptake_flux>` is read by two parsers in two different units** —
  `defineKinetics.hh` wants mol/mol cells/s, FBA wants mmol/gDW/h. A microbe
  using a combined reaction type would have to satisfy both with one number.
  There is now a dedicated `<fba_maximum_uptake_flux>` that takes precedence,
  and a warning when the shared tag is used by a combined-type microbe.
- **The joint clamp is per processor, not per voxel.** GLPK and surrogate
  microbes may both be active; each clamps independently against the same
  undecremented `avail`. Two organisms of *different solver types* sharing a
  substrate can still overdraw it. Fixing this properly means a shared
  per-voxel budget across processors, which is a larger change than this port.

## What I verified, and what I could not

**Verified here (after the fixes above)**

- Every fix re-verified by an independent pass: 12 of 12 PASS, including ASAN
  harnesses for the two crashes, a 600-trial randomised fuzz over the reaction
  processors (0 negative concentrations, 0 NaN), a 2625-state sweep of the
  surrogate limiter, 12,000 COBRApy solves with stable refcounts and flat RSS,
  and the surrogate network still bit-for-bit identical to the 2D original.

**Also verified**

- All four build configurations (`none`, `GLPK`, `COBRApy`, both) compile clean
  under `g++ -std=c++11 -Wall -Wextra` against a minimal Palabos stub, with every
  new template explicitly instantiated. Also clean with `-DFBA_BULK_ONLY`.
- `complab3d_glpkcpp.hh` compiles against GLPK 4.39, 4.65 and 5.0 stubs, and with
  no `glpk.h` present at all.
- The COBRApy bridge was run end-to-end against real CPython 3.11 + cobrapy
  0.31.1: 300 solves, refcount stable, RSS flat, feasible and infeasible paths
  both correct.
- `python3 src/complab3d_cobrapy.py` self-test passes.
- `extractMM.py` exercised on five paths including the self-overwrite refusal;
  round-trip PASS.
- The ported surrogate network reproduces the 2D original **bit-for-bit** on six
  test inputs.
- Both XML files are well-formed (the template originally was not — XML forbids
  `--` inside comments, which is easy to hit when drawing separator lines).
- Brace, paren and bracket balance across every edited file.

**Not verified — you will need to do this on Tahoma**

- **Nothing was compiled against real Palabos.** The stub covers the API surface
  the new code uses, but a full build is the real test. Expect the first
  compile to surface small signature mismatches.
- **Nothing was run.** No numerical result here has been checked against
  anything.
- **The `<biomass_molar_mass>` default of 24.6 g/mol is a placeholder.** Set it
  per organism before trusting any growth rate.
- **`fba_concentration_basis` defaults to `free`.** If you run FBA together with
  speciation, decide deliberately whether it should be `total`.
- **The GPL / AGPL licence question on `complab3d_glpkcpp.hh` is unresolved.**
  That file is derived from `glpkcc.cpp` (GLPKMEX / Octave-Forge, also shipped in
  the COBRA Toolbox), which is GPL. CompLaB carries AGPL headers. The provenance
  is now stated prominently at the top of the file with three resolution options,
  but this needs you and a person, not a code change.

**Suggested first run**: build with no options, point `CompLaB.xml` at a small
geometry, and confirm the flow-plus-transport template runs. Then turn on GLPK
with one microbe and a small model before scaling up.

---

## Things I deliberately did not change

- `initialize_complab()` still takes 60+ arguments by reference. The metabolic
  settings are parsed by a **separate** `initialize_metabolic()` reading the same
  file, so the existing parser is untouched and this whole feature can be removed
  by deleting three lines from `complab.cpp` and three includes.
- The FBA processors default to `bulkAndEnvelope`, matching `run_kinetics`
  exactly. That solves a linear program for every envelope cell, which is waste —
  `-DFBA_BULK_ONLY=ON` switches it off for a typical 20–40% saving. It is opt-in
  because it changes a correctness-critical setting and deserves a side-by-side
  check on your own case first.
- GLPK and COBRApy microbes still cannot be mixed in one run, as in 2D.
- The reaction-type numbering keeps `0 = none` and `1 = kinetics` from the
  existing 3D code, so existing input files are unaffected. The new solvers are
  appended as 2–7. Never renumber 0 or 1.

---

# Mineral dissolution (added after the metabolic layer)

The reverse of the precipitation feature: a mineral can now be eaten away until
its voxel opens back up to flow. Independent on/off switch, so you can run
precipitation only (exactly the old behaviour), dissolution only, both, or
neither.

## The idea: a "solid phase"

CompLB3D has two kinds of solid voxel and could not tell them apart — original
grains from `geometry.dat`, and voxels that *became* solid because a mineral
precipitated. That matters the moment dissolution exists: if the rule were "a
solid voxel with no mineral left becomes pore", an inert quartz grain is already
below any threshold and the whole grain pack would dissolve on step one.

So dissolution works on **declared phases**. A phase is four numbers:

| | precipitated FeS | calcite grain | inert quartz |
|---|---|---|---|
| `material_number` | its own code | 0 | 1 |
| `substrate` (holds the inventory) | 9 | 13 | — |
| `full_density` mol/L | 48.9 | 27.1 | — |
| `initial_fill` | 0, grows by reaction | 27.1, starts full | — |

Precipitate and rock become the *same* mechanism, differing only in how the
inventory is seeded. Anything not declared can never dissolve — quartz is safe by
construction, not by a special case.

## How a voxel is tracked

A new `phaseLattice`, alongside the existing `maskLattice` / `ageLattice` /
`distLattice`. It records which phase occupies each voxel, 0 for none: seeded
from geometry at start-up, stamped by `precipNodeConversion3D` when a voxel
seals, cleared by `dissolNodeConversion3D` when one reopens.

The **mask value of a sealed voxel stays `solid`**. That was deliberate — every
existing processor tests `mask == solid` to decide what to skip, and a new mask
code would have meant auditing all of them. The phase lattice is purely additive
and invisible to code that doesn't ask for it.

## Where the chemistry lives

Not in the new file. Exactly as precipitation takes its rate from whatever you
wrote in `defineAbioticKinetics.hh`, dissolution calls one hook in that same file:

```cpp
void defineDissolutionRate(plint phaseId, const std::vector<double>& C,
                           double mineral, std::vector<double>& subsR,
                           double& mineralR, plint mask);
```

You get the water *in contact with* the mineral (the open neighbours, averaged —
a solid voxel has no water of its own that moves) and how much mineral is left.
You return how fast it dissolves and what that releases. A worked FeS + H⁺ →
Fe²⁺ + HS⁻ example is in the file; the shipped default does nothing, so turning
dissolution on without editing it is harmless.

## Where the products go

Into the **open face neighbours**, split evenly. Anything released inside the
solid voxel would sit in water that never moves. A voxel with no open neighbour is
fully enclosed, unreachable by water, and does not dissolve at all.

This is the one piece with no precipitation counterpart, and it's where mass
conservation lives.

## Hysteresis

A voxel seals at `inventory >= full_density` and reopens below
`reopen_fraction × full_density`, default 0.9. With a single threshold a voxel
sitting on it would flip every update interval — and every flip triggers a full
Navier-Stokes re-solve. The gap prevents that. The code refuses to start if
`reopen_fraction` is not strictly between 0 and 1.

## Stability

One step may never remove more than half the mineral present, and when the
limiter bites it scales **everything you returned** by the same factor — so your
element balance survives it exactly. You don't need your own cap.

## Verified

A lattice-backed driver runs the real processors through a complete cycle: a
voxel sealed at 48.9 mol/L FeS, dissolving in acid, over 40,000 steps.

```
start: FeS=48.900000  mask=0 (solid)  phase=1
  FeS left      = 44.008146  (reopen threshold 44.0100)
  mask          = 2   (REOPENED)
  phase         = 0   (cleared)
  reopened at step 2107
  worst relative mass-balance error over all steps = 9.073e-13
  Fe(mineral) + Fe(dissolved) = 48.900000000   (started at 48.900000000)
  conservation error          = 1.104e-14
PASS
```

Note the 2108 dissolving steps out of 40,000: once the voxel reopened and its
phase was cleared, dissolution correctly stopped. Total iron is conserved to
1.1e-14 — the residual is floating-point representation, not a leak. Clean under
AddressSanitizer and UndefinedBehaviorSanitizer.

## Files

New: `src/dissolutionVOP.hh` (~470 lines) — `SolidPhase`, `DissolutionConfig`,
`initPhaseLattice3D`, `surfaceDissolutionKinetics3D`, `dissolNodeConversion3D`.

Changed: `defineAbioticKinetics.hh` (the new hook, with the worked example),
`src/precipitationVOP.hh` (`precipNodeConversion3D` gained an optional phase-id
argument, defaulted to 0 so the old two-lattice call still compiles),
`src/complab.cpp` (parse, phase lattice, two loop calls), `CompLaB.xml`.

## Not done

- **No saturation-state driver.** Real mineral dissolution follows
  `R = k·A·(1 − IAP/Ksp)`, but the equilibrium solver computes no solubility
  products, so it cannot supply the saturation state. First-order in the reactant
  is what the code can support today. Adding `Ksp` and an ion activity product to
  the tableau would be a separate piece of work.
- **No reactive surface area model.** The rate sees how much mineral is left,
  which is a crude proxy. A proper `A = A₀(V/V₀)^(2/3)` form would be a small
  change to the hook if you want it.
