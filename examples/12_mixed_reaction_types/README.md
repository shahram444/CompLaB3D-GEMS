# 12 — one organism, two rate paths: flux balance plus a compiled law

## 1. What this example does

A metabolic model describes what an organism's network *can* do. It does not
describe everything the organism actually does. The classic omission is
**maintenance**: the donor a cell consumes simply to stay alive — to hold its
membrane potential, turn over proteins, run its pumps — which yields no biomass
and appears in no flux-balance objective. Leave it out and the model
over-predicts growth wherever substrate is scarce, which in a pore space is
most of the domain.

This case runs **one organism through two rate paths at once**. Flux balance
analysis supplies the growth, exactly as in example 09; a compiled rate law in
`kinetics/defineKinetics.hh` supplies the maintenance draw the linear program
knows nothing about. The two rates are added.

That is what the combined reaction types are for:

```
    glpk_and_kinetics        cobrapy_and_kinetics        surrogate_and_kinetics
```

Each names a metabolic back end **plus** the compiled kinetics header, per
organism. This case uses `glpk_and_kinetics`.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 3 of them: `S`, `O`, `P` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Toybug` — flux balance analysis (GLPK) **plus** compiled kinetics |
| **Biomass** | `Toybug` — attached biofilm (seeded on its own material number), finite difference — biomass spreads by diffusion on the same grid |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 500 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The two rates, and how they are kept from colliding

**From the linear program** — the same four-reaction toy model examples 09 and
10 use, with the uptake bounds rebuilt each step from the local concentrations:

```
    maximise   c' v          the biomass reaction's flux
    subject to S v  = 0
               lb <= v <= ub  with  V = Vmax * C / ( Kc + C )
```

**From the compiled kinetics** — a maintenance drain, first order in biomass and
zero order in substrate as long as any is present:

```
    r_maint  =  m_S * B                     mol of donor per litre per second

    m_S = 2.0e-6    mol donor per gDW per second, yielding no biomass
```

capped at what the voxel actually holds, so a step cannot draw a substrate
negative.

**The sign discipline is the whole difficulty of a combined type.** The kinetics
header in this folder writes `subsR` — the substrate rates — and deliberately
leaves `bioR` at zero:

```cpp
    subsR[0] = -std::min(drain, S);

    // bioR deliberately left at zero: growth is the FBA path's job.
```

Growth belongs to the flux-balance path. Writing it in both places double-counts
it, and nothing complains: the run completes, the numbers look plausible, and
the organism grows at twice the rate the model says. That is the classic mistake
with combined reaction types, and the reason the comment is in the header rather
than only here.

**Two Vmax tags, because the two paths read different units.** This trips people
up often enough to be worth stating on its own:

| tag | read by | units |
|---|---|---|
| `<maximum_uptake_flux>` | the compiled kinetics path | mol per mol cells per second |
| `<fba_maximum_uptake_flux>` | the metabolic path | mmol gDW⁻¹ h⁻¹ |

One number cannot mean both things, so a combined microbe carries both tags.
Here `<maximum_uptake_flux>` is `0. 0. 0.` — the maintenance law does not use it,
it uses its own `m_S` — while `<fba_maximum_uptake_flux>` is `10. 4. 0.`, the
same bounds as example 09.

| tag | value | what it does |
|---|---|---|
| `<enable_kinetics>` | `true` | compiles and calls `defineKinetics.hh` |
| `<enable_fba_glpk>` | `true` | builds the in-process linear-program path |
| `<reaction_type>` | `glpk_and_kinetics` | routes this organism through **both** |
| `<decay_coefficient>` | `0.01` | first-order biomass decay, separate from maintenance |

Maintenance and decay are not the same thing and both are present here.
Maintenance is a **substrate** draw that produces no biomass; decay is a
**biomass** loss that returns nothing to the solutes. A cell pays the first while
alive and suffers the second regardless.

## 4. What to check

1. **growth is lower than example 09's**, on the same geometry and the same
   bounds, because the donor is now being spent on maintenance as well;
2. **the gap widens where substrate is scarce.** Deep in the domain, where
   diffusion delivers little, maintenance takes a larger share of a smaller
   supply — that is the whole reason to include it;
3. **the donor drops faster than the acceptor.** Only `S` carries a maintenance
   term, so its balance against biomass no longer follows the objective
   reaction's stoichiometry, while `O` still does;
4. **set `m_S` to zero and example 09 comes back.** That is the cleanest test
   that the two paths are being added rather than one silently replacing the
   other;
5. **no `[NEG!]` warnings.** The maintenance draw is capped at what the voxel
   holds and the solver clamps again downstream.

## 5. What this case does not do

No flow, no abiotic reaction, no geometry evolution, no thermodynamic control.
The maintenance law here is the simplest one that makes the point — constant per
biomass, on one substrate. A real maintenance term usually carries its own
half-saturation and often a temperature dependence; both go in the same header,
in the same place.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

**[v1.3] One prerequisite does live outside the folder: GLPK.** This case runs
its flux balance path through GLPK, so the GLPK development package must be
installed (INSTALL.md names it `libglpk-dev`) and the build must be configured
with `-DENABLE_GLPK=ON`. `pipeline.sh` already passes that flag; a hand build
does not, because the CMake option is OFF by default. If it is missed, the
failure is not a compile error: the build succeeds, and the solver then
terminates immediately at start-up with `<enable_fba_glpk> is true but this
executable was built without GLPK`.

```bash
./scripts/setup_case.sh 12_mixed_reaction_types run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and reports porosity, refusing to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Exports the metabolic model into the flat format the solver reads. Run before the first build. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's checks. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `input/toy_model.xml` | The metabolic model, in the flat `<Metabolic_Model>` format. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

### The rate law this case ships

`kinetics/defineKinetics.hh` **replaces** the shared header for this case:
`setup_case.sh` copies it to `defineKinetics.hh` at the case root, over the
default it laid down there, so the maintenance term
of section 3 is what gets compiled in. That is the only file in this folder that
changes what the solver computes.

### A note on the domain

`x = 0` and `x = nx-1` carry the boundary conditions named per substrate. The
solver gives the other four faces nothing — not a wall, not a symmetry plane, not
periodicity — so `preprocess.py` draws an inert wall there, and `NY` and `NZ` are
two larger than the pore space they hold. The layers are **added**, not taken out
of the pore space, so porosity and every count are what the case declares.

> **Why shared code is copied here rather than referenced.** So that this folder
> *is* the procedure. The cost is real and worth stating: a fix to a shared tool
> has to be applied to every case that carries it, and `tests/check_repo.sh`
> fails if a copy drifts from `tools/`.
