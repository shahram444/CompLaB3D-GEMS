# 09 — growth predicted by a metabolic model instead of prescribed

## 1. What this example does

Examples 05 to 08 tell the organism how fast to grow: a Monod law with a
`mu_max` you chose. This case does not. It gives the organism a **metabolic
model** — a stoichiometric matrix of its internal reactions — and asks a linear
program, in every voxel and at every step, how fast that network *can* run given
what diffusion has delivered there.

```
    maximise   c' v          the biomass reaction's flux
    subject to S v  = 0      every internal metabolite at steady state
               lb <= v <= ub  uptake bounded by the LOCAL concentrations
```

Growth is the answer to that program, not an input to it. Change the substrate
supply and the growth rate changes because the optimum moves, which is the point
of flux balance analysis and the reason it is worth the cost.

The linear program is solved by **GLPK, in process** — no interpreter, no file
round-trip. Example 10 is the same case through COBRApy instead, as a
cross-check.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 3 of them: `S`, `O`, `P` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Toybug` — flux balance analysis, GLPK in process |
| **Biomass** | `Toybug` — attached biofilm (seeded on its own material number), finite difference — biomass spreads by diffusion on the same grid |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 500 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The metabolic model

`input/toy_model.xml`, deliberately small enough to check by hand:

```
    metabolites      S_c    O_c    P_c
    reactions
      0  EX_S        S_c ->                    negative flux = import
      1  EX_O        O_c ->                    negative flux = import
      2  EX_P        P_c ->                    positive flux = export
      3  BIO         S_c + 0.5 O_c -> P_c      the objective

    Steady state forces  v0 = -v3  and  v1 = -0.5 v3, so with uptake bounds
    Vs and Vo the optimum is

        growth  =  min( Vs , 2 Vo )
```

Three metabolites, four reactions, and an optimum you can work out on paper.
That is the point: when the solver reports a growth rate you can check it by
hand, which you cannot do with a genome-scale model of 2000 reactions.

The model format is the flat `<Metabolic_Model>` one that `training/extractMM.py`
writes from SBML. Real models go through that exporter first; the repository root's
[`../../models/`](../../models/) carries *E. coli* core and iJO1366 for when you
want one.

## 4. How the linear program is coupled to transport

```
    per voxel, per step:

      1.  read the local concentrations
      2.  turn each into an uptake bound          Michaelis-Menten on the
                                                  local concentration, capped
                                                  by <fba_maximum_uptake_flux>
      3.  solve  max c'v   s.t.  S v = 0,  lb <= v <= ub
      4.  growth = the objective flux; substrate rates = the exchange fluxes
      5.  scale by the local biomass and the time step, and write the increment
```

Step 2 is what couples the linear program to transport: the bounds are not
constants, they are rebuilt from what diffusion has delivered to that voxel.
Step 5 converts mmol gDW-1 h-1 — the unit published metabolic models use — into
the mol L-1 s-1 the transport solver carries, through `<biomass_molar_mass>`.

| tag | value | what it does |
|---|---|---|
| `<model_filename>` | `toy_model` | which model this organism carries |
| `<exchange_reaction_indices>` | `0 1 2` | which reaction each substrate maps to |
| `<fba_maximum_uptake_flux>` | `10. 4. 0.` | Vmax per substrate, mmol gDW⁻¹ h⁻¹ |
| `<substrate_lower_bounds>` | `-10. -4. 0.` | the LP's lower bounds; negative = import allowed |
| `<objective_direction>` | `maximize` | reaches GLPK as the sense of the program |
| `<biomass_molar_mass>` | `24.6` | g mol⁻¹, the unit conversion of step 5 |

## 5. What to check

1. **growth responds to supply.** Raise `<fba_maximum_uptake_flux>` on the donor
   and growth should rise until the acceptor becomes binding — the toy model's
   `min(Vs, 2 Vo)` says exactly where that happens;
2. **substrates are consumed in the ratio the stoichiometry gives.** The
   objective reaction takes S and O in 1 : 0.5, so the two draws must hold that
   ratio wherever the LP is unconstrained;
3. **no substrate goes negative.** The uptake bound is built from what is
   locally available, and the increment is clamped again before it is applied;
4. **the LP solves in every voxel.** A failed solve is reported, not silently
   treated as zero growth.

## 6. What this costs

A linear program per voxel per step is the most expensive rate path in this
repository by some margin. Example 11 replaces it with a network fitted to its
answers, at roughly 1/400 of the cost; example 17 with a symbolic law you can
read. Both are trained from sweeps of exactly this solver.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

**[v1.3] One prerequisite does live outside the folder: GLPK.** This case runs
flux balance analysis through GLPK, so the GLPK development package must be
installed (INSTALL.md names it `libglpk-dev`) and the build must be configured
with `-DENABLE_GLPK=ON`. `pipeline.sh` already passes that flag; a hand build
does not, because the CMake option is OFF by default. If it is missed, the
failure is not a compile error: the build succeeds, and the solver then
terminates immediately at start-up with `<enable_fba_glpk> is true but this
executable was built without GLPK`.

```bash
./scripts/setup_case.sh 09_fba_glpk run/mycase
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
| `input/toy_model.xml` | The metabolic model of section 3, in the flat `<Metabolic_Model>` format. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

### The offline step

`training/extractMM.py` turns an SBML model into the flat format the solver
reads. It is carried here rather than referenced so this folder is the whole
procedure.

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
