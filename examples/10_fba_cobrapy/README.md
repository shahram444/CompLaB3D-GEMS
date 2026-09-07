# 10 — the same linear program, solved by COBRApy instead

## 1. What this example does

This is example 09's case — the same pore space, the same three solutes, the
same four-reaction metabolic model — run through a **different linear-program
solver**. Example 09 calls GLPK from inside the C++ solver; this one hands the
stoichiometry to COBRApy in an embedded Python interpreter and reads the fluxes
back.

The two must agree. They share the geometry, the transport step and the
uptake-bound construction, and they are handed the same linear program; the only
thing that differs is who solves it. So a disagreement between 09 and 10 is not
a modelling choice, it is a bug in one of the two couplings, and that is the
strongest check available on the flux-balance path in this code.

Why keep both. GLPK is one to two orders of magnitude faster and has no
interpreter to configure, so it is what you run. COBRApy is what the metabolic
modelling community reads models with, so it is what you check against — and it
is the path that reaches everything COBRApy can do that a bare linear program
cannot.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 3 of them: `S`, `O`, `P` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Toybug` — flux balance analysis, COBRApy through an embedded interpreter |
| **Biomass** | `Toybug` — attached biofilm (seeded on its own material number), finite difference — biomass spreads by diffusion on the same grid |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 500 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The metabolic model

`input/toy_model.xml`, the same file example 09 uses and small enough to check
by hand:

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

The uptake bounds are not constants. Each is rebuilt every step from the
concentration diffusion has delivered to that voxel:

```
    V  =  Vmax * C / ( Kc + C )
```

With `<fba_maximum_uptake_flux>` 10 and 4, half-saturation constants of 0.05,
and both concentrations well above them, `Vs` approaches 10 and `2 Vo`
approaches 8 — so the electron acceptor is what limits, and the growth rate
settles near 8 mmol gDW⁻¹ h⁻¹. **That number comes from the stoichiometry, not
from this code**, which is exactly what makes it a test rather than a
demonstration.

## 4. How COBRApy is reached from C++

```
    C++ side                                Python side
    --------                                -----------
    read input/toy_model.xml once     -->   build a cobra.Model from the arrays
    per voxel, per step:
      build lb/ub from local C        -->   set Reaction.bounds
      call run_FBA()                  -->   model.optimize()
      read growth + exchange fluxes   <--   solution.fluxes
```

The model file is opened by the C++ side, not by COBRApy: the stoichiometric
matrix, the bound vectors and the objective position cross the boundary as plain
arrays, and the Python half assembles a `cobra.Model` from them. That is what
makes the two back ends comparable — they are handed the same numbers, in the
same order, from the same parser.

The contract for that argument vector is written out in full at the top of
`training/complab3d_cobrapy.py`, and duplicated verbatim in
`src/complab3d_pythonAPI.hh`. If you change one, change the other.

| tag | value | what it does |
|---|---|---|
| `<enable_fba_cobrapy>` | `true` | builds the embedded-interpreter path |
| `<reaction_type>` | `cobrapy` | routes *this organism* through it |
| `<model_filename>` | `toy_model` | which model this organism carries |
| `<exchange_reaction_indices>` | `0 1 2` | which reaction each substrate maps to |
| `<fba_maximum_uptake_flux>` | `10. 4. 0.` | Vmax per substrate, mmol gDW⁻¹ h⁻¹ |
| `<substrate_lower_bounds>` | `-10. -4. 0.` | the LP's lower bounds; negative = import allowed |
| `<biomass_molar_mass>` | `24.6` | g mol⁻¹, converting the model's gDW basis to the solver's moles |

## 5. What to check

1. **09 and 10 agree.** Run both and compare `output/summary.csv` field by
   field. They should agree to solver tolerance; anything larger is a bug in one
   coupling, not a property of either solver;
2. **the growth rate approaches 8 mmol gDW⁻¹ h⁻¹**, because `2 Vo` is the
   binding constraint — raise `<fba_maximum_uptake_flux>` on the acceptor and it
   should rise until `Vs` becomes binding instead;
3. **P accumulates at the rate biomass is made.** The objective reaction exports
   exactly one P per unit growth, so the two curves are the same curve;
4. **S and O are drawn in a 1 : 0.5 ratio** wherever the program is
   unconstrained, which is the stoichiometry of the objective reaction;
5. **no `[NEG!]` warnings.** Uptake is bounded by what is locally present and
   the increment is clamped again before it is applied.

## 6. Two things that cost people time

**The interpreter.** COBRApy is imported by name at run time from
`<src_path>`; `src/complab3d_cobrapy.py` has to be there, and the embedded
interpreter has to be the same Python you installed `cobra` into. The solver
prints the interpreter path in its first few lines and `offline.sh` prints the
shell's, so a mismatch shows up before the run rather than during it.

**Positional exchange indices.** `<exchange_reaction_indices>` is positional:
`0 1 2` is correct only for this exact model file. Insert one reaction into the
model and every index after it points somewhere else, silently.
`<exchange_reaction_names>` names the reactions instead and resolves them
against the model at start-up, so a wrong name stops the run and prints the near
matches. It needs an SBML model — the flat matrix format `toy_model.xml` uses
carries no reaction names — which is why this case still uses indices.
**Example 16 shows the named form**, on a genome-scale model.

## 7. What this case does not do

No flow (`<Peclet>0</Peclet>`), no abiotic reaction, no geometry evolution, and
no thermodynamic control on the rate. Example 11 replaces the linear program
with a network fitted to its answers; example 12 runs one organism on both a
flux balance path and a compiled kinetics path in the same voxels. ([v1.3] One
organism, not two competing populations as this line used to say: case 12
declares `<number_of_microbes>1</number_of_microbes>`, so its `summary.csv`
carries one biomass column rather than two.)


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 10_fba_cobrapy run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and reports porosity, refusing to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Prints which Python `cobra` imports from, and checks the model file parses. The embedded interpreter finding a different Python from the one you installed cobra into is the usual failure here. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's checks. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `input/toy_model.xml` | The metabolic model of section 3, in the flat `<Metabolic_Model>` format. Read by the C++ side, not by COBRApy. |
| `src/complab3d_cobrapy.py` | The Python half of the bridge, imported by name at run time. `setup_case.sh` copies it from `tools/` into `src/`. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

### The offline code this case ships

| File | What it is |
|---|---|
| `training/complab3d_cobrapy.py` | The Python half of the COBRApy back end — what the embedded interpreter actually calls, and where the argument-vector contract is written down. |
| `training/extractMM.py` | Exports an SBML or BiGG model to the flat format above, and answers the first question worth asking of any new model: can it grow at all when nothing constrains it? |

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
