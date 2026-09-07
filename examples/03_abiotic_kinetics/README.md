# 03 — a chemical reaction with no organisms in it

## 1. What this example does

This is the first case with a reaction. Three species diffuse through the pore
space and two of them combine to make the third:

**A + B → C**

The rate is second order, one order in each reactant, coded in
`kinetics/defineAbioticKinetics.hh`:

```
    R   =   k · [A] · [B]              k = 5.0e-2  L mol-1 s-1

    dA/dt = -R        dB/dt = -R        dC/dt = +R
```

There are no organisms. `<biotic_mode>` is false, so no biomass field is
allocated at all and no rate path involving a microbe is compiled into the step.
This is the abiotic route: `<enable_abiotic_kinetics>` and one header you write
yourself, for mineral reactions, redox couples, sorption — anything that happens
in water whether or not something is alive.

**What makes it checkable.** The stoichiometry is 1:1:1, so every mole of A
consumed by the reaction appears as a mole of C, and the reaction takes A and B
away in equal measure. **[v1.3]** That constraint applies to the per step
*reaction increments*, the amounts the rate law adds to and removes from each
field in one time step, and not to the reported totals. A and B are each held at
a Dirichlet boundary, meaning a boundary that holds the concentration at a fixed
value, so both are fed from outside the domain and all three totals rise. This
case also leaves the checking to `postprocess.py`: its `<diagnostics>` block has
`<enabled>`, `<summary_csv>`, `<interval>` and `<tolerance>` but no `<conserve>`
entry, and the solver only checks a sum that a `<conserve>` tag creates, so
nothing is checked every interval here. Earlier text claimed `<diagnostics>`
checked these sums every interval; it does not.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 3 of them: `A`, `B`, `C` |
| **Abiotic reaction** | `defineAbioticKinetics.hh`, compiled in |
| **Biotic reaction** | none — `<biotic_mode>false</biotic_mode>` |
| **Biomass** | none; no biomass field is allocated |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 1500 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. Where the reaction is allowed to happen

The rate law's first line is a guard, and it is the one thing to copy when you
write your own:

```cpp
    if (C.size() < 3 || subsR.size() < 3) return;   // a short substrate list
    if (mask < 2) return;                           // solid or wall voxel
```

The first stops the law reading past the end of a concentration vector shorter
than it expects — `tests/check_kinetics_bounds.py` compiles every shipped rate
law under AddressSanitizer and calls it with substrate lists from length 0 to 5
to prove this holds. The second keeps chemistry out of walls and grains.

The concentrations are also floored at zero before use. The advection–diffusion
solver can leave a value fractionally below zero, and a rate law that squares or
multiplies such a value turns a rounding artefact into a reaction.

## 4. What to check

1. **[v1.3] all three totals rise, and C rises.** Earlier text told you to check
   that "A and B fall together" and that "C rises by the same amount A fell".
   Neither can happen in this case. `CompLaB.xml` gives A and B an
   `<initial_concentration>` of 0 and holds each at 1.0 on one face with a
   Dirichlet condition, so both are supplied from outside the domain and both
   totals *rise*. The 1:1:1 stoichiometry constrains the per step reaction
   increments, not the reported totals. What a reader can actually check is that
   C rises, since C is the only closed species here, and that the increments
   balance;
2. **the balance is reported by `postprocess.py`, not by the solver.** This case
   declares no `<conserve>` line, so it is `postprocess.py` that reports the
   change in the closed species; the tag exists and is documented in
   [`../../config/CompLaB.everything.xml`](../../config/CompLaB.everything.xml)
   if you want the solver to check a closed sum too. **[v1.3] Do not name A or B
   in a `<conserve>` sum.** A species fed from a Dirichlet boundary is not
   conserved and never will be, so such a sum would guarantee a mass-balance FAIL
   that means nothing;
3. **[v1.3] no `[NEG!]` warnings.** Earlier text said every increment is clamped
   as ΔC = max(RΔt, −C) before it is applied. There is no such clamp on the path
   this case uses. `run_abiotic_kinetics` in `src/complab3d_processors_part1.hh`
   computes `dC = subs_rate[iS] * dt` and applies it as computed, with no
   comparison against the local concentration; a positivity clamp exists only on
   the flux-balance and learned-law paths, neither of which is active here. So a
   time step or a rate constant large enough to consume more than a voxel holds
   *will* drive that voxel negative. The solver reports it as a `[NEG!]` warning,
   and the fix is to reduce `<ade_dt>` or the rate constant `k`. Checking for
   `[NEG!]` is still the right check; only the reason given for it was wrong;
4. **the reaction stops where the reactants run out**, not at a fixed time.

## 5. What this case does not do

No flow, no organisms, no geometry change. Example 13 uses the same abiotic route
to precipitate a mineral and seal the pore space; example 14 to dissolve one.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 03_abiotic_kinetics run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and reports porosity, refusing to continue if it does not percolate. Standard library only. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's checks. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `kinetics/defineAbioticKinetics.hh` | This case's own chemistry: the A + B -> C rate law above, compiled in. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

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
