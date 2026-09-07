# 04 — aqueous speciation, solved in every voxel at every step

## 1. What this example does

The three cases before this one move species around and react them on a
timescale. This one adds the reactions that are **already finished**: acid–base
speciation, which equilibrates far faster than anything transport can do, and is
therefore not a rate at all but a constraint the composition must satisfy
everywhere, at all times.

The system is the carbonate one, in a slab of water with acid entering from the
left:

```
    HCO3-              <->   CO3(2-)  +  H+          log K = -10.33
    HCO3-  +  H+       <->   H2CO3                   log K =  +6.35
```

Two components — total carbonate and protons — and four species. The solver
carries the two components as transported fields and reconstructs all four
species in every open voxel, every step, from the tableau in `CompLaB.xml`.

**The measured result.** Acid enters and drives the carbonate system: over 600
steps `HCO3` falls by 1.175 while `H2CO3` accumulates 0.0400 and `CO3` 0.0150,
each formed in place from components that were fed in. The speciation solve
converged in **1 900 801 of 1 900 801 voxel-steps (100.0 %)**, and the initial
composition was feasible in 7 iterations to a residual of 8.9e-16.

That convergence count is the number to watch. A speciation solver that fails
in even a small fraction of voxels leaves those voxels holding a composition
that satisfies nothing, and the error propagates through every later step.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 4 of them: `HCO3`, `H`, `CO3`, `H2CO3` |
| **Abiotic reaction** | none |
| **Biotic reaction** | none — `<biotic_mode>false</biotic_mode>` |
| **Biomass** | none; no biomass field is allocated |
| **Geometry evolution** | none — the pore space is fixed |
| **Aqueous speciation** | on — continued-fraction solver with Anderson acceleration |
| **Run length** | 600 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. How the speciation is solved

Not by Newton–Raphson on the full system. The solver uses **continued-fraction
iteration with Anderson acceleration**, which is chosen for a specific reason:
Newton needs a Jacobian and a good starting guess, and in a pore-scale run it is
being restarted in a million voxels a step with whatever the transport step left
behind. The continued-fraction form is derived directly from the mass-action
expressions, cannot produce a negative concentration, and does not need a
Jacobian; Anderson acceleration recovers most of the convergence rate Newton
would have given.

```
    tolerance   1e-10
    max iters   200
    reported    converged / did not converge, summed over the whole run and
                printed once, after the last step
```

## 4. The tableau

Written in `CompLaB.xml` as a stoichiometric matrix with a log K per row:

```
    components          HCO3    H       log K
    HCO3  (component)     1     0        0
    H     (component)     0     1        0
    CO3                   1    -1      -10.33
    H2CO3                 1     1       +6.35
```

A component's own row is the identity with log K 0. Everything else is written
as a combination of components, which is what lets the solver carry two
transported fields instead of four and reconstruct the rest.

## 5. What to check

1. **the convergence line reads 100 %.** Anything less means some voxels hold a
   composition that satisfies no equilibrium;
2. **the initial composition is feasible** — reported separately at start-up,
   before the first step. An infeasible start is a tableau error, not a
   numerical one;
3. **the complexes are formed in place.** `CO3` and `H2CO3` have closed
   boundaries: they are neither fed nor drained, so their totals change only by
   speciation. `postprocess.py` reports them separately for that reason;
4. **`CO3` need not be monotone.** It rises, then partly reverses as the acid
   front arrives and pushes carbonate toward `H2CO3`. That is the chemistry, not
   an instability.

## 6. What this case does not do

Speciation here is **ideal and isothermal**. There is no Debye–Hückel or Davies
activity correction, no temperature dependence, no gas phase, no redox couple
and no mineral saturation index. The log K values you supply are conditional
constants at your own ionic strength and 25 °C, and it is your responsibility
that they are.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 04_equilibrium run/mycase
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
