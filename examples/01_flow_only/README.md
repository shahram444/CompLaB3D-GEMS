# 01 — flow through a pore space, and a tracer carried by it

## 1. What this example does

This case solves Navier–Stokes flow through a 24 × 26 × 8 pore space — 3168 open
voxels, porosity 0.6346 — and then sweeps a conservative tracer through it. There
is no chemistry and no biology. It is the smoke test: if this fails, the build or
the geometry is wrong and nothing further is worth debugging.

What makes it more than a smoke test is **how the flow is driven**. You do not
supply a pressure drop. You supply a target Péclet number, and the solver works
backwards to the pressure drop that achieves it *on your geometry*:

```
    1.  solve the flow at a seed pressure drop        deltaP = 1e-6
    2.  measure the permeability of this pore space   k      = 1.083443e-01
    3.  compute the deltaP that hits the target Pe    deltaP = 6.922374e-02
    4.  solve again                                   Pe     = 0.99973
```

The seed value is discarded. Permeability depends on the pore space you happened
to build, so a pressure drop that gives Pe = 1 in one geometry gives something
else in another; asking for the Péclet number directly is what makes two
different geometries comparable.

The flow converges at iteration 904 against a tolerance of 1e-6, then the tracer
is released. Held at 1 on the inlet face and 0 at the outlet, its total rises
from 208 to 488 over 400 steps as it fills the pore space, and the profile stays
monotone in x — concentration falls with distance, never rises.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet of 1 |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 1 of them: `tracer` |
| **Abiotic reaction** | none |
| **Biotic reaction** | none — `<biotic_mode>false</biotic_mode>` |
| **Biomass** | none; no biomass field is allocated |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 400 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. How each piece is solved

**Flow** — D3Q19 lattice Boltzmann, BGK collision, run to a convergence tolerance
rather than a fixed step count. `<tau>` sets the viscosity:

```
    nu   =  ( tau - 1/2 ) / 3        tau = 0.8  ->  nu = 0.1, omega = 1.25
```

**The tracer** — its own D3Q7 advection–diffusion lattice, advected by the
converged velocity field and diffusing with its own coefficient. It reacts with
nothing, so its only sources are the two boundaries.

| setting | value | what it controls |
|---|---|---|
| `<Peclet>` | 1.0 | the target; 0 skips the flow solver entirely |
| `<tau>` | 0.8 | viscosity, through the relation above |
| `<delta_P>` | 1e-6 | a seed only, replaced at step 3 above |
| `<biotic_mode>` | false | no biomass fields are allocated at all |

## 4. Expected output

```
  [NS] tau=8.000000e-01, omega=1.250000e+00, nu=1.000000e-01
  [NS] Converged at iter=904
  [NS] Permeability k=1.083443e-01 (lattice)
  [NS] Re-running with corrected deltaP=6.922374e-02
  [NS] Pe achieved=9.997300e-01 (target=1.000000e+00)
```

Four things to check, in this order:

1. **the flow converges** — an iteration count, not the cap. Hitting
   `<ns_max_iT1>` instead means the tolerance was never met;
2. **the achieved Péclet is within a fraction of a percent of the target** —
   0.99973 against 1.0 here;
3. **no `[NEG!]` warnings** — a conservative tracer cannot go negative, so one
   would mean the transport step is wrong;
4. **the tracer profile is monotone in x** — it falls with distance from the
   inlet and never rises.

## 5. What this case does not do

No reaction of any kind, so nothing here exercises a rate path, the mass-balance
checks or the geometry evolution. Example 02 removes the flow and keeps the
transport; example 03 adds the first reaction.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 01_flow_only run/mycase
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
