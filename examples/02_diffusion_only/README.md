# 02 — diffusion with no flow, and the profile you can check by hand

## 1. What this example does

This case takes example 01's tracer and switches the flow solver off entirely.
`<Peclet>0</Peclet>` means the Navier–Stokes solve never runs, no velocity field
is built, and transport is pure diffusion. The tracer is held at 1 on the inlet
face and 0 at the outlet, and reacts with nothing.

That leaves a problem with a known answer. Diffusion with no reaction and no
flow, between two fixed values, has an analytic steady state: **a straight line**.

```
    d2C
    ───  =  0        ->      C(x)  =  1 - x / L
    dx2
```

So this is the one case in the repository whose correct answer can be written
down without running it, which is what makes it worth running. Curvature in the
converged profile means something is firing that should not be — a reaction, a
boundary that is not being held, or a transport step that is not conservative.
There is nothing else happening here to hide it.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 1 of them: `tracer` |
| **Abiotic reaction** | none |
| **Biotic reaction** | none — `<biotic_mode>false</biotic_mode>` |
| **Biomass** | none; no biomass field is allocated |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 2000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. What to check

1. **the profile is straight** between the two boundary values, once the run has
   settled. This is the whole point of the case;
2. **the tracer total has stopped changing between the last rows of
   `output/summary.csv`.** The advection-diffusion loop has no convergence test:
   `<ade_converge_iT>` is read from `CompLaB.xml` and then never used, so the loop
   always runs the full `<ade_max_iT>` steps, 2000 here, and never stops early.
   Steadiness is something you confirm from the record rather than something the
   solver reports, and a total still moving at the last row means the run is too
   short;
3. **no `[NEG!]` warnings**. A conservative tracer between 0 and 1 cannot leave
   that range;
4. **no dependence on z**. The geometry is uniform in z inside its walls, so the
   profile must be too. A gradient across z means the boundary treatment differs
   between the two z faces. **[v1.3]** Earlier text asked for no dependence on y
   either; that was wrong. This case's `preprocess.py` builds a *staggered slot*,
   meaning grain blocks that alternate between the two y walls every four voxels
   along x (blocks at `x % 4 == 0` against one wall and at `x % 4 == 3` against
   the other, 288 grain voxels in total), and its own docstring says "the field is
   not trivially symmetric". Diffusion has to go around those blocks, so the
   concentration necessarily varies with y near every block. That y variation is a
   property of the staggered geometry, not a sign of a boundary problem.

## 4. What this case does not do

No flow, no reaction, no biology, no geometry change. Example 01 is the same
tracer with the flow on; example 03 adds the first reaction.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 02_diffusion_only run/mycase
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
