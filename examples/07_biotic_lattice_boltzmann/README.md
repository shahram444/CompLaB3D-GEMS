# 07 — biomass transported on its own lattice, and the relaxation time that constrains it

## 1. What this example does

Example 05 with `<solver_type>LBM</solver_type>`. Biomass is given its own D3Q7
advection–diffusion lattice and transported exactly like a dissolved species —
the same machinery, the same boundary conditions, and it will advect with the
flow if there is any.

Over 1000 steps the population grows **126.00 → 126.23**, +0.185 %, in line with
examples 05 and 06 whose chemistry is the same file.

**But the console total looks like nothing happened**, and understanding why is
the point of this case. The visible total *falls* from 108.0 to 101.9 while the
biomass is in fact growing. Two effects stack:

```
    what computeDensity() reports          108.00  ->  101.92
    what is in transit at bounce-back walls  18.00  ->   24.31
    the conserved total                     126.00  ->  126.23
```

Palabos's `BounceBack::computeDensity` answers from a stored number and ignores
the populations the cell is holding. As an LBM biomass field spreads, more of it
is in flight at a wall at any instant — here 6 units of it — so the reported
total falls while the real one rises. The scalar record's `_held` column is that
difference, measured from the populations, and **`_total + _held` is the
conserved quantity**.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 2 of them: `donor`, `product` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Bug` — compiled kinetics (`defineKinetics.hh`) |
| **Biomass** | `Bug` — attached biofilm (seeded on its own material number), lattice Boltzmann — biomass gets its own D3Q7 lattice, like a solute |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 1000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The relaxation time, which is this solver's real constraint

`LBM` carries the biomass diffusivity as a lattice relaxation time:

```
                   D
    tau   =   ───────────  ·  ( tau_ref - 1/2 )   +   1/2
               D_ref
```

Example 06 uses a biomass diffusivity of 3e-13 m2 s-1 against a reference
solute of 5e-10. That ratio puts tau at **0.50018**, and BGK at a relaxation time
that close to 0.5 does not diffuse slowly — it **rings**. Measured on this case
with every reaction switched off and a fully closed box, the worst negative
biomass as a fraction of the peak:

| tau | worst negative | negative voxels at step 1000 |
|---|---|---|
| 0.50018 | 25 % | 842 |
| 0.510 | 1.6 % | 520 |
| 0.520 | 0.33 % | 102 |
| 0.550 | 0.02 % | 0 |
| 0.800 | 0 | 0 |

**So this case cannot use the same diffusivity as 06.** It uses 1e-10,
which puts tau at 0.56. The solver now audits every relaxation time that will
actually relax — skipping immobile species and `CA`/`FD` biomass, whose tau is
never used — and refuses below 0.51, so this cannot be walked into by accident.

That constraint is the reason slow biomass usually belongs on `CA` or `FD`:
neither has a relaxation time, `CA` needs no diffusivity at all, and `FD` carries
3e-13 without difficulty.

## 4. What LBM is for, and what it costs

It costs a full extra lattice per population. It is the right choice when biomass
genuinely moves with the flow — a planktonic population in a flowing pore space —
and overkill when it does not.

## 5. What to check

1. **the conserved total, not the visible one.** `Bug_total + Bug_held` in
   `output/summary.csv`. With no reaction at all this case conserves biomass to
   1.3e-14 over 1000 steps;
2. **no negative biomass anywhere.** `Bug_min` in the record. A negative value
   means the relaxation time is too close to 0.5, not that the chemistry is
   wrong;
3. **no `[NEG!]` warnings on the donor**;
4. **compare the total against 05 and 06.** Same chemistry, so the growth should
   be of the same size. It is not the same number: the LBM total differs from the
   CA and FD totals in the first significant figure, because the biomass lattice
   streams, so this is a sanity check on the same order of magnitude rather than
   a digit-for-digit match.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 07_biotic_lattice_boltzmann run/mycase
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
| `kinetics/defineKinetics.hh` | This case's own chemistry: the Monod rate law of example 05. Byte-identical to 05, 06 and 08. |

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
