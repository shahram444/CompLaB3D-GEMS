# 06 — the same organism, with biomass that spreads by diffusion

## 1. What this example does

Example 05 with one line changed: `<solver_type>FD</solver_type>`. The chemistry
is byte-identical — the same rate law, the same parameters, the same geometry —
so everything that differs between the two runs is the biomass solver and nothing
else.

Instead of holding biomass in place until a voxel fills, `FD` spreads it
continuously by diffusion on the same grid, with its own coefficient. The
measured consequence over 1000 steps:

```
                            total biomass        voxels holding biomass
    05, CA                  +0.23 %               108  ->  108
    06, FD                  +0.23 %               108  ->  876
```

**The two solvers grow the same amount and put it in different places.** That is
the whole comparison, and it is also a warning about how to read the output: the
*peak* density falls 4.05 % in this case, because the same biomass is spread over
about eight times as many voxels. A run that reported only the peak would say this
population shrank. It did not. The closing report prints the peak and the total
separately for exactly this reason.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 2 of them: `donor`, `product` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Bug` — compiled kinetics (`defineKinetics.hh`) |
| **Biomass** | `Bug` — attached biofilm (seeded on its own material number), finite difference — biomass spreads by diffusion on the same grid |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 1000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The rate law

Identical to examples 05, 07 and 08 — the same `kinetics/defineKinetics.hh`, so
any difference between those four runs is the solver and not the chemistry:

```
                              [S]
    growth  =  mu_max  ·  ───────────  ·  B                gDW L-1 s-1
                           Ks + [S]

    donor consumed   =   growth / Y            product made  =  growth · fP
    decay            =   kd · B                (first order, always on)
```

| symbol | value | unit | meaning |
|---|---|---|---|
| mu_max | 2.0e-4 | s⁻¹ | maximum specific growth rate |
| Ks | 5.0e-2 | mol L⁻¹ | half-saturation constant on the donor |
| Y | 0.4 | gDW mol⁻¹ | yield: biomass made per mole of donor |
| fP | 2.0 | mol gDW⁻¹ | product released per unit biomass made |
| kd | 1.0e-6 | s⁻¹ | first-order decay |

## 4. When to choose FD

`FD` is the middle option of three, and the choice is about what biomass
physically does in your system rather than about numerics.

- **`CA`** — an attached biofilm that holds its shape and only pushes outward
  when it is full. The classical rule, and the cheapest.
- **`FD`** — biomass that spreads down its own gradient continuously. Right when
  the population is motile, or when the aggregate is loose enough that growth
  redistributes it smoothly rather than in jumps.
- **`LBM`** — biomass carried like a solute on its own lattice, advecting with
  the flow. Right for a genuinely planktonic population; overkill otherwise, and
  it carries a numerical constraint the other two do not (see example 07).

`FD` is the only solver that **requires** `<biomass_diffusion_coefficients>`, and
it is rejected outright for a planktonic microbe: a population with no
`<microbeN>` entry under `<material_numbers>` cannot use it.

Unlike `LBM`, `FD` has no relaxation time, so it carries the case's biomass
diffusivity of 3e-13 m2 s-1 without difficulty. The same number on a
lattice-Boltzmann solver is unusable — example 07 explains why.

## 5. What to check

1. **the total biomass rises** — this is the growth. The peak falls, and that is
   the biomass spreading, not dying;
2. **it spreads from the seeded patch outward**, smoothly, with no jumps;
3. **the donor is drawn down where the biomass is**, and replenished from the
   left boundary;
4. **no `[NEG!]` warnings**;
5. **run 05, 06 and 07 and compare the totals.** 05 and 06 agree closely, because
   the chemistry is the same file and neither solver streams biomass. 07 does
   stream it, so its total is of the same size rather than the same number.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 06_biotic_finite_difference run/mycase
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
| `kinetics/defineKinetics.hh` | This case's own chemistry: the Monod rate law above. Byte-identical to examples 05, 07 and 08. |

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
