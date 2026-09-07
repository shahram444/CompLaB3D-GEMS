# 08 — two populations competing for one substrate

## 1. What this example does

Two populations share one pore space and one donor, and they are not the same
organism. `Bug1_fast` grows quickly but needs a lot of donor to do it;
`Bug2_efficient` grows half as fast but keeps going at concentrations where the
first has stalled. They are seeded on separate material numbers, so each starts
in its own patch, and from there they compete.

**This is also the case that proves the per-microbe settings are independent.**
The two organisms differ in three ways at once:

```
                     Bug1_fast          Bug2_efficient
    mu_max            2.0e-4 s-1         1.0e-4 s-1        twice as fast
    Ks                5.0e-2 mol/L       5.0e-3 mol/L      ten times hungrier
    solver_type       CA                 LBM               different transport
```

Only the last of those three is a tag in `CompLaB.xml`: `<solver_type>` is
declared once per microbe. There is no `<mu_max>` tag anywhere in the repository,
and both `mu_max` and the `Ks` this rate law uses come from the per-microbe table
in `kinetics/defineKinetics.hh`. `<half_saturation_constants>` is a per-microbe
tag and carries the same two numbers, but it feeds the flux-balance and surrogate
paths, which build their own uptake bounds, and a compiled rate law does not read
it.

Nothing in the solver assumes populations are alike. `<reaction_type>` and
`<solver_type>` are per microbe, and a run may mix a cellular automaton with a
lattice-Boltzmann field, or flux balance analysis with a hand-written rate law
(example 12), in the same step.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet of 0.6 |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 2 of them: `donor`, `product` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Bug1_fast` — compiled kinetics (`defineKinetics.hh`); `Bug2_efficient` — compiled kinetics (`defineKinetics.hh`) |
| **Biomass** | `Bug1_fast` — attached biofilm (seeded on its own material number), cellular automaton — biomass stays put until a voxel fills, then spills into neighbours; `Bug2_efficient` — attached biofilm (seeded on its own material number), lattice Boltzmann — biomass gets its own D3Q7 lattice, like a solute |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 1000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The rate law, and how it handles two organisms

The same `defineKinetics.hh` as examples 05 to 07. It is written for **any**
number of populations: it loops over `B.size()` and takes each organism's
parameters from a table indexed in `CompLaB.xml` order.

```
                                 [S]
    growth[m]  =  mu_max[m] · ───────────── · B[m]           gDW L-1 s-1
                               Ks[m] + [S]

    donor consumed   =   SUM over m of  growth[m] / Y[m]
    product made     =   SUM over m of  growth[m] · fP[m]
```

| symbol | `Bug1_fast` | `Bug2_efficient` | unit |
|---|---|---|---|
| mu_max | 2.0e-4 | 1.0e-4 | s⁻¹ |
| Ks | 5.0e-2 | 5.0e-3 | mol L⁻¹ |
| Y | 0.4 | 0.4 | gDW mol⁻¹ |
| fP | 2.0 | 2.0 | mol gDW⁻¹ |
| kd | 1.0e-6 | 1.0e-6 | s⁻¹ |

Both draw from the **same** `donor` field, and their draws are summed into one
increment before it is applied. Neither can take the donor below zero: the
increment is clamped as ΔC = max(RΔt, −C) after both contributions are in, so
competition is resolved on the shared budget rather than by whichever organism
happened to be evaluated first.

## 4. A note on the Péclet number

This case runs at `<Peclet>0.6</Peclet>`, not 1. At Pe = 1 this geometry gives a
**grid** Péclet of 2.6, above the value of 2 at which advection overshoots, and
the run produced product concentrations of −4.3e-4 against a peak of 4.5e-3 —
flagged `[NEG!]` at every interval. The solver warns about both the grid Péclet
and the Mach number at start-up; this case is that warning acted on. At 0.6 the
grid Péclet is 1.6 and the Mach number 0.20, both inside their limits.

If you raise the Péclet here, read the `[STABILITY]` block before believing the
output.

## 5. What to check

1. **which population wins, and where.** The fast one should dominate near the
   inlet where donor is plentiful; the efficient one should hold on further in;
2. **the donor is drawn down by both** and never goes negative;
3. **the two totals separately.** `output/summary.csv` carries a column per
   population, so the competition can be read as two curves rather than one sum;
4. **no `[NEG!]` and no `[STABILITY]` warning.**


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 08_two_microbes run/mycase
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
| `kinetics/defineKinetics.hh` | This case's own chemistry: the two-organism Monod law above. Byte-identical to examples 05, 06 and 07. |

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
