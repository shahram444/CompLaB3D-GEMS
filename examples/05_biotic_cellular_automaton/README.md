# 05 — the first organism, and biomass that spills when a voxel fills

## 1. What this example does

This is the first case with something alive in it. One population sits in a
seeded patch, consumes a dissolved donor, grows, and releases a product:

**donor → biomass + product**

The rate is Monod in the donor, so growth saturates as the donor becomes
plentiful and falls off as it runs out. Decay is first order and always on, so
the population approaches a steady value rather than growing without bound.

**What is specific to this case is how the biomass moves.** With
`<solver_type>CA</solver_type>` it does not diffuse at all. A voxel accumulates
biomass until it reaches `<maximum_biomass_density>`, and only then does the
excess spill into its neighbours — the classical cellular-automaton biofilm rule.
That is why the measured result is what it is: over 1000 steps the population
grows **+0.23 %** and stays in exactly the **108 voxels** it was seeded in, with
**0 CA triggers**. Nothing filled, so nothing spilled.

Compare that with example 06, which is the same case with `FD` instead: the same
+0.23 % of growth, spread over **876 voxels**.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 2 of them: `donor`, `product` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Bug` — compiled kinetics (`defineKinetics.hh`) |
| **Biomass** | `Bug` — attached biofilm (seeded on its own material number), cellular automaton — biomass stays put until a voxel fills, then spills into neighbours |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 1000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The rate law

Coded in `kinetics/defineKinetics.hh`, the same file examples 06, 07 and 08 use:

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

The parameters live in the header, next to the equation that uses them, rather
than in `CompLaB.xml`. A hand-written rate law is not obliged to read the XML's
`<half_saturation_constants>` — those feed the flux-balance and surrogate paths,
which build their own uptake bounds — and this one does not.

## 4. The three biomass solvers, and why the choice matters

Examples 05, 06 and 07 are the same case with `<solver_type>` changed. The
chemistry is byte-identical — the same `defineKinetics.hh`, the same substrate
concentrations, the same geometry. Only the way biomass moves differs.

| | how biomass moves | needs |
|---|---|---|
| `CA` | stays put until a voxel reaches `<maximum_biomass_density>`, then spills the excess into neighbours | nothing extra |
| `FD` | spreads continuously by diffusion on the same grid | `<biomass_diffusion_coefficients>` |
| `LBM` | gets its own D3Q7 lattice and is transported like a solute, advecting with the flow if there is any | `<biomass_diffusion_coefficients>`, and a relaxation time the lattice can carry — see example 07 |

**Run all three and compare the totals, not the peaks.** CA and FD both grow the
same amount, because the chemistry is the same; what differs is where the biomass
ends up. CA holds it in the 108 seeded voxels; FD spreads it to 876. That makes
the *peak* density fall in the FD case while the *total* rises, which is why the
run reports the two separately.

## 5. What to check

1. **biomass rises and then levels off** as growth and decay balance;
2. **the donor is drawn down around the colony** and replenished from the left
   boundary;
3. **a product plume spreads from the colony** — it is made only where biomass
   is, so its shape is a map of where the reaction ran;
4. **no `[NEG!]` warnings.** The donor must never go negative. **[v1.3]** Earlier
   text said every increment is clamped as ΔC = max(RΔt, −C) before it is
   applied. There is no such clamp on the path this case uses. `run_kinetics` in
   `src/complab3d_processors_part1.hh`, the compiled kinetics path used by
   examples 05 to 07, computes `dC = subs_rate[iS] * dt` and applies it as
   computed, with no comparison against the local concentration; a positivity
   clamp exists only on the flux-balance and learned-law paths, neither of which
   is active here. So a time step or a rate constant large enough to consume more
   than a voxel holds *will* drive that voxel negative. The solver reports it as a
   `[NEG!]` warning, and the fix is to reduce `<ade_dt>` or the rate constant.
   Checking for `[NEG!]` is still the right check; only the reason given for it
   was wrong;
5. **the total biomass, not the peak.** The run prints both, labelled. For this
   solver they move together; for `FD` and `LBM` they do not.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 05_biotic_cellular_automaton run/mycase
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
| `kinetics/defineKinetics.hh` | This case's own chemistry: the Monod rate law above, compiled in. Shared with examples 06, 07 and 08. |

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
