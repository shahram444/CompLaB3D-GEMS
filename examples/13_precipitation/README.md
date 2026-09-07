# 13 — a mineral that forms, fills its voxels, and blocks the flow

## 1. What this example does

Everything before this case runs on a pore space that never changes. This one
does not: a mineral forms where two solutes meet, accumulates in the voxels
where it formed, and once a voxel is full it becomes **solid**. The flow is then
re-solved around the new geometry, so the precipitate changes the transport that
produced it.

The reaction is iron monosulphide precipitation, the mineral that forms wherever
dissolved iron meets sulphide in anoxic sediment:

```
    Fe(2+)  +  HS(-)   ->   FeS(s)  +  H(+)
```

(The proton is not carried as a substrate in this case: the three fields are
`Fe2`, `HS` and `FeS`, and the rate law in `kinetics/defineAbioticKinetics.hh`
writes no `H`. Do not go looking for it in the output.)

Iron enters from the left face, sulphide from the right. Neither is present at
the start. They diffuse and advect toward each other and the mineral forms in a
band where they overlap — not at either boundary, and not uniformly, but where
transport puts them together. That band is the result: it is set by the two
supplies and by the pore geometry, and nothing in the configuration says where
it should be.

The feedback is the point of the case. Mineral fills a voxel, the voxel seals,
the flow finds another path, the reactants meet somewhere else. Porosity falls
as the run proceeds and `output/summary.csv` records it.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346), plus 72 voxels of a declared solid phase |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet of 1 |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 3 of them: `Fe2`, `HS`, `FeS`. `FeS` is immobile: never streams, only accumulates its reaction term |
| **Abiotic reaction** | `defineAbioticKinetics.hh`, compiled in |
| **Biotic reaction** | none — `<biotic_mode>false</biotic_mode>` |
| **Biomass** | none; no biomass field is allocated |
| **Geometry evolution** | precipitation: a mineral fills pore voxels and seals them |
| **Run length** | 2000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

The geometry is a **staggered slot**: grain blocks alternate between the two
walls every four voxels along x, so the flow has to weave and the field is not
trivially symmetric. The one-voxel layer on each block's inner face — the face
that touches water — is declared as a reactive solid phase rather than inert
wall. Mineral buried inside a grain can never be reached, so declaring it would
cost solve time and change nothing.

## 3. The chemistry

Second order, one term, in `kinetics/defineAbioticKinetics.hh`:

```
    R  =  k_FeS * [Fe2+] * [HS-]                  mol L-1 s-1

    k_FeS = 1.0    L mol-1 s-1
```

and the stoichiometry that rate is applied with:

```
    d[Fe2+]/dt = -R          d[HS-]/dt = -R          d[FeS]/dt = +R
```

This is a **kinetic** precipitation law, not an equilibrium one: the rate is
proportional to the product of the two concentrations at all times, with no
solubility product and no supersaturation threshold. That is the simplest law
that produces the behaviour this case exists to show. A real mineral would carry
a saturation-state term, `R = k A (Omega - 1)^n`, with `Omega` the saturation
index and `A` the reactive surface area; it goes in the same function, in the
same place.

**FeS is declared immobile.** `<immobile>true</immobile>` and a diffusion
coefficient of zero mean the mineral's lattice never streams: it accumulates
exactly in the voxel where it formed and stays there. A solid that diffused
would smear the band away, and the mass-balance diagnostics account for immobile
fields separately for the same reason.

## 4. From concentration to geometry

A voxel is converted to solid when the mineral it holds reaches the density of
the mineral itself. That number comes from the mineral, not from the model:

```
    max_precipRho  =  1000 / Vmolar          mol L-1, Vmolar in cm3 mol-1

    mackinawite    87.91 g mol-1
                    4.30 g cm-3     ->   20.44 cm3 mol-1   ->   48.9 mol L-1
```

`48.9` is what `<max_precipRho>` is set to, and if you change minerals you
recompute it the same way.

| tag | value | what it controls |
|---|---|---|
| `<solid_substrate>` | `2` | which substrate index is the mineral |
| `<max_precipRho>` | `48.9` | mol L⁻¹ in a completely full voxel, from the arithmetic above |
| `<surface_only>` | `1` | heterogeneous nucleation: mineral grows only on wall-adjacent voxels |
| `<perm_ratio>` | `0` | a filled voxel is an impermeable wall, not a low-permeability one |
| `<update_interval>` | `100` | steps between geometry re-checks and flow re-solves |

**`surface_only` is a physical choice, not a numerical one.** Set to 1, the
mineral can only grow in a voxel that touches an existing solid — heterogeneous
nucleation on a grain surface, which is how minerals overwhelmingly precipitate
in porous media. Set to 0 it nucleates anywhere in the water, which seals the
pore throat far more readily and far earlier than a real system would.

`<update_interval>` is a cost trade. The geometry check is cheap; the flow
re-solve that follows a change is not. Every 100 steps means the flow field lags
the geometry by at most 100 steps, which is acceptable while sealing is slow and
is not once it accelerates.

## 5. What to check

1. **porosity falls monotonically.** It is the headline result; read it from the
   `porosity` column of `output/summary.csv` rather than by eye in the log. This
   case only seals — it cannot reopen anything — so any rise is a bug;
2. **the mineral forms in a band, not at a boundary.** Iron enters left and
   sulphide right; a mineral layer plated onto either inlet face means one
   supply is overwhelming the other and the case is not showing what it
   intends to;
3. **Fe and HS are consumed one for one, and FeS gains what they lose.** The
   stoichiometry is 1 : 1 : 1, so the three totals are one balance;
4. **the flow re-solves after a change.** Look for the flow solver reporting
   again mid-run; if it never does, no voxel ever reached `max_precipRho`;
5. **no `[NEG!]` warnings.** A second-order rate can overshoot in a step where
   one reactant is nearly exhausted, which is what the clamp exists for.

## 6. What this case does not do

**The conversion is one way.** Nothing here reopens a sealed voxel: the mineral
can only grow. Example 14 runs the reverse — a declared mineral phase dissolving
and its voxels returning to pore — and example 15 runs both at once, which is
where the two paths have to agree about who owns a voxel.

There is also no aqueous speciation: `Fe2` and `HS` are the total dissolved
concentrations, not the free ions an activity model would give. Example 04 adds
the speciation solver.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 13_precipitation run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the staggered slot, lines the pore-facing grain surfaces with the reactive phase, and refuses to continue if the slot does not percolate. Standard library only. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and checks that porosity only falls. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

### The rate law this case ships

`kinetics/defineAbioticKinetics.hh` **replaces** the shared header for this
case: `setup_case.sh` copies it to `defineAbioticKinetics.hh` at the case root,
over the default it laid down there, so the
second-order law of section 3 is what gets compiled in. Its
`defineDissolutionRate()` is present but deliberately returns zero — this case
does not dissolve anything, and an unused entry point still has to link.

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
