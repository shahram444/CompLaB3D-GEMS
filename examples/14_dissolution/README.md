# 14 — acid eating the grains away, and the pore space reopening

## 1. What this example does

Example 13 seals pore space. This one opens it. Acid is delivered from the left
face, eats into the mineral coating on the grain surfaces, and the voxels it
empties are converted back to pore, after which the flow is re-solved through
the wider channel.

```
    CaCO3(s)  +  H(+)   ->   Ca(2+)  +  HCO3(-)
```

(The bicarbonate is not carried as a substrate in this case: the three fields are
`H`, `Ca` and `calcite`, and the rate law in `kinetics/defineAbioticKinetics.hh`
writes only `Ca` and `H`. Do not go looking for it in the output.)

Calcite dissolution in an acidified pore space is the textbook case for this:
it is fast, it is strongly transport-limited, and the feedback runs the other
way from precipitation. A dissolving channel widens, which lets more acid
through, which widens it faster. That instability is why dissolution fronts
finger rather than advance flat, and why a pore-scale model is worth running for
it at all.

The mineral does not have to have precipitated first. `<initial_fill>` starts
the declared grain surfaces completely solid, at full mineral density, and they
are eaten away from there.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346), plus 72 voxels of a declared solid phase |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet of 1 |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 3 of them: `H`, `Ca`, `calcite`. `calcite` is immobile: never streams, only accumulates its reaction term |
| **Abiotic reaction** | `defineAbioticKinetics.hh`, compiled in |
| **Biotic reaction** | none — `<biotic_mode>false</biotic_mode>` |
| **Biomass** | none; no biomass field is allocated |
| **Geometry evolution** | dissolution: a mineral is consumed and its voxels reopen |
| **Run length** | 2000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

The geometry is the same staggered slot example 13 uses. The 72 voxels of
declared phase are the one-voxel layer on each grain's **pore-facing** side;
the rest of each block stays inert wall. Mineral buried inside a grain can never
be reached by water, so declaring it would cost solve time and change nothing —
and would make the reactive surface area of the case a fiction.

## 3. What is allowed to dissolve, and why that is declared

Only material inside a `<phaseN>` block can dissolve. That restriction is
deliberate and it is not a convenience.

Consider the alternative rule — *a solid voxel with no mineral left becomes
pore*. An inert quartz grain has never held any mineral. It is already below any
threshold you could write. On the first step the entire grain pack would
dissolve. Declaring phases makes inert material safe **by construction** rather
than by a guard that someone has to remember to write.

```xml
    <phase0>
        <name>calcite_grain</name>
        <material_number>0</material_number>   which voxels in geometry.dat
        <substrate>2</substrate>               which substrate holds its inventory
        <full_density>27.1</full_density>      mol/L when completely solid
        <initial_fill>27.1</initial_fill>      starts solid
    </phase0>
```

`<full_density>` comes from the mineral, the same way `<max_precipRho>` does in
example 13:

```
    full_density  =  1000 / Vmolar          mol L-1, Vmolar in cm3 mol-1

    calcite       100.09 g mol-1
                    2.71 g cm-3   ->   36.9 cm3 mol-1   ->   27.1 mol L-1
```

## 4. The rate law, and the water it is given

`defineDissolutionRate()` is a **separate hook** from the precipitation one, in
the same header. First order in acid and first order in the mineral remaining:

```
    R  =  k_calcite * [H+] * m                mol L-1 s-1

    k_calcite = 1.0e-2   L mol-1 s-1 per mol L-1 of remaining mineral
    m                    the mineral inventory still in that voxel
```

```
    dm/dt = -R           d[Ca2+]/dt = +R           d[H+]/dt = -R
```

The `m` factor is what makes the rate fall as a voxel empties, instead of
running at full speed until the inventory hits zero and then stopping in one
step.

**Which water the rate sees.** A solid voxel has no water of its own that moves.
So the rate law is handed the concentration in the water **touching** the
mineral — averaged over that voxel's open neighbours — not a concentration
stored in the solid voxel itself. A declared voxel with no open neighbour is
unreachable and does not dissolve at all, which is the same physical statement
as the geometry note above.

| tag | value | what it controls |
|---|---|---|
| `<reopen_fraction>` | `0.9` | a voxel reopens once it falls below 0.9 × full density |
| `<surface_only>` | `1` | only wetted voxels dissolve |
| `<update_interval>` | `100` | steps between geometry re-checks and flow re-solves |

**`reopen_fraction` is a hysteresis knob.** Reopening at exactly full density
would flicker a voxel between solid and pore on the arithmetic noise of a single
step, and each flip costs a flow re-solve. 0.9 means a voxel has to have
genuinely lost a tenth of its inventory before the geometry changes.

## 5. Mass balance, and why this case is the one that catches boundary bugs

Calcium is created only by dissolution and destroyed by nothing. It enters
through no boundary. That makes it the tightest conservation statement available
in this repository:

```
    total Ca in the domain   ==   mineral removed from the declared phase
```

with no source and no sink to explain a discrepancy. Any gap is transport losing
or manufacturing mass, and this case is sensitive enough to see it. Both of
calcium's faces are `closed` — genuine bounce-back, exactly no-flux — rather
than zero-gradient, because a zero-gradient face is an **open outflow** that
copies whatever is inside it outward, and using one as a wall quietly
manufactures mass in exactly this kind of balance. This case declares no
`<conserve>` line: `postprocess.py` runs the comparison on the rows of
`output/summary.csv`, and it accounts separately for the mass held in bounce-back
voxels, which `computeDensity` cannot see. The tag exists and is documented in
[`../../config/CompLaB.everything.xml`](../../config/CompLaB.everything.xml) if
you want the solver to run the check as well.

## 6. What to check

1. **porosity rises monotonically.** Only dissolution is enabled, so nothing can
   close a voxel; any fall is a bug;
2. **the calcium balance closes** to solver tolerance — see section 5. This is
   the check to trust before any of the others;
3. **the front is not flat.** Dissolution should reach further along the
   high-flow paths, because those deliver more acid. A flat front means the flow
   field is not being re-solved, or transport is diffusion-dominated and the case
   is not showing what it is for;
4. **acid is consumed one for one with calcium released.** The stoichiometry is
   1 : 1;
5. **the mineral inventory only falls**, and no voxel goes below zero. The rate's
   `m` factor should make that impossible before the clamp is needed.

## 7. What this case does not do

Nothing precipitates: a voxel that reopens stays open. Example 15 enables both
paths at once, which is where the two have to agree about who owns a voxel and
what happens when a voxel is being filled and emptied in the same interval.

The acid is a bare `H` field with no carbonate system behind it — no
bicarbonate, no CO2, no pH buffering — so the front advances faster than a real
carbonate system would allow. Example 04 adds the aqueous speciation solver that
would supply the buffering.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 14_dissolution run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the staggered slot and coats the pore-facing grain surfaces with the declared phase, refusing to continue if the slot does not percolate. Standard library only. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and checks the calcium balance of section 5. Standard library only. |
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
first-order law of section 4 is what gets compiled in. Its
`defineAbioticRxnKinetics()` is present but returns zero — this case runs no
homogeneous reaction, and an unused entry point still has to link.

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
