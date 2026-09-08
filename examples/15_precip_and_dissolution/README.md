# 15 — a mineral that seals a voxel and can be dissolved back out of it

## 1. What this example does

Examples 13 and 14 each run the pore space one way: 13 only seals voxels, 14
only reopens them. This one runs **both directions on the same mineral**, so a
voxel can precipitate shut, sit sealed while the flow reroutes around it, and
later be dissolved back open by acid arriving down a different path.

```
    precipitation    Fe(2+)  +  HS(-)      ->  FeS(s)
    dissolution      FeS(s)  +  H(+)       ->  Fe(2+)  +  HS(-)
```

**[v1.3]** The precipitation line produces no proton. The shipped
`kinetics/defineAbioticKinetics.hh` writes only the iron, the sulphide and the
mineral in its precipitation hook, and never touches the proton entry; the
reaction in `CompLaB.xml` is written the same way. This is deliberate, and it is
why the two rate laws are not inverses of one another and why acid only ever
falls, as sections 3 and 5 describe.

Nothing is mineral at the start. Iron comes in from the left, sulphide from the
right, acid from the left; FeS forms where iron and sulphide overlap, and is
attacked wherever acid reaches it. Porosity is no longer a one-way trend, which
is the whole reason to run this case rather than the two before it.

**The two rate laws are not inverses of one another**, and this is the design
point worth understanding before you read the output. Precipitation here is
second order in the dissolved ions; dissolution is first order in acid and in
the mineral remaining. They are two separate laws that happen to move the same
mineral. That is how this code is meant to be used: there is no solubility
product anywhere in it and therefore no saturation state, so there is nothing to
drive a single reversible rate. If you need detailed balance — the forward and
back rates equal at equilibrium — you have to impose it yourself, in the two
laws you write.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346), plus 72 voxels of a declared solid phase |
| **Flow** | D3Q19 lattice Boltzmann, BGK collision, driven to a target Peclet of 1 |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 4 of them: `Fe2`, `HS`, `H`, `FeS`. `FeS` is immobile: never streams, only accumulates its reaction term |
| **Abiotic reaction** | `defineAbioticKinetics.hh`, compiled in |
| **Biotic reaction** | none — `<biotic_mode>false</biotic_mode>` |
| **Biomass** | none; no biomass field is allocated |
| **Geometry evolution** | precipitation: a mineral fills pore voxels and seals them; dissolution: a mineral is consumed and its voxels reopen |
| **Run length** | 3000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The two rate laws

Both live in `kinetics/defineAbioticKinetics.hh`, in the two separate hooks the
solver calls.

**Precipitation** — `defineAbioticRxnKinetics()`, a homogeneous reaction between
two dissolved species:

```
    Rp  =  kp * [Fe2+] * [HS-]                mol L-1 s-1        kp = 1.0 L mol-1 s-1

    d[Fe2+]/dt = -Rp      d[HS-]/dt = -Rp      d[FeS]/dt = +Rp
```

**Dissolution** — `defineDissolutionRate()`, a surface reaction, handed the
concentration in the water touching the mineral because a solid voxel has no
water of its own:

```
    Rd  =  kd * [H+] * m                      mol L-1 s-1        kd = 5.0e-2 L mol-1 s-1

    dm/dt = -Rd     d[Fe2+]/dt = +Rd     d[HS-]/dt = +Rd     d[H+]/dt = -Rd
```

`m` is the mineral inventory left in that voxel, so the rate falls as the voxel
empties rather than running flat out and stopping in one step.

Because a proton is consumed by dissolution and none is supplied by
precipitation in this header, acid is a genuine reactant here, delivered from the
left face and consumed as it works. Where it runs out, dissolution stops and
whatever has sealed stays sealed.

## 4. How the two paths share a voxel

The two configuration blocks are independent, which is why both can be on at
once. What connects them is `<phase0>`:

```xml
    <phase0>
        <name>FeS_precipitate</name>
        <material_number>0</material_number>
        <substrate>3</substrate>
        <full_density>48.9</full_density>
        <initial_fill>0</initial_fill>          nothing there at the start
        <is_precipitate>true</is_precipitate>   voxels sealed by precipitation
    </phase0>
```

`<is_precipitate>true</is_precipitate>` is the link. A voxel sealed by the
precipitation path is **stamped** as belonging to this phase, which is what makes
it eligible for the dissolution path later. Without that stamp the two halves
would not know they were operating on the same mineral, and a sealed voxel would
be permanent, as it is in example 13.

`48.9 mol L⁻¹` is mackinawite's density, and it appears in both blocks for the
same reason: `<max_precipRho>` is the density a voxel must reach to seal, and
`<full_density>` is the density it is measured against when reopening. They are
the same physical number and must be written the same in both places.

```
    max_precipRho = full_density = 1000 / Vmolar

    mackinawite   87.91 g mol-1 / 4.30 g cm-3  =  20.44 cm3 mol-1  ->  48.9 mol L-1
```

**`<reopen_fraction>` earns its keep here and nowhere else.** With both
directions active a voxel can sit right at the sealing threshold. On a single
threshold it would flip open and shut every update interval, and every flip
re-solves the flow field. `0.9` puts a gap between sealing at full density and
reopening below 0.9 × full, so a marginal voxel settles instead of chattering.
Watch porosity for this specifically: it should move in steps and then hold, not
oscillate at the update interval.

| tag | value | what it controls |
|---|---|---|
| `<solid_substrate>` | `3` | which substrate index holds the mineral |
| `<max_precipRho>` / `<full_density>` | `48.9` | the same mineral density, in both blocks |
| `<surface_only>` | `1` (both) | growth and attack only on wetted, wall-adjacent voxels |
| `<perm_ratio>` | `0` | a sealed voxel is an impermeable wall |
| `<reopen_fraction>` | `0.9` | the hysteresis gap described above |
| `<update_interval>` | `100` (both) | steps between geometry re-checks and flow re-solves |

## 5. What to check

1. **porosity moves both ways, and settles.** It should fall as the mineral
   forms and rise where acid reaches it, in steps at the update interval, and
   then hold. Oscillation with a period of exactly `<update_interval>` is the
   chattering `reopen_fraction` exists to prevent — that is the failure this case
   is built to expose;
2. **iron and sulphide balance against the mineral.** Both are 1 : 1 with FeS in
   both directions, so the sum of what is dissolved and what is mineral changes
   only through the boundaries;
3. **acid only ever falls.** It is consumed by dissolution and produced by
   nothing in this header;
4. **mineral appears before any of it dissolves.** `<initial_fill>` is 0, so at
   step 0 there is nothing to attack — a dissolution rate reported before any
   precipitation has occurred means a phase is being stamped that should not be;
5. **the flow re-solves more than once.** Both a seal and a reopening trigger
   one;
6. **no `[NEG!]` warnings.** Two rate laws competing for the same mineral in the
   same interval is the situation most likely to overdraw it.

## 6. What this case does not do

No biology and no aqueous speciation. `Fe2`, `HS` and `H` are total dissolved
concentrations, not free ions: there is no FeHS⁺ complex, no sulphide
protonation equilibrium, and no ionic-strength correction, all of which control
where real mackinawite forms. Example 04 adds the speciation solver.

And, as section 1 says, the forward and back laws are independent. They do not
share a solubility product, so the mineral in this case does not have a
well-defined equilibrium — it has whichever steady state the two rate constants
and the transport happen to produce.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 15_precip_and_dissolution run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the staggered slot and reports porosity, refusing to continue if it does not percolate. Standard library only. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and reports the **net** change in porosity only. **[v1.3]** It does not test for the chattering of section 5: it compares the first and last rows and prints FELL, ROSE or unchanged, so a porosity that oscillates and returns to its starting value prints "unchanged". To check for chattering yourself, plot the porosity column of `output/summary.csv` against iteration and look for oscillation with a period equal to `<update_interval>`. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
