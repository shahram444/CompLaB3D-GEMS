# 19 — a reaction that stops where the energy runs out

## 1. What this example does

Every rate path in this repository answers one question: **how fast can this
organism go, given what is here.** A Monod law, a linear program, a fitted
network — all of them take the substrates and return a rate.

None of them asks the second question: **is there enough energy in this
reaction, at these local concentrations, for it to run at all.**

In a pore-scale domain those two questions have different answers in different
voxels, and this case is built to make the difference visible. A single microbial
aggregate, 120 µm across, sits in open water and oxidises methane with sulphate:

```
    CH4  +  SO4(2-)   ->   HS(-)  +  HCO3(-)  +  H2O
```

At the rim, methane and sulphate arrive from the surrounding water and sulphide
is carried away, so the reaction has around 29 kJ mol⁻¹ to give and runs freely.
Forty micrometres further in, the same reaction has raised its own products until
only about 12 kJ mol⁻¹ is left — less than the 12.5 kJ mol⁻¹ the organisms need
to make ATP. **A dual-Monod law keeps the reaction running there anyway**, because
it never looks at a product: in that same voxel its Monod factor is still well
above a half.

`<thermodynamics>` supplies the factor the rate law is missing. The aggregate
stops being a uniform reactor with a slow middle and becomes what it is in
reality: **a reacting shell around a dead core**, with an internal front whose
position is set by the local energetics rather than by the geometry.

The case runs twice — once gated, once with the gate switched off and nothing
else changed — because a gated run on its own says nothing about what the gate
did.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 32 x 32 x 10 voxels at 10 um, 7680 open (porosity 0.7500) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 4 of them: `CH4`, `SO4`, `HS`, `HCO3` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `ANME` — compiled kinetics (`defineKinetics.hh`) |
| **Biomass** | `ANME` — attached biofilm (seeded on its own material number), cellular automaton — biomass stays put until a voxel fills, then spills into neighbours |
| **Geometry evolution** | none — the pore space is fixed |
| **Thermodynamic control** | on — every rate multiplied by F_T from `input/*.thm` |
| **Upscaling** | on — reports the effectiveness factor and Thiele modulus |
| **Run length** | 12000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The intrinsic rate law

`kinetics/defineKinetics.hh`, an ordinary dual Monod in the electron donor and
the terminal electron acceptor, with **no thermodynamics in it whatsoever**:

```
                              [CH4]                [SO4]
    growth  =   mu_max  ·  ───────────────  ·  ───────────────  ·  B     gDW L-1 s-1
                            K_CH4 + [CH4]       K_SO4 + [SO4]

                growth
    r       =   ──────                                                  mol L-1 s-1
                  Y
```

| symbol | value | unit | meaning | source |
|---|---|---|---|---|
| K_CH4 | 1.0 × 10⁻³ | mol L⁻¹ | half-saturation constant, methane | Nauhaus et al. (2002) |
| K_SO4 | 5.0 × 10⁻⁴ | mol L⁻¹ | half-saturation constant, sulphate | Nauhaus et al. (2002) |
| Y | 0.25 | gDW mol⁻¹ | growth yield, ≈1 % of the methane carbon assimilated | Nauhaus et al. (2007) |
| μ_max | 1.0 × 10⁻¹ | s⁻¹ | maximum specific growth rate — **not the measured value, see §11** | — |
| kd | 0 | s⁻¹ | decay off, so the gate is the only thing acting on the biomass | — |

Sulphide and bicarbonate do not appear. That is the omission the rest of this
case is about.

## 4. The thermodynamic factor

At every voxel, for every organism that has a block in the `.thm` file:

```
                    nu_i
    (1)   Q   =   PROD  a_i                          reaction quotient
                    i

    (2)   dG  =   dG0  +  R T ln Q                   Craig & Meile, Eq. 3.6

    (3)   f   =   dG  +  m dG_ATP  +  dG_loss        Craig & Meile, Eq. 3.4

                                    (    f    )
    (4)   F_T =   max[ 0 ,  1 - exp ( ─────── ) ]    Jin & Bethke (2003);
                                    (  X R T  )      Craig & Meile, Eq. 3.3


    realised rate   =   r  x  F_T
```

`a_i` are activities, `nu_i` the stoichiometric coefficients of §1 (negative for
reactants), and `X` the average stoichiometric number, 1 here.

Optionally, an energy-based growth yield replaces the fixed `Y` of §3:

```
                                       1.8              [           2  0.16                  ]
    (5)   dG_dis  =  200 + 18 (6 - C)      +   exp      [ ((3.8 - g) )      (3.6 + 0.4 C)     ]

                            dG
    (6)   Y       =   ───────────────────
                       -( dG_ana + dG_dis )
```

Equation 5 is Heijnen & van Dijken (1992), Craig & Meile Eq. 3.9; equation 6 is
Craig & Meile Eq. 3.10. `C` is the carbon number of the substrate and `g` its
degree of reduction — 1 and 8 for methane.

**On the sign convention.** `dG` here is negative when the reaction releases
energy. Craig & Meile write Eq. 3.4 as `f = -dG - m dG_ATP`, a positive "energy
available", and carry the minus sign inside the exponential instead. The two are
the same expression written twice, and `F_T` is identical either way.

**Why the core stops rather than merely slowing.** Substitute this stoichiometry
into equation 1:

```
              [HS-] [HCO3-]
    Q   =   ─────────────────
             [CH4] [SO4(2-)]
```

Products in the numerator, reactants in the denominator. Inward, the products
have accumulated and the reactants are depleted, so `Q` rises, `dG` rises with
`ln Q`, and at some depth `f` crosses zero and `F_T` vanishes identically. Not a
small number — zero.

**The transition is sharp, and that is the equation, not a defect.** With
`chi = 1` the whole fall from open to shut happens over about 5 kJ mol⁻¹, a few
per cent of `dG0`; the half-width is exactly `chi · R T · ln 2`. In an aggregate
that appears as a front rather than a gradual fade. A larger `chi` widens it.

## 5. Why this is not a seventh rate path

It is a **factor**, not a rate. It multiplies whichever path each organism
already uses, so it composes with all six — `kinetics`, `glpk`, `cobrapy`,
`surrogate`, `symbolic`, `graphnet` — and none of them had to be rewritten to
accept it. This case puts it on top of compiled kinetics because that is the
plainest demonstration: the rate law of §3 is the same shape as example 05's.

**One restriction, on the kinetics path only.** `defineKinetics.hh` returns one
combined rate vector for every organism in a voxel at once, so there is no way to
give two organisms two different gates there. A `.thm` file with more than one
reaction block is refused at start-up when any organism is on that path. The
other five paths compute rates per organism and have no such restriction.

## 6. The energetics live in a file

```
    temperature 277.15
    energy_units kJ
    concentration_scale 1

    reaction AOM
      microbe   ANME
      stoich    CH4 -1  SO4 -1  HS 1  HCO3 1
      dG0       -33.12
      atp       0.25
      dGatp     50.0
      chi       1.0
      cmin      1e-9
      dGloss    25.7
      yield     on
```

Same reasoning as the `.sym` file of example 17 and the `.gnn` of example 18:
`dG0` and the ATP threshold are **cited claims** that get revised and quoted in a
paper. They belong in something a reader can open, diff and cite, not inside a
binary. The whole file is echoed into the log, so the result carries the
energetics that produced it.

| symbol | value | unit | meaning | source |
|---|---|---|---|---|
| dG0 | −33.12 | kJ mol⁻¹ | standard free energy at 4 °C, the seep temperature | He et al. (2021) |
| m · dG_ATP | 0.25 × 50 = 12.5 | kJ mol⁻¹ | minimum energy quantum for ATP synthesis | Hoehler et al. (2001) |
| chi | 1.0 | — | average stoichiometric number | Jin & Bethke (2003) |
| cmin | 1 × 10⁻⁹ | mol L⁻¹ | activity floor, so a species driven to zero cannot send ln Q to −∞ | — |
| dG_loss | 25.7 | kJ mol⁻¹ | interspecies electron-transfer loss, ANME to SRB (Craig & Meile Eq. 3.5) | calibrated, see §10 |

`CompLaB.xml` carries two lines:

```xml
    <thermodynamics>
        <enabled>true</enabled>
        <energetics_file>input/aom.thm</energetics_file>
    </thermodynamics>
```

Absent or false, the run is bit-identical to one built before this block existed.

## 7. Four things to get right in a file of your own

**1. Quote `dG0` at the temperature you are running at.** −33.12 kJ mol⁻¹ is the
value at 4 °C. The 25 °C number is different, and using it puts the threshold in
the wrong place by more than the entire width of the transition. There is no way
for the solver to catch this: both numbers are plausible.

**2. Write every product into the `stoich` line.** The gate exists precisely
because products accumulate. Leaving one out gives a gate that never closes,
which is indistinguishable in the output from having no gate at all.

**3. State the energy threshold, do not guess it.** `atp × dGatp` is what the
reaction must release before the organism can conserve any of it. 10–20 kJ mol⁻¹
is the band reported for anaerobic syntrophs (Schink 1997; Hoehler et al. 2001).
Copying 40 kJ mol⁻¹ from an aerobic paper shuts the gate everywhere.

**4. Check `concentration_scale`.** `ln Q` is not scale-free. A case posed in
mmol L⁻¹ and one posed in mol L⁻¹ give different energies from the same numbers.
This case is in mol L⁻¹, so the scale is 1; the log says so either way.

## 8. Look at the gate before running the solver

`./offline.sh` reads `input/aom.thm`, sweeps the composition path this case
follows from the aggregate rim to its centre, and prints `dG`, `F_T` and the
yield along it, together with the depths at which `F_T` falls through 0.5 and
reaches zero and the width of the transition in kJ mol⁻¹.

This is the counterpart of fitting a `.sym` law or training a `.gnn` network: the
step where you find out whether the thing you are about to run says anything at
all. **A gate that never closes along the path the run takes produces a finished
simulation that means nothing**, and that is trivial to see beforehand and
impossible to see afterwards.

## 9. The gate is self-limiting, so it settles on the threshold

The offline sweep says where `F_T` would reach zero **if** the composition kept
moving along that path. In the run it does not, and the reason is worth
understanding before reading the output.

The gate is inside the feedback loop. As `F_T` falls the reaction slows, so it
stops making the very products that were closing it, and the deepest voxel
settles asymptotically **on** the threshold rather than passing through it — `dG`
there sitting at the ATP cost exactly, with the rate suppressed by a large
factor. So the run reports a small minimum `F_T`, not a zero one, and that is
correct behaviour rather than a shortfall.

`F_T` reaches exactly zero only where **transport**, rather than the local
reaction, pushes the composition past the threshold: a product plume arriving
from a neighbouring aggregate, or a reactant removed.

## 10. The transfer loss, and where its one number comes from

Craig & Meile compute `dG_net = dG_cat - dG_act - dG_ohm` (Eq. 3.5): the
activation and ohmic losses along the conductive network between the two
syntrophic partners. This code has no solver for the cytochrome redox state, so
it cannot compute those terms. The `.thm` format carries them as a constant
instead — `dGloss`, added to `dG`.

**That number is not measured, and it is the only tuned quantity in the case.**
It is fixed by the one observational constraint available. The growth efficiency
of these consortia has been measured at 1.0–2.75 % of the methane carbon
(Nauhaus et al. 2007; Orphan et al. 2009). Equation 6 with no transfer loss
returns 5.18 % at the aggregate rim; a loss of 25.7 kJ mol⁻¹ brings it to
2.75 %, the top of the measured band.

**[v1.3]** The unloaded figure is 5.18 %, not the 4.7 % this section used to
quote. Recomputed from the shipped constants (`dG0` -33.12, `T` 277.15, `dGana`
-30.0, `carbon` 1, `reduction` 8.0, `yieldmax` 0.05) at the rim composition of
`offline/thermo_curve.py`: RT is 2.304225 kJ mol⁻¹, RT ln Q is -21.6760
kJ mol⁻¹, the Heijnen dissipation for one carbon at degree of reduction 8 is
1088.09 kJ per C-mol so the denominator is 1058.09, and dG is -54.796, giving
5.18 %. Note also that `yieldmax` 0.05 then clamps what the solver reports, so a
reader who follows the instruction below to set `dGloss` to 0 will see 5.00 %,
not 5.18 % and not 4.7 %. The 2.75 % figure is unclamped and stands as written.

That single number then determines where the gate closes, with nothing further
adjusted. It is the honest form of a fitted parameter: **one degree of freedom,
fixed against a measurement that is not the one being predicted.** Set `dGloss`
to 0 to recover the direct-contact assumption and watch the yield and the gate
move together — which is the correct coupling, since a loss that makes a reaction
less favourable must make it both slower and less productive.

## 11. On the rate constant

`mu_max` in `defineKinetics.hh` is about six orders of magnitude faster than a
real ANME–SRB consortium, which doubles on a timescale of months (Orphan et al.
2009). At the measured rate this case would need a simulated year to build the
product gradient it builds in a minute.

What is actually being matched is not the rate constant but the **Thiele
modulus**, `r L² / (D C)` — the dimensionless group that decides whether a
reaction inside a diffusive aggregate is supply-limited, and therefore whether an
interior gradient forms at all. In a real aggregate over months that group is
large; raising `mu_max` is how a short run reproduces the same regime. **Every
other number is the measured one. Do not quote that rate constant.**

## 12. The control run

`pipeline.sh` runs the case twice — once as shipped, once with
`<enabled>false</enabled>` and nothing else changed — so every difference between
the two field sets is the factor and only the factor. The ungated run oxidises
substantially more methane, because it keeps the reaction going in the interior
where the gated run has shut it down.

That comparison is what exposed an old defect worth recording: the solver used to
accept a file name on the command line and ignore it, so `./complab variant.xml`
ran `CompLaB.xml`. The control run therefore came back byte-identical to the
gated one, for a reason that had nothing to do with the gate. It is fixed, and
`postprocess.py` now checks that the two logs actually differ, so the failure
cannot recur silently.

## 13. From this aggregate to a continuum rate

`<upscaling>` is enabled in this case, and it reduces the resolved aggregate to
the one number a column- or reservoir-scale model needs — the **effectiveness
factor**:

```
              < r >_aggregate
    eta  =  ───────────────────
                r( C_bulk )
```

what the aggregate actually does, over what it would do if every point inside it
saw the bulk composition. The numerator is not recomputed: the rate processors
have already written their result into the increment lattices, and the diagnostic
samples those after the rate processors run and before the increments are
applied, so it averages the number the solver is about to use rather than a
second model of it.

**The gate is what makes this more than a textbook exercise.** For a sphere with
first-order kinetics the classical result is

```
                    3    [       1        1  ]                       ___
    eta_class  =  ───── · [  phi ─────── - 1  ]      phi  =  R · / k
                  phi^2   [      tanh phi     ]                    \/  D
```

and it is printed beside the measured value. Both the modulus and the formula
are written for the sphere RADIUS. That pairing matters: the same solution
written for the generalized modulus, `(R/3)sqrt(k/D)`, looks quite different, and
mixing the two conventions evaluates the classical curve at three times the right
argument. Versions before v1.3 did exactly that, so every classical value they
printed was too low by about 0.15. Under Monod kinetics alone the two
track each other. With the gate on they do not, because the core does not merely
run slowly — it **stops**, at a moving internal boundary no first-order Thiele
analysis contains. So a second group is reported beside `phi`: `F_T(C_bulk)`, the
gate at the bulk composition.

This case measures `eta` while the aggregate is still growing and while the
surrounding water is being flushed, which is enough to demonstrate the diagnostic
but not enough to quote. **Example 20 is the same aggregate posed so that `eta`
is a steady-state quantity** — every species held Dirichlet at both ends and the
biomass field frozen — and it sweeps radius and bulk sulphide to give
`eta(phi, F_T)` rather than a single point. Read that case before using any
number from this one.

## 14. What to check

1. **the `.thm` file is echoed at start-up**, and the reaction, the threshold and
   the temperature in it are the ones you meant;
2. **the offline curve closes the gate somewhere along the path** — §8. If it
   does not, stop before running the solver;
3. **the gated and ungated runs differ.** `postprocess.py` checks that the two
   logs are not identical, and that check exists because they once were, for the
   wrong reason;
4. **the minimum `F_T` is small but not zero.** §9 explains why that is the
   expected result, and a run reporting exactly zero everywhere in the core means
   transport, not the local reaction, is pushing the composition past the
   threshold;
5. **the reaction front is sharp.** With `chi = 1` the gate closes over a few
   kJ mol⁻¹; a gradual fade means the threshold or `chi` is not what you think;
6. **no `[NEG!]` warnings**, and the four species stay in 1 : 1 : 1 : 1 balance —
   the gate scales the whole rate vector, so it cannot break the stoichiometry.

## 15. What this case does not do

No flow (`<Peclet>0</Peclet>`), no abiotic reaction, no geometry evolution, and
no aqueous speciation — the activities in equation 1 are concentrations, with no
ion-pairing, no ionic-strength correction and no carbonate equilibrium behind
`HCO3`. Example 04 adds the speciation solver.

The consortium is one biomass pool: nothing here resolves the electron transfer
between ANME and SRB, which is exactly why that loss enters as the constant
`dGloss` of §10 rather than being computed.

And the effectiveness factor measured here is not a steady-state one. Example 20
is the case posed for that.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 19_thermodynamic_gate run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the aggregate — a sphere of biomass in open water, in a box walled on four sides — and reports the equivalent-sphere radius the solver will use. Standard library only. |
| **offline** | `offline.sh` | Runs `offline/thermo_curve.py`: the gate these energetics describe, along the composition path the run will take, before the solver starts. §8. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). `pipeline.sh` runs it twice — gated and ungated. |
| **post** | `postprocess.py` | Asks, in order: did the gate act, do the two runs differ, and is the run worth believing. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order, including the ungated control of §12. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The aggregate and the water around it. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `input/aom.thm` | The energetics of §6 — dG0, the ATP threshold, the transfer loss, the yield. Echoed into the log. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

### The rate law and the offline code this case ships

| File | What it is |
|---|---|
| `kinetics/defineKinetics.hh` | The **ungated** dual-Monod law of §3. `setup_case.sh` copies it to `defineKinetics.hh` at the case root, over the default it laid down there. It contains no thermodynamics at all — that is the point of §5. |
| `offline/thermo_curve.py` | The gate these energetics describe, evaluated along the rim-to-centre composition path, before any solver runs. |
| `offline/upscale.py` | Sweeps aggregate radius and bulk sulphide, one solver run per point, and fits the correction the classical curve needs. Not part of `pipeline.sh`. See example 20, which is posed for this measurement. |

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

## References

- Craig, K. and Meile, C. (2024). Ch. 3, Eqs. 3.1–3.10 — the thermodynamic rate
  formulation this case implements, and the electron-transfer loss term of
  equation 3. PhD thesis, University of Georgia.
- Heijnen, J. J. and van Dijken, J. P. (1992). In search of a thermodynamic
  description of biomass yields for the chemotrophic growth of microorganisms.
  *Biotechnology and Bioengineering* **39**, 833–858.
- He, X. et al. (2021). Thermodynamics of anaerobic methane oxidation.
  *Environmental Science & Technology* **55**, 6739–6749.
- Hoehler, T. M. et al. (2001). Apparent minimum free energy requirements for
  methanogenic archaea and sulfate-reducing bacteria. *FEMS Microbiology Ecology*
  **38**, 33–41.
- Jin, Q. and Bethke, C. M. (2003). A new rate law describing microbial
  respiration. *Applied and Environmental Microbiology* **69**, 2340–2348.
- Jin, Q. and Bethke, C. M. (2005). Predicting the rate of microbial respiration
  in geochemical environments. *Geochimica et Cosmochimica Acta* **69**,
  1133–1143.
- Nauhaus, K. et al. (2002). In vitro demonstration of anaerobic oxidation of
  methane coupled to sulphate reduction. *Environmental Microbiology* **4**,
  296–305.
- Nauhaus, K. et al. (2007). In vitro cell growth of marine archaeal–bacterial
  consortia during anaerobic oxidation of methane with sulfate. *Environmental
  Microbiology* **9**, 187–196.
- Orphan, V. J. et al. (2009). Geobiological signatures of anaerobic methane
  oxidation in marine sediments. *Geobiology* **7**, 128–140.
- Schink, B. (1997). Energetics of syntrophic cooperation in methanogenic
  degradation. *Microbiology and Molecular Biology Reviews* **61**, 262–280.
