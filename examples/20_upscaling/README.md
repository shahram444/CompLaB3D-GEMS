# 20 — from one aggregate to a continuum rate

## 1. What this example does

This case simulates one microbial aggregate — 120 µm across, 896 voxels of it —
held in seawater of fixed composition, and measures one number: **how much of that
aggregate is actually working.**

The aggregate oxidises methane with sulfate. Methane and sulfate arrive from the
surrounding water and diffuse inward; sulfide and bicarbonate are produced inside
and diffuse outward. Both gradients run the same way — deeper into the aggregate
means less reactant and more product — and two things follow from that.

First, **the rate falls with depth**, because the Monod terms fall with the
reactants. Second, **beyond a certain depth the rate is exactly zero**, because
the accumulated products have pushed ΔG above the 12.5 kJ mol⁻¹ the organisms need
to make ATP. The aggregate is not a uniform reactor with a slow middle. It is a
reacting shell around a dead core.

A sediment-column model cannot resolve any of that. It has one grid block, one
methane concentration in it, and it needs one rate. Put that single concentration
through the same rate law and you get **0.000349 mol L⁻¹ s⁻¹**. Sum what the
resolved aggregate is really doing and you get **0.000288**. The first number is
21 % too high, and refining the column-scale grid does not help, because the whole
error lives inside the aggregate.

The ratio of those two numbers is the **effectiveness factor**:

```
    eta  =  0.000288 / 0.000349  =  0.823
```

That is the correction. Multiply the single-concentration estimate by it and the
column-scale model agrees with the resolved one again.

So this run does three things: it integrates the resolved aggregate until the
internal profile stops changing, it measures `eta`, and it reports `eta` against
the two dimensionless groups that determine it — so the result can be applied at
sizes and compositions this run did not simulate.

The second of those groups is the point of the case. One is the Thiele modulus,
which is a century old and covers the falling rate. The other is the energy
threshold, which covers the dead core, and which no classical treatment contains.

### What is simulated

| | |
|---|---|
| **Geometry** | 32 x 34 x 10 voxels at 10 um, 8192 open (porosity 0.7529) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 4 of them: `CH4`, `SO4`, `HS`, `HCO3` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `ANME` — compiled kinetics (`defineKinetics.hh`) |
| **Biomass** | `ANME` — attached biofilm (seeded on its own material number), cellular automaton — biomass stays put until a voxel fills, then spills into neighbours |
| **Geometry evolution** | none — the pore space is fixed |
| **Thermodynamic control** | on — every rate multiplied by F_T from `input/*.thm` |
| **Upscaling** | on — reports the effectiveness factor and Thiele modulus |
| **Run length** | 24000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 2. The reaction

Sulfate-dependent anaerobic oxidation of methane (AOM), catalysed by a syntrophic
consortium of ANME archaea and sulfate-reducing bacteria in a dense aggregate. It
is the dominant methane sink in marine sediments, consuming most of the upward
methane flux before it reaches the water column.

**CH₄ + SO₄²⁻ → HS⁻ + HCO₃⁻ + H₂O**

The stoichiometry is 1:1:1:1, and two consequences follow directly from it:

- the **reactants are supplied externally** and must diffuse inward;
- the **products are generated internally** and must diffuse outward.

The interior therefore sits at both ends of two opposing gradients: reactant
concentrations decrease inward, product concentrations increase inward.

## 3. The intrinsic rate law

The intrinsic rate — the volumetric rate a voxel achieves at its own local
concentrations, in the absence of any transport limitation — is a dual Monod
expression in the electron donor and the terminal electron acceptor, coded in
`kinetics/defineKinetics.hh`:

```
                            [CH4]                [SO4]              B
    r   =   mu_max  ·  ───────────────  ·  ───────────────  ·  ─────────
                        K_CH4 + [CH4]       K_SO4 + [SO4]           Y


    r    volumetric reaction rate,  mol L-1 s-1
    B    biomass concentration,     gDW L-1
```

| symbol | value | unit | meaning | source |
|---|---|---|---|---|
| K_CH4 | 1.0 × 10⁻³ | mol L⁻¹ | half-saturation constant, methane | Nauhaus et al. (2002) |
| K_SO4 | 5.0 × 10⁻⁴ | mol L⁻¹ | half-saturation constant, sulfate | Nauhaus et al. (2002) |
| Y | 0.25 | gDW mol⁻¹ | growth yield, ≈1 % of methane carbon assimilated | Nauhaus et al. (2007) |
| μ_max | 1.0 × 10⁻¹ | s⁻¹ | maximum specific growth rate — **not the measured value, see §11** | — |

**The rate law is a function of the reactants alone.** Sulfide and bicarbonate do
not appear in it. Unmodified, it sustains the intrinsic rate everywhere the local
reactant concentrations permit, irrespective of how far the products have
accumulated — which is precisely the regime in which a Monod formulation ceases
to be adequate.

## 4. The thermodynamic constraint, and the dead core

A reaction proceeds only where it is exergonic, and a microbially catalysed one
only where the energy yield exceeds the organism's ATP requirement. The
`<thermodynamics>` block supplies that condition as a dimensionless factor on the
intrinsic rate, from `input/aom.thm`:

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

`a_i` are activities, `nu_i` the stoichiometric coefficients of §2 (negative for
reactants), and `X` is the average stoichiometric number, 1 here.

The formulation follows Craig & Meile (2024), Ch. 3. One sign convention is worth
noting: Craig & Meile write Eq. 3.4 as `f = -dG - m dG_ATP`, a positive driving
force, and this implementation carries the negative of it, so that `f < 0` means
the reaction can proceed. The resulting `F_T` is identical either way.

Every parameter is cited in `input/aom.thm`; those that set the threshold are:

| symbol | value | unit | meaning | source |
|---|---|---|---|---|
| dG0 | −33.12 | kJ mol⁻¹ | standard free energy at 4 °C, the seep temperature | He et al. (2021) |
| m · dG_ATP | 0.25 × 50 = 12.5 | kJ mol⁻¹ | minimum energy quantum for ATP synthesis | Hoehler et al. (2001) |
| dG_loss | 25.7 | kJ mol⁻¹ | interspecies electron-transfer loss, ANME to SRB (Craig & Meile, Eq. 3.5) | calibrated against measured growth efficiency |

**Why the core goes to zero rate rather than to a low one.** Substitute this
stoichiometry into equation 1:

```
              [HS-] [HCO3-]
    Q   =   ─────────────────
             [CH4] [SO4(2-)]
```

Products in the numerator, reactants in the denominator. Toward the centre of the
aggregate the products have accumulated and the reactants are depleted, so `Q`
increases, `dG` increases with `ln Q`, and at some depth the driving force `f`
crosses zero and `F_T` vanishes identically.

The aggregate therefore develops a **reacting shell and a dead core**, separated
by an internal front whose position is set by the local energetics rather than by
the geometry. That front is the reason this case exists.

## 5. The upscaling problem

A pore-scale simulation resolves everything above. A column-, basin- or
reservoir-scale reactive-transport model cannot: it carries a single methane
concentration per grid block, possibly a metre across, and requires a single
volumetric rate to close the mass balance.

The obvious closure is to evaluate the intrinsic rate law at the block-averaged
concentration. It is wrong, because **the volume average of a nonlinear rate is
not the rate evaluated at the volume-averaged concentration**. The standard
closure is the effectiveness factor:

```
              < r >_aggregate
    eta  =  ───────────────────
                r( C_bulk )
```

the volume-averaged realised rate over the rate the same law returns at the bulk
composition. A continuum model that applies `eta` to `r(C_bulk)` is closed
consistently; one that omits it carries an error equal to `1 - eta`.

## 6. The two dimensionless groups that set eta

**The Thiele modulus** — the classical group, the ratio of the intrinsic reaction
rate to the rate of diffusive supply over the aggregate's characteristic length
(Thiele 1939):

```
                    ___                        r( C_bulk )
    phi   =   R  · / k        where    k   =  ─────────────      s-1
                 \/  D                            C_bulk

    R    aggregate radius, m        D    diffusivity inside the aggregate, m2 s-1
```

For a sphere with first-order kinetics there is an exact solution (Aris 1975):

**[v1.3]** A caveat this document did not previously carry: the classical
solution is posed on a **sphere**, while what `preprocess.py` actually draws is a
circular **cylinder** spanning the full depth in z. Its membership test is on x
and y only, and a sphere is not available at this size anyway, the drawn diameter
being 12 voxels with only 8 z layers of water. Part of any departure between the
measured curve and the classical one is therefore a shape effect, and not only
the kinetics and the energetics.

```
                    3     [        1         ]
    eta_class  =  ─────  · [  phi ─────── - 1 ]
                  phi^2   [      tanh phi     ]
```

[v1.3] That box used to carry the OTHER spelling of the same solution,
`(1/phi)[1/tanh(3 phi) - 1/(3 phi)]`, which belongs to the generalized (Aris)
modulus `(R/3)sqrt(k/D)` and not to the radius modulus defined just above it.
Pairing the two evaluates the curve at three times the right argument: at this
case's own `phi = 0.924008` the old box gives 0.700 where the solver reports
0.947345. That mispairing was the v1.3 bug described in section 9; it was fixed
in `src/complab3d_thiele.hh` and in every number this document quotes, and this
formula box was the last copy of the old convention left standing.

The solver reports it alongside the measured value. Posing the problem on the
geometry for which the classical solution exists is what makes any departure
interpretable: it is then attributable to the kinetics and the energetics rather
than to the shape.

**The thermodynamic factor at the bulk** — `F_T(C_bulk)`, equation 4 of §4
evaluated at the exterior composition. This group appears in no classical
treatment. Under Monod kinetics alone the interior is merely rate-limited; with
the energy constraint imposed it is inactive, and no first-order analysis
predicts a dead core. `eta` is therefore reported against both groups, and the
sweep in §10 has two axes.

## 7. Conditions required for a steady measurement

`eta` is defined at steady state for a fixed bulk composition. Three settings
differ from example 19 to satisfy that definition; each is necessary.

**The bulk composition is imposed.** All four species take Dirichlet conditions at
both ends of the domain, at their initial values. This is the boundary condition
every effectiveness-factor definition assumes.

Imposing it on the *products* as well does not predetermine the result. The
thermodynamic factor falls because sulfide accumulates **within** the aggregate
relative to the bulk, and that internal gradient remains free. Only the exterior
composition is fixed, which is precisely the single value a continuum grid block
carries.

**The biomass field is held fixed.** `<freeze_biomass>true</freeze_biomass>`
discards the biomass increments while applying the solute increments: consumption
and production continue, but the catalyst distribution does not evolve.

Without it `eta` does not converge. Biomass accumulates faster at the rim than in
the interior, because the rim is better supplied, so the measured value continues
to drift on the growth timescale for reasons unrelated to transport. Every run in
the sweep reported `NOT STEADY` at five diffusion times until this was imposed.

**The integration is long enough.** The intra-aggregate concentration profile
relaxes on the diffusive timescale `R^2/D`. At `R` = 60 µm and
`D` = 1.5e-10 m2/s that is 24 s; the run covers five of them — 24 000 steps at
5 ms each.

## 8. Running it

```bash
./scripts/setup_case.sh 20_upscaling run/mycase
cd run/mycase
./pipeline.sh
```

Four steps, about four minutes on one core:

| | what happens |
|---|---|
| 1 | `preprocess.py` draws the aggregate and the box around it |
| 2 | `cmake` builds the solver against your Palabos |
| 3 | `./complab CompLaB.xml` runs, writing `output/upscaling.csv` |
| 4 | `postprocess.py` says whether the number is worth quoting |

## 9. Expected output

The closing report carries the whole result. These are the numbers from the
shipped configuration, so you have something to compare against:

```
  [UPSCALE] one aggregate, reduced to the two numbers a continuum model needs
  [UPSCALE]   aggregate 896 voxel(s), bulk 6784 voxel(s)
  [UPSCALE]   bulk concentration      0.00975972 mol/L
  [UPSCALE]   rate at bulk            0.000349462 mol/L/s
  [UPSCALE]   measured mean rate      0.000287729 mol/L/s
  [UPSCALE]   effectiveness factor    0.823348
  [UPSCALE]   Thiele modulus          0.924008   (classical eta 0.947345)
  [UPSCALE]   gate at bulk            0.980596
  [UPSCALE]   steady: the last three intervals agree to 0.00011%
```

**Read `eta` = 0.823 two ways, both correct.** 17.7 % of the aggregate's intrinsic
capacity is lost to internal transport and thermodynamic limitation; equivalently,
a continuum model applying `r(C_bulk)` unchanged over-predicts the methane
oxidation rate by 21.5 %, since it divides rather than multiplies by `eta`.

Two things to check before anything else.

**At iteration 0, `eta` must be exactly 1.** Every voxel is still at the bulk
composition, so the mean rate inside the aggregate is the rate at the bulk by
construction. The shipped case reports `1.0000000002`. That is the measurement
checking itself, and if it is not 1 then nothing else in the record means
anything. `postprocess.py` tests it first for that reason.

**`NOT STEADY` means the value is not converged.** The solver compares the last
three diagnostic intervals and reports it if they differ by more than 1 %. An
effectiveness factor evaluated during a transient is not a steady-state quantity.

**Why measured `eta` falls short of the classical value here.** 0.823 against
0.947, a ratio of 0.87. Diffusion alone, acting on a first-order law, would cost
this aggregate only about 5 %. It loses 18 %, so something beyond transport is
holding the interior back, and the something is the thermodynamic gate: the
products the reaction makes accumulate towards the middle faster than diffusion
carries them out, and the free energy available falls with them. At this bulk
composition the gate is barely engaged in the pore water (`F_T` = 0.98) and
still costs three times what diffusion does inside the aggregate, which is the
whole reason the gate cannot be evaluated at the bulk composition and applied
as a constant.

Monod saturation pushes the other way, and the two do not cancel: a saturating
law sustains a near-maximal rate at a concentration where a first-order law
would already have declined, which on its own would lift the ratio above one.
The measured ratio of 0.87 is what is left after the gate has more than paid
that back.

**[v1.3] This conclusion is the reverse of what earlier versions reported.**
`classicalEta` was written for the generalized Thiele modulus, `(R/3)sqrt(k/D)`,
and was being handed the radius-based modulus this case prints, so every
classical value was evaluated at three times the right argument and came back
about 0.15 too low. The reported ratio was 1.18 rather than 0.87 and the
measured factor appeared to sit above the classical curve. Every `eta_classical`
number in this README, in `offline/expected_sweep.csv` and in the `[UPSCALE]`
line above has been recomputed.

## 10. Mapping the curve: the parameter sweep

A single run yields `eta` at one aggregate radius and one bulk composition. A
continuum model will be evaluated over a range of both, so it requires
`eta(phi, F_T)` rather than a single value:

```bash
python3 offline/upscale.py --radii 3,4,5 --bulk 1e-5,3e-4,1e-3
```

Nine solver runs. **Sulfide is the second axis** because it is the product that
drives `F_T`: raising the bulk sulfide moves the entire aggregate closer to the
energy threshold. Each run is sized by its own diffusive timescale, so doubling
the radius quadruples the integration — start with `--radii 3,4` while checking
the setup.

`offline/expected_sweep.csv` holds six points produced this way:

| R (vox) | bulk HS⁻ | phi | F_T(bulk) | eta | eta_class | ratio |
|---|---|---|---|---|---|---|
| 3 | 1 × 10⁻⁵ | 0.608 | 0.993 | 0.955 | 0.976 | 0.979 |
| 3 | 3 × 10⁻⁴ | 0.601 | 0.970 | 0.949 | 0.977 | 0.972 |
| 3 | 1 × 10⁻³ | 0.584 | 0.916 | 0.936 | 0.978 | 0.957 |
| 4 | 1 × 10⁻⁵ | 0.715 | 0.989 | 0.921 | 0.967 | 0.952 |
| 4 | 3 × 10⁻⁴ | 0.707 | 0.965 | 0.913 | 0.968 | 0.943 |
| 4 | 1 × 10⁻³ | 0.686 | 0.910 | 0.894 | 0.970 | 0.922 |

The script fits

```
    eta
  ─────────   =   0.602  +  0.368 · F_T( C_bulk )          R2 = 0.42
  eta_class
```

**How to read the table.** Along the sulfide axis at fixed radius the ratio
decreases monotonically: 0.979 to 0.972 to 0.957 at `R` = 3, and 0.952 to 0.943
to 0.922 at `R` = 4. This is the thermodynamic constraint. As the bulk
composition approaches the threshold the dead core grows and the measured factor
falls further below the classical curve. The ratio also falls with radius at
fixed composition, 0.979 down to 0.952 at the low-sulfide end, which the Thiele
modulus alone would not produce: it is already in `eta_class`. Both trends say
the same thing, that the sulfide axis carries a dependence the Thiele modulus
does not contain.

**Six points is not a constitutive relation.** R² = 0.42 over this range, which
is weak, and the span of the ratio is narrow, 0.922 to 0.979. Every point lies
below the classical curve, so the gate more than pays back what Monod saturation
lends. Extend the sulfide axis far enough for the ratio to move by more than its
own scatter before quoting a fit, and do not read the intercept as a physical
quantity.

**On evaluating `F_T` at the bulk composition as a substitute.** It is not one.
`r(C_bulk)` already contains `F_T(C_bulk)` — it is in the denominator of `eta` —
so whatever `eta` departs from 1 is exactly the error that shortcut incurs. In
the shipped run, 17.7 %.

## 11. On the rate constant

`mu_max` in `defineKinetics.hh` is about six orders of magnitude faster
than a real ANME–SRB consortium, which doubles on a timescale of months (Orphan
et al. 2009). At the measured rate this case would need a simulated year to build
the product gradient it builds in a minute.

What is being matched is not the rate constant but the **Thiele modulus** — the
dimensionless group that decides whether an interior gradient forms at all. In a
real aggregate over months that group is large; raising `mu_max` is how a
short run reproduces the same regime. **Every other number is the measured one. Do
not quote that rate constant.**

## 12. Files

Everything this case needs is in this folder. Nothing to fetch.

| file | what it is |
|---|---|
| `preprocess.py` | builds the aggregate: a 120 µm circular cylinder spanning the full depth in z, in open water, in a box walled on four sides. Reports the equivalent-sphere radius the solver will use. |
| `CompLaB.xml` | the case. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml) |
| `kinetics/defineKinetics.hh` | the ungated dual-Monod rate law of §3 |
| `input/aom.thm` | the energetics of §4: Δ*G*°, the ATP threshold, the yield |
| `input/geometry.dat` | the pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything |
| `postprocess.py` | asks three questions in order: did the measurement chain work, has it stopped moving, what did the energy limit do |
| `offline/upscale.py` | the sweep of §10. Not part of `pipeline.sh` — it is one solver run per point |
| `offline/expected_sweep.csv` | the six reference points, to check your own against |
| `pipeline.sh` | runs steps 1 to 4 of §8 |

### A note on the domain

`x = 0` and `x = nx-1` carry the boundary conditions. The solver gives the other
four faces nothing — not a wall, not a symmetry plane, not periodicity — so
`preprocess.py` draws an inert wall there, and `NY` and `NZ` are two larger than
the water they hold. The layers are **added**, not taken out of the water, so the
aggregate and the bulk are exactly the volumes this case declares.

## 13. References

- Aris, R. (1975). *The Mathematical Theory of Diffusion and Reaction in Permeable
  Catalysts.* Oxford University Press.
- Craig, K. and Meile, C. (2024). Ch. 3, Eqs. 3.1–3.10 — the thermodynamic rate
  formulation this case implements, and the electron-transfer loss term of
  equation 3. PhD thesis, University of Georgia.
- He, X. et al. (2021). Thermodynamics of anaerobic methane oxidation.
- Hoehler, T. M. et al. (2001). Apparent minimum free energy requirements for
  methanogenic archaea and sulfate-reducing bacteria. *FEMS Microbiology Ecology*
  **38**, 33–41.
- Jin, Q. and Bethke, C. M. (2003). A new rate law describing microbial
  respiration. *Applied and Environmental Microbiology* **69**, 2340–2348.
- Nauhaus, K. et al. (2002). In vitro demonstration of anaerobic oxidation of
  methane coupled to sulphate reduction. *Environmental Microbiology* **4**,
  296–305.
- Nauhaus, K. et al. (2007). In vitro cell growth of marine archaeal–bacterial
  consortia during anaerobic oxidation of methane with sulfate. *Environmental
  Microbiology* **9**, 187–196.
- Orphan, V. J. et al. (2009). Geobiological signatures of anaerobic methane
  oxidation in marine sediments.
- Thiele, E. W. (1939). Relation between catalytic activity and size of particle.
  *Industrial & Engineering Chemistry* **31**, 916–920.
