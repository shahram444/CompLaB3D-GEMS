# 21 — pinning the by-products that plain FBA leaves free

## 1. What this example does

Example 09 solves one linear program per voxel and takes what the simplex
returns. For the growth rate that is enough: the optimum is unique. For the
**by-products it is not**.

Give the E. coli core model glucose and oxygen, both bounded at 10 mmol gDW⁻¹
h⁻¹, and ask flux variability analysis what the exchange fluxes may be *at
exactly the optimal growth rate of 0.5591 h⁻¹*:

| exchange | minimum | maximum |
|---|---|---|
| acetate | 9.9057 | 11.5033 |
| formate | 5.1129 | 11.5033 |

Every value in those ranges is optimal. Acetate and formate can absorb the same
electrons and the objective does not distinguish them, so the solver returns
whichever its pivoting rule happens to reach. Nothing about the answer is
*wrong* — it is one member of an optimal set — but it is arbitrary, it need not
survive a change of solver version, and it cannot be compared against a measured
yield.

This example replaces the single program with a **chain of programs**, each
optimising the next quantity while holding every earlier one at its optimum:

```
    1.  max  biomass                                  ->  mu*
    2.  max  acetate export,   biomass >= a1 * mu*
    3.  max  formate export,   biomass >= a1 * mu*
                               acetate >= a2 * (step 2)
```

After step 3 nothing is free. The same model at the same bounds returns the same
numbers, on any solver, every time.

This is the **multi-step FBA** of Song et al. (2025); the same device appears as
lexicographic optimisation in DFBAlab (Gomez, Höffner & Barton 2014). The
difference is only that the retain fractions `a_k` are free parameters here
rather than fixed at 1.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 × 26 × 8 voxels at 10 µm, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 4 of them: `glc`, `o2`, `ac`, `for` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Ecoli` — flux balance analysis, GLPK in process, **three-stage lexicographic chain** |
| **Metabolic model** | `e_coli_core.xml`, read natively as SBML — 72 metabolites, 95 reactions |
| **Biomass** | `Ecoli` — attached biofilm, finite difference — biomass spreads by diffusion on the same grid |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 500 advection-diffusion steps |

Glucose and oxygen are held at the left face and free to leave at the right.
Acetate and formate start at zero, are produced by the organism, and leave
through the right face — so the two stages that pin them are pinning quantities
that actually reach the transport step.

## 3. What the retain fractions do

`<retain_fraction>` is the fraction of a stage's optimum that later stages must
preserve. At 1.0 the chain is strict: nothing may give up any of an earlier
objective.

That is mathematically clean, and for real organisms often wrong. A cell that
maximises growth to the last decimal excretes nothing. Song et al. could not
reproduce the measured pyruvate and acetate yields of *Shewanella oneidensis*
until they let biomass fall to about 67 % of its theoretical maximum. Lowering
the first fraction here does the same thing:

| `a₁` | growth (h⁻¹) | acetate | formate |
|---|---|---|---|
| 1.0 | 0.5591 | 11.5033 | 5.1129 |
| 0.9 | 0.5031 | 12.3529 | 10.0423 |
| 0.6 | 0.3354 | 14.9020 | 16.5210 |

Growth falls, excretion rises. That is the trade, and it is a **calibration**,
not a mechanism: the value comes from fitting predicted to measured yields, not
from the network. Anything below 1.0 should be reported wherever the results
are.

This case ships at 1.0, so the run carries no fitted parameter at all while
still removing the ambiguity. That is the recommended starting point.

## 4. Why the reactions are named and not numbered

`<multi_step>` requires an SBML model, because it addresses reactions by name:

```xml
<stage_reactions>  Biomass_Ecoli_core   EX_ac_e   EX_for_e  </stage_reactions>
```

A column index stops meaning anything the moment the model file is regenerated,
and a wrong index is not caught at all — the chain would pin the wrong reaction
and the run would look fine. A wrong name stops the run at start-up and prints
the candidates it looked at. The SBML `R_` prefix is stripped by the reader, so
it is not written here.

Two stages resolving to the same reaction is also refused at start-up: it would
pin a quantity against itself and make every later stage vacuous.

## 5. What to check when it has run

**Growth and fluxes must come from the same vector.** The growth rate reported
is read out of the *final* stage's solution, not from the first stage's
objective value. With every fraction at 1.0 those are equal; at `a₁ = 0.9` they
differ by 10 %, and at 0.6 by 40 %. Using the first would have the solver add
biomass at one rate while removing substrate at a rate belonging to a different
solution, and mass would stop balancing with no sign of it in the output. If you
lower a fraction, confirm `summary.csv` still closes.

**The answer should be reproducible.** Run the case twice. Every flux should
agree to the last digit. If they do not, a stage is missing.

**Cross-check against COBRApy.** Change `<reaction_type>` to `cobrapy` and
`<enable_fba_glpk>`/`<enable_fba_cobrapy>` to match. The two back ends walk the
same chain through independent code, so a disagreement is a defect in one of the
couplings rather than a property of either solver. This is the strongest check
available on this path.

**Cost.** Three stages is three linear programs instead of one, but the stages
after the first change one bound and restart from the previous basis, so they
finish in a few pivots. Expect about 1.6× a single solve, not 3×.
`<track_performance>true</track_performance>` is on in this case so the closing
report gives the measured figure.

## 6. Files

| | |
|---|---|
| `CompLaB.xml` | the case |
| `input/geometry.dat` | 24 × 26 × 8 packed-sphere pore space, shared with example 09 |
| `input/e_coli_core.xml` | the metabolic model, SBML, uncompressed from `models/` |

## 7. References

- Song, H.-S., Ahamed, F., Lee, J.-Y., Henry, C. S., Edirisinghe, J. N., Nelson,
  W. C., Chen, X., Moulton, J. D., Scheibe, T. D.: Coupling flux balance analysis
  with reactive transport modeling through machine learning for rapid and stable
  simulation of microbial metabolic switching. *Scientific Reports* **15**, 6042
  (2025). https://doi.org/10.1038/s41598-025-89997-9
- Gomez, J. A., Höffner, K., Barton, P. I.: DFBAlab: a fast and reliable MATLAB
  code for dynamic flux balance analysis. *BMC Bioinformatics* **15**, 409
  (2014). https://doi.org/10.1186/s12859-014-0409-8
- Orth, J. D., Thiele, I., Palsson, B. O.: What is flux balance analysis?
  *Nature Biotechnology* **28**, 245–248 (2010).
  https://doi.org/10.1038/nbt.1614
