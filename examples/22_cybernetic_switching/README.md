# 22 — an organism that switches carbon source when the first runs out

## 1. What this example does

Example 21 fixes *which optimum* the solver returns. This one fixes a different
problem: *which substrate the organism is living on*, when more than one will do.

E. coli offered glucose and acetate together does not eat both. It eats the
glucose, **excretes acetate while doing so**, and only turns round and consumes
that acetate once the glucose is gone.

Flux balance analysis on its own cannot produce this. Hand the model both carbon
sources at once and it takes whichever combination maximises growth, together,
from the first step. A Monod term cannot produce it either, and for a sharper
reason: acetate is a *product* in the first regime and a *substrate* in the
second, and a concentration alone does not say which.

`<cybernetic>` treats the organism as a competition between separate growth
options, one per carbon source:

**1. Solve the model once per source**, each time with the other sources shut.
Two sources, two flux vectors — two different ways the cell could be making a
living right now.

**2. Score each option by the carbon it would bring in:**

```
    r_k  =  k_k C_k / (K_k + C_k)             the unregulated kinetic rate
    u_k  =  n_k r_k / SUM_j n_j r_j           the cybernetic variable
```

`n_k` is the carbon count — glucose 6, acetate 2 — which is what makes the
competition about carbon rather than about molecules. The `u_k` are
non-negative and sum to one.

**3. Blend the flux vectors with those weights:**

```
    v  =  SUM_k  u_k v_k                      growth and every exchange
```

This is the **cybernetic approach** of Ramkrishna & Song (2012), in the form
Song et al. (2025) used to couple it to reactive transport. Carbon uptake rate
as the objective of the competition is theirs.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 × 26 × 8 voxels at 10 µm, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 3 of them: `glc`, `o2`, `ac` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Ecoli` — flux balance analysis, GLPK in process, **two competing growth options** |
| **Metabolic model** | `e_coli_core.xml`, read natively as SBML — 72 metabolites, 95 reactions |
| **Biomass** | `Ecoli` — attached biofilm, finite difference — biomass spreads by diffusion on the same grid |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 500 advection-diffusion steps |

Glucose and oxygen are held at the left face. Acetate starts at 0.4 everywhere
so both options are live from the first step, and is free to leave through the
right face. Because glucose is supplied and acetate is not, the pore space
develops a gradient in *which regime the organism is in* — glucose-fed near the
inlet, acetate-fed deep in the biofilm — which is the pore-scale version of the
switch and the reason it is worth doing in 3-D at all.

## 3. What the two options look like

Solved separately on this model, at these bounds:

| option | growth (h⁻¹) | glucose | oxygen | acetate |
|---|---|---|---|---|
| glucose only | 0.5591 | −10.000 | −10.000 | **+9.906** |
| acetate only | 0.0524 | 0.000 | −6.510 | **−4.400** |

Negative is uptake, positive is release. The acetate column is the whole point:
one option makes it, the other eats it.

Blending those two with the cybernetic weights as glucose depletes:

| [glc] | [ac] | u_glc | u_ac | µ | glucose flux | acetate flux |
|---|---|---|---|---|---|---|
| 1.000 | 0.00 | 1.000 | 0.000 | 0.5591 | −10.000 | **+9.906** |
| 0.500 | 0.20 | 0.886 | 0.114 | 0.5011 | −8.857 | +8.270 |
| 0.100 | 0.40 | 0.836 | 0.164 | 0.4762 | −8.364 | +7.566 |
| 0.020 | 0.50 | 0.682 | 0.318 | 0.3978 | −6.818 | +5.354 |
| 0.002 | 0.50 | 0.224 | 0.776 | 0.1658 | −2.239 | **−1.197** |
| 0.000 | 0.50 | 0.000 | 1.000 | 0.0524 | 0.000 | **−4.400** |

The acetate flux crosses zero on its own. Nothing in the configuration says
where that happens; it falls out of the weights. That sign change is the switch.

## 4. Two things the implementation gets right, and why they matter

**Shutting a source forbids uptake, never release.** An organism growing on
glucose must stay free to excrete the acetate that the *other* option will later
consume. Closing the exchange outright would remove the very thing being
switched to, and the second regime would never arrive. In the code this is
`lb → 0`, with the upper bound untouched.

**The blend conserves mass exactly.** Each `v_k` satisfies `S v_k = 0`, and a
convex combination of vectors in the null space of `S` is itself in the null
space of `S`. So the blended vector is stoichiometrically consistent exactly,
not approximately. Averaging model outputs is not usually safe; here it is, and
that is the reason.

## 5. What it costs

Two sources is two linear programs per organism per voxel per step. That is the
honest cost — the method rests on having two separate flux vectors to blend.

Two things make it cheaper than 2× suggests. The programs differ only in which
bounds are open, so each warm-starts from the basis the previous left. And any
source whose weight falls below `<weight_floor>` is not solved at all — in a run
where one substrate dominates most of the domain most of the time, that skip is
most of the cost. Set `<weight_floor>0</weight_floor>` to solve every source
every time; the shipped 1e-3 moves the answer in the fourth decimal.

## 6. What this does not do

The weights depend on concentration only. No enzyme state is carried between
time steps, so the organism switches as fast as the concentrations move: there
is no lag and no diauxic plateau of the kind a full cybernetic model with
inducible enzyme levels produces.

That is the same simplification Song et al. made, and it is adequate exactly
when the switch is fast compared with transport — which is the pore-scale case
this code is for. It would not be adequate for a well-mixed batch reactor, where
the lag is most of what you are trying to reproduce.

## 7. What to check when it has run

**Look for the sign change.** Plot the acetate flux across the biofilm. Near the
inlet it should be positive, deep inside it should be negative. If it is
positive everywhere, glucose is reaching everywhere and the domain is too small
or the biomass too sparse to develop the gradient.

**Weights must sum to one.** Where any carbon is present the two weights sum to
1 exactly. Where there is none, both are zero and the organism is skipped
entirely — which is the right answer, not an infeasibility, and should not
appear in the infeasible count.

**Mass should close.** `summary.csv` should show the same closure as example 09.
If it does not, the blend is being applied to fluxes and growth from different
solves.

**Combining with example 21.** `<multi_step>` and `<cybernetic>` compose: each
of the two options is then solved as a lexicographic chain rather than a single
program, which is exactly the arrangement Song et al. ran. Add a `<multi_step>`
block to this case to see it. The cost multiplies — two sources times three
stages is six linear programs — so try it on a small domain first.

## 8. Files

| | |
|---|---|
| `CompLaB.xml` | the case |
| `input/geometry.dat` | 24 × 26 × 8 packed-sphere pore space, shared with example 09 |
| `input/e_coli_core.xml` | the metabolic model, SBML, uncompressed from `models/` |

## 9. References

- Song, H.-S., Ahamed, F., Lee, J.-Y., Henry, C. S., Edirisinghe, J. N., Nelson,
  W. C., Chen, X., Moulton, J. D., Scheibe, T. D.: Coupling flux balance analysis
  with reactive transport modeling through machine learning for rapid and stable
  simulation of microbial metabolic switching. *Scientific Reports* **15**, 6042
  (2025). https://doi.org/10.1038/s41598-025-89997-9
- Ramkrishna, D., Song, H.-S.: Dynamic models of metabolism: review of the
  cybernetic approach. *AIChE Journal* **58**, 986–997 (2012).
  https://doi.org/10.1002/aic.13734
- Ramkrishna, D., Song, H.-S.: *Cybernetic Modeling for Bioreaction
  Engineering*. Cambridge University Press (2018).
