# 17 — a rate law you can read, in a file rather than the binary

## 1. What this example does

This is example 05's colony — the same aerobic organism, the same cellular
automaton spreading its biomass — with one thing changed: **the rate law is not
compiled into the program**. It is four lines of algebra in a text file, read at
start-up, printed into the log, and evaluated as an expression tree in every
voxel at every step.

The reaction is aerobic respiration of acetate:

```
    CH3COO(-)  +  2 O2   ->   2 HCO3(-)  +  H(+)          +  biomass
```

and the law itself, `input/growth.sym`:

```
    vars    acetate o2 Bug

    rate    growth  = 0.35 * acetate / (0.05 + acetate) * o2 / (0.01 + o2)
    rate    acetate = -2.5 * growth * Bug
    rate    o2      = -5.0 * growth * Bug

    range   acetate 1e-6 5.0e-3
    range   o2      1e-6 2.0e-3
```

That is a **dual-Monod** law: growth is limited by acetate and by oxygen
independently, and the product of the two saturation terms means whichever is
scarcer controls the rate. `0.35 h⁻¹` is the maximum specific growth rate;
`0.05` and `0.01 mol L⁻¹` are the half-saturation constants for the donor and
the acceptor.

The two substrate lines are written as multiples of `growth`, not fitted
separately. `-2.5` is one over a yield of 0.4 mol biomass per mol acetate;
`-5.0` is two oxygen per acetate, straight from the reaction above. Writing them
this way means the **stoichiometry cannot drift**: change the growth expression
and the substrate draws follow it exactly.

Why it matters that this is a file rather than a header. In example 05 a
different rate law is a rebuild. Here a different rate law is a different file —
which makes it practical to run the same case against several candidate laws, to
ship a law alongside a dataset, and to have a law come out of a fitting procedure
rather than out of a person's judgement about which form to assume.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 2 of them: `acetate`, `o2` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Bug` — a symbolic rate law read from a `.sym` file |
| **Biomass** | `Bug` — attached biofilm (seeded on its own material number), cellular automaton — biomass stays put until a voxel fills, then spills into neighbours |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 1000 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. How the file is bound to the simulation

**The variable names are the binding.** Every name on the `vars` line must match
a `<name_of_substrates>` or the organism's `<name_of_microbes>` exactly:

```xml
    <name_of_substrates>acetate</name_of_substrates>       vars acetate
    <name_of_substrates>o2</name_of_substrates>            vars o2
    <name_of_microbes>Bug</name_of_microbes>               vars Bug
```

A mismatch stops the run at start-up and names the offending variable. It does
not bind to the wrong lattice and it does not quietly evaluate to zero — which
is the failure the positional index conventions elsewhere in this code are
vulnerable to, and the reason names are worth the parser.

`<half_saturation_constants>` in `CompLaB.xml` is **unused** on this path and is
set to `0 0` to say so. The `.sym` file carries its own constants; there is one
place a half-saturation constant lives and it is the expression.

**The ranges are enforced, not advisory.** Every evaluation clamps each input to
its `range` line, and the closing report says how many evaluations were clamped.
A law fitted on a range means nothing outside it, and this is the one path in
the repository where that limit is written down beside the law itself and
checked at run time. The inlet concentrations here — `4.0e-3` acetate,
`1.5e-3` oxygen — sit inside the fitted ranges deliberately.

| tag | value | what it does |
|---|---|---|
| `<symbolic><enabled>` | `true` | builds and enables the expression-tree path |
| `<expressions_file>` | `input/growth.sym` | which law to read |
| `<reaction_type>` | `symbolic` | routes *this organism* through it |
| `<solver_type>` | `CA` | biomass spreads by cellular automaton, as in example 05 |
| `<decay_coefficient>` | `0.0` | no decay, so growth is the only biomass term |

## 4. Where a law like this comes from

The case ships `input/growth.sym`, so `offline.sh` is optional. It is how you
would find a law from data of your own.

`training/fit_symbolic.py` is a **symbolic regression**: it searches over the
*form* of the expression, not over the constants in a form you chose. Ordinary
fitting starts from an assumption — you decide the law is Monod and the computer
finds the three numbers. If the truth is not Monod you get the best Monod there
is, and no hint that you asked the wrong question. This builds expressions out
of `+ - * /` and your variables, breeds the ones that fit, and returns the
expression. Nobody tells it about Monod; if the data is Monod it finds Monod.

**What comes back is a Pareto set, not an answer.** One expression per node
count, from a crude two-term form up to a long one that fits better. The choice
among them is yours, and the shortest expression whose accuracy you can live
with is almost always right: a longer expression that fits marginally better is
usually fitting noise and will not survive extrapolation.

```
    1.  fit        search over expression forms      -> a Pareto set
    2.  choose     take the elbow, not the best fit
    3.  finish     write the substrate lines by hand as multiples of the
                   fitted rate, so stoichiometry cannot drift
    4.  bound      add one range line per variable
```

Steps 3 and 4 are yours because the search cannot do them: it fits **one output
column** and knows nothing about your stoichiometry, your units, or which
variable is a concentration and which a biomass.

`training/growth_samples.csv` ships with the case — 400 points drawn from the
same dual-Monod law with 2% noise on the growth column, so the search has
something realistic to work on. The search writes
`input/growth_discovered.sym` and does **not** install it; compare it against the
shipped law before copying it over.

For production work PySR searches harder and is a published tool people know.
`fit_symbolic.py` is here so the loop is complete without a Julia install.

## 5. What to check

1. **the law is echoed at start-up.** The parsed expressions are printed; read
   them and confirm they are the law you meant. This is the cheapest check in
   the repository;
2. **growth is co-limited.** Lower the acetate inlet and growth falls; lower the
   oxygen inlet and it falls too. A law responding to only one input has a
   mis-bound variable or a saturation term that never engages;
3. **the clamp count stays low.** A large number of range clamps means the
   simulation is asking the law about conditions it was never fitted for;
4. **acetate and oxygen are drawn 1 : 2.** That ratio is written into the two
   substrate lines and nothing in the run can change it;
5. **no `[NEG!]` warnings.** The expression has no stoichiometry enforcing
   non-negativity, so the increment is clamped against what is locally present
   before it is applied.

## 6. What this case does not do

No flow, no abiotic reaction, no geometry evolution, no thermodynamic control on
the rate. The law is exactly what the file says — there is no pH dependence, no
temperature dependence, and no inhibition term, because none was written.

Example 05 is the same colony with the law compiled in, which is faster per
evaluation; example 18 is the same idea with a graph network in place of the
expression; example 19 multiplies a rate law of this shape by a thermodynamic
factor.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 17_symbolic_law run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and reports porosity, refusing to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Optional. Searches for a rate law from a data table and writes a candidate `.sym`. The case ships a law, so nothing has to be fitted before running it. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's checks. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `input/growth.sym` | The rate law of section 1. Read at start-up, echoed into the log, evaluated per voxel. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

### The offline code this case ships

| File | What it is |
|---|---|
| `training/fit_symbolic.py` | The symbolic regression of section 4. Searches over expression forms and returns a Pareto set. |
| `training/growth_samples.csv` | 400 sample points with 2% noise, so the search has something to work on out of the box. |

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
