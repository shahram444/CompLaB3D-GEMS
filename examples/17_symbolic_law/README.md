# 17 - Symbolic law: a rate law read from a text file instead of compiled in

## 1. The scenario

This is example 05's colony, the same aerobic organism, the same cellular
automaton spreading its biomass, with one thing changed: the rate law is not
compiled into the program. It is four lines of algebra in a text file, read
at start-up, printed into the log, and evaluated as an expression tree in
every voxel at every step.

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

That is a dual-Monod law: growth is limited by acetate and by oxygen
independently (Monod is the microbial-growth form of Michaelis-Menten
saturation kinetics), and the product of the two saturation terms means
whichever substrate is scarcer controls the rate. `0.35` per hour is the
maximum specific growth rate; `0.05` and `0.01` mol/L are the half-saturation
constants for the donor and the acceptor, the concentration at which growth
runs at half its maximum. The two substrate lines are written as multiples of
`growth` rather than fitted separately: `-2.5` is one over a yield of 0.4 mol
biomass per mol acetate, and `-5.0` is two moles of oxygen per mole of
acetate, straight from the reaction above. Writing them this way means the
stoichiometry cannot drift; change the growth expression and the substrate
draws follow it exactly.

Why it matters that this is a file rather than a header: in example 05, a
different rate law is a rebuild. Here, a different rate law is a different
file, which makes it practical to run the same case against several
candidate laws, to ship a law alongside a dataset, or to use a law that came
out of a fitting procedure rather than out of a modeller's judgement about
which functional form to assume.

## 2. The picture

```
   24 x 26 x 8 voxels at 10 um per voxel  =  240 x 260 x 80 um

   x=0                                                      x=23
    |                                                         |
    |....................  ####  ..........................  |
acet|...................  ####  ..........................  |zero
4e-3|@@@@@@@@@@..........  ########  ......................  |gradient
o2  |@@@@@@@@@@..........                                     |
1.5e-3                                                        |
held|                                                         |
    +---------------------------------------------------------+
      no flow: Peclet = 0. acetate and o2 diffuse in from x=0
      @@@@@ = the seeded biomass patch, 108 of the 3168 open voxels,
              against the wall blocks at x % 4 == 0

      #### = solid grains        porosity 0.6346, 3168 open voxels
```

Only the two end faces `x = 0` and `x = 23` carry a boundary condition. The
other four faces are an inert wall, drawn by `preprocess.py`. The biomass
patch sits on material number 3, at every `x` where `x % 4 == 0`, in the
three voxel rows next to the wall block there, through the full depth in z:
108 voxels in total, carved out of open pore rather than added to it, so the
porosity is unchanged from the abiotic cases.

## 3. What goes in

| | |
|---|---|
| **Domain** | 24 x 26 x 8 voxels, dx = 10 um, 3168 open voxels, porosity 0.6346 |
| **Flow** | none. `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Transport** | D3Q7 advection-diffusion lattice Boltzmann, two species: `acetate`, `o2` |
| **Diffusivity** | 5e-10 m2/s for `acetate`; 2e-9 m2/s in pore and 1e-9 m2/s in biofilm for `o2` |
| **Chemistry** | none compiled in; a symbolic expression read from `input/growth.sym` |
| **Biology** | `Bug`, one population, seeded on material 3 at initial density 1.0e-5, cellular automaton biomass movement, rate routed through the symbolic path |
| **Geometry change** | none. The pore space is fixed for the whole run |
| **Run length** | 1000 advection-diffusion steps, output every 200 |

Boundary conditions:

| Species | Left face (x = 0) | Right face (x = 23) | Initial |
|---|---|---|---|
| **acetate** | Dirichlet, held at 4.0e-3 | Neumann, zero gradient | 2.0e-3 |
| **o2** | Dirichlet, held at 1.5e-3 | Neumann, zero gradient | 1.0e-3 |
| **Bug** (biomass) | closed | closed | 1.0e-5 in the seeded patch, 0 elsewhere |

*Dirichlet* means the concentration at that face is pinned to a value no
matter what the interior does, a reservoir that keeps the colony fed.
*Neumann with a zero gradient* means nothing flows across that face, an open
outlet neither substrate reaches. *Closed* means the same thing on both
ends: biomass, which this solver moves only voxel to voxel, never crosses a
domain face at all.

## 4. The reaction

**The variable names are the binding.** Every name on the `vars` line in
`growth.sym` must match a `<name_of_substrates>` or the organism's
`<name_of_microbes>` exactly:

```xml
    <name_of_substrates>acetate</name_of_substrates>       vars acetate
    <name_of_substrates>o2</name_of_substrates>            vars o2
    <name_of_microbes>Bug</name_of_microbes>                vars Bug
```

A mismatch stops the run at start-up and names the offending variable. It
does not bind to the wrong lattice and it does not quietly evaluate to zero,
which is the failure a positional convention would be vulnerable to and the
reason names are worth the parser.

`<half_saturation_constants>` in `CompLaB.xml` is unused on this path and is
set to `0 0` to say so: the `.sym` file carries its own constants, so there
is exactly one place a half-saturation constant lives, the expression.

**The ranges are enforced, not advisory.** Every evaluation clamps each
input to its `range` line, and the closing report says how many evaluations
were clamped. A law fitted on a range means nothing outside it, and this is
the one path in the repository where that limit is written down beside the
law itself and checked at run time. The inlet concentrations in section 3,
`4.0e-3` acetate and `1.5e-3` oxygen, sit inside the fitted ranges
deliberately.

| tag | value | what it does |
|---|---|---|
| `<symbolic><enabled>` | `true` | builds and enables the expression-tree path |
| `<expressions_file>` | `input/growth.sym` | which law to read |
| `<reaction_type>` | `symbolic` | routes this organism through it |
| `<solver_type>` | `CA` | biomass spreads by cellular automaton, as in example 05 |
| `<decay_coefficient>` | `0.0` | no decay, so growth is the only biomass term |

## 5. What happens each step

The solver repeats this 1000 times:

1. **Stream and collide** `acetate` and `o2` on their D3Q7 lattices, which
   advances diffusion by one step. There is no velocity field to advect them
   with.
2. **Apply the boundary conditions**, re-pinning `acetate` to 4.0e-3 and
   `o2` to 1.5e-3 at the left face.
3. **Evaluate the symbolic expression tree** in every open voxel that holds
   biomass: clamp each input to its declared range, evaluate `growth` and
   the two substrate lines, and write the increments for `acetate`, `o2`
   and `Bug`.
4. **Apply the biomass increment** in place.
5. **Run the cellular-automaton rule.** Wherever a voxel's biomass now
   exceeds `<maximum_biomass_density>`, spill the excess into an open
   neighbouring voxel.

The stability caveat is on step 3: the expression has no stoichiometry
enforcing non-negativity, so the increment is applied as evaluated, with no
positivity clamp. A time step too large for the local concentrations
produces a `[NEG!]` warning exactly as the compiled kinetics paths do.

## 6. What comes out

```
output/
  acetate_0000200.vti  ...  acetate_0001000.vti     concentration of acetate
  o2_*.vti                                            concentration of oxygen
  Bug_*.vti                                            biomass density of the colony
  rate_acetate_*.vti  rate_o2_*.vti                    the reaction rate as a field, mol/L/s
  summary.csv                                          one row every 100 steps
  run.log                                              the whole run, including the [SYM] block
```

Numbers to expect, as orders of magnitude and directions rather than exact
values:

| Quantity | What it should do |
|---|---|
| Total biomass | Rises from the seeded patch and levels off as growth and the maximum density balance |
| Total acetate | Falls near the colony, replenished from the left boundary |
| Total o2 | Falls near the colony, replenished from the left boundary |
| `rate_acetate` | Negative wherever the colony is active; zero elsewhere |
| Clamp count in the `[SYM]` block | Low, since the inlet values sit inside the ranges deliberately |
| Minimum of either substrate | Zero or above. A negative value is a real failure |
| Wall clock | Seconds on one core. This is a small case on purpose |

## 7. What to check

1. **The law is echoed at start-up.** The parsed expressions are printed in
   the `[SYM]` block; read them and confirm they are the law you meant. This
   is the cheapest check available.
2. **Growth is co-limited.** Lower the acetate inlet and growth falls; lower
   the oxygen inlet and it falls too. A law responding to only one input has
   a mis-bound variable or a saturation term that never engages.
3. **The clamp count stays low.** A large number of range clamps means the
   simulation is asking the law about conditions it was never fitted for.
4. **Acetate and oxygen are drawn 1 : 2.** That ratio is written into the two
   substrate lines and nothing in the run can change it.
5. **No `[NEG!]` warnings.**

`postprocess.py` reports the change in every field's total and flags any
negative minimum, covering checks 4 and 5 indirectly through the totals, and
5 directly. It also greps `run.log` for the `[SYM]` block and prints it, so
checks 1 and 3 are read straight off that output. Check 2 needs two runs
with different inlet values, which is outside what a single `postprocess.py`
pass can do. This case declares no `<conserve>` tag: both substrates are fed
from a Dirichlet boundary, so no total here is closed, and naming either in
a mass-balance check would guarantee a failure that means nothing.

## 8. What this case demonstrates

**A rate law read from a file at start-up, rather than compiled in, and
produced by symbolic regression rather than assumed by hand.** The variable
names are the binding between the file and the simulation, the ranges are
enforced rather than advisory, and the same organism, geometry and boundary
conditions as example 05 now run on a law that lives outside the executable
entirely.

Where a law like this comes from: the case ships `input/growth.sym`, so
`offline.sh` is optional. It shows how you would find a law from data of
your own.

`training/fit_symbolic.py` performs a **symbolic regression**: it searches
over the *form* of the expression, not over the constants in a form you
chose. Ordinary curve fitting starts from an assumption, you decide the law
is Monod and the computer finds the three numbers; if the truth is not
Monod, you get the best Monod there is and no hint that you asked the wrong
question. This search instead builds expressions out of `+ - * /` and your
variables, breeds the ones that fit the data best, and returns the
expression it converged on. Nobody tells it about Monod; if the data is
Monod-shaped, it finds Monod.

**What comes back is a Pareto set, not one answer.** One expression per node
count, from a crude two-term form up to a long one that fits better. The
choice among them is yours, and the shortest expression whose accuracy you
can live with is almost always right: a longer expression that fits only
marginally better is usually fitting noise and will not survive
extrapolation.

```
    1.  fit        search over expression forms      -> a Pareto set
    2.  choose     take the elbow, not the best fit
    3.  finish     write the other substrate lines by hand as multiples of
                   the fitted rate, so stoichiometry cannot drift
    4.  bound      add one range line per variable
```

Steps 3 and 4 are yours because the search cannot do them: it fits one
output column and knows nothing about your stoichiometry, your units, or
which variable is a concentration and which a biomass.

### Letting the solver write the substrate lines instead

Step 3 above is the dangerous one, because those two lines carry a ratio that
nothing checks:

```
rate    acetate = -2.5 * growth * Bug
rate    o2      = -5.0 * growth * Bug
```

5.0 over 2.5 is two oxygen per acetate, which **is** the reaction. Mistype one
digit and the run still finishes, the fields still look smooth, and the model
creates or destroys oxygen every step for twelve thousand steps. Nothing in the
output can tell you, because a rate law is entitled to whatever numbers it was
given.

So the file can state the chemistry instead and let the solver do the
arithmetic. `input/growth_stoich.sym` ships beside `input/growth.sym` and is the
same law written this way:

```
reaction  acetate -1   o2 -2      the balanced reaction, negative consumed
yield     acetate 0.4             biomass made per unit of acetate
biomass   Bug                     which variable is the organism itself

rate      growth = 0.35 * acetate / (0.05 + acetate) * o2 / (0.01 + o2)
```

At start-up the solver writes one substrate line per species in the reaction:

```
    coefficient  =  ( stoichiometry / |stoichiometry of the yield species| ) / yield

    acetate:  -(1 / 1) / 0.4  =  -2.5
    o2:       -(2 / 1) / 0.4  =  -5.0
```

and prints both in the log marked `<- from the reaction, not the file`. Their
ratio is now arithmetic rather than typing, so it cannot disagree with the
chemistry. `tests/test_sym_stoich.cpp` proves the two files give bit-identical
rates at 40 different compositions, and that every way of writing the block
wrongly stops the run instead of being guessed at.

To use it, point the configuration at the other file:

```xml
<expressions_file>input/growth_stoich.sym</expressions_file>
```

Three rules keep it honest. A species may have a `rate` line **or** a place in
the `reaction` line, never both, since two sources of truth for one number is
what this removes. Species outside the reaction keep their hand-written lines,
so a law can be partly derived. And omitting all three keywords changes nothing,
so every `.sym` file written before this existed still loads untouched.

The same reasoning already governs example 18: its graph network predicts one
extent per reaction and forms the species rates from the stoichiometry, which is
why that file records a stoichiometric residual of exactly zero.

**The rule underneath all of it:** fit what you do not know, declare what you do.
You do not know the growth law, so let the search find it. You do know the
stoichiometry, so state it and never let a fit near it.

`training/growth_samples.csv` ships with the case, 400 points drawn from the
same dual-Monod law with 2 percent noise on the growth column, so the search
has something realistic to work on. Running `offline.sh` writes
`input/growth_discovered.sym` and does not install it: compare it against
the shipped law before copying it over.

This tool is here so the loop is complete without a separate install; a
published symbolic-regression tool would generally search harder for
production work.

### What the search actually finds here, and what it takes

The defaults in `offline.sh`, `--pop 200 --gens 15`, finish in under a minute
and return

```
    growth = 647.6 x acetate x o2                       about 4% error
```

a plain product with no saturation in it. That is not a failure of the search.
The samples span acetate up to 5e-3 mol/L against a half-saturation constant of
0.05, and oxygen up to 2e-3 against 0.01, so every point sits in the linear part
of both Monod terms and the product is the correct answer for the data it was
shown. The saturation shoulder is simply not in the training range.

A longer search does reach it:

```bash
POP=600 GENS=60 ./offline.sh
```

and the list then contains expressions carrying **0.010** and **0.050** inside
them, which are the true oxygen and acetate half-saturation constants, recovered
without the search ever being told the law is Monod. Rearranged into the Monod
normal form, one such expression came out as

```
    found   7.03 x a x o / (0.01005 + o + 0.187 a + ...)
    true    7.00 x a x o / (0.01000 + o + 0.200 a + ...)
```

The automatic pick will not choose it. It takes the steepest gain per node,
which lands on the crude bilinear form, so read the list and use `--pick`.

### Reproducibility

A seed alone was not enough, and the reason is worth stating because it is not
the obvious one.

`fit_symbolic.py` ranks candidate expressions by how well their constants fit.
Those constants come from a least-squares fit, and the last digits of a
least-squares fit are not repeatable: a threaded BLAS adds its terms in whatever
order the threads finish in, and even on one thread a generated expression is
often over-parameterised, so the fit has a whole valley of equally good answers
rather than one. A difference of one part in 10^10 then decided which candidate
survived a tournament, and by the third generation two runs of the same command
were exploring different regions entirely. The same command with the same
`--seed` returned a different rate law.

Four changes fix it, and they are worth knowing about because the same failure
can appear in any fitting code:

1. **The thread pools are pinned to one** before numpy loads, so reductions add
   in a fixed order. `COMPLAB_THREADS` overrides it and says so on stderr.
2. **Each expression gets its own random stream**, derived from the expression
   itself and the run's seed, so fitting an expression no longer depends on when
   during the search it was first seen.
3. **Every comparison is quantised** to six significant digits and every tie is
   broken by the expression itself, so noise below the noise floor cannot
   reorder two candidates.
4. **A tiny penalty on the size of the constants** makes an over-parameterised
   fit well posed: among all the constant vectors that fit equally well it
   prefers the smallest, which is one point rather than a valley.

**What was measured, rather than assumed.** With those four changes, three
separate runs of

```bash
python3 training/fit_symbolic.py --data training/growth_samples.csv \
        --target growth --inputs acetate,o2 --pop 150 --gens 10 --seed 1 --out a.sym
```

produced byte-identical files. The graph-network trainer of example 18 and the
surrogate trainer of example 11 were checked the same way and were already
reproducible once the threads were pinned.

**What is still not guaranteed, and why.** A much longer search, `--pop 600
--gens 60`, did not always agree with itself. The reason is arithmetic rather
than logic: `scipy.optimize.least_squares` is not bit-reproducible on a
rank-deficient problem even on a single thread, a long search makes tens of
thousands of comparisons, and it only takes one of them landing on a rounding
boundary to send two runs down different paths. No environment variable fixes
that, and neither does a fifth change to this script.

So the claim this case makes is a checked one rather than a promised one.
`offline.sh` runs the search twice and compares, and the comparison is on what
the expressions **compute**, not on how they are spelled. That distinction
matters: the search regularly returns the same law written two ways, `647.606 x
acetate x o2` in one run and `acetate x o2 / 0.00154415` in the other, and
calling that a reproducibility failure would be wrong twice over. It is the same
function, and crying wolf teaches a reader to ignore the warning that counts.

Three outcomes, and the check says which one you got:

```
  the two runs agree, expression for expression and character for character.
```

```
  the two runs found the SAME LAW at every complexity: the expressions compute
  identical values on every sample. Some are written differently.
```

```
  PARTLY REPRODUCIBLE.
    same law in both runs at 1, 3, 5 nodes.
    DIFFERENT law at 7, 9, 11, 13, 15, 17, 19, 21 nodes.
```

The third is the common one on a long search, and the pattern in it is not
random: the short expressions repeat and the long ones do not. That is the same
rank-deficiency again. A long expression carries more constants than 400 samples
can pin down, so its fit has no unique answer, so it is exactly the part of the
list that cannot repeat. It is also the part you should not be quoting. The
elbow is the answer, and when the elbow sits in the agreeing set, the expression
you would actually use is reproducible.

Set `VERIFY=0` to skip the check. Every `.sym` the search writes carries the
command, the seed, the thread setting and the library versions in its header,
because a fitted law whose provenance has been lost is not a result, while one
that cannot be re-derived on somebody else's machine still is, provided the file
itself is kept and cited.

**What it deliberately leaves out.** No flow, no abiotic reaction, no
geometry evolution, and no thermodynamic control on the rate: the law is
exactly what the file says, with no pH dependence, no temperature
dependence and no inhibition term, because none was written. Example 05 is
the same colony with the law compiled in, which is faster per evaluation;
example 18 is the same idea with a graph network in place of the
expression; example 19 multiplies a rate law of this shape by a
thermodynamic factor.

## 9. How to build and run

Two files do this, and they are deliberately separate because compiling and
running belong in different places on a cluster.

| File | What it is | Where you run it |
|---|---|---|
| `COMPILE.txt` | The interactive session that builds `./complab`, step by step | Once, by hand, on an interactive node |
| | `run.sh` | The SLURM batch job: modules, pre-processing, the solver, the checks. Submit with `sbatch run.sh`. |
| | `COMPILE.txt` | The interactive session that builds `./complab`, one step at a time. |

**First time in this case folder,** open `COMPILE.txt` and follow it. It asks
for an interactive node, loads the modules, runs cmake with the flags this case
actually needs, and compiles. It also says exactly when you have to come back
and recompile, which for most cases is almost never.

**Every run after that** is one command:

```bash
sbatch run.sh
squeue -u $USER                 # watch it
tail -f output/run.log          # read it while it runs
```

`run.sh` carries a comment on every line: the SLURM header, the module loads
that must match what you compiled with, the pre-processing, the solver call and
the post-processing checks. It runs on one rank on purpose, and the note at the
bottom of the file says why.

To work in a scratch copy instead of dirtying this folder:

```bash
./scripts/setup_case.sh 17_symbolic_law run/mycase
cd run/mycase
```

### On a laptop instead of the cluster


**Build requirements for this case:** a C++11 compiler, CMake 3.5 or newer,
and Palabos v2.3.0. No optional solver is needed. Cases 09, 12, 21 and 22
additionally need GLPK, and case 10 needs Python with COBRApy; this one
needs neither, so the plain cmake line is enough:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DPALABOS_ROOT=/path/to/palabos-v2.3.0
cmake --build build -j
```

**When you must recompile.** This case has no chemistry compiled in: the law
lives in `input/growth.sym`, read at start-up. Editing that file, or
`CompLaB.xml`, or `input/geometry.dat`, never needs a rebuild. A rebuild is
only needed if you change the solver source itself, such as the expression
parser in `src/`.

### The files

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and the seeded biomass patch, reports porosity, refuses to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Optional. Searches for a rate law from a data table and writes a candidate `.sym`. The case ships a law, so nothing has to be fitted before running it. Runs once, before the build. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| | `input/growth.sym` | The rate law of section 1. Read at start-up, echoed into the log, evaluated per voxel. |
| | `input/geometry.dat` | The pore space and the patch. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, runs this case's checks. Standard library only. |
| | `run.sh` | All of the above in order, one comment per line. |
| | `pipeline.sh` | The older, terser version of the same chain. |
| **train** | `training/fit_symbolic.py` | The symbolic regression of section 8. Searches over expression forms and returns a Pareto set. |
| | `training/growth_samples.csv` | 400 sample points with 2 percent noise, so the search has something to work on out of the box. |

Only two things are inherited: the shared kinetics defaults in
[`../../config/kinetics/`](../../config/kinetics/) and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top; this
case does not ship a `kinetics/` header of its own, since its rate law comes
from `input/growth.sym` instead.
