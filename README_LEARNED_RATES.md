# CompLaB3D, Python workflow, with learned rate laws

This is version 2's lineage, unchanged in how you work with it: `extractMM.py` still builds
the metabolic model file, COBRApy still runs through the embedded interpreter, and the fitting
tools are Python scripts you run yourself. Nothing was moved into C++ that was not already
there. Version 3 is the one where that machinery is native.

What is added here is two more ways to supply a reaction rate.

Two ways to get a rate law into a simulation without recompiling it, added as
separate reaction types with their own XML switches. Neither replaces anything;
both sit beside `defineKinetics.hh`, GLPK, COBRApy and the surrogate, and
nothing changes for a run that does not ask for them.

| | `<symbolic>` | `<graphnet>` |
|---|---|---|
| what it holds | a formula | a fitted network |
| file | `.sym`, human-writable | `.gnn`, produced by the trainer |
| output | the rates you write down | one rate per species, plus growth |
| structure | whatever you wrote | the stoichiometric matrix, as a graph |
| use it when | symbolic regression handed you an expression, or you want to compare four candidate mechanisms | you need every species' rate and you already know the reaction network |
| `<reaction_type>` | `symbolic` (8), `symbolic_and_kinetics` (9) | `graphnet` (10), `graphnet_and_kinetics` (11) |

They are two answers to the same question — *how does a rate law get from a
fitting run into a simulation without a rebuild?* — and they are deliberately
independent. You can use either, both, or neither.

## What is in here

```
src/complab3d_symbolic.hh                NEW  the expression parser and evaluator
src/complab3d_processors_symbolic.hh     NEW  its per-voxel driver
src/complab3d_graphnet.hh                NEW  the .gnn loader and forward pass
src/complab3d_processors_graphnet.hh     NEW  its per-voxel driver
src/complab3d_processors_abiotic_learned.hh  NEW  both of the above, swept abiotically
src/complab3d_metabolic.hh           PATCHED  four new reaction types, 8..11
src/complab.cpp                      PATCHED  the wiring for both

tools/fit_symbolic.py                    discovers a rate law and writes the .sym
tools/train_graphnet.py                  fits a network and writes the .gnn
examples/ecoli.sym                       a worked expression file, written by hand
examples/ecoli_sweep.csv                 a metabolic-model sweep to fit against
examples/ecoli_discovered.sym            what fit_symbolic.py found in that sweep
examples/aom.gnn                         a worked network: AOM in the SMTZ
examples/aom_stoich.csv, aom_samples.csv what it was fitted from
examples/symbolic_block.xml              the XML for <symbolic>
examples/graphnet_block.xml              the XML for <graphnet>, and for the abiotic sweep
examples/abiotic.sym                     a worked abiotic law
tests/                                   see "Testing" below
run_tests.sh                             runs every check
```

Copy `src/` over version 2's `src/`. The build system needs no edit: all four
new headers are header-only, and neither `complab3d_symbolic.hh` nor
`complab3d_graphnet.hh` depends on Palabos at all.

## `<symbolic>` — a formula in a text file

```
# fitted by PySR from the e_coli_core sweep, hold-out R2 = 0.9987
units  per_hour
vars   glucose o2 Ecoli

rate   growth  = 0.0412 * glucose * o2 / (0.083 + o2)
rate   acetate = 0.31 * growth
rate   glucose = -1.0 * growth / 0.0412

range  glucose 0.1 10.0
range  o2      0.1 20.0
```

A `rate` names either a substrate or the reserved word `growth`. Later rates may
use earlier ones, so a stoichiometric ratio is one line. A variable is either a
substrate name or the organism's own biomass.

### Where the expression comes from

You can write it yourself if you already have kinetic constants. Otherwise
`tools/fit_symbolic.py` will find one for you, and it is worth being clear about
what that means: it does **not** fit constants into a form you chose. It searches
over the form itself, building expressions out of your variables and `+ - * /`,
and hands back the formula. Nobody tells it about Monod; if the data is Monod it
will find something Monod-shaped, and if it is not, it will find something else.

```bash
python tools/fit_symbolic.py --data examples/ecoli_sweep.csv \
       --target growth --out input/rates.sym --units per_hour
```

The input is a table: run the mechanistic model over a grid of concentrations,
record the growth rate, save it as CSV. The output is not one answer but a short
list, one expression per level of complexity, with its typical and worst-case
error. **Read that list yourself.** The tool marks the steepest gain per node as
a default, but that is a rule of thumb and often not the one you want; the
shortest expression whose accuracy you can live with usually is. `--pick <nodes>`
takes a different one.

Two things it deliberately does not do. It fits **one** column, so fit `growth`
and then write the substrate lines by hand as multiples of it — that keeps them
in exact stoichiometric ratio, which four independent fits would not. And it
judges the fit on **relative** error by default, because a rate law is used
across orders of magnitude and plain least squares is excellent at the top of the
range while being wrong by hundreds of percent at the bottom — which is exactly
where a reaction front sits. `--loss absolute` switches that off if you really
want it.

For production work PySR is the better tool: faster, searches harder, and
published. `fit_symbolic.py` exists so the loop closes without a Julia install,
and so a `.sym` file can come out of a table in one command.

### Does it really find the form?

`examples/ecoli_sweep.csv` is 300 samples generated from

    growth = 0.8 * glucose/(5e-4 + glucose) * O2/(2e-4 + O2)

so there is a right answer to compare against. The search was given the three
columns and nothing else. What it returned, at 17 nodes:

```
(O2 / ((0.0002 + O2) / 0.151699)) * (glucose / (((0.00557898/11.1633) + glucose) / 5.27156))
```

Rearranged, the first bracket is `0.151699 * O2/(2e-4 + O2)` and the second is
`5.27156 * glucose/(4.9976e-4 + glucose)`. The two prefactors multiply to
**0.79969**, against a true 0.8; the half-saturation constants come out
**4.9976e-4** and **2.0000e-4**, against 5e-4 and 2e-4. Worst deviation over the
300 samples: **0.03%**.

It was never told about Monod. It was given numbers and found the law.

`examples/ecoli_discovered.sym` is that result with the three substrate lines
added by hand from the stoichiometry, which is the workflow in miniature.

```xml
<symbolic>
    <enabled>true</enabled>
    <expressions_file>input/rates.sym</expressions_file>
</symbolic>
```

## `<graphnet>` — a network that already knows the chemistry

A dense network from concentrations to rates has to learn from data which
species affect which. That information is not missing: it is the stoichiometric
matrix. Supplying it as *structure* rather than making the network infer it
means fewer weights, less training data, and a model that cannot invent a
coupling the chemistry does not have.

The graph is bipartite — one node per species, one per reaction, an edge
wherever `s[i][r]` is non-zero carrying that coefficient — and one round of
message passing is

```
reactions   hR_r = tanh( Wsr . SUM_i s[i][r] hS_i  +  Wda da_r  +  bR )
species     hS_i = tanh( Wss . hS_i  +  Wrs . SUM_r s[i][r] hR_r  +  bS )
```

The sums are stoichiometry-weighted, so a species consuming two moles per
reaction speaks twice as loudly into that reaction as one consuming a single
mole. `da_r` is a per-reaction Damköhler number carried on the edges, which is
what lets one trained network cover a range of transport regimes instead of one
— the individual-Da idea from the workplan, put where it belongs.

Fitting it:

```bash
python tools/train_graphnet.py \
    --stoich examples/aom_stoich.csv \
    --data   examples/aom_samples.csv \
    --out    input/aom.gnn \
    --width 6 --rounds 2 --epochs 250 --units per_hour --da 1.0
```

`--stoich` is a CSV, one row per species, one column per reaction. `--data` is a
CSV with a header: one column per species concentration, then one `<species>_rate`
column each, then `growth` if there is one. numpy is the only dependency — no
PyTorch — for the same reason the `.srg` trainer has none.

```xml
<graphnet>
    <enabled>true</enabled>
    <network_file>input/aom.gnn</network_file>
</graphnet>
```

The network's `species` line must match the simulation's substrate names
exactly. A name that does not resolve stops the run and says which name, the
same choice `<exchange_reaction_names>` made and for the same reason: an index
that silently points at the wrong species produces a well-posed simulation with
a wrong answer in every voxel.

## The worked example: AOM in the sulfate–methane transition zone

`examples/aom.gnn` is a real fit, not a hand-written illustration. The reaction
is the one the SMTZ is defined by,

```
CH4 + SO4(2-)  ->  HS(-) + HCO3(-) + H2O
```

and `aom_samples.csv` was generated from a dual-Monod law, so the network can be
checked against the law it came from. The command that produced it is in the
file's own `provenance` lines, which is the point of having them.

## Biotic and abiotic are two different sweeps

`defineKinetics.hh` and `defineAbioticKinetics.hh` are separate for a reason, and so are these.

A law bound to an organism, through `<expressions_file>` or `<network_file>`, is run by a sweep
that skips any cell holding no biomass. That is right for a microbe: a cell that is not there
cannot eat. It is wrong for abiotic chemistry, which happens in open water whether anything is
alive nearby or not — put an abiotic reaction on that path and it would fire inside the biofilm
and stop dead at its edge, with no warning.

So an abiotic law is pointed at separately and swept over every fluid cell:

```xml
<symbolic>  <abiotic_file>input/abiotic.sym</abiotic_file>  </symbolic>
<graphnet>  <abiotic_file>input/abiotic.gnn</abiotic_file>  </graphnet>
```

It belongs to no organism, so two things follow, and both are refused at load rather than quietly
ignored: **every variable must be a substrate** (a biomass name means the file was meant for the
other key) and **there is no growth rate** (a `growth` line, or a network with a growth output,
is a stop). Increments land in the same `dC[]` the hand-written abiotic kinetics uses, so the two
can run side by side and are applied together in one pass.

Between the four keys, both kinetics headers are now replaceable — biotic and abiotic — either
one at a time or together.

## Two things both paths share, deliberately

**The training box is enforced, not advised.** A formula fitted between 0.1 and
10 mM is no more valid at 50 mM than a neural network is, and neither one fails
out there — both return a confident number. So both carry their valid range in
the file, both clamp to it, both count the clamping, and the end of a run prints
how often it happened with a warning above 5%. That count is the only thing
standing between an extrapolated run and a result.

**Production is not clamped away.** `runFBA_cobrapy3D` and `run_surrogate3D` both
end their budget loop with `if (draw > T()) draw = T();`, which discards any net
release — so an organism on those two paths can consume but cannot excrete.
`runFBA_glpk3D` has no such line. Both new processors follow GLPK: the draw is
limited by what the voxel holds, and production passes through. A rate law that
makes a product has to be able to write it.

(That difference is also what the GLPK-versus-COBRApy comparison in the GMD
paper turned up, and it is worth fixing on those two paths as well. Not done
here — it changes the behaviour of existing runs, so it should be its own
change.)

## The two loops, end to end

Both methods now have the same shape, and that is the point:

| | symbolic | graph network |
|---|---|---|
| 1. generate data | run the mechanistic model over a grid | same |
| 2. fit | `tools/fit_symbolic.py` | `tools/train_graphnet.py` |
| 3. result | `.sym`, a formula you can read | `.gnn`, a fitted network |
| 4. run | `<symbolic>` in the XML | `<graphnet>` in the XML |
| 5. check | the stoichiometric ratios in the output | same |

Neither step 2 happens inside the solver. Fitting is slow and seed-dependent, and
neither belongs in a simulation. Fit offline, write the file, run.

## Testing

```bash
./run_tests.sh
```

All eight pass. None needs Palabos — only g++ and python3 with numpy.

- `test_sym`, `test_file` — the expression language: precedence, right-associative
  `^`, unary minus, every function, the four cases a rate law must survive without
  producing an infinity (divide by zero, `log(0)`, `sqrt(-1)`, fractional powers of
  a negative), and eight malformed expressions each checked to be refused.
- `test_gnn` — a network small enough to work out on paper, checked against the
  message-passing rule written as a single expression; the clamp and its counter;
  the `per_hour` factor of 3600; and seven malformed files each refused rather than
  half-loaded.
- `test_abiotic` — the abiotic sweep, and above all the one thing it exists for: a reaction
  firing in a cell with **no biomass anywhere in the lattice vector**. Plus the same budget,
  deposit, mask and `useTotals` checks the biotic tests make, so the two sweeps cannot drift.
- `test_symbolic_processor`, `test_graphnet_processor` — the processors **run**,
  on a stub lattice (`tests/palabos_stub.hh`, used by nothing else). These check
  the five things that would otherwise be wrong in a way a simulation would not
  announce: that a rate reaches the substrate it names, that the budget cannot
  take a concentration below zero, that a product *is* written, that growth scales
  with biomass, that the D3Q7 deposit's seven weights sum to the whole, and that
  `useTotals` charges a draw once rather than twice.
- `xval_gnn.py` — the C++ forward pass against the Python one, over the same
  network and the same inputs, half of them deliberately outside the training box.
  They agree to about 1e-15 of the largest output, so the file format is not a
  convention anyone has to maintain by hand.
- `examples/check_aom.py` — the shipped example network against the dual-Monod law
  it was fitted from, on 300 samples it never saw: R = 0.95 to 0.96 per species,
  peak rate within 1%, and the signs right in 298 of 300 (the two are rates near
  zero, where the sign carries no information).

## What is still not tested

`complab.cpp` and `complab3d_metabolic.hh` have not been compiled, because that
needs Palabos. Both processors now have been, against the real
`complab3d_metabolic.hh`, so `MetabolicConfig` is genuinely checked. What was
verified on `complab.cpp` instead: braces balance, the constructor takes nine
parameters and the call passes nine for both processors, every new variable is
declared at an enclosing scope, and the lattice-vector length assertion covers
both new pointer lists. Treat the first `make` as part of the job.

## First thing to run

For `<symbolic>`: a closed box, one organism, `rate substrate = -k*substrate`,
checked against `A(t) = A0*exp(-k*t)`. That exercises the parser, the binding,
the budget and the deposit in one case whose answer is known in closed form.

For `<graphnet>`: the same box with `examples/aom.gnn` and a large excess of
sulfate, so the rate collapses to first order in methane and the same closed form
applies.
