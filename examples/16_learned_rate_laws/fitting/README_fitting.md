# What is in this folder

Two worked systems, one per method, plus the data each was fitted from.

## The E. coli system — for `<symbolic>`

    glucose + 6 O2 = 6 CO2 + 6 H2O

| file | what it is |
|---|---|
| `ecoli.sym` | a rate law written by hand, to show the format |
| `ecoli_sweep.csv` | 300 samples from a metabolic-model sweep: glucose, O2, growth |
| `ecoli_discovered.sym` | what `fit_symbolic.py` found in that sweep, plus the substrate lines |
| `symbolic_block.xml` | the XML to paste into CompLaB.xml |

`ecoli_sweep.csv` was generated from a known law,

    growth = 0.8 * glucose/(5e-4 + glucose) * O2/(2e-4 + O2)

so there is a right answer to check against. The search was given the three columns and
nothing else, and returned an expression that rearranges to vmax **0.79969** with
half-saturation constants **4.9976e-4** and **2.0000e-4** — against a true 0.8, 5e-4 and
2e-4, worst deviation **0.03%**. It was never told the form. Reproduce it with:

```bash
python ../../tools/fit_symbolic.py --data ecoli_sweep.csv --target growth \
       --out mine.sym --units per_hour --pop 300 --gens 40 --seed 3
```

It takes a few minutes and prints the whole list of candidate expressions, not just one.

## The AOM system — for `<graphnet>`

    CH4 + SO4(2-) = HS(-) + HCO3(-)

| file | what it is |
|---|---|
| `aom_stoich.csv` | the stoichiometry, one row per species, one column per reaction |
| `aom_samples.csv` | 400 samples: four concentrations, four rates, and growth |
| `aom.gnn` | the trained network |
| `graphnet_block.xml` | the XML, including the abiotic keys |
| `check_aom.py` | holds `aom.gnn` up against the law it was fitted from |

`check_aom.py` is part of `run_tests.sh`. On 300 unseen samples the network tracks its
source law at R = 0.95 to 0.96 per species with the peak rate within 1%.

## The abiotic example

`abiotic.sym` — iron sulfide formation, **Fe(2+) + HS(-) = FeS + H(+)**. Pointed at by
`<abiotic_file>`, which is what makes it sweep every fluid cell rather than only the cells
holding biomass. No organism, so no `growth` line and every variable must be a substrate;
both are refused at load rather than ignored.

## Also here, unused on purpose

`ecoli_stoich.csv` and `ecoli_gnn_samples.csv` are the E. coli system laid out for the
graph network — the same reaction, with a rate column per species. There is deliberately
**no** trained `ecoli.gnn` shipped: the one I first fitted reached only R = 0.69 on CO2 and
put the O2:glucose ratio at 5.67 instead of 6.00, which is not an example worth copying.
Train your own and watch that ratio; it is the honest measure of whether it converged:

```bash
python ../../tools/train_graphnet.py --stoich ecoli_stoich.csv \
       --data ecoli_gnn_samples.csv --out ecoli.gnn \
       --width 6 --rounds 2 --epochs 600 --lr 0.03 --units per_hour --da 1.0
```
