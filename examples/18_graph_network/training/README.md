# Training the graph network

**In:** a reaction, and a table of concentrations with the rates they produced.
**Out:** `input/aom_retrained.gnn`, a text file the solver reads at start-up.
**Run:** `./offline.sh` from the case directory.
**Needed for:** this case only if you want your own network. The case ships
`input/aom.gnn` and runs without ever coming here.

---

## What is in this directory

| | |
|---|---|
| `aom_stoich.csv` | the reaction: one line per species, its stoichiometric number |
| `aom_samples.csv` | 400 rows: four concentrations, the four species rates they produced, and the growth rate |
| `make_training_data.py` | writes a table like that one from a law you can read |
| `train_graphnet.py` | fits a network to such a table |
| `xval_gnn.py` | holds a trained network up against data it was not fitted to |

---

## The reaction

Anaerobic oxidation of methane coupled to sulfate reduction:

```
    CH4  +  SO4(2-)   ->   HS(-)  +  HCO3(-)
```

One turn consumes one methane and one sulfate and makes one sulfide and one
bicarbonate. That is exactly what `aom_stoich.csv` says, and the network is
built on it: it predicts **one number per reaction**, how fast that reaction
turns, and every species rate is its stoichiometric number times that number.

That is why this case reports a stoichiometric residual of exactly zero. The
network is never asked to predict four rates and hope they stay in proportion.
It predicts the turn, and the proportions are arithmetic.

---

## Where the shipped table came from

`aom_samples.csv` was drawn from a law written down in `make_training_data.py`:

```
    extent = 4.0e-3 * CH4/(1.0e-3 + CH4) * SO4/(5.0e-4 + SO4)      mol/L/h
    growth = 5.0 * extent
```

So this case has an answer key. `tests/check_aom.py` holds the trained network
up against that same law on 300 samples it never saw, which is a stronger claim
than a training error.

**Two things about the shipped table worth knowing before you read any accuracy
figure.**

It carries **no noise at all**. Whatever error the network shows is its own
approximation error, not scatter it could never have fitted. Pass `--noise` to
`make_training_data.py` if you want to see how it copes with a realistic table.

Its **sampling reaches past both half-saturation constants**: methane to five
times its constant, sulfate to sixteen times its. That matters more than it
sounds. Below its half-saturation constant the law is a straight line in that
variable, and every value of the constant fits a straight line equally well, so
a table that stops short cannot identify the constant no matter how many rows
it has. If you build your own table, check this first.

The two products, sulfide and bicarbonate, are sampled but do not enter the law.
They are there so the network has to learn that they do not matter, rather than
being told.

---

## Making your own table

```bash
python3 make_training_data.py                        # an equivalent table
python3 make_training_data.py --n 2000 --noise 0.02  # more rows, with noise
python3 make_training_data.py -o wider.csv --ch4-hi 2e-2
```

It writes a companion `<name>_truth.json` recording the law, the constants, the
noise and the seed, so the answer key cannot drift away from the data it goes
with.

**Regenerating means retraining.** `aom_samples.csv` and the shipped
`aom.gnn` were made as a pair. Overwrite the table and the network no longer
belongs with it, so run `./offline.sh` afterwards. The script defaults to a
different filename to make that hard to do by accident.

---

## Training

```bash
./offline.sh
```

which runs

```bash
python3 training/train_graphnet.py \
        --stoich training/aom_stoich.csv \
        --data   training/aom_samples.csv \
        --da 1.0 --rounds 2 --width 6 --epochs 250 --units per_hour \
        --out input/aom_retrained.gnn
```

`EPOCHS=800 ./offline.sh` if you want a longer fit.

| flag | what it does |
|---|---|
| `--rounds` | how many times a species passes a message to its reaction and back |
| `--width` | how many numbers each message carries |
| `--da` | one Damkohler number per reaction, which sets the scale the rate is measured in |
| `--units` | `per_hour` or `per_second`; the solver works in seconds and the file says which it holds |
| `--readout` | `extent` predicts the turn of the reaction, `species` predicts each rate separately and gives up the exact stoichiometry |

---

## Reading what comes out

The trainer writes its own provenance into the file:

```bash
grep '^provenance' input/aom_retrained.gnn
```

Compare it with the shipped network's:

```bash
grep '^provenance' input/aom.gnn
```

Two numbers matter. **R** is how much of the variation the network accounts for,
and should sit very close to 1. **Scaled mse** is the error in the units of the
rate itself. A fit that looks good on R and poor on scaled mse is usually a fit
that got the large rates right and the small ones wrong, and the small ones are
where a substrate is running out, which is the half of the domain that decides
where the reaction stops.

Then check it against data it never saw:

```bash
python3 training/xval_gnn.py
python3 ../../tests/check_aom.py
```

---

## Installing it

The retrained file is **not** installed. `CompLaB.xml` still points at
`input/aom.gnn`. Copy it over only once you are satisfied:

```bash
cp input/aom_retrained.gnn input/aom.gnn
```

---

## The one that catches people

**The species names in the file must match `<name_of_substrates>` exactly.** The
names are the binding between the network and the simulation, and a mismatch
stops the run at start-up naming the offending species. That is deliberate: the
alternative is binding to the wrong lattice and producing a plausible, wrong
answer.
