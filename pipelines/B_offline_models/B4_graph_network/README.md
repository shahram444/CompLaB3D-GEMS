# B4 — Training the graph network

**In:** a stoichiometric matrix, a sample table, and one Damkohler number per
reaction.
**Out:** a `.gnn` text file the solver reads at start-up.
**Run:** `./run.sh`, which calls `tools/train_graphnet.py`.
**Needed for:** the graph-network path only.

## The three inputs

**`stoich.csv`** — one row per species, one column per reaction. This is the
same matrix your `CompLaB.xml` already describes, handed to the trainer as
structure rather than left to be inferred. It is what makes the network small:
457 parameters against 1280 for a dense network on the same four-species system.

```csv
-1,0
-4.8,0
4.8,-1
0,-1
```

**`samples.csv`** — one concentration column and one rate column per species,
plus a growth column if the model has one.

**The Damkohler numbers** — one per reaction, entering at the edge level so that
a single trained network spans a range of transport regimes rather than being
tied to the one it was fitted in.

## Running it

```bash
python ../../../tools/train_graphnet.py \
       --stoich stoich.csv --samples samples.csv \
       --da 1.0,1.0 --rounds 1 --width 8 --epochs 4000 \
       --out network.gnn
```

`--rounds` is worth understanding. One round lets each species hear from the
reactions it touches. Two rounds let it hear from species two edges away, which
matters when a species influences another only through an intermediate.

## Check what came out

```bash
python ../../../tests/xval_gnn.py --net network.gnn --samples samples.csv
```

The check that matters is not the correlation, it is the **sign and the ratio**.
A network that predicts the right magnitude but occasionally produces a reactant
instead of consuming it has violated the stoichiometry it was built from. On the
shipped example the median ratio of the Fe(III) rate to the acetate rate comes
out at 4.772 against a true 4.800, within 0.6 per cent — and that ratio was
never a training target. It comes out right because the stoichiometry is
structure.

## The file carries its own range

`trainmin` and `trainmax` are written into the `.gnn` file by the trainer, from
the observed range of the samples. They travel with the model rather than beside
it, and they are enforced at every evaluation and counted. At the end of a run
the solver reports how many evaluations were clamped, and says so plainly when
the fraction is large.
