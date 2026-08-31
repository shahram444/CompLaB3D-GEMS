# Example 08: two microbes competing, two different solvers

Two organisms on opposite walls sharing one donor. Bug1 uses the
cellular automaton and grows fast; Bug2 uses lattice Boltzmann and is
more efficient at low substrate.

## What to expect

- both populations growing, then their **ratio going flat**
- Bug1 ahead, because it sits nearer the inlet where the donor is
  plentiful and its higher mu_max pays off
- a donor shadow between the two colonies

## The switches this case exercises

Per microbe `<solver_type>` and `<reaction_type>`, two `<microbeN>`
material numbers, `bfilm_count` of 2.

## Notes

Swap the two `<solver_type>` lines and rerun. The ratio should barely
move: the solver decides how biomass spreads, not how much of it
there is. If it moves a lot, that is worth understanding before you
trust either solver on a real problem.

## Everything this case needs is in this folder

No cross-referencing, nothing to fetch. The metabolic model, the training
code, the rate-law file and the pore space are all here, beside the
configuration that uses them.

```bash
./scripts/setup_case.sh 08_two_microbes run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space: a staggered slot, 24 × 24 × 6, with the biomass patches seeded. Reports porosity and refuses to continue if it does not percolate. Standard library only. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's balance check. Standard library only. |
| | `pipeline.sh` | Runs the four in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |

### This case's own chemistry

| File | What it is |
|---|---|
| `kinetics/defineKinetics.hh` | This case's own chemistry, compiled in on top of the shared default. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

> **Why the training code is copied here rather than shared.** So that this
> folder *is* the procedure. The cost is real and worth stating: a fix to a
> trainer has to be applied to every case that carries it, and
> `tests/check_repo.sh` fails if a copy drifts from `tools/`.
