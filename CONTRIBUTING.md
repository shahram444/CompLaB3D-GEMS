# Contributing

## The rule the tree is built on

**Outside `examples/`, nothing exists twice.** One reference configuration, one
copy of each shared kinetics header, one copy of each tool.
`tests/check_repo.sh` fails on any two identical files over 200 bytes outside
`examples/`.

**Inside `examples/`, duplication is deliberate.** A case folder holds
everything that case needs — its pore space, its metabolic model, the training
code that produced what it was trained on — so that the folder is the whole
procedure rather than a set of pointers into the rest of the tree. Several cases
therefore carry the same geometry, and three carry the same model exporter.

The cost of that is real and there is no way around it: **a fix to a shared tool
has to be applied to every example that carries a copy.** `check_repo.sh`
enforces this rather than hoping — it compares every copy under `examples/*/`
against its original and fails if one has drifted. The list of pairs it checks
is in the script; add a line to it whenever you copy a new tool into a case.

So the order when you fix a tool is:

1. fix it in `tools/` (or `tests/`),
2. copy it over every `examples/*/training/` that carries it,
3. run `./tests/check_repo.sh` — the drift check tells you if you missed one.

## Before you add an XML tag

Two files in `config/` document the configuration, and both have to gain a line:

- `config/CompLaB.everything.xml` — one line for every ability, with its unit
  and its default. This is the index, and `tests/check_repo.sh` fails if the
  solver reads a tag this file does not list, so you will find out immediately.
- `config/CompLaB.reference.xml` — the same tag with the reasoning: why it
  exists, and what goes wrong if someone sets it wrong.

A tag documented in neither is a feature nobody can find.

## Before you change the chemistry interface

If you change the signature of `defineKinetics` or `defineAbioticKinetics`, you
have to update:

- `config/kinetics/defineKinetics.default.hh`
- `config/kinetics/defineAbioticKinetics.default.hh`
- `config/kinetics/defineKinetics.biotic.hh` — shared by examples 05–08
- the six files under `examples/*/kinetics/` that override them

and nothing else. That is the point of the arrangement — before it, the same
change meant editing thirty-two files, twenty-three of which were identical.

## Adding an example

An example carries its own `preprocess.py`, `postprocess.py` and `pipeline.sh`,
plus an `offline.sh` if it needs offline work. Write those four **for that
case** rather than copying a neighbour's — a generic script repeated twenty
times teaches nothing, and the whole point is that reading one folder tells you
what that case does. They use the standard library only, so the folder is
readable without following an import.

Shared *tools* are different: copy those in verbatim, into `training/`, and add
the pair to the drift check in `check_repo.sh`.

`check_repo.sh` also fails if an example reaches back into `../../tools` or
`../../models`. If a case needs something, it gets its own copy.

## The geometry an example runs on

Each example draws its own, in its own `preprocess.py`, and that generator is the
authority: `check_repo.sh` fails if a shipped `input/geometry.dat` is not what
its generator writes. The four files in `config/geometry/` are unpadded starting
points, not shared inputs — no example reads them. If you add one there, document
it in `config/geometry/README.md` alongside the others.

Whatever shape you draw, **add** the inert wall layer on the four faces the
solver does not condition rather than converting pore into wall, so porosity and
every voxel count stay what the case declares.

## The README an example ships

Every case README has the same shape, described in
[`examples/README.md`](examples/README.md): what the case is for, a **What is
simulated** table, the physics with every constant carrying its unit and its
source, **What to check**, and what the case does not do. Generate the table from
the case's own `CompLaB.xml` and `input/geometry.dat` rather than typing it —
a hand-written one is wrong the first time a switch changes and nothing
notices.

## Before you change a configuration tag

`config/CompLaB.reference.xml` is the only place a tag is documented. Add the
tag there, with its units and its default, in the same style as its neighbours.
Do not copy the documentation into an example: the examples are working subsets,
not references.

## Tests

```bash
./tests/check_repo.sh     # the tree — fast, no dependencies
./tests/run_tests.sh      # the solver
```

Run `check_repo.sh` after any change to the tree: it fails on a duplicated file,
a `case.files` line pointing at nothing, an example that no longer assembles, a
shell script that does not parse, or a dead link in a README.

Everything must pass. New behaviour needs a test that fails without it. The
suite runs against a Palabos stub rather than the real library, so it is fast
and has no external dependency — keep it that way.

Two properties the suite exists to protect, which are easy to break and hard to
notice:

- **the mass budget.** No species may be drawn below zero, at any step length.
  A fitted rate law is free to be wrong; it is not free to create negative mass.
- **trainer and solver agree.** A model is written by Python and evaluated by
  C++, and nothing about a text file by itself stops the two drifting apart.
  `tests/xval_gnn.py` compares them across the training box and outside it.

## Style

C++ follows the surrounding file. Python is PEP 8 with four-space indents.
Comments explain why, not what.

## Reporting a problem

Open an issue with the `CompLaB.xml`, the geometry dimensions, the first twenty
lines of the log, and what you expected. The start-up lines say what the solver
decided, which is usually enough to see the problem.

## Publishing

[`docs/PUBLISHING.md`](docs/PUBLISHING.md) has the whole path from a folder on
disk to a tagged, DOI-bearing release, and the day-to-day loop after that.
