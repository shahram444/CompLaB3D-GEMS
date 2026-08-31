# The CompLB3D User Guide

`CompLB3D_User_Guide.pdf` is built from the LaTeX sources here.

```bash
make            # produces CompLB3D_User_Guide.pdf
```

## The structure

| file | contents |
|---|---|
| `main.tex` | preamble, styling, title page, "how to read this guide" |
| `part1.tex` | Getting started: what it does, installing, first run, the input file |
| `part2.tex` | Tutorial: one chapter per capability, built on the 15 shipped examples |
| `part3.tex` | Going further: your own kinetics, surrogate training, clusters, output |
| `part4.tex` | Reference: every XML tag, every error message, known limits, file map |

## Keeping it honest

The reference chapter and the examples describe the same code, so they drift
apart the moment one changes without the other. Two habits keep that from
happening:

- when you add or rename an XML tag, update Part IV in the same commit
- when you add an example, add its chapter to Part II and its row to the table
  in `examples/README.md`

`examples/validateExamples.py` extracts the tag whitelist mechanically from the
source, so a renamed tag surfaces there as an unknown tag. That catches the
code-side half of the drift; the manual-side half is on you.

## Things deliberately not in the manual

- **Simulation results.** Nothing has been run against real Palabos yet, so
  there are no result figures. Add them after the first Tahoma runs; the
  diagrams are all schematic and none of them will need changing.
- **A derivation of the lattice-Boltzmann method.** The Palabos user guide does
  that better, and this guide links to it.
