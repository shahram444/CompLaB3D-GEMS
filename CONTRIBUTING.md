# Contributing to CompLB3D

Bug reports, questions and patches are all welcome. This file says where things
live and what a change needs to clear before it goes in.

## Reporting a problem

Open an issue with:

- what you ran: the `CompLaB.xml`, the geometry dimensions, and the cmake line
- what happened, including the terminal output. CompLaB is verbose at start-up
  and that banner usually contains the answer
- what you expected instead

If it is a crash, the smallest input that still reproduces it is worth more
than a large one.

## Before opening a pull request

**Run the example suite.** `examples/` has fifteen cases, one per capability,
each of which runs in seconds:

```bash
python3 examples/validateExamples.py .      # structural check, one second
./examples/runAllExamples.sh .              # build and run all fifteen
```

A change that breaks a case needs either a fix or an explanation in the pull
request of why the expected result changed.

Four cases have answers that do not come from this code, and those are the ones
worth watching: example 02 has an analytic steady state, examples 03 and 13 to
15 must conserve mass, and examples 09 and 10 must agree with each other and
with `min(Vs, 2*Vo)`.

## Where things live

| | |
|---|---|
| `src/complab.cpp` | the driver: parse, set up lattices, time loop |
| `src/complab_functions.hh` | XML parsing and lattice setup |
| `src/complab3d_processors*.hh` | the data processors, one file per concern |
| `src/complab3d_metabolic.hh` | the optional metabolic layer and its config |
| `defineKinetics.hh` | **your** biotic rate laws, compiled in |
| `defineAbioticKinetics.hh` | **your** abiotic rate laws and the dissolution hook |
| `surrogateModel.hh` | the trained surrogate network |
| `surrogate_training/` | how to fit a new one |
| `examples/` | the test suite |

## Things to know before changing the processors

**The argument layout is the API.** Palabos data processors receive a flat
vector of lattices, and the meaning of each slot is positional:
`[C..., B..., dC..., dB..., mask]`. Adding a lattice means updating every
offset that indexes past it. `dCloc`, `dBloc` and `maskLloc` in the processor
headers are those offsets; keep them as named constants rather than open-coding
the arithmetic.

**Rate laws are per second.** A rate constant published per day must be divided
by 86400. That single mistake causes more dead runs than anything else here.

**Never call `exit()` inside a data processor.** Under MPI it leaves the other
ranks hanging and destroys the run. Print a warning and return a safe value.

## Style

Follow the file you are editing. New parsing should degrade the way the
existing code does: an absent optional tag falls back to a documented default,
and a present but malformed one stops the run with a message naming the tag.

## Licence

CompLB3D is AGPL-3.0-or-later. Contributions are accepted under the same terms.
