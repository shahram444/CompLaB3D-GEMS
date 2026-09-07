# A — Pre-processing

Two things may need making before a run. The first always does.

| | What it makes | Needed |
|---|---|---|
| [A1](A1_geometry/) | `geometry.dat`, the pore space | always |
| [A2](A2_initial_fields/) | starting concentrations and biomass | only when you want something other than uniform |

Nothing here knows anything about chemistry. Both stages produce plain arrays
the solver reads at start-up, and both are the same whichever rate path you go
on to use.
