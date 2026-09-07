# Method guides

One document per rate path, plus one for the thermodynamic gate that multiplies
all of them and one for the two mineral directions that edit the geometry
underneath. Each is written to be read straight through by someone who has never
opened the source, and each ends with a worked example that puts real numbers
into every equation, one substitution at a time.

| Guide | Covers | The example it works through |
|---|---|---|
| `Method_1_FBA_guide.docx` | the linear program, in process through GLPK and across the boundary through COBRApy, and its two refinements `<multi_step>` and `<cybernetic>` | one voxel of example 09, solved by hand: bounds, program, increments |
| `Method_2_Surrogate_guide.docx` | the network fitted to the program's answers | one voxel of example 11, one neuron of layer 1 written out in full |
| `Method_3_Symbolic_guide.docx` | rate laws read from a `.sym` file at start-up | one voxel of example 17, from the clamp to the increments |
| `Method_4_GraphNetwork_guide.docx` | the bipartite graph network over a `.gnn` file | a three-species network small enough to check with a pencil |
| `Method_5_Thermodynamics_guide.docx` | the `F_T` gate and the `.thm` file | two voxels of the same aggregate, differing in one concentration |
| `Appendix_D-precipitation_and_dissolution_NEW.docx` | the two mineral directions and the geometry they edit | one voxel filling, one voxel being eaten away |

Each guide's Figure 2 is a flowchart of that method end to end: what is done offline
before any simulation exists, what the solver does once at start-up, and what runs in
every voxel at every step. The boxes are lettered and the caption explains each letter.
Appendix D numbers its figures D1 to D11 rather than 1 to *n*, and its Figure D1 is
the same kind of end-to-end chart for the two mineral directions.

[v1.3] Method 1 carries two more of these, one per refinement of the linear
program that v1.3 added. **Figure 11** draws multi-step flux balance analysis
(`<multi_step>`) end to end, from choosing the stage order offline to the one
implementation detail that matters: growth is read from the last stage's flux
vector, never from the first stage's objective value. **Figure 14** draws
cybernetic switching (`<cybernetic>`), including why shutting a source must
forbid its uptake and never its release. Both captions walk every lettered box,
and both are followed by the XML block that turns the feature on, quoting what
examples 21 and 22 actually ship.

Methods 3, 4 and 5 have a matching `-figures.pptx`. Where one exists, every
figure in that document is a render of one slide of the deck, so the editable
original of a figure is the slide and a figure is changed by editing the deck and
re-exporting rather than by redrawing.

## Where these sit relative to the rest of the documentation

The manual in the directory above (`docs/main.tex` and its parts) is the
reference: it covers every tag, every solver, and the numerical scheme.
A reader who wants to know what `<readout>extent</readout>` does looks it up in
the manual; a reader who wants to know why the extent readout exists, and what
it buys, reads Method 4. All six are guides rather than reference, and they
overlap the manual on purpose.

[v1.3] Rate paths 1 and 2 now have their guides here too; earlier text said they
lived with the paper rather than in the repository.

## What is in each one

**Method 1, the linear program.** What a genome-scale model is as a matrix, and
the three statements the program makes about it, one of which is local and is the
only place the pore-scale condition enters; where the solve sits in the loop; how
GLPK runs it in process and how COBRApy runs it across a Python boundary; then
the two cases where one program is not enough, the chain (`<multi_step>`) and the
competition (`<cybernetic>`), each with its own end-to-end figure and its
configuration block. The worked example is one voxel of a three-metabolite,
four-reaction model small enough to solve by hand, carried from concentrations to
uptake bounds to fluxes to the increments the lattice receives, and then repeated
in the same voxel once the donor has run down.

**Method 2, the surrogate network.** What the network is fed and why the inputs
are the uptake bounds rather than the concentrations; the rescaling and what it
records; the network layer by layer, and the one place the growth-only and
multi-output forms differ; how the training points are placed and what the fit
minimises; the range check and what it costs when it bites. The worked example
writes out one neuron of the first hidden layer in full, then the remaining three
layers and the output, and ends at the increments the transport solver receives.

**Method 3, symbolic rate laws.** What symbolic regression returns and why it is
a list rather than an answer; the expression language; the tree the solver walks;
the `.sym` file line by line; the two guard rails, units and range, and why both
failures are silent. The worked example carries 1.2 mM acetate and 0.6 mM oxygen
through the tree, the unit scale and the substrate budget, then repeats it in a
voxel where the range clamp bites and costs a factor of 1.38 in growth rate.

**Method 4, the graph network.** The reaction network drawn as a graph beside
the matrix it is; one round of message passing, drawn and then written; where the
fitted numbers actually are and why they do not belong to any one species; the
two readouts, and the single number that separates them. The worked example runs
three species, one reaction, two numbers per node and one round, so every value
in it can be checked by hand, and it shows the same voxel giving a sulfide per
methane ratio of exactly 1 through one readout and 1.120 through the other.

**Method 5, the thermodynamic gate.** Where the factor lives, which is the
question the rest of the document depends on: one value per biomass-bearing
voxel per reaction step, from that voxel's own concentrations, discarded and
recomputed at the next step. Then the four equations, the energetics file, the
optional growth yield, and what the gate does not do. The worked example takes
two voxels of one aggregate that differ only in bicarbonate, 9.6 against 20.0 mM,
and shows the factor falling from 0.547 to 0.056 on that alone.

**Appendix D, precipitation and dissolution.** The two mineral directions as one
mechanism run in opposite signs: the surface-controlled rate that acts only in
the wetted rows, the fill and reopen thresholds, what happens to a voxel at the
moment it becomes solid or becomes pore again, and the flow feedback that follows
from it. Its figures are numbered D1 to D11 and step through one pore
cross-section as it seals, then through a calcite-lined pore as it opens.
