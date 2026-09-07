# Changelog

## v1.3: a correctness pass over the whole tree

No new capability in this pass. What follows is the result of reading every
source file and every example against what it claims, and it exists because
three of the things that pass looked exactly like the things that fail: a run
that finishes, fields that look smooth, and a number that is wrong. The release
also carries the two refinements of the linear program described in the next
section, `<multi_step>` and `<cybernetic>`, with examples 21 and 22; those are
capability, and they are kept separate here so that the correctness pass can be
read on its own.

**The one that changes a published conclusion.** `classicalEta()` was the sphere
solution written for the generalized Thiele modulus, `(R/3)sqrt(k/D)`, while
`complab.cpp` computed and every README printed the modulus on the radius,
`R sqrt(k/D)`. Feeding one into the other evaluates the curve at three times the
right argument. Every `eta_classical` the solver has ever reported was `eta(3
phi)`, about 0.15 too low over the range example 20 covers, and the headline
finding of that case inverts: the measured effectiveness factor sits **below**
the classical curve, at 0.823 against 0.947, not above it at 0.823 against 0.700.
Below is also the direction the physics says: the gate closes the core, and
`offline/upscale.py`'s own header comment predicted exactly that while the code
contradicted it.

The formula now lives in `src/complab3d_thiele.hh`, which states the convention
three times, and `tests/test_upscale.cpp` pins it three ways a single shared
mistake cannot satisfy at once: against the closed form written out
independently, against the small-phi series (`1 - phi^2/15` here against
`1 - 0.6 phi^2` for the other convention, a factor of nine in the first
correction term), and against the diffusion-limited asymptote (`3/phi` against
`1/phi`). `offline/expected_sweep.csv` and both READMEs carry the recomputed
numbers.

### Fixed

- **A voxel that seals or reopens now stops or starts conducting solute.** Every
  substrate lattice took its dynamics at start-up from the static `geometry`
  field, which precipitation and dissolution never touch: they change the pore
  space by writing the mask. So a fully clogged throat went on diffusing at the
  full pore diffusivity through a voxel the flow solver had already turned into
  a wall, and a voxel the mineral had vacated kept its bounce-back, so acid
  could never reach the fresh surface. `updateSoluteDynamics3D` could not do this
  job: it skips `mask == solid` outright and only ever calls `setOmega`, which is
  a no-op on bounce-back. New `updateSoluteSolidDynamics3D`, called from both
  conversion blocks and only when the mask actually moved. Examples 13, 14 and 15
  were all affected; 13 seals 114 voxels in its shipped configuration.
- **The CA and FD mask update no longer walks past the end of its lattice list.**
  `updateLocalMaskNtotalLattices3D` loops over the `bio` vector it is handed and
  reads `lattices[iM]` for each row. It was handed `bio_dynamics`, which has one
  row per BIOFILM microbe, while `ptr_ca_lattices` and `ptr_fd_mask` carry only
  the microbes on that solver. Shipped example 08 is exactly such a run, microbe0
  on the CA and microbe1 on the LBM: the loop read into the copy lattices, so
  in-transit biomass was counted twice and the second microbe's biomass was never
  summed into `totalbFilmLattice` at all. With enough microbes it would have
  indexed off the end. Two correctly shaped vectors, `bio_ca` and `bio_fd`, are
  built beside the lattice lists they belong to.
- **A `.thm` file cannot silently discard one organism's energetics any more.**
  The overlap check compared resolved microbe ids, so a catch-all block (`microbe
  all`, id -1) never collided with a named one, and `registerModel()` then let the
  catch-all overwrite every slot in file order. A block written for one organism
  was discarded whenever a catch-all appeared after it, and that organism was
  gated by the wrong reaction's `dG0`, stoichiometry and threshold in every voxel
  of every step, with `describe()` still printing both blocks as though both
  applied. Mixing the two is ambiguous however it is resolved, so it is refused
  at start-up with both block names.
- **A surrogate `.srg` with no training range is refused rather than loaded.** An
  absent `trainmin`/`trainmax`/`xoffset`/`xgain` line left the vector empty,
  `valid()` returned true, the error string was empty, and `evalBound()` then read
  `trainMin[0]` off a zero-length vector at the first voxel. `complab3d_graphnet.hh`
  already refused a `.gnn` with no box on the grounds that a fitted model whose
  valid range is unknown cannot be trusted anywhere; the surrogate path, whose own
  header argues hardest that the box must travel with the weights, now does too,
  and names the missing line.
- **The Heijnen dissipation correlation no longer returns NaN above six carbons.**
  `pow(6 - NoC, 1.8)` on a negative base with a fractional exponent is NaN, so
  glucose was the last carbon source that worked: sucrose or a fatty acid printed
  `dG_dis nan` and returned a silent zero yield. The correlation is a distance
  from the six-carbon reference and Heijnen tabulates it on both sides, so the
  term is `|6 - NoC|`. The reduction term two lines down already squared its
  argument for exactly this reason.
- **A zero pore diffusivity on substrate 0 is caught instead of producing an
  all-NaN field.** Every solute's relaxation time is scaled against it, so
  `ade_dt` was `inf` and every `rate * dt` increment `inf` or `NaN`. The tau
  screen below it tests `< TAU_REJECT` and `> 2.0`, both false for NaN, so
  nothing caught it. Easy to hit: an immobile species is declared with
  `<in_pore>0.</in_pore>`, and putting one first in `<name_of_substrates>` was
  enough. The message says to move it down the list.
- **The flow solver no longer reports a run that hit the iteration cap as
  converged.** `[NS] Converged at iter=` printed unconditionally after a loop that
  may have exhausted `ns_max_iT1`, so the first line a reader checks could not be
  believed. It now says which happened, and the dissolution reopen re-solve, which
  used to exit silently, warns the same way before coupling its velocity field
  into every mobile solute.
- **The equilibrium solver's species list keeps water.** `setSpeciesNames()`
  dropped `H2O`, while the caller passes the full substrate list and reads the
  answer back by substrate index and the stoichiometry matrix keeps all rows. A
  run naming a substrate `H2O` read one past the end of the returned vector and
  matched every species after it against the wrong row. The list is kept whole
  and the solvent is skipped where it actually matters.

### Fixed, in the second pass

The first pass read the driver, the mineral paths and the twenty cases. This one
read everything else: all thirty-three files in `src/`, all ten programs in
`tools/`, the whole test suite, and every file in every example. It found more,
and two of these change results.

- **Biomass inside a biofilm diffused at the wrong rate.** The unit system is
  fixed in one line, `ade_dt = refNu * dx^2 / vec_solute_poreD[0]`, so the lattice
  viscosity of any physical diffusivity `D` is `refNu * D / vec_solute_poreD[0]`.
  That is what the substrate lines do and what `bioNUinPore` does. The
  biofilm biomass line divided by `vec_bMass_poreD[iM]` instead, applying the
  conversion a second time as a ratio to itself, so it was right only by accident
  when a microbe's pore biomass diffusivity happened to equal substrate 0's.
  Example 07 asks for 1e-10 in both places against a substrate 0 of 5e-10 and got
  tau 0.56 in pore but tau 0.8 in biofilm: biomass spreading five times too fast
  through exactly the voxels where biofilm biomass lives. Example 08's lattice
  Boltzmann microbe had the same factor. The tau screen could not catch it,
  because it recomputes the number from the same wrong expression.
- **`tools/geometry.py` wrote every geometry transposed.** `readGeometry()` takes
  one `ny*nz` x-slice per iteration, and all twenty shipped `preprocess.py`
  scripts write `for x: for y: for z:` with a comment saying so. `write_dat()`
  wrote z slowest and x fastest, which is the natural order for the `(nz, ny, nx)`
  array the module carries but the reverse of the file contract, and `read_dat()`
  made the same assumption. A geometry built or imported with this tool had the
  right number of values and no error anywhere, and the flow ran across what the
  user had drawn as depth; `geometry.py inspect` on a correct file reported the
  wrong porosity structure and tested the wrong axis for percolation. Both
  directions are corrected, and reading a shipped `geometry.dat` and writing it
  back now reproduces the file byte for byte.
- **A guard that excluded nothing.** `updateSoluteNbiomassLatticesDynamics3D`
  tested `mask != bb || mask != solid`. Those are distinct material numbers, so
  one inequality always holds and the condition is a tautology. All three sibling
  processors in the same file write the `&&` that was meant.
- **The cellular automaton indexed an empty vector.** When an over-capacity
  biofilm voxel has no non-wall neighbour, which is the sealed-throat case its own
  `[CA-ROBUST]` comment describes, the neighbour list is empty and the fallback
  branch read `nbrsLocMask[0..2]` from a default-constructed vector, whose data
  pointer is null. The comment said the voxel is left at the cap and the sweep
  moves on; now it is.
- **Two heap overflows on every single run**, before any branch:
  `strcat(strdup(str_inputDir.c_str()), ns_filename)` appends past an allocation
  sized for the directory alone, and the checkpoint path then appends `".chk"`
  past the end of that. Plus eight one-past-the-end NUL writes,
  `buf[item.size()+1] = 0` on a `calloc(item.size()+1)`. Silent with short paths,
  a `malloc` abort somewhere unrelated with long ones.
- **`half_saturation_constants` of the wrong length** printed "Terminating the
  simulation" and did not terminate, and it is the one branch that also pushes no
  row, so `vec_Kc` was left short and later indexed per microbe with no bound
  check.
- **`maximum_biomass_density` was only required with exactly one CA microbe**
  (`ca_count==1`), so a two-CA run that omitted it got a cap of 1e9 kg/m3, the
  redistribution never fired, and biomass piled up past capacity in one voxel for
  the whole run. The paired `thrd_biofilm_fraction` check ten lines up has the
  correct predicate.
- **`thrd_bFilmFrac` had no initialiser** and is only assigned when the tag is
  present, so an abiotic run, or a finite-difference or lattice Boltzmann biofilm
  run that omits it, read an indeterminate double and used it as the pore/biofilm
  reclassification threshold, which decides voxel identity, solute relaxation and
  the flow geometry.
- **A duplicated term in a guard**, `bMass_poreD[iT]>0 && bMass_poreD[iT]>0`,
  where the second was meant to be `bMass_bFilmD`. The absent sentinel is -99, so
  a microbe given one diffusivity and not the other passed the guard.
- **`ny` was computed with `getNx()`** in all four VTI writers, including the one
  that writes every substrate, microbe and mask volume. Harmless on the shipped
  cases, which are all at least as wide in x as in y; any domain with `ny > nx+2`
  silently wrote volumes truncated in y.
- **The COBRApy positivity guard was inverted.** It added the running increment in
  the `free` branch, where that increment is provably still zero, and added
  nothing in the `total` branch, which is the one its own comment is about. The
  stated protection was a no-op.
- **An explicit `<objective_direction>maximize</objective_direction>` could be
  silently flipped** to minimize by an SBML `fbc:type="minimize"`, because "the
  user asked for maximize" and "nobody said" were stored as the same value and the
  override could not tell them apart.
- **Four counters in the closing report were never reduced across ranks.** The
  symbolic, graph network, surrogate and thermodynamic gate reports are built
  inside data processors, so each rank has its own; `pcout` prints rank 0's. The
  gate's two verdicts, "closed everywhere, every time" and "never closed
  anywhere", are global claims that were being drawn from one rank, and each
  report returns the empty string at zero evaluations, so on a decomposition where
  rank 0 holds no biomass the report vanished entirely. The dissolution counters
  beside them have been reduced since v1.2 for exactly this reason.
- **The flow re-solve could run out of iterations in silence.** Both in-loop
  re-solves cleared the saturation flag when the outlet was still flowing and
  carried on, coupling a velocity field that never reached steady state into every
  mobile solute. That is the failure the warning on the initial solve was added to
  eliminate; the two in-loop paths now warn once per run.
- **`[DIAG] verdict: PASS` could be printed having compared nothing.** The first
  recorded row only establishes a baseline, so a run producing fewer than two
  diagnostic rows left the worst drift at zero and reported PASS. Checks now count
  their comparisons and a run with none reports NOT CHECKED.
- **The four outer-face boxes overlapped** on the four x-parallel edges, so both
  the open count and the total in the "N of M voxels have no boundary condition"
  message double counted `4*nx` voxels, and unequally: the total always, the count
  only where those edges were open.
- **Two ways a truncated `.srg` loaded and then misbehaved.** `W` and `B` are
  resized in pairs, so a file ending before its last bias line left that bias an
  empty vector the evaluator indexed; and nothing compared the layer count against
  the declared architecture, so a file cut at a layer boundary loaded as a shorter
  network in which a hidden layer silently became the output layer. Separately,
  the "no complete line" check for `trainmin`, `trainmax`, `xoffset` and `xgain`
  could only ever detect an absent line, never a short one, because the vector was
  pre-sized before parsing. A short `trainmin` meant the clamp report came out
  reassuring and wrong; a short `xgain` made that input's predictions meaningless.
- **`<walls>` was unreachable on the `import_raw` path**, being read inside the
  `<generate>` block, while `importRaw()` still acted on its `"y"` default. Two
  whole planes of an imported CT volume were turned into wall with no log line,
  and `<walls>none</walls>` could not switch it off.
- **A one-block `.thm` naming an organism other than microbe0** was accepted on the
  compiled kinetics path, which asks for the gate by hard-coded index 0. The gate
  then never applied, the report returned nothing, and the run was bit-identical to
  the thermodynamics being off, with the start-up echo of the unused block as the
  only trace. It is now refused with an explanation.
- **`makeKinetics.py` dropped fractional reactants.** `int(round(-c))` turns a
  coefficient of 0.5 into a reaction order of zero, so `0.5 A + B -> C` emitted a
  rate with no `A` factor at all while the stoichiometry stayed right, and the
  generated self-test passed because it checks conserved moieties, which is a
  property of the stoichiometry and not of the rate law.
- **The generated self-test reported PASS having checked nothing.** With no abiotic
  block the left null space was returned as the identity, so the "nothing to check"
  note never printed, every trial was skipped for having an all-zero rate, and the
  test closed with "PASS: every conserved moiety is conserved". It now says
  NOTHING CHECKED, and counts the trials that actually contributed.
- **`postprocess.py` never matched the filenames the solver writes.** It split each
  name into a stem and a trailing field index, on the belief that snapshots are
  called `subsLattice0_0000500.vti`. They are named after the field:
  `A_0000500.vti`. Every column came out as `A0_total` where the docstring, the
  `--conserve` help and every case README promise `A_total`, so `--conserve "A+C"`
  reported `SKIPPED: no column` and exited 0. A mass-balance check that silently
  checked nothing.
- **Four pipeline stages could not run at all.** `A1_geometry`, `A2_initial_fields`,
  `D1_fields` and `D2_balances` invoked `geometry.py` and `postprocess.py` with
  flags neither program defines, and without the arguments both require, so each
  exited 2 before doing anything. All four are rewritten against the real
  interfaces and were run end to end.

### Corrected in the documentation

Fifty-one verified factual errors, each checked against the code or the input
files before it was touched. The ones worth naming: the seeded biomass patch in
examples 05 and 06 is 108 voxels rather than 72, which also corrects "twelve
times as many voxels" to about eight; example 07 attributed a biomass
diffusivity to example 05, which declares none; example 09 pointed at a `models/`
directory it does not have; example 11 described the shipped surrogate as fitted
on *E. coli* core when it came from an iAF987 sweep that is not bundled; three
places claimed a `<conserve>` check runs in the solver when no example declares
one and `postprocess.py` does the balance; three referenced a `runAllExamples.sh`
that does not exist; nineteen `postprocess.py` files printed the binary at
`./build/complab` when it is written to the case root; and seventeen places
described a `closed` boundary as "Neumann both ends" when closed is bounce-back
and example 14's own README explains that a Neumann face is an open outflow that
would break the very balance it was being credited with.

From the second pass, the ones a reader would have been misled by:
`examples/08_two_microbes/CompLaB.xml` was not well-formed XML, carrying a `--`
inside a comment, and was the only file in the tree that `postprocess.py` could
not parse; example 20's section 6 still printed the pre-v1.3 generalized
effectiveness-factor formula, which disagrees with its own section 9 by 0.25 at
the case's own modulus and was the last copy of the old convention left standing;
example 20's aggregate is a circular cylinder spanning the full depth in z, not
the sphere it was called in four places, which matters because the classical
solution being compared against is posed on a sphere; `expected_sweep.csv` named a
command that cannot produce it, the default being 4 diffusion times where the file
was recorded at 5, and its bulk voxel counts are whole-domain counts while the
sampler excludes the inlet and outlet planes; example 19's calibration of its one
tuned parameter quoted a no-transfer-loss yield of 4.7 % where the shipped
constants give 5.18 %, reported as 5.00 % after the `yieldmax` cap; example 07 was
described as planktonic and advected by three of its four scripts, one of which
told the reader to delete the `<material_numbers>` entry that makes it the
attached biofilm the XML actually configures; example 03 asserted in four places a
total mass balance that its Dirichlet inlets make impossible, and example 02 asked
the reader to verify a y-independence its own staggered grain blocks rule out;
example 15's section 1 released a proton that its rate law never produces, and its
pipeline table promised a chattering check `postprocess.py` does not perform;
example 18's section 3 showed a `.gnn` `stoich` line in a syntax the parser
refuses; two READMEs claimed to be self-contained while needing GLPK at build
time; example 10 said the two flux-balance back ends cannot be mixed, a
restriction removed in v1.2; example 12 was called two competing populations in
three places when it is one organism on two rate paths; example 11's install step
named the copy of `surrogateModel.hh` the build does not compile, so a retrained
network silently had no effect; eighteen `preprocess.py` docstrings gave the
pre-padding grid as the size of the file they write; and two `<diagnostics>`
comments named an interval their own tag contradicts.

### Added

- `docs/methods/`: three reader's guides, one per rate path that is not a linear
  program plus one for the thermodynamic gate, each ending in a worked example
  that puts real numbers into every equation one substitution at a time. Every
  figure is a render of a slide in the matching `-figures.pptx`, so a figure is
  changed by editing the deck rather than redrawing it.
- `tests/test_upscale.cpp`, `src/complab3d_thiele.hh`.

## v1.3 — the two refinements of the linear program

Two things the published coupling literature needs that CompLaB3D could not do.
Both are off by default, both are additive, and a run that does not ask for
them solves exactly the program v1.2 solved.

### Multi-step (lexicographic) FBA — `<multi_step>`

Plain flux balance analysis maximises one objective and stops. The growth rate
is unique; the exchange fluxes very often are not. On the E. coli core model at
10 mmol gDW⁻¹ h⁻¹ glucose and oxygen, growth comes out at 0.5591 h⁻¹ every time
while acetate export may be anything between 9.9057 and 11.5033 — all of it
optimal, and whichever the simplex reached first.

`<multi_step>` solves a chain of programs instead, each holding the earlier ones
at their optimum, until nothing is free. The answer is then reproducible across
solvers, and the retain fractions let an organism trade growth for excretion,
which is what makes measured by-product yields reachable at all.

Implemented in both back ends — `run_glpk_lex()` in `complab3d_lexicographic.hh`
and `solve_multistep()` in `tools/complab3d_cobrapy.py` — so the two can be
cross-checked against each other, which is the only independent check either
has. Stages are named, not numbered, and resolved against the model at start-up.

Method: Song et al. (2025), *Sci Rep* **15**:6042. Example 21.

### Metabolic switching — `<cybernetic>`

An organism offered several carbon sources eats the best one and switches when
it runs out. FBA cannot produce that on its own, and neither can a Monod term,
because the same compound is a product in one regime and a substrate in the
next.

`<cybernetic>` solves the model once per source with the others shut, weights
the options by the carbon each would bring in, and blends the flux vectors. The
switch is continuous, needs no threshold, and adds no parameter the Monod rates
did not already carry. The blend conserves mass exactly, because a convex
combination of vectors in the null space of S is in the null space of S.

Written once, above the solver, so it works through GLPK, COBRApy and — in
principle — a trained surrogate, which is the combination Song et al. ran.

Method: Ramkrishna & Song (2012); Song et al. (2025). Example 22.

### One defect found while building the above

`run_glpk_lex()` first reported the growth rate from the chain's **first** stage
while returning the **last** stage's flux vector. With every retain fraction at
1.0 those agree, so it passed the obvious test. At 0.9 they differ by 10 % and
at 0.6 by 40 %, which would have had the solver add biomass at one rate while
removing substrate at a rate belonging to a different solution — mass silently
not conserved, in proportion to how far the calibration sat below 1. Growth and
fluxes are now both read out of the final vector, in both back ends.

### Not done

No enzyme state is carried between time steps, so `<cybernetic>` produces no lag
and no diauxic plateau. That is the same simplification the source paper makes,
and it is adequate when the switch is fast compared with transport.

## v1.2 — thermodynamic control, upscaling, and three passes through Known limitations

Two releases' worth of change: the thermodynamic factor described below, then
three passes through every entry in Known limitations to fix what could be fixed
and measure what could not, and finally the pore-to-aggregate upscaling those
passes had left as the last open item.

Eight of the original entries are gone, three are replaced by something narrower,
and seven defects nobody had documented were found on the way — a biomass solver
that produced negative biomass, four faces of the domain that never had a boundary
condition at all, and a boundary type that was being used to mean "closed" while
behaving as an open outflow.

**The one number that says whether this release did what it claims:** a fully
closed box with every reaction switched off used to lose 7 % of its protons. It
now conserves to 3.9 × 10⁻¹³.

### Fixed, in the limitations pass

- **The compiled surrogate now enforces its training range.** A network pasted
  into `surrogateModel.hh` used to extrapolate silently outside the box it was
  fitted over — the one failure a fitted model has that a hand-written rate law
  does not, and the only path where it was invisible. It now clamps to the box
  and counts, exactly as the `.srg`, `.sym` and `.gnn` paths do, and the closing
  report gives the percentage per path. `tests/test_surrogate_range.cpp` checks
  that the answer outside the box equals the answer at it, that one evaluation
  clamping two inputs counts once, and that the two paths keep separate totals.
- **The graph network is now stoichiometrically exact.** The readout ran on each
  species node, so the stoichiometry shaped the message passing but did not
  constrain the answer: on the shipped AOM network, HS produced over CH₄ consumed
  ran between 0.76 and 1.11 over the upper 80 % of the rate range, and between
  −1.57 and 12.2 once near-zero rates were included. `readout extent` runs it on
  each *reaction* node instead, giving one extent per reaction, and forms the
  species rates as *r = S ξ*. Worst deviation from exact is now **zero**, at
  every sample. It also fits better with fewer numbers — R rose from 0.9679 to
  0.9720 — because the network no longer spends capacity learning four rates that
  were always one rate times four coefficients. A file with no `readout` line
  still means the old behaviour, so every network written before this evaluates
  to exactly what it always did, and `tests/xval_gnn.py` cross-checks both modes
  against the Python trainer.
- **COBRApy and the surrogate can excrete.** Both discarded every net release
  (`if (draw > T()) draw = T();`), so a COBRApy organism could consume and never
  secrete, and a multi-output surrogate — trained specifically to return the
  fluxes the model excretes — had those fluxes thrown away. The GLPK path never
  had that line. It is gone from both.
- **COBRApy gained the repair loop.** Several organisms in one voxel were each
  given an uptake bound built from the full local supply and then scaled back by
  the joint budget, which keeps concentrations positive but splits the substrate
  by clamping rather than by re-optimising. It now marks the over-drawn
  substrates and re-solves against the pooled biomass, up to three passes, as the
  GLPK path has since the port. With one organism it never fires.
- **An immobile species inside a solid voxel now reports its own contents.**
  `NoDynamics::computeDensity` returns its stored density and ignores the
  populations, so every mineral inventory read back as exactly 1.0 mol/L whatever
  was seeded into it. On example 14 that meant 27.1 mol/L of calcite reading as
  1.0, the reopening test finding every grain below threshold, and the entire
  geometry converting to pore on **iteration zero** — after which the case
  dissolved nothing for the rest of the run. Solid voxels on an immobile
  species' lattice now take a dynamics that sums its populations. The case
  retains its calcite and performs 106 252 product deposits where it previously
  performed none.
- **The dissolution deposit now obeys the positivity clamp.** Every other
  increment in the solver obeys ΔC = max(RΔt, −C); this one did not, so a
  dissolution reaction that *consumes* a species drove its neighbours negative
  without limit.
- **A growth-only surrogate says so at start-up**, along with what follows from
  it: consumption falls back to a Monod term, and the organism cannot release a
  product at all.

### Fixed, in the second limitations pass

- **GLPK and COBRApy organisms can share a run.** The ban's stated reason — that
  the two back ends write into the same increment lattices with different flux
  index conventions — was no longer true: both address substrates through
  `cfg->subsLoc[globalMicrobe][substrate]` and convert through the same
  `fluxToDeltaC()`. What actually broke was narrower. Both processors were handed
  one list of FBA microbes, built over `usesFBA()`, so each tried to solve for
  the other's organisms, and `cfg->vec_lp[gM]` is null for a COBRApy microbe.
  `complab.cpp` now builds `glpk_globalId` and `cpy_globalId` separately, each
  with its own lattice vector, exactly as the surrogate, symbolic and graphnet
  paths already coexist. Verified on a two-organism run — the same toy model, one
  organism on each back end, in the same voxels: the two biomass fields are
  **bit-identical at every output step**. The one honest caveat, printed at
  start-up, is that shared-substrate exhaustion is repaired within a back end and
  not across the two.
- **A lattice-Boltzmann relaxation time near 0.5 is now rejected.** One diffusion
  coefficient is the reference; every other is carried onto the lattice as a
  ratio to it, and each ratio becomes its own relaxation time. Only the reference
  was ever checked. Example 07 shipped with a biomass diffusivity 1700× below the
  reference solute, which puts τ at **0.50018**, and BGK there does not diffuse
  slowly — it rings. Worst negative biomass as a fraction of the peak, measured
  on that case with every reaction off and a fully closed box: 25 % at τ =
  0.50018, 1.6 % at 0.510, 0.33 % at 0.520, 0.02 % at 0.550, zero at 0.80. The
  solver now audits every relaxation time that will actually relax — skipping
  immobile species and `CA`/`FD` biomass, whose τ is never used, so examples 13,
  14 and 15 are not rejected for a number nothing reads — warns below 0.55 and
  refuses below 0.51, naming the three ways out. Examples 07 and 08 were given a
  diffusivity the lattice can carry. After this, `CA` and `FD` agree on total
  biomass to five digits and **no field holds a negative value anywhere**.
- **The `FD` biomass solver is exercised.** It was listed as covered by no
  shipped case. Examples 05, 06 and 07 are now run side by side as a three-way
  comparison: `CA` and `FD` both give +0.229 % total biomass over 1000 steps and
  agree closely, with `FD` spreading the patch from 108 voxels to 876 while
  `CA` holds it in place — which is the difference between the two solvers, shown
  rather than asserted.
- **Every shipped rate law is proven safe on a short substrate list.** All
  thirteen already guarded their own lengths on entry; nothing checked that,
  which is how the next one comes to be written without a guard.
  `tests/check_kinetics_bounds.py` compiles each header on its own under
  AddressSanitizer and calls it with every substrate count from 0 to 5 and every
  microbe count from 0 to 3. Confirmed to fire by removing one guard. The README
  entry about the shipped `defineKinetics.hh` indexing `C[0]`–`C[94]` is gone:
  the shipped default writes zeros and says why.

### Fixed, in the third limitations pass

- **`closed` is now a boundary type, and Neumann means what it says.** There were
  two kinds: Dirichlet, which holds a value, and Neumann, which sets the plane to
  its neighbour every step. Neumann is an OUTFLOW, and every case that wanted a
  closed species had been written with it. Where the interior concentration rises,
  `FlatAdiabaticBoundaryFunctional3D` tops the plane up from nothing and streams
  that back in: on example 14 it manufactured **32 % of the calcium**. `closed` is
  bounce-back and exactly conservative, and every species in every shipped case
  that had Neumann on both sides has been switched to it. One subtlety that cost a
  measurement: the closed planes must be installed AFTER `initializeAtEquilibrium`,
  because it cannot reach a bounce-back cell — installed before, they start at
  density 1 and push it inward, and protons rose by 1.34 mol/L in a closed box
  where the reaction consumed 0.23.
- **Example 14's calcium balance closes.** Three defects, each found by
  measurement. The `released` counter summed every substrate's share, so it
  reported exactly twice the mineral dissolved (8.6052 against 4.3026); it now
  breaks out per species and prints the mineral removed beside them. The Neumann
  boundary above was the second. The third was the accounting, below. The case's
  own check now reads **residual −3.2 × 10⁻¹², relative 7.3 × 10⁻¹³, PASS**, from
  26 %.
- **Mass held where `computeDensity()` cannot see it is now measured.** Palabos's
  `BounceBack` and `NoDynamics` both answer from a stored number and ignore their
  populations, so mass in flight at a wall, mass resting in a grain, and mass in a
  closed boundary plane were absent from every total. Each field in the scalar
  record now carries a `_held` column — the difference between a full-block
  population sum and the reported total, so it cannot miss a category — and
  `_total + _held` is what `<conserve>` is checked against. Biomass is in the
  record too, which it never was. **In a fully closed box with every reaction off,
  the conserved sum is now constant to 3.9 × 10⁻¹³ for solutes and 1.3 × 10⁻¹⁴ for
  biomass; `_total` alone drifted by 7 %.**
- **The four unconditioned faces are closed by every shipped geometry.** Each
  `preprocess.py` now PADS rather than converts: one added layer of inert wall
  outside any face still open, with `NY`/`NZ` grown by two to hold it, so the pore
  space, the grains and the patches keep exactly the volumes the case declares.
  Walling existing voxels instead would have spent a third of a six-deep slab on
  the boundary condition. All 20 geometries regenerated, all XML dimensions synced.
- **Every shipped geometry is the one its generator writes.** They had drifted:
  examples 13, 14 and 15 shipped a pore space with 696 mineral voxels while their
  own `preprocess.py` wrote 72. Both ran. Which one a published number came from
  was unknowable. `tests/check_repo.sh` now regenerates each geometry, diffs it
  against the shipped file, and cross-checks `<nx>/<ny>/<nz>`.

### Added

- **`<upscaling>`: one resolved aggregate, reduced to a continuum rate.** Reports
  the effectiveness factor *η* = ⟨r⟩/r(C_bulk), the Thiele modulus, the classical
  sphere result beside it, and the thermodynamic gate at the bulk composition. The
  numerator is read from the increment lattices after the rate processors run and
  before the increments are applied, so it averages the number the solver is about
  to use rather than a second model of it. At iteration zero, when every voxel is
  still at the bulk composition, it reports **η = 1.0000000002** — the measurement
  checking itself. `<freeze_biomass>` holds the catalyst while the concentration
  profile relaxes; without it *η* never settles, because biomass grows faster at
  the rim than in the core, and every run reported NOT STEADY at five diffusion
  times. The solver checks steadiness itself and refuses to call a transient a
  result.
- **Example 20, `upscaling`**, and `tools/upscale_sweep.py`, which sweeps aggregate
  radius against bulk sulfide — the product that shuts the gate, so the second axis
  is how close the bulk sits to the energy threshold — and fits the correction the
  classical Thiele curve needs once the reaction has an energy limit.

### Documentation

- **Every example README rewritten to one structure.** Each of the twenty now
  opens with what the case is for in plain language and the reaction it runs,
  then a **What is simulated** table — geometry, flow, solute transport, abiotic
  and biotic reaction, how each organism's biomass is solved, geometry evolution,
  run length — then the physics of that case with every constant carrying its
  unit and its source, then **What to check**, then what the case does *not* do
  and which case does it instead. The table is generated from that case's own
  `CompLaB.xml` and `input/geometry.dat`, so it cannot describe a switch the case
  does not set.
- **`examples/README.md`, `config/geometry/README.md` and `CONTRIBUTING.md`
  brought back into line with the tree.** They still described eighteen examples,
  four shared geometries serving all of them, and a domain size of 24 × 12 × 12.
  There are twenty examples; every one generates its own geometry and **no
  example reads the four shared files**, which are unpadded starting points; and
  the slot cases run at 24 × 26 × 8, because `preprocess.py` adds the wall layer
  the solver does not supply on four faces.
- **Example 10's offline check no longer fails on its own model file.** It called
  `cobra.io.read_sbml_model()` on `input/toy_model.xml`, which is the flat
  `<Metabolic_Model>` format the C++ side parses, not SBML — COBRApy never opens
  that file, it is handed the stoichiometry as arrays. The check now confirms
  `cobra` imports, prints which interpreter it came from, and validates the flat
  model's dimensions.

### Measured rather than fixed

- **The four unconditioned faces remain the geometry's job, not the solver's.**
  Every shipped case now walls them, and every run counts any that are open and
  says what follows from them. But a domain written by hand is still open there
  unless you wall it. `<outer_faces>sealed</outer_faces>` will do it at the cost
  of deleting that layer of pore. Two attempts at a zero-gradient closure — which
  is what a thin slab actually wants, and costs no volume — were tried and are not
  shipped; `src/complab3d_outerfaces.hh` records exactly how each failed, both on
  the same Palabos convention that `NoDynamics` returns its own stored density and
  ignores the populations.
- **Periodicity is still not offered.** Palabos would wrap the lattices, but the
  cellular automaton, the finite-difference biomass step and the dissolution
  gather all walk their neighbours with hand-written index arithmetic that stops
  at the edge. A periodic option would be true for the lattice-Boltzmann fields
  and false for everything else in the same run.

### Resolved by the bounce-back fix

- **Examples 17 and 18 no longer clamp.** Both used to report about 99.8 % of
  evaluations outside the fitted box. Re-run against the fix, both report
  **0 clamped, 0.00 %**. The rate paths were never at fault.

## v1.2 — thermodynamic control

Every rate path answered how fast an organism *can* run its reaction. None asked
whether the reaction releases enough energy, at the local concentrations, to be
run at all. This release adds the factor that answers it, and fixes two things
found while testing it that were quietly wrong for every case in the repository.

### Added

- **`<thermodynamics>`: a per-voxel thermodynamic factor** (`src/complab3d_thermo.hh`).
  ΔG = ΔG° + RT ln Q, then F_T = max(0, 1 − exp((ΔG + m ΔG_ATP)/(χRT))) after
  Jin & Bethke (2003) and Craig (2024, Eqs. 3.3–3.6). It is not a seventh rate
  path: it multiplies whichever path each organism already uses, so all six
  inherit it and none was rewritten. Energetics live in a `.thm` file read at
  start-up and echoed into the log, for the same reason a `.sym` rate law does.
- **Energy-based growth yield**, Heijnen & van Dijken (1992) dissipation and
  Craig Eq. 3.10, optional per reaction block.
- **`complab_thermo::fbaGrowthEfficiency`**, Craig Eq. 3.20: turns a biomass
  flux into a number comparable with a measured growth efficiency.
- **`examples/19_thermodynamic_gate`** — an ANME–SRB aggregate whose dual-Monod
  rate law never reads its own products. Ships an offline step that reports
  where the gate closes *before* the solver runs, and a pipeline that runs the
  case twice, gated and ungated, because a gated run alone says nothing about
  what the gate did.
- **`tests/test_thermo.cpp`** — the four equations checked against the same
  arithmetic written out longhand, and thirteen malformed files each checked to
  be refused rather than half-loaded.
- **`config/CompLaB.everything.xml`** — one line for every ability the solver
  has, each with its unit and default, grouped into the twelve blocks the XML
  actually has. It is an index, not a case; `CompLaB.reference.xml` remains the
  manual. `tests/check_repo.sh` walks the source for every XML path the solver
  opens and fails on any the index does not list, so the claim stays true by
  construction rather than by diligence.
- **`tests/check_repo.sh`** now fails if any biotic rate processor stops
  consulting the gate.

### Fixed

- **Bounce-back voxels started at C = 1 mol/L on every advection–diffusion
  lattice.** Palabos stores an AD population as a deviation from equilibrium at
  density one, and `BounceBack::computeEquilibrium` returns zero whatever
  density it is handed, so `initializeAtEquilibrium` could not reach those
  nodes. Bounce-back conserves mass, so nothing removed it: every wall voxel
  streamed C = 1 fluid into its neighbours from the first step. In a millimolar
  case that is a source three orders of magnitude larger than the chemistry. It
  presented as every field rising several-fold over the first hundred steps and
  then relaxing, which looks like a transient and is not. Measured on a closed
  domain with no reaction: methane total 10.75 → 65.13 with walls, held to five
  digits without them. This is the cause of the field divergence noted in
  examples 17 and 18.
- **The solver ignored its command-line argument.** All six readers opened the
  literal string `CompLaB.xml`, so `./complab variant.xml` ran the wrong file
  and said nothing. A control run and its experiment came back byte-identical,
  and the natural reading of that is that the setting under test does nothing.
  One path, `complab_input::configPath()`, set once from argv.
- **`pipeline.sh` in every example, and the run instructions in `README.md`,
  `INSTALL.md` and `docs/`, pointed at `./build/complab`.** The binary is
  written to the case root, not into `build/`.

## v1.1 — the learned rate paths are wired

Two rate paths shipped in v1.0 implemented, unit-tested and never called. This
release connects them, and adds the checks that would have caught it.

### Fixed

- **The symbolic and graph-network rate paths are now reachable from the
  solver.** `reaction_type` accepts `symbolic`, `graphnet` and their
  `_and_kinetics` combinations; `<symbolic><enabled>` and
  `<graphnet><enabled>` are read; `<expressions_file>`, `<network_file>` and the
  two `<abiotic_file>` tags are loaded at start-up, bound to this run's
  substrate and microbe names by name rather than by index, and refused with a
  message that says which name did not resolve; and `run_symbolic3D`,
  `run_graphnet3D`, `run_symbolic_abiotic3D` and `run_graphnet_abiotic3D` are
  constructed in `src/complab.cpp`. Verified by building examples 17 and 18
  against Palabos v2.3.0 and running both end to end.
- **`scripts/setup_case.sh` did not copy `surrogateModel.hh`**, so *no*
  assembled case compiled — the solver includes that header unconditionally.
  It is now laid down with the other two user-editable defaults.
- **The out-of-range counters counted values, not evaluations.** A two-variable
  rate law could report 199 % of its evaluations as clamped and a four-species
  network 400 %, which makes the one number a reader uses to judge a result
  unreadable. Both now count one per evaluation that clamped anything.

### Added

- `tests/test_wiring.cpp` — every spelling parses, the predicates do not
  overlap, a program and a network bind to this run's names and refuse a
  mismatch, and the runtime registry works.
- Two structural checks in `tests/check_repo.sh`: every rate processor is
  constructed in `src/complab.cpp`, and every header in `src/` is reachable
  from it. Either would have failed on v1.0.

### Known, and stated in the README rather than discovered later

- Examples 17 and 18 run, but their fields leave the fitted box within a
  hundred steps. Setting `reaction_type` to `none` gives the identical
  trajectory, so the cause is in those cases' transport setup, not the rate
  path.
- The graph network holds the stoichiometric ratio approximately, not exactly.

## v1.0 — first release of CompLaB3D-GEMS

The six capabilities, the eighteen examples and the pipelines, gathered into one
repository.

### Added

- **Mineral precipitation** (`src/precipitationVOP.hh`). Surface-controlled
  reaction, node conversion at a declared fill density, flow re-solve on the
  reduced pore space, and clogging detection.
- **Mineral dissolution** (`src/dissolutionVOP.hh`). Declared solid phases,
  averaging over open faces, mass-conserving release into the neighbouring
  pore, and hysteresis between the sealing and reopening thresholds.
- **Flux balance analysis** through GLPK in process (`complab3d_glpkcpp.hh`)
  and through COBRApy in an embedded interpreter (`complab3d_pythonAPI.hh`),
  with warm-started simplex between voxels.
- **Surrogate growth network** (`src/surrogateModel.hh`), fitted offline to
  flux balance output.
- **Symbolic regression** (`src/complab3d_symbolic.hh`), a rate law read from a
  `.sym` text file at start-up, with ranges enforced and clamps counted.
- **Graph network** (`src/complab3d_graphnet.hh`), message passing over the
  bipartite species-reaction graph, returning the whole rate vector at once.
- **`pipelines/`** — what to run before and after the solver, with a runnable
  `run.sh` at every stage.
- **`scripts/setup_case.sh`** — assembles a runnable working directory from the
  shared defaults, the case's `case.files` manifest, and the case's own files.
- **Corrected several offline commands that had never been run.** `extractMM.py`
  takes a positional model and `-o`, not `--model/--out`, and has no `--check`;
  `fit_symbolic.py` takes `--inputs`, not `--vars`; `train_graphnet.py` takes
  `--data`, not `--samples`; `xval_gnn.py` is a parity self-test that takes no
  arguments; neither the model converter nor the sweep reads `.gz`, so a gunzip
  step was missing; and the *E. coli* core objective in the shipped model is
  `Biomass_Ecoli_core`, not the `BIOMASS_Ecoli_core_w_GAM` quoted in the
  pipeline stages. Every offline step in the repository has now been executed
  end to end rather than only written.
- **An example folder is now the whole procedure.** Each of the eighteen holds
  its pore space, the metabolic model it reads, the training code that produced
  whatever it was trained on, the data that code was run on, its rate-law file,
  and its chemistry — beside the configuration that uses them. Nothing is
  fetched from elsewhere in the tree, and `case.files` is gone.
- **Duplication inside `examples/` is now deliberate and checked.**
  `tests/check_repo.sh` still forbids two identical files outside `examples/`,
  and additionally compares every copy under `examples/*/training/` against its
  original, so a tool fixed in one place and not the others is caught rather
  than shipped.
- **Every example carries its own pipeline.** Each of the eighteen now holds a
  `preprocess.py` that builds that case's pore space and checks it percolates, a
  `postprocess.py` that knows which of that case's totals are closed and runs its
  balance check, an `offline.sh` where there is offline work, and a `pipeline.sh`
  that runs the chain. The pre- and post-processing scripts import nothing beyond
  the standard library, so an example directory is the whole procedure rather than
  a pointer into a shared tools directory.
- **Examples 17 and 18** — a symbolic rate law and a graph network, each with
  the file it reads, so both paths have a case of their own rather than only a
  shipped artefact. Neither needs offline preparation to run.
- **`tests/test_surrogate_parity.cpp`** — the surrogate forward pass exists
  twice in this codebase: once generically in `complab3d_surrogate.hh`, once
  written out with literal loop bounds in the generated header, because that is
  the shape MATLAB's `genFunction` emits and what every previously exported
  network looks like. Two implementations of one calculation are a liability
  unless something checks they agree; this is that check.
- **`tests/test_surrogate_multi.cpp`** — round-trips a multi-output network
  through the file, checks the outputs are scaled independently, checks a
  single-output file still loads and still means the same thing, and checks that
  a file whose declared output count disagrees with its architecture is refused
  rather than accepted and mis-wired.
- **`tests/check_repo.sh`** — structural checks on the tree itself: no
  duplicated file, every `case.files` entry resolving, every example assembling,
  every internal link alive. Runs in a second, needs nothing installed.

### Changed

- **The surrogate can now return the exchange fluxes, not only growth.** A
  growth-only network leaves the solver to guess consumption with a Monod term,
  which is exact only where the swept uptake bound was the binding constraint
  and cannot release a product at all. The sweep
  (`tools/surrogate/generateTrainingData.py`) now records the whole flux vector
  rather than the objective alone; `trainSurrogate.py` fits every column and
  emits a header whose last layer is a matrix; `verifyExport.py` checks each
  output separately; `inspectSurrogate.py` reports each one. `--growth-only`
  reproduces the previous behaviour exactly.
- **`trainSurrogate.py --srg` writes the run-time format too.** Before this the
  only producer of a `.srg` was the in-run trainer inside the solver, which fits
  growth alone — so a multi-output network could not be loaded at run time at
  all, only compiled in. Both files now come from one fit.
- **The `.srg` run-time format carries multiple outputs**, with per-output
  scaling, an `outputs` count and an `outputnames` line, all validated against
  the architecture on load. Switching between a growth-only and a
  growth-and-fluxes network is therefore one line of XML — `<weights_file>` —
  with no rebuild. Files written before this change load unchanged.
- **No change was needed in the solver.** `defineSurrogateModel` already
  received `Fout` beside `bioR`, and the per-voxel processor already read both;
  the slot was there and nothing was filling it.

- **One reference configuration.** `config/CompLaB.reference.xml` documents
  every tag; it now covers the symbolic and graph-network paths, which were
  previously described only in loose example fragments.
- **The examples deduplicated, entirely.** Thirty-nine byte-identical files
  removed: twenty-three copies of the two kinetics headers, twelve copies of
  four pore geometries, and two copies of the FBA toy model. What is shared now
  lives once — `config/kinetics/`, `config/geometry/`, `models/toy_model.xml` —
  and each example names what it borrows in a `case.files` manifest that
  `scripts/setup_case.sh` resolves. Six genuine chemistry overrides kept.
- **`scripts/setup_case.sh` assembles in three passes** — shared defaults, then
  the case's `case.files`, then the case's own files — so an example can
  override anything it inherits.

### Known limitations

- **The symbolic and graph-network processors are not yet called from
  `src/complab.cpp`.** Both are implemented and covered by the regression
  suite, but `<expressions_file>` and `<network_file>` are not read by the
  shipped solver, so examples 17 and 18 are correct configuration waiting on a
  call site. The `.sym` and `.gnn` loaders, the expression evaluator and both
  per-voxel drivers are exercised by `tests/run_tests.sh`.
- Dissolution products are lost at MPI block boundaries, so dissolution results
  depend on the processor count. Run on one process where a quantitative mass
  balance matters.
- The COBRApy path clamps positive exchange draws to zero and does not run the
  repair loop the GLPK path uses.
- The surrogate path does not enforce its training range; the symbolic and
  graph-network paths do.
