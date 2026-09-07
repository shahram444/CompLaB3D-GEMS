# 11 — a fitted network in place of the linear program

## 1. What this example does

Flux balance analysis answers the right question but pays a linear program for
it, in every voxel, at every step. That is the most expensive rate path in this
repository. This case replaces the linear program with a **small feed-forward
network fitted to its answers offline** — a few hundred multiplications instead
of a solve, two to three orders of magnitude faster — and runs the same kind of
pore-scale case on top of it.

The organism is *Geobacter metallireducens* respiring **acetate with Fe(III) as
the terminal electron acceptor**, which is the reaction that dominates iron
reduction in anoxic sediments:

```
    CH3COO(-)  +  8 Fe(III)  +  4 H2O   ->   2 HCO3(-)  +  8 Fe(II)  +  9 H(+)
```

Neither the growth rate nor the uptake fluxes are written down anywhere in this
case. They are what the network returns when it is handed the two local
concentrations, and the network returns them because it was fitted to a sweep of
the linear program over exactly those two bounds.

**The single thing to understand before running this.** A fitted network is
valid only inside the box it was trained on, and outside that box it does not
fail — it extrapolates, confidently, and returns a number that looks like every
other number. The shipped network was trained over

```
    acetate    0.00089   to  9.998    mmol gDW-1 h-1
    Fe(III)    0.000035  to  0.4997   mmol gDW-1 h-1
```

so this case sets `<fba_maximum_uptake_flux>` to `8.` and `0.4`, safely inside.
Raise the second above 0.5 and the run leaves the training box. It will still
print results. They will not mean anything. `training/inspectSurrogate.py`
recovers those bounds from any network by inverting its `mapminmax` scaling, and
the evaluator clamps to the box and counts how often it had to.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 2 of them: `acetate`, `Fe3` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Geobacter` — a fitted surrogate network in place of the linear program |
| **Biomass** | `Geobacter` — attached biofilm (seeded on its own material number), finite difference — biomass spreads by diffusion on the same grid |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 500 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. The network, and how it is evaluated

The shipped `src/surrogateModel.hh` is a MATLAB-style network written out as
plain C++:

```
    2 inputs                 the two uptake bounds, mmol gDW-1 h-1
      -> mapminmax           scale each input to -1 .. +1 using the
                             training minimum and maximum
      -> 4 hidden layers     10 tansig neurons each
      -> 1 linear output     no squashing, so the output is unbounded
      -> mapminmax reverse   back to a growth rate
```

Every layer is one matrix multiply, one bias add and one `tanh`. The header
writes them out one line per layer, weights and all, so the whole forward pass
is readable rather than a call into a library:

```
    a1 = tansig( b1 + IW1_1 * x )
    a2 = tansig( b2 + LW2_1 * a1 )
    a3 = tansig( b3 + LW3_2 * a2 )
    a4 = tansig( b4 + LW4_3 * a3 )
    y  =          b5 + LW5_4 * a4        the last layer is LINEAR
```

The uptake bounds handed in are built exactly as they are on the flux-balance
path — Michaelis–Menten on the local concentration, capped by
`<fba_maximum_uptake_flux>`:

```
    V  =  Vmax * C / ( Kc + C )
```

so the only difference between this case and example 09 is who turns those two
bounds into a growth rate.

| tag | value | what it does |
|---|---|---|
| `<enable_surrogate>` | `true` | builds the surrogate path |
| `<reaction_type>` | `surrogate` | routes *this organism* through it |
| `<fba_maximum_uptake_flux>` | `8. 0.4` | Vmax per substrate; **both inside the training box** |
| `<half_saturation_constants>` | `0.05 0.01` | Kc in the expression above, mol L⁻¹ |
| `<biomass_molar_mass>` | `24.6` | g mol⁻¹, converting the network's gDW basis to the solver's moles |
| `<weights_file>` | *(unset)* | set it to a `.srg` file to swap networks without rebuilding |

## 4. Where the network comes from

`offline.sh` is the longest preparation in this repository and the only one that
ends in a recompile. The repository ships a fitted header, so **this case runs
without any of it**; the script demonstrates the procedure on the bundled
*E. coli* core model, and is how you would make a header for a model of your own.
The shipped `src/surrogateModel.hh` came from a sweep of the *Geobacter
metallireducens* model iAF987, which is not bundled, so running `offline.sh` as it
stands fits a new network rather than reproducing the one in the repository.

```
    1.  sweep      solve the linear program on a GRID x GRID mesh of the two
                   uptake bounds                        -> training_data.csv
    2.  fit        train the network on that sweep      -> surrogate_weights.hh
                                                           surrogate_weights.srg
    3.  verify     check the exported header reproduces the trainer,
                   output by output
    4.  inspect    look at the response surface before trusting it
```

Step 1 is the cost: `GRID^2` linear programs. The default grid of 41 is 1681
solves and a couple of minutes; the shipped header used 141, which is 19 881.

Step 3 is not optional. A transcription error between the fitted weights and the
C++ header is invisible in the numbers and fatal in the results, so `verifyExport.py`
re-evaluates the header against the trainer for every output.

The sweep records more than growth. A linear program returns a whole flux
vector, and the exchange fluxes in it are the ones the solver would otherwise
have to guess with a Monod term — a guess that is exact only where the swept
bound was the binding constraint and wrong everywhere else. So every swept
exchange is fitted alongside growth. Those flux columns are the harder fit:
growth is a smooth surface, while an exchange flux is piecewise linear with
kinks where the binding constraint changes, and a network smooths exactly those
kinks. Look at them in step 4 before believing them.

Two files come out of one fit and hold the same numbers:
`surrogate_weights.hh` is compiled in (fastest, needs a rebuild) and
`surrogate_weights.srg` is read at start-up through `<weights_file>` (swap
networks by editing `CompLaB.xml`). `tests/test_surrogate_parity.cpp` checks the
two evaluators agree.

## 5. What to check

1. **the run stays inside the training box.** The evaluator clamps to it and
   counts how often; a large count means the case has wandered outside the
   region the network knows and the numbers should not be read;
2. **growth responds to both substrates.** Lower the acetate boundary and growth
   should fall; lower Fe(III) and it should fall too — a network that responds
   to only one input was fitted on a sweep where the other never bound;
3. **acetate and Fe(III) are drawn together**, in something near the 1 : 8 ratio
   of the reaction above, wherever neither bound is limiting;
4. **no `[NEG!]` warnings.** A surrogate has no stoichiometry enforcing
   non-negativity, so the increment is clamped against what is locally present
   before it is applied — that clamp firing often is itself a signal;
5. **against example 09.** The surrogate is only worth its speed if it agrees
   with the program it replaced, on the region it was fitted over.

## 6. What this case does not do

No flow, no abiotic reaction, no geometry evolution, no thermodynamic control.
The network predicts growth from two numbers and nothing else — no pH, no
temperature, no inhibition that was not in the sweep. Example 17 fits a
**symbolic** law to the same kind of sweep instead, which you can read and
extrapolate from; example 18 fits a graph network over the species–reaction
graph.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 11_surrogate run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and reports porosity, refusing to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Sweeps the linear program, fits the network, verifies the exported header against the trainer, and plots the response surface. Not needed to run the case — a fitted header ships with the repository. |
| **run** | `CompLaB.xml` | What the solver reads. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's checks. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `surrogateModel.hh` | The fitted network of section 3, compiled in, at the case root beside `CMakeLists.txt`. Ships with the repository, fitted on iAF987; `offline.sh` is the same procedure, run on the bundled *E. coli* core model. **[v1.3]** A newly fitted header must be installed over THIS copy. `setup_case.sh` also lays one down in `src/`, but that one is never compiled: `src/complab3d_processors_surrogate.hh` includes the header as `"../surrogateModel.hh"`. Copying a new fit into `src/` produces bit identical results from the old network, with no error. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

### The offline code this case ships

| File | What it is |
|---|---|
| `training/generateTrainingData.py` | Sweeps the linear program over the two uptake bounds and records growth plus every swept exchange flux. |
| `training/trainSurrogate.py` | Fits the network and writes both the compiled header and the run-time `.srg`. |
| `training/trainSurrogate.m`, `exportSurrogateHeader.m`, `testExporterParity.m` | The MATLAB route to the same two files, for anyone whose networks are fitted there. |
| `training/verifyExport.py` | Re-evaluates the exported header against the trainer, output by output. Run it every time. |
| `training/inspectSurrogate.py` | Recovers a network's training range from its scaling blocks and evaluates it at a point. |
| `training/extractMM.py` | Exports an SBML or BiGG model to the flat format the sweep reads. |
| `training/README_SURROGATE_TRAINING.md` | The long form of section 4. |
| `models/e_coli_core.xml.gz` | The model `offline.sh` sweeps by default, so the whole procedure can be run here. The shipped header came from an iAF987 sweep, which is not bundled. `models/NOTICE.md` carries this model's provenance and licence. |

### A note on the domain

`x = 0` and `x = nx-1` carry the boundary conditions named per substrate. The
solver gives the other four faces nothing — not a wall, not a symmetry plane, not
periodicity — so `preprocess.py` draws an inert wall there, and `NY` and `NZ` are
two larger than the pore space they hold. The layers are **added**, not taken out
of the pore space, so porosity and every count are what the case declares.

> **Why shared code is copied here rather than referenced.** So that this folder
> *is* the procedure. The cost is real and worth stating: a fix to a shared tool
> has to be applied to every case that carries it, and `tests/check_repo.sh`
> fails if a copy drifts from `tools/`.
