# 16 — model fetched, network trained and verified, all from one config file

## 1. What this example does

Every earlier metabolic case needs something prepared before the solver starts:
a model exported by a Python script, a network fitted by an offline sweep, a
header recompiled into the binary. This case needs none of it. The metabolic
model is unpacked from the bundle in this folder, the exchange reactions are
matched **by name** against that model's own reaction list, the surrogate network
is **trained during start-up** from a sweep of the linear program, and the fit is
checked before the first transport step. Copy the folder, build, run.

The organism is *E. coli* on the `e_coli_core` model — glucose in, oxygen in,
acetate out — growing as attached biofilm in a diffusive pore space:

```
    glucose  +  O2   ->   biomass  +  acetate
```

Acetate is the interesting output. It is not a waste product the configuration
declares; it appears because *E. coli* overflows to fermentation when the carbon
supply outruns the oxygen supply, and the linear program finds that on its own.
In a pore space with both delivered by diffusion from the same face, that
happens deep in the domain where oxygen has been drawn down but glucose has not.

## 2. What is simulated

| | |
|---|---|
| **Geometry** | 24 x 26 x 8 voxels at 10 um, 3168 open (porosity 0.6346) |
| **Flow** | **not solved** — `<Peclet>0</Peclet>`, so transport is pure diffusion |
| **Solute transport** | D3Q7 advection-diffusion, one lattice per species — 3 of them: `glucose`, `o2`, `acetate` |
| **Abiotic reaction** | none |
| **Biotic reaction** | `Ecoli` — a fitted surrogate network in place of the linear program |
| **Biomass** | `Ecoli` — attached biofilm (seeded on its own material number), finite difference — biomass spreads by diffusion on the same grid |
| **Geometry evolution** | none — the pore space is fixed |
| **Run length** | 400 advection-diffusion steps |

Everything above is read from `CompLaB.xml` and `input/geometry.dat`; nothing in
that table is a description that can drift from the case.

## 3. Where the model comes from

```xml
    <model_source>bigg:e_coli_core</model_source>
    <model_bundle>models</model_bundle>
    <model_cache>input</model_cache>
    <allow_download>false</allow_download>
```

`e_coli_core` ships gzipped in `models/`. At start-up it is unpacked into
`<model_cache>` and checked against `models/manifest.txt`, so a model that has
been revised upstream is visible in the log rather than silently different.

`<allow_download>false</allow_download>` is the default, and it is why this case
runs on a compute node with no route to the internet. Nothing is fetched.

## 4. Named exchange reactions

```xml
    <exchange_reaction_names>EX_glc__D_e EX_o2_e EX_ac_e</exchange_reaction_names>
```

One reaction name per substrate, in `<name_of_substrates>` order, resolved
against the model's own reaction list at start-up.

This is the safe form, and the reason is worth being blunt about. The
positional `<exchange_reaction_indices>` that examples 09, 10 and 12 use **cannot be
checked**: a wrong index produces a well-posed linear program with the wrong
biology, which runs to completion and reports plausible numbers. A wrong *name*
stops the run and prints the near matches. Try it — change `EX_o2_e` to
`EX_oxygen_e` and rerun.

There is no `<model_filename>` in this case; `<model_source>` supplied it.

## 5. The surrogate, fitted at start-up

```xml
    <surrogate>
        <enabled>true</enabled>
        <weights_file>output/ecoli.srg</weights_file>
        <train_if_missing>true</train_if_missing>
        <train>
            <samples>300</samples>          points swept from the linear program
            <restarts>2</restarts>          independent fits, best kept
            <epochs>400</epochs>
            <verify_points>64</verify_points>   LP re-solves the fit never saw
            <layers>8</layers>
            <seed>20260814</seed>
            <input0>glucose 0.1 10.0</input0>
            <input1>o2      0.1 20.0</input1>
        </train>
    </surrogate>
```

No weights file exists the first time, so one is made: the linear program is
swept over the two ranges, a network is fitted to what it returned, and the fit
is checked by re-solving the linear program at 64 points it never saw. The
result is written to `<weights_file>`, and every later run of this folder loads
it in milliseconds instead of retraining. Delete `output/ecoli.srg` to force a
retrain.

**The training ranges and `<fba_maximum_uptake_flux>` are the same numbers, and
that is not a coincidence.** Michaelis–Menten bounding means no voxel can ever
request more than Vmax, so the top of each training range *is* Vmax:

```
    input0   glucose  0.1 .. 10.0        fba_maximum_uptake_flux  10.
    input1   o2       0.1 .. 20.0                                 20.
```

Raise one without retraining and the run leaves the box the network was fitted
in. It will say so, and say how often — but it will not stop, and outside the box
a network does not fail, it extrapolates. Example 11 is the longer discussion of
that failure mode.

## 6. Two settings that are not optional for a biofilm organism

```xml
    <maximum_biomass_density>50</maximum_biomass_density>
    <thrd_biofilm_fraction>0.01</thrd_biofilm_fraction>
```

Leaving them out is not harmless, and the comment in `CompLaB.xml` records what
happens because it happened.

`<thrd_biofilm_fraction>` defaults to 0, which makes the *biomass at or above
threshold* test true in **every** pore voxel, empty ones included. An empty voxel
has no species to label, so the mask update reaches its "should never happen"
branch and prints an error — once per empty voxel per sweep. The first run of
this case wrote a 28 MB log of `Error: Updating mask failed.` and the results
were still fine, which is the worst combination: loud, and not actually wrong.

`<maximum_biomass_density>` defaults to 1e9 kg m⁻³, which is not a density any
biofilm reaches, so the cap never engages and biomass grows without bound in the
voxels that are best supplied.

## 7. What to check

1. **the model resolves.** The start-up lines name the model, its reaction
   count, and each exchange reaction matched by name. If any of that is not what
   you meant, stop there;
2. **the fit is verified before the run starts.** The 64 held-out points are
   re-solved with the linear program and compared; a poor agreement reported at
   start-up means the network is not a stand-in for the program and nothing
   after it is worth reading;
3. **the run stays inside the training box**, and the out-of-box count stays
   small;
4. **acetate appears where oxygen is drawn down.** That is the overflow
   metabolism of section 1, and it is the result that would be hard to get from
   a hand-written rate law;
5. **the second run is fast.** It should load `output/ecoli.srg` rather than
   retrain — if it retrains, the weights file was not written where
   `<weights_file>` says.

## 8. What this case does not do

No flow (`<Peclet>0</Peclet>`), no abiotic reaction, no geometry evolution, no
thermodynamic control. And the surrogate is fitted over two inputs only:
everything else the model could respond to is held at whatever the sweep held it
at. A network is a stand-in for the linear program **on the axes it was swept
over**, and on no others.


## Everything this case needs is in this folder

No cross-referencing, nothing to fetch.

```bash
./scripts/setup_case.sh 16_complete_pipeline run/mycase
cd run/mycase
./pipeline.sh
```

### The pipeline

| | File | What it does |
|---|---|---|
| **pre** | `preprocess.py` | Builds the pore space and reports porosity, refusing to continue if it does not percolate. Standard library only. |
| **offline** | `offline.sh` | Optional. The model is unpacked and the network trained by the solver itself at start-up; this script is here for inspecting the bundled model by hand. |
| **run** | `CompLaB.xml` | What the solver reads — including the model source, the named exchange reactions and the surrogate training block. Every tag is documented once, in [`../../config/CompLaB.reference.xml`](../../config/CompLaB.reference.xml). |
| **post** | `postprocess.py` | Reads `output/summary.csv` and the log, reports every field's change, flags negatives, and runs this case's checks. Standard library only. |
|  | `pipeline.sh` | Runs the steps above in order. |

### What the solver reads

| File | What it is |
|---|---|
| `input/geometry.dat` | The pore space. `preprocess.py` rebuilds it; this copy is here so the case runs before you have run anything. |
| `models/e_coli_core.xml.gz` | The metabolic model, unpacked into `input/` at start-up and checked against `models/manifest.txt`. |
| `output/ecoli.srg` | The fitted network. Written on the first run, loaded on every later one. Delete it to force a retrain. |

### What it inherits

Only the two shared kinetics headers, from
[`../../config/kinetics/`](../../config/kinetics/), and the solver sources.
`setup_case.sh` lays those down and copies this folder whole on top.

### The offline code this case ships

| File | What it is |
|---|---|
| `training/extractMM.py` | Exports an SBML or BiGG model to the flat format the GLPK path reads. Not needed here — the solver reads the bundled model itself — but carried so this folder is the whole procedure. |
| `models/NOTICE.md`, `models/manifest.txt` | Provenance, licence and checksum of the bundled model. |

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
