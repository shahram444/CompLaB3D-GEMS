# The .h5 viewer

A window for looking at the two kinds of `.h5` this project produces, on
Windows, macOS or Linux. It is written for someone who wants to see what is in a
file without writing any code.

```
Windows            double click  h5_viewer.bat
macOS, Linux       ./h5_viewer.sh
any system         python3 prt_h5_viewer.py            then File, Open
                   python3 prt_h5_viewer.py file.h5    open that file at once
```

The first run installs `numpy`, `h5py` and `matplotlib` if they are missing.
Nothing else is needed, and nothing is installed system wide.

---

## The two kinds of file

| | Written by | What is in it |
|---|---|---|
| **campaign dataset** | `collect_to_h5.py`, on the CompLaB3D side | the pore spaces, and for every simulation the concentration, reaction rate and velocity volumes at every stored time, plus the exact input record of each run |
| **results** | `make_results_h5.py`, on the PRT3D side | how accurate every fitted model is on the training, validation and test sets, the simulated and predicted volumes of the held out snapshots, the design of the campaign, the settings each model was fitted with, the fitting history, the weights, and the notes |

The viewer works out which one it has from the file itself, so there is nothing
to choose. Anything else still opens, because the Browse tab is a plain HDF5
tree and works on any file.

## The panels

Every panel is offered for both kinds of file. A panel whose data a particular
file does not carry says so in place of its picture rather than vanishing, so
the window looks the same whatever was opened.

| Panel | What it shows |
|---|---|
| **Overview** | what the file is, how big it is, and every heading it carries |
| **Pore spaces** | each pore structure in slices or in 3D: the grains, the geodesic distance, the straight line distance, with porosity and tortuosity |
| **Simulations** | the simulated chemistry. For a dataset, any run at any stored time, as concentration, reaction rate, velocity magnitude or a time series. For a results file, the simulated volume of any held out snapshot |
| **Run inputs** | for a dataset, the literal input record of a chosen run. For a results file, every simulation in the campaign with the two numbers that defined it and which side of the split its pore space fell on |
| **Predictions** | prediction beside simulation for any held out snapshot, as a comparison, in 3D, or as the difference between them |
| **Settings** | every input value the run or the fit was given, with a sentence per value saying what it does. See below |
| **Notes** | a plain text file of your own, read from beside the `.h5` |
| **Browse** | the raw HDF5 tree, for anything the other panels do not cover |

## The Settings panel

This is the panel for a reader who wants to know what a run was told to do
without opening any code. Each row is one setting, the value this file was made
with, and a sentence saying what changing it would change about the science.
Settings are grouped by what they affect rather than by the order the program
happens to declare them in.

For a results file it lists what each model was fitted with, grouped into what
is being predicted, how the data was divided, how the fitting was done, and the
optional features that were switched off. Underneath it adds what the fitting
actually did: how many passes were made, which pass was kept, and why the run
stopped.

For a campaign dataset it lists what one simulation was given, grouped into the
two numbers that define a case, the block, the flow, the chemistry, time, and
housekeeping.

Nothing is dropped. Anything in the file the glossary does not cover is printed
at the end, so the table is always the complete record.

## The Notes panel

The Notes panel shows no notes of its own. It looks beside the `.h5` for a text
file and prints whatever is in it, so what appears there is written by whoever
is using the file rather than by whoever made it. Any of these work:

```
<the same name as the .h5>.txt
<the same name as the .h5>_notes.txt
notes.txt        or  NOTES.txt       beside the .h5
notes.txt        or  NOTES.txt       beside this viewer
```

## Nothing is loaded until it is asked for

A 586 MB dataset opens instantly, because only the attributes and the small
index arrays are read when a file is opened. A volume is read one snapshot at a
time, 32 kB at a time, when a control changes. That is why the viewer is usable
on a laptop on files far larger than its memory.

## If something goes wrong

| What you see | What it means |
|---|---|
| `h5py is required` | the install step was skipped. Run `python3 -m pip install -r requirements.txt` |
| `No module named tkinter` | this Python was built without it. On Debian or Ubuntu `sudo apt install python3-tk`; on macOS with Homebrew `brew install python-tk` |
| the window opens but a panel says "not in this file" | that file does not carry that data. The Overview panel lists what it does carry |
| a tab is missing | only Predictions is ever absent, and only for a campaign dataset, which has no model predictions in it |

---

## Where each file is made

This repository writes the **campaign dataset**. `collect_to_h5.py` under
`ABC_sweep_package/tools/` gathers the output of a whole campaign into one
`.h5`, and `verify_dataset.py` beside it checks that what came out is complete
and self consistent before anything is done with it.

The **results** file is written on the other side, by `make_results_h5.py` in
the PRT-DeepONet-3D repository. The viewer here opens both, so a results file
handed over by a collaborator can be read without installing anything else.
