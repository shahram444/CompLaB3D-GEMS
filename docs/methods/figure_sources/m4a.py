#!/usr/bin/env python3
"""Method 4, the graph network path as code."""
from flow import *          # noqa: F401,F403
from flow import Fig, YEL, BLU, GRN, RED, PINK, EDGE, REDEDGE, CX, BW, PANEL_X, PANEL_W

F = Fig(1900, 2630,
        "The graph network path as code: from a table of states to a change in one voxel",
        "Read straight down. The yellow band runs once, before any simulation exists. The blue band runs once when the solver starts. The green band runs\n"
        "in every wet voxel that holds biomass, at every reaction step, and is the box marked in red in Figure 1.")

F.legend([(YEL, "offline, once, outside CompLaB3D", EDGE),
          (BLU, "start-up, once per run", EDGE),
          (GRN, "per voxel, per reaction step", EDGE),
          (RED, "a guard rail: it stops the run, or changes the number", REDEDGE),
          ("white", "a question the code asks itself", EDGE)])

# =============================================================== OFFLINE =====
a1 = F.box(300, 84, YEL,
           "A stoichiometry file, and a table of states: one row per\n"
           "composition, the rate of every species in it, and the growth rate",
           sub="training/make_training_data.py writes the example table from a law you can read,\n"
               "so the fit can be asked whether it found the right answer and not merely a close one")

a2 = F.box(474, 96, YEL,
           "Put every concentration on one scale, and project the\n"
           "species rates onto one extent per reaction",
           sub="the projection is what lets the readout be exact later: a rate vector that is not a set\n"
               "of reaction rates cannot be represented at all, which is the point")

a3 = F.box(650, 96, YEL,
           "Fit the weights: plain mean squared error, Adam,\n"
           "gradients by finite difference on 24 weights per array")

a4 = F.dia(796, 92, "epochs left?", w=330)

a5 = F.box(950, 110, YEL,
           "Write the .gnn file: the stoichiometry, both scalings, the\n"
           "training box, the readout mode, and every fitted weight")

F.band(a1["cy"] - 42 - 46, a5["cy"] + 55 + 10, YEL, "offline",
       "training/make_training_data.py   then   tools/method_4_graphnet/train_graphnet.py", a1["cy"] - 42)
F.down(a1, a2); F.down(a2, a3); F.down(a3, a4); F.down(a4, a5, "no")
F.lane([(CX - a4["w"] / 2, a4["cy"]), (CX - BW / 2 - 90, a4["cy"]),
        (CX - BW / 2 - 90, a3["cy"]), (CX - BW / 2, a3["cy"])])
F.t(CX - a4["w"] / 2 - 14, a4["cy"] - 16, "yes", fs=12.5, c="#7F7F7F", ha="right")

hy = F.handover(a5, "the .gnn file", "is the only thing that crosses this line")

# =============================================================== START-UP ===
b1 = F.box(hy + 72, 72, BLU,
           "Read the file keyword by keyword, and check every array is\n"
           "the shape its header line promised")
b2 = F.dia(b1["cy"] + 150, 124,
           "is every species name in\n<name_of_substrates>, and is\nevery array the right size?", fs=13.5, w=580)
b3 = F.box(b2["cy"] + 160, 110, BLU,
           "Bind each species to its lattice by name, not by position.\n"
           "Store the training box, the unit scale and the readout mode")

F.band(b1["cy"] - 36 - 46, b3["cy"] + 55 + 10, BLU, "start-up",
       "src/complab3d_graphnet.hh", b1["cy"] - 36)
F.down(b1, b2); F.down(b2, b3, "yes")

stop = F.box(b2["cy"], 76, RED,
             "Stop before the first step,\nnaming the species and listing\nwhat this run does have",
             ec=REDEDGE, cx=1330, w=480, fs=14)
F.arrow(CX + b2["w"] / 2, b2["cy"], 1330 - 240, b2["cy"])
F.t(CX + b2["w"] / 2 + 14, b2["cy"] - 18, "no", fs=12.5, c="#7F7F7F")
F.t(1330, b2["cy"] + 58, "a network trained on CH4 cannot be bound to a run that calls it methane,\n"
                         "and binding it to the wrong lattice would be worse than not binding it",
    fs=12.5, c="#7F7F7F", it=True, ha="center", va="top")

# ============================================================== PER VOXEL ===
c1 = F.box(b3["cy"] + 130, 66, GRN,
           "Read this voxel's concentrations and its biomass")
c2 = F.box(c1["cy"] + 120, 72, RED,
           "Pull each concentration back inside the training box, and\n"
           "count the evaluation if any of them moved", ec=REDEDGE)
c3 = F.box(c2["cy"] + 148, 72, GRN,
           "Encode: scale each concentration, and give its species node\n"
           "a short vector of its own")
c4 = F.box(c3["cy"] + 128, 72, GRN,
           "Pass messages: the species speak to the reactions,\n"
           "the reactions speak back. Twice.")
c5 = F.box(c4["cy"] + 128, 72, GRN,
           "Read one extent off each reaction node, then build every\n"
           "species rate from the stoichiometry")
c6 = F.box(c5["cy"] + 136, 96, GRN,
           "Apply the unit scale and the thermodynamic factor. Growth\n"
           "multiplies biomass, and the increments go to transport,\n"
           "floored at what the voxel actually holds")
c7 = F.dia(c6["cy"] + 130, 92, "more voxels?", w=330)
end = F.box(c7["cy"] + 112, 58, PINK, "back to the solve loop of Figure 1",
            ec=REDEDGE, w=520, fs=14.5, r=29)

F.band(c1["cy"] - 33 - 46, end["cy"] + 29 + 10, GRN, "every step",
       "src/complab3d_processors_graphnet.hh", c1["cy"] - 33)
for a, b in ((c1, c2), (c2, c3), (c3, c4), (c4, c5), (c5, c6), (c6, c7)):
    F.down(a, b)
F.down(c7, end, "no")
F.loop(c7, c1["cy"], "the next voxel")

F.note(CX, c2["cy"] + 36 + 16,
       "the count is printed at the end of the run, and a run with more than five per cent of its\n"
       "evaluations clamped is told to widen the ranges and fit again", color=REDEDGE)

# ---------------------------------------------------------------- panels ----
F.panel(PANEL_X, 570, PANEL_W,
        "inside one epoch of the fit",
        "There is no held-out set. The best training loss ever seen is what is kept,\nso the honest check is fresh samples afterwards, not the fit itself.",
        ["take two dozen weights at random from each array",
         "nudge each one up and then down by a ten-thousandth",
         "the difference between the two losses is its gradient",
         "Adam moves every weight a little way downhill",
         "five epochs with no gain: halve the step and rewind"],
        anchor=(CX + BW / 2, a3["cy"]))

F.panel(PANEL_X, c3["cy"] - 60, PANEL_W,
        "inside one voxel, when the path is a graph network",
        "No linear program. The stoichiometry is not a penalty term in a loss: it is\nthe wiring, and it decides which messages exist at all.",
        ["each species node starts as its own scaled concentration",
         "a reaction adds up the species nodes that feed it",
         "a species adds up the reactions it takes part in",
         "the same weights at every node and in every voxel",
         "one extent per reaction comes out, not one rate per species"],
        anchor=(CX + BW / 2, c4["cy"]))

F.save("/home/claude/fig3/m4a")
