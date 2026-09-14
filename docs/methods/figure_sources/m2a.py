#!/usr/bin/env python3
"""Method 2, the surrogate network path as code."""
from flow import Fig, YEL, BLU, GRN, RED, PINK, EDGE, REDEDGE, CX, BW, PANEL_X, PANEL_W

F = Fig(1900, 2900,
        "The surrogate path as code: from a sweep of linear programs to a change in one voxel",
        "Read straight down. The yellow band runs once, before any simulation exists, and it is the only place a linear program is ever solved. The blue band\n"
        "runs once when the solver starts. The green band runs in every wet voxel that holds biomass, at every reaction step.")

F.legend([(YEL, "offline, once, outside CompLaB3D", EDGE),
          (BLU, "start-up, once per run", EDGE),
          (GRN, "per voxel, per reaction step", EDGE),
          (RED, "a guard rail: it stops the run, or changes the number", REDEDGE),
          ("white", "a question the code asks itself", EDGE)])

# =============================================================== OFFLINE =====
a1 = F.box(300, 84, YEL,
           "Sweep the uptake bounds: one linear program per sample, on\n"
           "a grid or a Latin hypercube over one to three exchanges",
           sub="generateTrainingData.py. An infeasible solve is recorded as zero growth and kept rather\n"
               "than dropped, so the network learns where the feasible region ends instead of guessing")

a2 = F.box(474, 84, YEL,
           "Put every column on its own scale, then split the rows\n"
           "seventy, fifteen and fifteen",
           sub="the two held-back fifteens are the only thing standing in for regularisation: there is no\n"
               "weight decay, no dropout and no per-output weighting anywhere in the fit")

a3 = F.box(648, 84, YEL,
           "Fit four hidden layers of ten, five times over, from five\n"
           "different random starting points")

a4 = F.dia(792, 92, "restarts left?", w=330)

a5 = F.box(946, 110, YEL,
           "Keep the restart with the smallest error on the first fifteen\n"
           "per cent. Write the .srg file: the architecture, both scalings,\n"
           "the training box, and every weight")

F.band(a1["cy"] - 42 - 46, a5["cy"] + 55 + 74, YEL, "offline",
       "tools/method_2_surrogate/generateTrainingData.py   then   trainSurrogate.py", a1["cy"] - 42)
F.down(a1, a2); F.down(a2, a3); F.down(a3, a4); F.down(a4, a5, "no")
F.lane([(CX - a4["w"] / 2, a4["cy"]), (CX - BW / 2 - 90, a4["cy"]),
        (CX - BW / 2 - 90, a3["cy"]), (CX - BW / 2, a3["cy"])])
F.t(CX - a4["w"] / 2 - 14, a4["cy"] - 16, "yes", fs=12.5, c="#7F7F7F", ha="right")
F.note(CX, a5["cy"] + 55 + 14,
       "and 512 checkpoints beside it, so the written file can be shown to agree with the fit that\n"
       "produced it to one part in ten billion before anybody runs anything")

hy = F.handover(a5, "the .srg file", "is the only thing that crosses this line", gap=210)

# =============================================================== START-UP ===
b1 = F.box(hy + 72, 72, BLU,
           "Read the file, and check every layer carries as many\n"
           "numbers as its shape line promised")
b2 = F.dia(b1["cy"] + 158, 132,
           "is the training box complete, and is\nevery input name in\n<name_of_substrates>?",
           fs=13.5, w=580)
b3 = F.box(b2["cy"] + 176, 110, BLU,
           "Bind each input to its substrate by name. Say in the log that\n"
           "a file with flux outputs will have them read, checked and\n"
           "then dropped, because nothing maps them to substrates")

F.band(b1["cy"] - 36 - 46, b3["cy"] + 55 + 10, BLU, "start-up",
       "src/complab3d_surrogate.hh", b1["cy"] - 36)
F.down(b1, b2); F.down(b2, b3, "yes")

stop = F.box(b2["cy"], 96, RED,
             "Refuse to load the file.\nA network evaluated outside the box\n"
             "it was fitted over does not fail:\nit returns a confident number",
             ec=REDEDGE, cx=1370, w=500, fs=13.5)
F.arrow(CX + b2["w"] / 2, b2["cy"], 1370 - 250, b2["cy"])
F.t(CX + b2["w"] / 2 + 14, b2["cy"] - 18, "no", fs=12.5, c="#7F7F7F")

# ============================================================== PER VOXEL ===
c1 = F.box(b3["cy"] + 160, 66, GRN,
           "Read this voxel's concentrations and its biomass")
c2 = F.box(c1["cy"] + 154, 96, GRN,
           "Work out how much of each substrate this organism could\n"
           "take up here, by Monod, capped by what the voxel holds.\n"
           "That, and not the concentration, is what the network is asked",
           sub="the sweep varied bounds on exchange reactions, so the network's inputs are bounds, in\n"
               "millimoles per gram dry weight per hour. Handing it a concentration would be a category error")
c3 = F.box(c2["cy"] + 180, 72, RED,
           "Pull each input back inside the training box, and count the\n"
           "evaluation if any of them moved", ec=REDEDGE)
c4 = F.box(c3["cy"] + 176, 96, GRN,
           "Scale, then four layers of ten and one linear layer. Undo the\n"
           "output scaling, raise ten to it if the file says so, and floor\n"
           "anything under a hundred-millionth at exactly zero")
c5 = F.box(c4["cy"] + 148, 84, GRN,
           "Apply the thermodynamic factor. Growth multiplies biomass;\n"
           "an organism whose growth landed on the floor decays instead,\n"
           "at its own decay coefficient")
c6 = F.box(c5["cy"] + 142, 84, GRN,
           "Take the Monod uptake as the draw on each substrate, add it\n"
           "up over every organism in this voxel, and floor the total at\n"
           "what the voxel actually holds")
c7 = F.dia(c6["cy"] + 128, 92, "more voxels?", w=330)
end = F.box(c7["cy"] + 112, 58, PINK, "back to the solve loop of Figure 1",
            ec=REDEDGE, w=520, fs=14.5, r=29)

F.band(c1["cy"] - 33 - 46, end["cy"] + 29 + 10, GRN, "every step",
       "src/complab3d_processors_surrogate.hh", c1["cy"] - 33)
for a, b in ((c1, c2), (c2, c3), (c3, c4), (c4, c5), (c5, c6), (c6, c7)):
    F.down(a, b)
F.down(c7, end, "no")
F.loop(c7, c1["cy"], "the next voxel")

F.note(CX, c3["cy"] + 36 + 16,
       "once per evaluation, not once per input: a two-input network counting each value would report\n"
       "two hundred per cent of its evaluations as out of range", color=REDEDGE)

# ---------------------------------------------------------------- panels ----
F.panel(PANEL_X, 560, PANEL_W,
        "what one row of the sweep is",
        "This is the only place a linear program is solved on this path. Everything the\nnetwork knows about the cell, it knows from these rows.",
        ["set the exchange bounds this sample asks for",
         "solve the program with every other bound as the model has it",
         "record the growth, and the flux on every swept exchange",
         "an infeasible solve is a real answer: zero growth, kept",
         "a few thousand of these, and the cell is never opened again"],
        anchor=(CX + BW / 2, a1["cy"]))

F.panel(PANEL_X, c4["cy"] - 150, PANEL_W,
        "what the run-time file actually returns",
        "A file can carry flux outputs as well as growth, and on this path they are read,\nchecked for shape, and then dropped. Only growth is used.",
        ["the format names its outputs by exchange reaction id",
         "nothing in the run maps a reaction id to a substrate",
         "so the substrate draw comes from the Monod term above",
         "which means this path can consume but cannot excrete",
         "the compiled path can, at the cost of a hand-written mapping"],
        anchor=(CX + BW / 2, c4["cy"]))

F.save("/home/claude/fig3/m2a")
