#!/usr/bin/env python3
"""Method 1, the flux balance path as code."""
from flow import (Fig, YEL, BLU, GRN, RED, PINK, PUR, EDGE, REDEDGE, CX, BW,
                  PANEL_X, PANEL_W, GREY)

F = Fig(1900, 2760,
        "The flux balance path as code: from a metabolic model on disk to a change in one voxel",
        "Read straight down. The yellow band happens before any simulation exists. The blue band happens once when the solver starts, and it is where the\n"
        "cost of this method is paid down. The green band happens in every wet voxel that holds biomass, at every reaction step, and it is the box marked in\n"
        "red in Figure 3. The two back ends differ in one step out of nine, and that step is marked.")

F.legend([(YEL, "before the run, outside CompLaB3D", EDGE),
          (BLU, "start-up, once per run", EDGE),
          (GRN, "per organism, per voxel, per reaction step", EDGE),
          (PUR, "the one step the two back ends do differently", EDGE),
          (RED, "a guard rail: it stops the run, or changes the number", REDEDGE),
          ("white", "a question the code asks itself", EDGE)])

# =============================================================== OFFLINE =====
a1 = F.box(322, 84, YEL,
           "Have a metabolic model: a list of reactions, the metabolites\n"
           "each one moves, and which reaction makes biomass",
           sub="an SBML file. The shipped example is three metabolites and four reactions so that it can\n"
               "be checked by hand; a genome-scale model is a few thousand of each and changes nothing here")

a2 = F.box(506, 84, YEL,
           "Say which exchange reaction stands for which of the run's\n"
           "substrates, and how fast the organism can take each one up",
           sub="two lines of CompLaB.xml: <name_of_substrates> in the order the model lists the exchanges,\n"
               "and <fba_maximum_uptake_flux>, an enzyme capacity in mmol per gram dry weight per hour")

F.band(a1["cy"] - 42 - 46, a2["cy"] + 42 + 78, YEL, "before the run",
       "an SBML model, and two lines of CompLaB.xml", a1["cy"] - 42)
F.down(a1, a2)

hy = F.handover(a2, "the model, and the names",
                "nothing else crosses this line; the file is never read again", gap=200)

# =============================================================== START-UP ===
b1 = F.box(hy + 76, 84, BLU,
           "Read the model once. Turn the stoichiometry into arrays, and\n"
           "find the column of every exchange the run names, by name",
           sub="src/complab3d_metabolic.hh. Positional indexing from here on, which is why the order of\n"
               "<name_of_substrates> has to match the order the model lists its exchanges")

b2 = F.dia(b1["cy"] + 188, 132,
           "does every named substrate have a\ncolumn? does the biomass reaction\nexist? are the capacities set?",
           fs=13.5, w=600)

b3 = F.box(b2["cy"] + 190, 96, BLU,
           "Build ONE problem object and keep it for the whole run:\n"
           "the matrix, the fixed bounds and the objective are the same\n"
           "in every voxel, so they are set up once and never again",
           sub="the COBRApy route builds its model from these same arrays rather than re-reading the file,\n"
               "so a disagreement between the two back ends is a coupling defect and not a modelling one")

F.band(b1["cy"] - 42 - 46, b3["cy"] + 48 + 78, BLU, "start-up",
       "src/complab3d_metabolic.hh", b1["cy"] - 42)
F.down(b1, b2); F.down(b2, b3, "yes")

stop = F.box(b2["cy"], 84, RED,
             "Stop before the first step.\nA substrate with no column would be\n"
             "consumed by nothing, silently, for\nthe whole run",
             ec=REDEDGE, cx=1300, w=560, fs=13.5)
F.arrow(CX + b2["w"] / 2, b2["cy"], 1300 - 280, b2["cy"])
F.t(CX + b2["w"] / 2 + 14, b2["cy"] - 18, "no", fs=12.5, c=GREY)

# ============================================================== PER VOXEL ===
c1 = F.box(b3["cy"] + 220, 66, GRN,
           "Read this voxel's concentrations and its biomass")

c2 = F.box(c1["cy"] + 152, 96, GRN,
           "Build the only thing that tells the program where it is: the\n"
           "lower bound on each uptake reaction, from Monod, held above\n"
           "the model's own floor and above what the voxel can give up",
           sub="the matrix, the capacities and the objective do not move. One number per substrate is all\n"
               "that changes from one voxel to the next, and it is the whole of the coupling. The bound is\n"
               "negative because an exchange is written as the metabolite leaving, so uptake is a negative flux")

c3 = F.box(c2["cy"] + 178, 72, PUR,
           "Solve. GLPK does it in this process, warm-started from the\n"
           "basis the last voxel left. COBRApy does it across a language\n"
           "boundary, from cold, about ten times slower", fs=14.5)

c4 = F.dia(c3["cy"] + 166, 116,
           "would the organisms here between\nthem draw more of some substrate\nthan the voxel holds?",
           fs=13.5, w=580)

c5 = F.box(c4["cy"], 84, RED,
           "Mark every organism that wanted that substrate as sharing it,\n"
           "give them one pooled supply between them, and solve again.\n"
           "Three passes, then clamp the draw and carry on", ec=REDEDGE,
           cx=1300, w=620, fs=13.5)

c6 = F.box(c4["cy"] + 176, 72, GRN,
           "Apply the thermodynamic factor to the solved fluxes and to\n"
           "the growth, not to the bounds that produced them")

c7 = F.box(c6["cy"] + 156, 84, GRN,
           "Turn each flux into a change in concentration, over every\n"
           "organism present, and never let the total take a substrate\n"
           "below zero")

c8 = F.box(c7["cy"] + 150, 72, GRN,
           "Growth multiplies biomass. An organism the program could not\n"
           "solve for, or that it gave no growth, decays instead")

c9 = F.dia(c8["cy"] + 140, 92, "more voxels?", w=330)
end = F.box(c9["cy"] + 112, 58, PINK, "back to the solve loop of Figure 3",
            ec=REDEDGE, w=520, fs=14.5, r=29)

F.band(c1["cy"] - 33 - 46, end["cy"] + 29 + 10, GRN, "every step",
       "src/complab3d_processors_fba.hh", c1["cy"] - 33)
for a, b in ((c1, c2), (c2, c3), (c3, c4), (c4, c6), (c6, c7), (c7, c8), (c8, c9)):
    F.down(a, b)
F.down(c9, end, "no")
F.loop(c9, c1["cy"], "the next voxel")

# the repair branch: out to the right, and back into the solve
F.arrow(CX + c4["w"] / 2, c4["cy"], 1300 - 310, c4["cy"])
F.t(CX + c4["w"] / 2 + 14, c4["cy"] - 18, "yes", fs=12.5, c=GREY)
F.lane([(1300, c5["cy"] - 42), (1300, c3["cy"]), (CX + BW / 2, c3["cy"])])
F.t(1312, (c3["cy"] + c5["cy"] - 42) / 2, "at most three times", fs=12.5,
    c=GREY, ha="left", rot=90)
F.t(CX, c4["cy"] + 58 + 8, "no", fs=12.5, c=GREY, ha="center")

# ---------------------------------------------------------------- panels ----
F.panel(PANEL_X, 600, PANEL_W,
        "the program, and the one block of it that moves",
        "Three statements about the cell, and only the third knows which voxel it is in.\nThat is what makes one problem object serve a whole run.",
        ["S v = 0, the cell is at steady state inside",
         "l <= v <= u, every reaction has a capacity",
         "max c'v, the organism grows as fast as it can",
         "only l on the uptake columns is rebuilt per voxel"],
        anchor=(CX + BW / 2, a1["cy"]))

F.panel(PANEL_X, c7["cy"] - 50, PANEL_W,
        "why the answer is checked and sometimes solved again",
        "Each organism solves as though the voxel were its own. Two of them in one voxel\ncan therefore each be granted the whole supply.",
        ["every organism solves against the same concentrations",
         "so the draws are added up afterwards, not before",
         "an over-draw marks the organisms that wanted it",
         "they are re-solved against one pooled supply",
         "after three passes the draw is clamped, never refused"],
        anchor=(CX + BW / 2, c7["cy"]))

F.save("/home/claude/fig3/m1a")
