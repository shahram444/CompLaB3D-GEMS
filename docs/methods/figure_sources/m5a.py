#!/usr/bin/env python3
"""Method 5, the thermodynamic gate as code."""
from flow import Fig, YEL, BLU, GRN, RED, PINK, EDGE, REDEDGE, CX, BW, PANEL_X, PANEL_W

F = Fig(1900, 2480,
        "The thermodynamic gate as code: from a file of energies to a rate that is allowed to run",
        "Read straight down. The gate makes no rate of its own. It multiplies the rate whichever path the organism already uses has just returned, which is\n"
        "why the green band sits inside the box marked in red in Figure 1 rather than beside it.")

F.legend([(YEL, "before the run, and none of it is required", EDGE),
          (BLU, "start-up, once per run", EDGE),
          (GRN, "per gated organism, per voxel, per reaction step", EDGE),
          (RED, "a guard rail: it stops the run, or changes the number", REDEDGE),
          ("white", "a question the code asks itself", EDGE)])

# =============================================================== OFFLINE =====
a1 = F.box(300, 96, YEL,
           "Write the .thm file by hand: the reaction and its signs, the\n"
           "standard free energy, the ATP the organism must pay for a\n"
           "turn of it, and which organism the block belongs to",
           sub="every number comes from the literature at the temperature of the case, which is why the\n"
               "temperature line is mandatory and why quoting a free energy measured at 25 C is a real error")

a2 = F.box(494, 84, YEL,
           "Optional: sweep the compositions this run will actually\n"
           "visit, and print the factor at each one",
           sub="offline/thermo_curve.py. Needed for nothing, and worth running anyway: a gate open\n"
               "everywhere and a gate shut everywhere are both invisible in the output afterwards")

F.band(a1["cy"] - 48 - 46, a2["cy"] + 42 + 74, YEL, "before the run",
       "offline/thermo_curve.py", a1["cy"] - 48)
F.down(a1, a2)

hy = F.handover(a2, "the .thm file", "is what the solver reads; the sweep only reads it too", gap=170)

# =============================================================== START-UP ===
b1 = F.box(hy + 72, 72, BLU,
           "Read the file, then scale every energy to kJ once the whole\n"
           "file has been read, never as the lines arrive")
b2 = F.dia(b1["cy"] + 164, 148,
           "is every species of the reaction in\n<name_of_substrates>? is chi above zero,\n"
           "and cmin above zero? does any organism\ncarry two blocks?", fs=13, w=640)
b3 = F.box(b2["cy"] + 172, 96, BLU,
           "Bind each species to its lattice by name. Build the map from\n"
           "organism to block, so an organism with no block costs one\n"
           "array lookup and returns a factor of one")

F.band(b1["cy"] - 36 - 46, b3["cy"] + 48 + 10, BLU, "start-up",
       "src/complab3d_thermo.hh", b1["cy"] - 36)
F.down(b1, b2); F.down(b2, b3, "all fine")

stop = F.box(b2["cy"], 96, RED,
             "Stop before the first step.\nTwo blocks on one organism would\n"
             "multiply, and a cmin of zero lets\nln Q reach minus infinity",
             ec=REDEDGE, cx=1370, w=500, fs=13.5)
F.arrow(CX + b2["w"] / 2, b2["cy"], 1370 - 250, b2["cy"])
F.t(CX + b2["w"] / 2 + 14, b2["cy"] - 18, "no", fs=12.5, c="#7F7F7F")

# ============================================================== PER VOXEL ===
c1 = F.box(b3["cy"] + 140, 84, GRN,
           "The organism's own rate path returns its potential rate:\n"
           "compiled kinetics, a linear program, a surrogate, a symbolic\n"
           "law or a graph network. The gate does not care which")
c2 = F.box(c1["cy"] + 144, 72, RED,
           "Read the same concentrations that path read, and floor\n"
           "each one at cmin before taking its logarithm", ec=REDEDGE)
c3 = F.box(c2["cy"] + 156, 96, GRN,
           "The reaction quotient, then the free energy actually\n"
           "available here, then what is left of it once the ATP\n"
           "the organism must make has been paid for")
c4 = F.box(c3["cy"] + 146, 72, GRN,
           "The factor: one minus e to the power of that, taken as zero\n"
           "at or past the threshold and as one far below it")
c5 = F.box(c4["cy"] + 130, 72, GRN,
           "Multiply every substrate rate and the growth rate by it.\n"
           "Nothing else about the metabolic answer is touched")
c6 = F.dia(c5["cy"] + 130, 100, "more gated organisms\nin this voxel?", w=430, fs=13.5)
c7 = F.dia(c6["cy"] + 134, 92, "more voxels?", w=330)
end = F.box(c7["cy"] + 112, 58, PINK, "back to the solve loop of Figure 1",
            ec=REDEDGE, w=520, fs=14.5, r=29)

F.band(c1["cy"] - 42 - 46, end["cy"] + 29 + 10, GRN, "every step",
       "inside each rate path's own processor", c1["cy"] - 42)
for a, b in ((c1, c2), (c2, c3), (c3, c4), (c4, c5), (c5, c6), (c6, c7)):
    F.down(a, b)
F.down(c7, end, "no")
F.loop(c6, c1["cy"], "the next organism")

F.note(CX, c2["cy"] + 36 + 16,
       "the floor is written so that it catches a value that is not a number as well, which is the only\n"
       "place in this path where that can arrive", color=REDEDGE)

# ---------------------------------------------------------------- panels ----
F.panel(PANEL_X, 524, PANEL_W,
        "what the file has to say, and why each line is there",
        "Only two of these have defaults worth trusting. The rest are measurements,\nand a wrong sign on the first of them shuts the gate everywhere.",
        ["temperature: the free energy must be quoted at it",
         "the reaction, with consumed negative and produced positive",
         "dG0, and the ATP cost, which sets the threshold",
         "chi and cmin: the width of the transition, and the floor"],
        anchor=(CX + BW / 2, a1["cy"]))

F.panel(PANEL_X, c3["cy"] - 40, PANEL_W,
        "why it is a smooth factor and not a switch",
        "The width of the transition is chi R T ln 2, which for the shipped case is\n1.6 kJ per mole. A reaction front therefore has a thickness.",
        ["at the threshold the factor is exactly zero",
         "one width past it, one half",
         "several widths past it, indistinguishable from one",
         "so the rate falls away over a distance, not at an edge",
         "and the distance is set by chemistry, not by a parameter"],
        anchor=(CX + BW / 2, c4["cy"]))

F.save("/home/claude/fig3/m5a")
