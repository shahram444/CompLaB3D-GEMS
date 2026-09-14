#!/usr/bin/env python3
"""Method 1, the two refinements of the solve step, as code."""
from flow import (Fig, YEL, BLU, GRN, RED, PINK, PUR, ORG, EDGE, REDEDGE, CX,
                  BW, PANEL_X, PANEL_W, GREY)

F = Fig(1900, 2380,
        "The two refinements, as code: what each one does to the solve step, and what it costs",
        "Neither of these is a different method. Both sit inside the one purple box of Figure 4, and everything around that box is unchanged. The first replaces\n"
        "the single program with a chain of them. The second wraps whichever solve is configured and runs it once per carbon source. They can be used together,\n"
        "in which case each option of the second is itself a chain.")

F.legend([(GRN, "<multi_step>: it replaces the single solve", EDGE),
          (ORG, "<cybernetic>: it wraps whatever solve is configured", EDGE),
          (RED, "a guard rail: it stops the run, or changes the number", REDEDGE),
          ("white", "a question the code asks itself", EDGE),
          (PINK, "back into the per-voxel loop of Figure 4", REDEDGE)])

# ================================================================= CHAIN =====
a1 = F.box(330, 72, GRN,
           "Solve the first program: maximise growth, at this voxel's\n"
           "bounds. Call what it reaches f1")

a2 = F.box(a1["cy"] + 164, 96, GRN,
           "Write a floor on the growth column at a fraction of f1.\n"
           "It is a bound on one column, not a new row, so the matrix\n"
           "never changes and the warm start survives",
           sub="a general constraint on a weighted sum would need a row, would resize the basis, and would\n"
               "throw the warm start away in every voxel. Every objective here is a single column, so it need not")

a3 = F.box(a2["cy"] + 186, 72, RED,
           "A floor on a column the caller re-bounds every voxel goes\n"
           "into the caller's bound vectors, not into the problem object",
           sub="the next stage resets every column it was given from those vectors on entry, so a floor\n"
               "written straight into the problem would be silently undone before it could bind",
           ec=REDEDGE)

a4 = F.box(a3["cy"] + 176, 84, GRN,
           "Solve the next program: maximise, or minimise, the next\n"
           "quantity, under every floor set so far. Each stage pins one\n"
           "more number, so after the last one nothing is free")

a5 = F.dia(a4["cy"] + 150, 92, "more stages?", w=330)

a6 = F.box(a5["cy"] + 136, 84, RED,
           "Read the growth out of the LAST stage's flux vector, never\n"
           "out of the first stage's objective value. Then undo every\n"
           "floor and restore the objective", ec=REDEDGE)

F.band(a1["cy"] - 36 - 46, a6["cy"] + 42 + 12, GRN, "the chain",
       "src/complab3d_lexicographic.hh", a1["cy"] - 36)
F.down(a1, a2); F.down(a2, a3); F.down(a3, a4); F.down(a4, a5)
F.down(a5, a6, "no")
F.lane([(CX - a5["w"] / 2, a5["cy"]), (CX - BW / 2 - 90, a5["cy"]),
        (CX - BW / 2 - 90, a4["cy"]), (CX - BW / 2, a4["cy"])])
F.t(CX - a5["w"] / 2 - 14, a5["cy"] - 16, "yes", fs=12.5, c=GREY, ha="right")
F.t(CX - BW / 2 - 79, (a4["cy"] + a5["cy"]) / 2, "the next stage", fs=12.5,
    c=GREY, ha="center", rot=90)

hy = F.handover(a6, "one determined flux vector",
                "the same model at the same bounds now returns the same numbers, on any solver", gap=210)

# =========================================================== COMPETITION =====
b1 = F.box(hy + 74, 84, ORG,
           "Score every carbon source: an ordinary Monod rate on that\n"
           "source alone, multiplied by the number of carbon atoms the\n"
           "molecule carries. Divide through, so the scores sum to one",
           sub="the weights depend on concentration and on nothing else. No enzyme level is carried between\n"
               "steps, so the organism switches as fast as the concentrations move and there is no lag")

b2 = F.dia(b1["cy"] + 176, 100, "is this source's weight\nabove the floor?",
           w=430, fs=13.5)

b3 = F.box(b2["cy"] + 166, 84, ORG,
           "Solve the whole model with only this source's uptake open.\n"
           "Whatever solve is configured does the work: one program, or\n"
           "the chain above, or a trained network")

b4 = F.box(b2["cy"], 100, RED,
           "Shutting a source forbids its UPTAKE\nand never its release. An organism\n"
           "living on one source has to stay free\nto excrete what a later option eats",
           ec=REDEDGE, cx=1300, w=560, fs=13.5)

b5 = F.dia(b3["cy"] + 156, 92, "more sources?", w=350)

b6 = F.box(b5["cy"] + 138, 84, ORG,
           "Blend the flux VECTORS with the weights: the growth and\n"
           "every exchange at once. Each one satisfies the steady state,\n"
           "so their weighted sum does too, exactly",
           sub="averaging model outputs is not usually safe. It is safe here because the thing being averaged\n"
               "is a set of vectors in the null space of one matrix, and that space is closed under the sum")

end = F.box(b6["cy"] + 168, 58, PINK, "on to the substrate budget of Figure 4",
            ec=REDEDGE, w=560, fs=14.5, r=29)

F.band(b1["cy"] - 42 - 46, end["cy"] + 29 + 10, ORG, "the competition",
       "src/complab3d_cybernetic.hh", b1["cy"] - 42)
F.down(b1, b2); F.down(b2, b3, "yes"); F.down(b3, b5); F.down(b5, b6, "no")
F.down(b6, end)
F.lane([(CX - b5["w"] / 2, b5["cy"]), (CX - BW / 2 - 90, b5["cy"]),
        (CX - BW / 2 - 90, b2["cy"]), (CX - BW / 2, b2["cy"])])
F.t(CX - b5["w"] / 2 - 14, b5["cy"] - 16, "yes", fs=12.5, c=GREY, ha="right")
F.t(CX - BW / 2 - 79, (b2["cy"] + b5["cy"]) / 2, "the next source", fs=12.5,
    c=GREY, ha="center", rot=90)
F.arrow(CX + b2["w"] / 2, b2["cy"], 1300 - 280, b2["cy"])
F.t(CX + b2["w"] / 2 + 14, b2["cy"] - 18, "no, skip it", fs=12.5, c=GREY)

# ---------------------------------------------------------------- panels ----
F.panel(PANEL_X, 560, PANEL_W,
        "the chain, written out",
        "Three stages on one organism. Each line is one linear program, and the part after\nthe comma is the floor the stages before it left behind.",
        ["LP 1   max growth                        -> f1",
         "LP 2   max product 1,   growth >= a1 f1  -> f2",
         "LP 3   max product 2,   growth >= a1 f1,",
         "                        product 1 >= a2 f2"],
        anchor=(CX + BW / 2, a2["cy"]))

F.panel(PANEL_X, 920, PANEL_W,
        "what the retain fraction buys, and what it costs",
        "At one, the chain gives nothing away and adds no parameter. Below one it trades\ngrowth for release, which is what real cells do and what the model will not.",
        ["a = 1: the first stage and the last agree exactly",
         "a = 0.8 on the small network: 3.6364 against 2.9091",
         "report the first and biomass outruns the fluxes by a quarter",
         "a below one is a fitted parameter; say so in the methods"],
        anchor=(CX + BW / 2, a6["cy"]))

F.panel(PANEL_X, b3["cy"] - 90, PANEL_W,
        "what the competition costs, and why it is less than it looks",
        "One program per carbon source per organism per voxel per step is the honest\nprice. Two things bring it down, and one removes it.",
        ["the options differ only in which bounds are open",
         "so each warm-starts from the basis the last one left",
         "an option below the weight floor is never solved",
         "a surrogate in place of the solver removes the cost",
         "the blend itself is m multiplies and an add"],
        anchor=(CX + BW / 2, b6["cy"]))

F.save("/home/claude/fig3/m1c")
