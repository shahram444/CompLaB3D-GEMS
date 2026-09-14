#!/usr/bin/env python3
"""
Figure 3 for Method 3 (symbolic rate laws), drawn in the same visual language as
Figure 1: the same five fills, the same thin grey box borders, the same legend
block top right, the same dashed side panels.

Where Figure 1 answers "where in the solve loop does this run", this answers
"what does the code actually do, in order, from a table of numbers to a change
in one voxel".
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from draw import Rec, to_pptx
from matplotlib.patches import FancyBboxPatch, Polygon, FancyArrowPatch, Rectangle
from matplotlib.path import Path
from matplotlib.patches import PathPatch

# ----------------------------------------------------------------- palette --
YEL, BLU, GRN, RED = "#FFE699", "#BDD7EE", "#C6E0B4", "#F2B2B2"
EDGE, REDEDGE = "#9AA5B1", "#B04A4A"
TXT, GREY, ARR = "#202020", "#7F7F7F", "#595959"
PBG, PED, ROWED = "#FBFBFD", "#C9C9D2", "#D9D9E0"
PINK = "#F8A5A5"

SANS = "Carlito"
MONO = "DejaVu Sans Mono"

W = 1900
H = 2320

fig = plt.figure(figsize=(W / 100.0, H / 100.0), dpi=100)
ax = Rec(fig.add_axes([0, 0, 1, 1]))
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.axis("off")
ax.add_patch(Rectangle((0, 0), W, H, fc="white", ec="none", zorder=0))


def box(cx, cy, w, h, fill, text, sub=None, ec=EDGE, lw=1.2, fs=15, bold=True,
        tc=TXT, z=3):
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=7",
                                fc=fill, ec=ec, lw=lw, zorder=z))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, color=tc,
            family=SANS, fontweight="bold" if bold else "normal",
            linespacing=1.35, zorder=z + 1)
    if sub:
        note(cx, cy + h / 2 + 16, sub)
    return dict(cx=cx, cy=cy, w=w, h=h)


def note(cx, y, text, color=GREY):
    ax.text(cx, y, text, ha="center", va="top", fontsize=12.5, color=color,
            family=SANS, style="italic", linespacing=1.35, zorder=6,
            bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="none"))


def tag(x, y, text):
    ax.text(x, y, text, ha="left", va="bottom", fontsize=12.5, color=GREY,
            family=MONO, zorder=6)


def diamond(cx, cy, w, h, text, fs=14):
    ax.add_patch(Polygon([(cx, cy - h / 2), (cx + w / 2, cy),
                          (cx, cy + h / 2), (cx - w / 2, cy)],
                         closed=True, fc="white", ec=EDGE, lw=1.2, zorder=3))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, color=TXT,
            family=SANS, linespacing=1.3, zorder=4)
    return dict(cx=cx, cy=cy, w=w, h=h)


def arrow(x0, y0, x1, y1, color=ARR, lw=1.6, z=2, style="-|>", ms=11):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=ms, color=color, lw=lw,
                                 shrinkA=0, shrinkB=0, zorder=z))


def elbow(pts, color=ARR, lw=1.6, z=2, ms=11):
    """Orthogonal polyline with an arrow head on the last segment."""
    for i in range(len(pts) - 2):
        ax.plot([pts[i][0], pts[i + 1][0]], [pts[i][1], pts[i + 1][1]],
                color=color, lw=lw, solid_capstyle="round", zorder=z)
    arrow(pts[-2][0], pts[-2][1], pts[-1][0], pts[-1][1], color, lw, z, ms=ms)


def down(a, b, label=None):
    y0 = a["cy"] + a["h"] / 2
    y1 = b["cy"] - b["h"] / 2
    arrow(a["cx"], y0, b["cx"], y1)
    if label:
        ax.text(a["cx"] + 12, (y0 + y1) / 2, label, ha="left", va="center",
                fontsize=12.5, color=GREY, family=SANS)


# ------------------------------------------------------------------ header --
ax.text(40, 52, "The symbolic path as code: from a table of numbers to a change in one voxel",
        fontsize=25, fontweight="bold", color=TXT, family=SANS, va="center")
ax.text(40, 100,
        "Read straight down. The yellow band runs once, before any simulation exists. The blue band runs once when the solver starts. "
        "The green band runs\nin every wet voxel at every reaction step, and is the box marked in red in Figure 1.",
        fontsize=15, color=GREY, family=SANS, style="italic", va="top",
        linespacing=1.45)

# ------------------------------------------------------------------ legend --
LX, LW = 1010, 840
ly = 200
ax.text(LX, ly - 26, "what the colours mean", fontsize=15, fontweight="bold",
        color=TXT, family=SANS, va="center")
for fill, lab, ec in (
        (YEL, "offline, once, outside CompLaB3D", EDGE),
        (BLU, "start-up, once per run", EDGE),
        (GRN, "per voxel, per reaction step", EDGE),
        (RED, "a guard rail: it stops the run, or changes the number", REDEDGE),
        ("white", "a question the code asks itself", EDGE)):
    ax.add_patch(FancyBboxPatch((LX, ly), LW, 46,
                                boxstyle="round,pad=0,rounding_size=6",
                                fc=fill, ec=ec, lw=1.2, zorder=3))
    ax.text(LX + LW / 2, ly + 23, lab, ha="center", va="center", fontsize=14.5,
            color=TXT, family=SANS, style="italic", zorder=4)
    ly += 58

# ------------------------------------------------------------- band labels --
def band(y0, y1, color, name):
    ax.add_patch(Rectangle((44, y0), 9, y1 - y0, fc=color, ec="none", zorder=2))
    ax.text(30, (y0 + y1) / 2, name, rotation=90, ha="center", va="center",
            fontsize=15, fontweight="bold", color=TXT, family=SANS)


CX = 560          # centre of the main column
BW = 700          # box width

# =============================================================== OFFLINE ====
y = 300
a1 = box(CX, y, BW, 72, YEL,
         "A table of data: one row per state, the concentrations\nin that state, and the one rate you want reproduced",
         sub="make_training_data.py writes the example table, and the law it drew from beside it,\nso this case can ask whether the search found the RIGHT law and not merely one that fits")

y += 176
a2 = box(CX, y, BW, 96, YEL,
         "Breed a population of expression trees out of\n+  −  ×  ÷  and your variable names, fit the numbers inside\neach one by least squares, and score it against the data")

y += 136
a3 = diamond(CX, y, 330, 96, "generations left?")

y += 140
a4 = box(CX, y, BW, 96, YEL,
         "Keep the best tree at each length. What comes back is that\nlist, not one answer, and you pick one length off it. --verify\nruns the same search twice to say which lengths are stable")

y += 152
a6 = box(CX, y, BW, 110, YEL,
         "Write the .sym file: the fitted rate line, the reaction,\nthe yield, the biomass name, one range line per variable,\nthe units line, and a header saying how it was made")

band(a1["cy"] - 36 - 46, a6["cy"] + 55 + 10, YEL, "offline")
tag(CX - BW / 2, a1["cy"] - 36 - 14, "training/make_training_data.py   then   tools/method_3_symbolic/fit_symbolic.py")

down(a1, a2)
down(a2, a3)
down(a3, a4, "no")
down(a4, a6)

# the generations loop, on the left
elbow([(CX - a3["w"] / 2, a3["cy"]), (CX - BW / 2 - 90, a3["cy"]),
       (CX - BW / 2 - 90, a2["cy"]), (CX - BW / 2, a2["cy"])])
ax.text(CX - a3["w"] / 2 - 14, a3["cy"] - 16, "yes", ha="right", va="center",
        fontsize=12.5, color=GREY, family=SANS)

# ------------------------------------------------------------- the handover -
hy0 = a6["cy"] + 55
hy1 = hy0 + 128
arrow(CX, hy0, CX, hy1, lw=2.4, ms=15)
ax.text(CX + 20, hy0 + 44, "the .sym file", ha="left", va="center",
        fontsize=15.5, fontweight="bold", color=TXT, family=SANS)
ax.text(CX + 20, hy0 + 76, "is the only thing that crosses this line",
        ha="left", va="center", fontsize=13.5, color=GREY, family=SANS,
        style="italic")
ax.plot([100, 1870], [hy0 + 60, hy0 + 60], color="#C9C9D2", lw=1.2,
        ls=(0, (6, 5)), zorder=1)

# =============================================================== START-UP ===
y = hy1 + 36 + 36
b1 = box(CX, y, BW, 72, BLU,
         "Read the file. Parse each rate line exactly once\ninto a small tree of nodes, held in one shared arena")

y += 140
b2 = diamond(CX, y, 560, 112,
             "every name known to this run,\nand no line using an output\ndefined below it?", fs=13.5)

y += 152
b3 = box(CX, y, BW, 110, BLU,
         "Bind each variable to its lattice by name, not by position.\nWork the substrate lines out from the reaction and the yield.\nStore the unit scale and the fitted range of each variable")

band(b1["cy"] - 36 - 46, b3["cy"] + 55 + 10, BLU, "start-up")
tag(CX - BW / 2, b1["cy"] - 36 - 14, "src/complab3d_symbolic.hh")

down(b1, b2)
down(b2, b3, "yes")

stop = box(1330, b2["cy"], 480, 76, RED,
           "Stop before the first step,\nnaming the line and the name",
           ec=REDEDGE)
arrow(b2["cx"] + b2["w"] / 2, b2["cy"], stop["cx"] - stop["w"] / 2, stop["cy"])
ax.text(b2["cx"] + b2["w"] / 2 + 14, b2["cy"] - 18, "no", ha="left",
        va="center", fontsize=12.5, color=GREY, family=SANS)
ax.text(1330, b2["cy"] + 56, "a name that does not exist is a typo, and a typo\n"
        "that reached the first step would be a silent zero",
        ha="center", va="top", fontsize=12.5, color=GREY, family=SANS,
        style="italic", linespacing=1.3)

# ============================================================== PER VOXEL ===
y = b3["cy"] + 55 + 74
c1 = box(CX, y, BW, 66, GRN,
         "Read this voxel's concentrations and its biomass")

y += 120
c2 = box(CX, y, BW, 72, RED,
         "Pull each input back inside the range the fit was shown,\nand count every time you have to", ec=REDEDGE)

y += 128
c3 = box(CX, y, BW, 72, GRN,
         "Walk each tree from the leaves up, in file order.\nA later line may read an earlier line's output")

y += 136
c4 = box(CX, y, BW, 96, GRN,
         "Apply the unit scale. Multiply the growth output by biomass and\nadd the substrate outputs up over organisms, then hand the\nincrements to transport, floored at what the voxel holds")

y += 130
c6 = diamond(CX, y, 330, 92, "more voxels?")

y += 112
ax.add_patch(FancyBboxPatch((CX - 260, y - 29), 520, 58,
                            boxstyle="round,pad=0,rounding_size=29",
                            fc=PINK, ec=REDEDGE, lw=1.2, zorder=3))
ax.text(CX, y, "back to the solve loop of Figure 1", ha="center", va="center",
        fontsize=14.5, fontweight="bold", color=TXT, family=SANS, zorder=4)
end = dict(cx=CX, cy=y, w=520, h=58)

band(c1["cy"] - 33 - 46, end["cy"] + 29 + 10, GRN, "every step")
tag(CX - BW / 2, c1["cy"] - 33 - 14, "src/complab3d_integration.hh")

down(c1, c2)
down(c2, c3)
down(c3, c4)
down(c4, c6)
down(c6, end, "no")

elbow([(CX - c6["w"] / 2, c6["cy"]), (CX - BW / 2 - 90, c6["cy"]),
       (CX - BW / 2 - 90, c1["cy"]), (CX - BW / 2, c1["cy"])])
ax.text(CX - c6["w"] / 2 - 14, c6["cy"] - 16, "yes", ha="right", va="center",
        fontsize=12.5, color=GREY, family=SANS)

# ------------------------------------------------------- the clamp reminder -
note(CX, c2["cy"] + 36 + 16,
     "the count is printed at the end of the run: a law used outside the range it was\n"
     "fitted over does not fail, it returns a confident number", color=REDEDGE)


# ---------------------------------------------------------------- panels ----
def panel(x, y0, w, title, italic, rows, anchor=None, accent=EDGE):
    line_h = 46
    hgt = 64 + (26 if italic else 0) + len(rows) * (line_h + 10) + 26
    ax.add_patch(FancyBboxPatch((x, y0), w, hgt,
                                boxstyle="round,pad=0,rounding_size=8",
                                fc=PBG, ec=PED, lw=1.3, ls=(0, (6, 5)),
                                zorder=2))
    yy = y0 + 34
    ax.text(x + 26, yy, title, ha="left", va="center", fontsize=15.5,
            fontweight="bold", color=TXT, family=SANS, zorder=4)
    yy += 30
    if italic:
        ax.text(x + 26, yy, italic, ha="left", va="top", fontsize=13,
                color=GREY, family=SANS, style="italic", linespacing=1.35,
                zorder=4)
        yy += 42
    for r in rows:
        ax.add_patch(FancyBboxPatch((x + 26, yy), w - 52, line_h,
                                    boxstyle="round,pad=0,rounding_size=5",
                                    fc="white", ec=ROWED, lw=1.0, zorder=3))
        ax.text(x + w / 2, yy + line_h / 2, r, ha="center", va="center",
                fontsize=13.5, color=TXT, family=SANS, zorder=4)
        yy += line_h + 10
    if anchor:
        ax.plot([anchor[0], x], [anchor[1], y0 + hgt / 2], color=PINK,
                lw=1.3, ls=(0, (6, 5)), zorder=1)
    return y0 + hgt


panel(
    1010, 570, 840,
    "inside one generation of the search",
    "Nobody tells it about Monod. The shape of the formula is what is\nbeing searched over, not just the numbers in a shape you chose.",
    ["take two trees that scored well",
     "swap a branch between them, or change one node",
     "fit that new tree's numbers by least squares",
     "score it: the error it leaves, plus a penalty for its length",
     "keep it only if it beats the best tree of its own length"],
    anchor=(CX + BW / 2, a2["cy"]))

panel(
    1010, c3["cy"] - 100, 840,
    "inside one voxel, when the path is a symbolic law",
    "No linear program, no network, no linear algebra. One pass over a\ntree that was already built before the first step ran.",
    ["the leaves are read: a number, or a concentration",
     "each node combines the two children below it",
     "growth is worked out first, because it is first in the file",
     "the substrate lines reuse that growth number, unscaled",
     "one growth rate and one source term per substrate come out"],
    anchor=(CX + BW / 2, c4["cy"]))

fig.savefig("/home/claude/fig3/figA.png", dpi=100, facecolor="white")
to_pptx(ax, W, H, fig=fig, path= "/home/claude/fig3/figA.pptx")
print("wrote fig3.png", W, H)
