#!/usr/bin/env python3
"""Method 4, companion: what one evaluation does to the graph, drawn.

Same rule as the symbolic companion. The flow down the left has no colour, so
that every coloured object in the figure is part of the network being operated
on. The example is the shipped anaerobic methane oxidation case: four species,
one reaction, width six, two rounds.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Rectangle

from draw import Rec, to_pptx

YEL, BLU, GRN, RED = "#FFE699", "#BDD7EE", "#C6E0B4", "#F2B2B2"
PINK = "#F8A5A5"
EDGE, REDEDGE = "#9AA5B1", "#B04A4A"
TXT, GREY, ARR = "#202020", "#7F7F7F", "#595959"
PBG, PED, ROWED = "#FBFBFD", "#C9C9D2", "#DCDCE3"
RING = "#B04A4A"
STRIPE = {"setup": "#8FAADC", "each": "#9CC183", "out": "#E8C56A"}
SANS, MONO = "Carlito", "DejaVu Sans Mono"

W, H = 1900, 2630
fig = plt.figure(figsize=(W / 100.0, H / 100.0), dpi=100)
ax = Rec(fig.add_axes([0, 0, 1, 1]))
ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
ax.add_patch(Rectangle((0, 0), W, H, fc="white", ec="none", zorder=0))


def t(x, y, s, fs=13, c=TXT, f=SANS, w="normal", ha="left", va="center",
      it=False, z=7, lsp=1.4, rot=0):
    ax.text(x, y, s, fontsize=fs, color=c, family=f, fontweight=w, ha=ha, va=va,
            style="italic" if it else "normal", zorder=z, linespacing=lsp,
            rotation=rot)


def arrow(x0, y0, x1, y1, color=ARR, lw=1.6, z=6, ms=11):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=ms, color=color, lw=lw,
                                 shrinkA=0, shrinkB=0, zorder=z))


def lane(pts, color=ARR, lw=1.6, z=2, ms=11):
    for i in range(len(pts) - 2):
        ax.plot([pts[i][0], pts[i + 1][0]], [pts[i][1], pts[i + 1][1]],
                color=color, lw=lw, solid_capstyle="round", zorder=z)
    arrow(pts[-2][0], pts[-2][1], pts[-1][0], pts[-1][1], color, lw, z, ms)


def cell(x, y, w, h, fill, text, fs=11.5, ec=EDGE, z=6, f=MONO, r=4, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=%g" % r,
                                fc=fill, ec=ec, lw=1.1, zorder=z))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=TXT, family=f, fontweight="bold" if bold else "normal",
            zorder=z + 1)


def vector(cx, cy, vals, fill=GRN, cw=44, ch=24, fs=10.5):
    """A hidden vector, drawn as its cells."""
    x = cx - len(vals) * cw / 2.0
    for v in vals:
        cell(x, cy - ch / 2, cw, ch, fill, v, fs=fs)
        x += cw


SPECIES = ["CH4", "SO4", "HS", "HCO3"]
NU = ["-1", "-1", "+1", "+1"]


def network(cx, top, edges=True, arrows=None, ringr=False, sw=96, gap=150):
    """The AOM graph: four species nodes over one reaction node."""
    xs = [cx + (i - 1.5) * gap for i in range(4)]
    ry = top + 132
    for i, s in enumerate(SPECIES):
        cell(xs[i] - sw / 2, top - 13, sw, 26, BLU, s, fs=11.5)
    cell(cx - 58, ry - 15, 116, 30, YEL, "R1", fs=12, bold=True)
    if ringr:
        ax.add_patch(FancyBboxPatch((cx - 70, ry - 27), 140, 54,
                                    boxstyle="round,pad=0,rounding_size=10",
                                    fc="none", ec=RING, lw=1.7, ls=(0, (5, 4)),
                                    zorder=9))
    if edges:
        for i in range(4):
            ax.plot([xs[i], cx], [top + 13, ry - 15], color=EDGE, lw=1.1,
                    zorder=4)
            t((xs[i] + cx) / 2, (top + 13 + ry - 15) / 2, NU[i], fs=11,
              c=GREY, f=MONO, ha="center")
    if arrows == "up":
        for i in range(4):
            arrow(xs[i], top + 15, cx + (xs[i] - cx) * 0.12, ry - 18,
                  color=RING, lw=1.4)
    if arrows == "down":
        for i in range(4):
            arrow(cx + (xs[i] - cx) * 0.12, ry - 18, xs[i], top + 15,
                  color=RING, lw=1.4)
    return xs, ry


# ------------------------------------------------------------------ nodes ---
FX, FW = 470, 580
PX, PW = 810, 1050
L1, L2 = 140, 62


def node(cy, h, text, phase=None, pill=False, fs=15, sub=None, num=None):
    x0 = FX - FW / 2
    r = h / 2 if pill else 8
    ax.add_patch(FancyBboxPatch((x0, cy - h / 2), FW, h,
                                boxstyle="round,pad=0,rounding_size=%g" % r,
                                fc=PINK if pill else "white",
                                ec=REDEDGE if pill else EDGE, lw=1.3, zorder=5))
    if not pill and phase:
        ax.add_patch(Rectangle((x0 + 2, cy - h / 2 + 2), 8, h - 4,
                               fc=STRIPE[phase], ec="none", zorder=6))
    if num is not None:
        ax.add_patch(FancyBboxPatch((x0 + 20, cy - 15), 30, 30,
                                    boxstyle="round,pad=0,rounding_size=15",
                                    fc="#4A4A55", ec="none", zorder=7))
        ax.text(x0 + 35, cy, str(num), ha="center", va="center", fontsize=13,
                color="white", family=SANS, fontweight="bold", zorder=8)
    ax.text(FX + 28 if num is not None else FX, cy, text, ha="center",
            va="center", fontsize=fs, color=TXT, family=SANS, fontweight="bold",
            linespacing=1.3, zorder=7)
    if sub:
        t(FX, cy + h / 2 + 12, sub, fs=12, c=GREY, it=True, ha="center", va="top")
    return dict(cx=FX, cy=cy, w=FW, h=h)


def dia(cy, h, text, w=520, fs=14):
    ax.add_patch(Polygon([(FX, cy - h / 2), (FX + w / 2, cy), (FX, cy + h / 2),
                          (FX - w / 2, cy)], closed=True, fc="white", ec=EDGE,
                         lw=1.3, zorder=5))
    ax.text(FX, cy, text, ha="center", va="center", fontsize=fs, color=TXT,
            family=SANS, linespacing=1.3, zorder=7)
    return dict(cx=FX, cy=cy, w=w, h=h)


def panel(y0, h, title, tie=None, x=PX, w=PW):
    ax.add_patch(FancyBboxPatch((x, y0), w, h,
                                boxstyle="round,pad=0,rounding_size=9",
                                fc=PBG, ec=PED, lw=1.2, zorder=1))
    t(x + 22, y0 + 24, title, fs=13, c=GREY, w="bold", z=3)
    if tie is not None:
        ax.plot([tie["cx"] + tie["w"] / 2, x], [tie["cy"], y0 + h / 2],
                color="#C0C0CA", lw=1.2, ls=(0, (5, 4)), zorder=1)
    return x + 22, y0 + 44


def down(a, b, label=None, side=1):
    y0, y1 = a["cy"] + a["h"] / 2, b["cy"] - b["h"] / 2
    arrow(FX, y0, FX, y1)
    if label:
        t(FX + 12 * side, (y0 + y1) / 2, label, fs=12.5, c=GREY,
          ha="left" if side > 0 else "right")


# ------------------------------------------------------------------ header --
t(40, 48, "Inside the graph network: what one evaluation does, drawn", fs=25, w="bold")
t(40, 88,
  "The flow on the left has no colour on purpose. All the colour belongs to the network on the right, because that is what the evaluation moves through.",
  fs=15, c=GREY, it=True, va="top")
t(40, 138, "the network is drawn as:", fs=13.5, w="bold")
kx = 280
for fl, lab in ((BLU, "a species node, one per substrate"),
                (YEL, "a reaction node, one per reaction"),
                (GRN, "a number the fit chose")):
    ax.add_patch(FancyBboxPatch((kx, 126), 32, 24,
                                boxstyle="round,pad=0,rounding_size=4",
                                fc=fl, ec=EDGE, lw=1.0, zorder=5))
    t(kx + 42, 138, lab, fs=13)
    kx += 360

# ============================================ the one idea: it is a graph ====
n0 = node(238, 54, "START:  this voxel's four concentrations", pill=True)
px, py = panel(180, 268, "the reaction is the wiring: an edge exists only where the stoichiometry has a number", tie=n0)
network(px + 330, py + 34)
t(px + 640, py + 24,
  "The shipped case is anaerobic methane\n"
  "oxidation: CH4 + SO4 → HS + HCO3.\n"
  "One reaction, four species.\n\n"
  "The numbers on the edges are the\n"
  "stoichiometry and they are not fitted.\n"
  "They decide which messages can exist\n"
  "at all, which is why a rate this network\n"
  "returns is a set of reaction rates by\n"
  "construction rather than by a penalty\n"
  "term in a loss.",
  fs=13, va="top")

# ============================================================== 1 clamp =====
n1 = node(516, 66, "Pull each concentration inside\nthe training box", "each", num=1)
px, py = panel(470, 166, "the box is the range the fit was shown, and it is enforced rather than advised", tie=n1)
for i, (s, c, lo, hi, out) in enumerate((
        ("CH4", "1.2e-3", "2.1e-5", "5.0e-3", False),
        ("SO4", "9.4e-3", "5.2e-5", "7.9e-3", True),
        ("HS", "4.0e-4", "1.2e-6", "2.0e-3", False),
        ("HCO3", "8.0e-4", "1.1e-4", "4.0e-3", False))):
    yy = py + 6 + i * 26
    t(px + 16, yy, s, fs=11.5, f=MONO)
    t(px + 96, yy, c, fs=11.5, f=MONO, c=RING if out else TXT)
    t(px + 190, yy, "in  " + lo + " .. " + hi, fs=11.5, f=MONO, c=GREY)
    if out:
        t(px + 440, yy, "→  pulled back to 7.9e-3", fs=11.5, f=MONO, c=RING)
t(px + 16, py + 120,
  "One value moved, so this evaluation is counted once. Not once per species: a four-species network\n"
  "counting each value would report four hundred per cent of its evaluations as out of range.",
  fs=12.5, c=GREY, it=True, va="top")

# ============================================================ 2 encode ======
n2 = node(742, 66, "Scale, and give every species node\na vector of its own", "each", num=2)
px, py = panel(690, 200, "one scaling per species, fitted offline, and then the same six weights at every node", tie=n2)
cell(px + 30, py + 30, 110, 28, "white", "1.2e-3", fs=11.5)
arrow(px + 152, py + 44, px + 212, py + 44)
t(px + 182, py + 30, "scale", fs=11, c=GREY, it=True, ha="center", va="bottom")
cell(px + 228, py + 30, 90, 28, "white", "-0.21", fs=11.5)
arrow(px + 330, py + 44, px + 390, py + 44)
t(px + 360, py + 30, "tanh", fs=11, c=GREY, it=True, ha="center", va="bottom")
vector(px + 600, py + 44, ["0.41", "-0.08", "0.77", "0.12", "-0.55", "0.30"])
t(px + 600, py + 76, "the CH4 node", fs=11.5, c=GREY, it=True, ha="center")
t(px + 30, py + 106,
  "Six numbers, because the file says width 6. Every species node gets six, every reaction node gets six,\n"
  "and the same encoding weights are used for all of them. What tells the nodes apart is what goes in.",
  fs=12.5, va="top")

# ===================================================== 3 and 4 messages =====
n3 = node(1020, 66, "The species speak to the reactions", "each", num=3)
px, py = panel(962, 254, "a reaction adds up the nodes of its own species, each weighted by its stoichiometric number", tie=n3)
network(px + 320, py + 34, arrows="up")
t(px + 640, py + 30,
  "The sum runs over the species that have\n"
  "an edge. A zero in the stoichiometry is\n"
  "not a multiply by zero: there is no edge,\n"
  "so there is no message at all.\n\n"
  "The reaction node then takes tanh of that\n"
  "sum, its own bias, and a term in the\n"
  "Damkohler number the file carries.",
  fs=13, va="top")

n4 = node(1296, 66, "The reactions speak back", "each", num=4)
px, py = panel(1238, 254, "a species adds up the reactions it takes part in, and keeps a little of what it already was", tie=n4)
network(px + 320, py + 34, arrows="down")
t(px + 640, py + 30,
  "Every species node becomes tanh of two\n"
  "sums: one over the reactions that touch\n"
  "it, one over what it already was.\n\n"
  "Steps 3 and 4 together are one round.\n"
  "The shipped file asks for two, which is\n"
  "what lets a species feel a species it\n"
  "shares no reaction with.",
  fs=13, va="top")

d5 = dia(1470, 96, "rounds left?", w=380)
down(n4, d5)
lane([(FX - d5["w"] / 2, d5["cy"]), (L1, d5["cy"]), (L1, n3["cy"]),
      (FX - FW / 2, n3["cy"])])
t(FX - d5["w"] / 2 - 12, d5["cy"] - 17, "yes", fs=12.5, c=GREY, ha="right")
ax.text(L1 + 12, (d5["cy"] + n3["cy"]) / 2, "the next round", rotation=90,
        ha="center", va="center", fontsize=12.5, color=GREY, family=SANS, zorder=6)

# ============================================================ 5 readout =====
n5 = node(1638, 66, "Read one extent off the reaction node", "out", num=5)
px, py = panel(1580, 300, "one number per reaction, not one per species, and this is the choice the whole method turns on", tie=n5)
network(px + 250, py + 6, ringr=True)
arrow(px + 250, py + 156, px + 250, py + 186)
cell(px + 180, py + 186, 140, 34, GRN, "3.1e-3", fs=12.5)
t(px + 250, py + 232, "the extent", fs=11.5, c=GREY, it=True, ha="center")
t(px + 560, py + 24,
  "The other readout takes a number off each species\n"
  "node instead. It is the default only because a file\n"
  "with no readout line means it.\n\n"
  "The difference is checkable. Chemistry fixes sulfide\n"
  "made per methane eaten at exactly one. Read off the\n"
  "reaction node it is one at every input; read off the\n"
  "species nodes it runs from 0.76 to 1.11.",
  fs=13, va="top")

# ============================================================ 6 rates =======
n6 = node(1970, 66, "Build the four rates from\nthe stoichiometry", "out", num=6)
px, py = panel(1912, 236, "this multiplication is the reason the ratios come out exact", tie=n6)
cell(px + 30, py + 40, 140, 34, GRN, "3.1e-3", fs=12.5)
t(px + 186, py + 57, "×", fs=16, ha="center")
for i, (s, nu) in enumerate(zip(SPECIES, NU)):
    cell(px + 220, py + 6 + i * 30, 84, 26, YEL, nu, fs=11.5)
    t(px + 316, py + 19 + i * 30, s, fs=11.5, f=MONO, c=GREY)
t(px + 400, py + 57, "=", fs=16, ha="center")
for i, (s, v) in enumerate(zip(SPECIES, ["-3.1e-3", "-3.1e-3", "+3.1e-3", "+3.1e-3"])):
    cell(px + 430, py + 6 + i * 30, 130, 26, "white", v, fs=11.5)
    t(px + 572, py + 19 + i * 30, s + "  per hour", fs=11.5, f=MONO, c=GREY)
t(px + 30, py + 130,
  "Growth is read separately, off the average of the four species nodes, and it keeps its own scaling.\n"
  "Then every one of these five numbers is divided by 3600, because the file says per hour and the\n"
  "solver works per second.",
  fs=12.5, va="top")

n7 = node(2222, 62, "Divide by 3600, apply the gate,\nmultiply growth by biomass", "each",
          sub="the same three lines every fitted rate path ends with")
n8 = node(2350, 58, "Hand the increments to transport, floored", "each")
d9 = dia(2462, 96, "more voxels?", w=380)
stop = node(2574, 54, "STOP:  back to the solve loop of Figure 1", pill=True)

for a, b in ((n0, n1), (n1, n2), (n2, n3), (n3, n4)):
    down(a, b)
down(d5, n5, "no")
down(n5, n6); down(n6, n7); down(n7, n8); down(n8, d9)
down(d9, stop, "no")
lane([(FX - d9["w"] / 2, d9["cy"]), (L2, d9["cy"]), (L2, n1["cy"]),
      (FX - FW / 2, n1["cy"])])
t(FX - d9["w"] / 2 - 12, d9["cy"] - 17, "yes", fs=12.5, c=GREY, ha="right")
ax.text(L2 - 12, (d9["cy"] + n1["cy"]) / 2, "the next voxel", rotation=90,
        ha="center", va="center", fontsize=12.5, color=GREY, family=SANS, zorder=6)

fig.savefig("/home/claude/fig3/m4b.png", dpi=100, facecolor="white")
to_pptx(ax, W, H, "/home/claude/fig3/m4b.pptx", fig=fig)
print("wrote m4b.png and m4b.pptx", W, H)
