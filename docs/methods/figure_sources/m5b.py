#!/usr/bin/env python3
"""Method 5, companion: the four equations on two real voxels, with the factor drawn.

Every number here is computed from the constants in the shipped aom.thm at the
two compositions offline/thermo_curve.py sweeps between, so a reader can check
any line with a calculator.
"""
import math

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
RIM, INNER = "#2E6DA4", "#B04A4A"
SANS, MONO = "Carlito", "DejaVu Sans Mono"

# ----------------------------------------------------------- the chemistry --
R, T = 8.314462618e-3, 277.15
RT = R * T
dG0, dGloss, ATP, dGATP, CHI = -33.12, 25.7, 0.25, 50.0, 1.0
THRESH = ATP * dGATP
NU = [("CH4", -1), ("SO4", -1), ("HS", +1), ("HCO3", +1)]
VOX = [("rim", {"CH4": 1.0e-2, "SO4": 2.8e-2, "HS": 1.0e-5, "HCO3": 2.3e-3}, RIM),
       ("inner", {"CH4": 2.0e-3, "SO4": 2.0e-2, "HS": 1.2e-2, "HCO3": 1.45e-2}, INNER)]


def solve(c):
    lnQ = sum(n * math.log(c[s]) for s, n in NU)
    dG = dG0 + RT * lnQ + dGloss
    f = dG + THRESH
    x = f / (CHI * RT)
    F = 0.0 if x >= 0 else min(1.0, max(0.0, 1.0 - math.exp(x)))
    return lnQ, dG, f, x, F


W, H = 1900, 2360
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


def cell(x, y, w, h, fill, text, fs=11.5, ec=EDGE, z=6, f=MONO, r=4, tc=TXT):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=%g" % r,
                                fc=fill, ec=ec, lw=1.1, zorder=z))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, family=f, zorder=z + 1)


FX, FW, PX, PW = 470, 580, 810, 1050
STRIPE = {"in": "#8FAADC", "calc": "#9CC183", "out": "#E8C56A"}


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


def dia(cy, h, text, w=420, fs=14):
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


def down(a, b, label=None):
    arrow(FX, a["cy"] + a["h"] / 2, FX, b["cy"] - b["h"] / 2)
    if label:
        t(FX + 12, (a["cy"] + a["h"] / 2 + b["cy"] - b["h"] / 2) / 2, label,
          fs=12.5, c=GREY)


# ------------------------------------------------------------------ header --
t(40, 48, "Inside the gate: the four equations, on two voxels of one aggregate", fs=25, w="bold")
t(40, 88,
  "Every number below is computed from the constants in the shipped aom.thm at the two compositions the offline sweep runs between, so any line can be\n"
  "checked with a calculator. The rim voxel is drawn in blue and the interior voxel in red throughout.",
  fs=15, c=GREY, it=True, va="top")
for lab, col, x in (("the rim of the aggregate", RIM, 40), ("deep inside it", INNER, 500)):
    ax.add_patch(FancyBboxPatch((x, 140), 26, 20,
                                boxstyle="round,pad=0,rounding_size=4",
                                fc=col, ec="none", zorder=5))
    t(x + 38, 150, lab, fs=13.5, w="bold")

# ============================================================= START ========
n0 = node(238, 54, "START:  the rate its path just returned", pill=True)
px, py = panel(186, 128, "the gate never computes a rate; it only says what fraction of one the energy balance allows", tie=n0)
t(px + 10, py + 4,
  "The compiled kinetics of example 19 hand it a growth rate of 3.06e-6 per second at the rim.\n"
  "Everything below decides what that becomes. The same four equations sit in front of a linear\n"
  "program, a surrogate, a symbolic law and a graph network without changing a line.",
  fs=13, va="top")

# ============================================================= 1 lnQ ========
n1 = node(430, 66, "The reaction quotient, as a sum\nof logarithms", "calc", num=1)
px, py = panel(352, 312, "one term per species, the sign taken straight from the reaction, and every activity floored at cmin", tie=n1)
t(px + 300, py + 2, "at the rim", fs=12.5, w="bold", ha="center", z=5)
t(px + 660, py + 2, "deep inside", fs=12.5, w="bold", ha="center", z=5)
t(px + 10, py + 2, "nu   species", fs=12, c=GREY, f=MONO)
for i, (s, n) in enumerate(NU):
    yy = py + 30 + i * 30
    t(px + 10, yy, "%+d   %-5s" % (n, s), fs=12, f=MONO)
    for k, (_, c, col) in enumerate(VOX):
        xx = px + 300 + k * 360
        t(xx, yy, "%-9.3g  ->  %+8.4f" % (c[s], n * math.log(c[s])), fs=12,
          f=MONO, ha="center", c=col)
ax.plot([px + 120, px + 840], [py + 148, py + 148], color=EDGE, lw=1.1, zorder=4)
for k, (_, c, col) in enumerate(VOX):
    lnQ = solve(c)[0]
    t(px + 300 + k * 360, py + 168, "ln Q = %+8.4f" % lnQ, fs=13, f=MONO,
      w="bold", ha="center", c=col)
t(px + 10, py + 196,
  "A negative ln Q means the products are scarce and the reaction is far from equilibrium. Deep in the\n"
  "aggregate the sulfide it has already made has built up, so the same reaction is close to its own\n"
  "equilibrium and ln Q has changed sign. Nothing else about the two voxels differs.",
  fs=12.5, va="top")

# ============================================================= 2 dG =========
n2 = node(722, 66, "The free energy actually\navailable here", "calc", num=2)
px, py = panel(676, 250, "the standard value, plus what the local composition does to it, plus the transfer loss", tie=n2)
for k, (nm, c, col) in enumerate(VOX):
    lnQ, dG, f, x, Ft = solve(c)
    yy = py + 16 + k * 78
    t(px + 10, yy, "%-7s" % nm, fs=13, w="bold", c=col)
    t(px + 110, yy, "dG  =  %.2f  %+.4f  %+.1f  =  %+.4f  kJ/mol"
      % (dG0, RT * lnQ, dGloss, dG), fs=13, f=MONO, c=col)
    t(px + 110, yy + 26, "        dG0        R T ln Q       loss", fs=11.5,
      f=MONO, c=GREY)
t(px + 10, py + 152,
  "R T is 2.304 kJ/mol at 4 C, which is why the temperature line is mandatory: a free energy quoted at\n"
  "25 C is a different number. The 25.7 is the one tuned quantity in the whole case, fixed by the\n"
  "measured growth efficiency of these consortia and not by anything in the energy accounting.",
  fs=12.5, va="top")

# ============================================================= 3 threshold ==
n3 = node(996, 66, "What is left once the ATP\nhas been paid for", "calc", num=3)
px, py = panel(950, 204, "the organism has to make a quarter of an ATP for every turn of the reaction", tie=n3)
for k, (nm, c, col) in enumerate(VOX):
    lnQ, dG, f, x, Ft = solve(c)
    yy = py + 18 + k * 36
    t(px + 10, yy, "%-7s" % nm, fs=13, w="bold", c=col)
    t(px + 110, yy, "f  =  %+8.4f  +  %.2f x %.1f  =  %+8.4f  kJ/mol"
      % (dG, ATP, dGATP, f), fs=13, f=MONO, c=col)
t(px + 10, py + 100,
  "Negative means there is energy left over. Positive means the reaction cannot pay for the ATP and\n"
  "nothing runs. The threshold, 12.5 kJ/mol, sits inside the 10 to 20 band measured for methanogens\n"
  "and sulfate reducers, which is the check that it is a physical number and not a fitted one.",
  fs=12.5, va="top")

# ============================================================= 4 factor =====
n4 = node(1252, 66, "The factor", "out", num=4)
px, py = panel(1170, 480, "one minus e to the power of f over chi R T, taken as zero at the threshold and one far below it", tie=n4)

GX, GY, GW, GH = px + 60, py + 40, 560, 300
ax.add_patch(Rectangle((GX, GY), GW, GH, fc="white", ec=EDGE, lw=1.1, zorder=3))
LO, HI = -36.0, 2.0


def gx(v):
    return GX + (v - LO) / (HI - LO) * GW


def gy(v):
    return GY + GH - v * GH


for v in (0.0, 0.5, 1.0):
    ax.plot([GX, GX + GW], [gy(v), gy(v)], color="#E4E4EA", lw=1.0, zorder=3)
    t(GX - 10, gy(v), "%.1f" % v, fs=11.5, f=MONO, c=GREY, ha="right")
for v in (-30, -20, -12.5, 0):
    ax.plot([gx(v), gx(v)], [GY, GY + GH], color="#E4E4EA", lw=1.0, zorder=3)
    t(gx(v), GY + GH + 14, "%g" % v, fs=11.5, f=MONO, c=GREY, ha="center")
N = 44
xs, ys = [], []
for i in range(N + 1):
    g = LO + (HI - LO) * i / float(N)
    ff = g + THRESH
    val = 0.0 if ff >= 0 else min(1.0, max(0.0, 1.0 - math.exp(ff / (CHI * RT))))
    xs.append(gx(g)); ys.append(gy(val))
for i in range(N):
    ax.plot([xs[i], xs[i + 1]], [ys[i], ys[i + 1]], color="#4A4A55", lw=2.0,
            zorder=6)
ax.plot([gx(-THRESH), gx(-THRESH)], [GY, GY + GH], color=GREY, lw=1.3,
        ls=(0, (5, 4)), zorder=5)
t(gx(-THRESH) + 8, GY + 18, "the threshold,\n−12.5 kJ/mol", fs=11.5, c=GREY,
  it=True, va="top")
for nm, c, col in VOX:
    lnQ, dG, f, x, Ft = solve(c)
    ax.add_patch(FancyBboxPatch((gx(dG) - 7, gy(Ft) - 7), 14, 14,
                                boxstyle="round,pad=0,rounding_size=7",
                                fc=col, ec="white", lw=1.4, zorder=8))
    t(gx(dG) + (46 if Ft > 0.5 else -6), gy(Ft) + (28 if Ft > 0.5 else -28), "%s\nF = %.3f" % (nm, Ft),
      fs=11.5, c=col, w="bold", ha="center", va="center")
t(GX + GW / 2, GY + GH + 44, "the free energy available,  dG,  in kJ per mole",
  fs=12.5, c=GREY, ha="center")
t(GX - 44, GY + GH / 2, "the factor", fs=12.5, c=GREY, ha="center", rot=90)

t(px + 680, py + 16,
  "The curve is the whole of the method.\n\n"
  "Left of the threshold the reaction can pay\n"
  "for its ATP and the factor climbs to one.\n"
  "Right of it nothing runs at all.\n\n"
  "The climb is not a step. Its width is chi\n"
  "R T ln 2, here 1.6 kJ per mole, so the two\n"
  "voxels are not on and off: they are 0.999\n"
  "and 0.000, and everything between them\n"
  "takes a real distance to cross.\n\n"
  "That distance is set by the chemistry and\n"
  "the temperature. There is no parameter\n"
  "in it to tune.",
  fs=13, va="top")

# ============================================================= 5 multiply ===
n5 = node(1730, 66, "Multiply the rate, and the growth,\nby the factor", "out", num=5)
px, py = panel(1672, 200, "the metabolic answer is kept and only its size is changed", tie=n5)
for k, (nm, c, col) in enumerate(VOX):
    Ft = solve(c)[4]
    yy = py + 20 + k * 42
    t(px + 10, yy, "%-7s" % nm, fs=13, w="bold", c=col)
    t(px + 110, yy, "3.06e-6  x  %.5f  =  %.3e  per second" % (Ft, 3.06e-6 * Ft),
      fs=13, f=MONO, c=col)
t(px + 10, py + 110,
  "Scaling the answer rather than tightening a bound is deliberate. A linear program asked for a\n"
  "smaller uptake would redistribute its whole flux pattern and could return different by-products,\n"
  "which is a modelling change. This asks only what fraction of the same answer the energy allows.",
  fs=12.5, va="top")

n6 = node(1932, 62, "Count it, and carry the smallest and\nlargest factor seen anywhere", "calc",
          sub="the end of the run says so plainly if the gate never closed, or closed everywhere")
d7 = dia(2066, 100, "more gated organisms,\nmore voxels?", w=450, fs=13.5)
stop = node(2196, 54, "STOP:  back to the solve loop of Figure 1", pill=True)

for a, b in ((n0, n1), (n1, n2), (n2, n3), (n3, n4), (n4, n5), (n5, n6), (n6, d7)):
    down(a, b)
down(d7, stop, "no")
ax.plot([80, 80], [d7["cy"], n0["cy"] + 40], color=ARR, lw=1.6, zorder=2)
ax.plot([80, FX - d7["w"] / 2], [d7["cy"], d7["cy"]], color=ARR, lw=1.6, zorder=2)
arrow(80, n0["cy"] + 44, 80, n0["cy"] + 8)
ax.plot([80, FX - FW / 2], [n0["cy"] + 8, n0["cy"] + 8], color=ARR, lw=1.6, zorder=2)
t(FX - d7["w"] / 2 - 12, d7["cy"] - 17, "yes", fs=12.5, c=GREY, ha="right")
ax.text(68, (d7["cy"] + n0["cy"]) / 2, "the next one", rotation=90, ha="center",
        va="center", fontsize=12.5, color=GREY, family=SANS, zorder=6)

fig.savefig("/home/claude/fig3/m5b.png", dpi=100, facecolor="white")
to_pptx(ax, W, H, "/home/claude/fig3/m5b.pptx", fig=fig)
print("wrote m5b.png and m5b.pptx", W, H)
