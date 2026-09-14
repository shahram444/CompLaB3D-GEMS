#!/usr/bin/env python3
"""Method 2, companion: one evaluation of the shipped network, on two real voxels.

Every number here is computed from the weights, the scalings and the training
box written into src/surrogateModel.hh, at two compositions example 11 actually
passes through, so a reader can check any line with a calculator.
"""
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Rectangle

from draw import Rec, to_pptx
from netcalc import IW, b1, LW2, b2, LW3, b3, LW4, b4, LW5, b5
from netcalc import xo, xg, xym, ygain, yym, yoff, LO, HI, layers, net, scale

YEL, BLU, GRN, RED = "#FFE699", "#BDD7EE", "#C6E0B4", "#F2B2B2"
PINK = "#F8A5A5"
EDGE, REDEDGE = "#9AA5B1", "#B04A4A"
TXT, GREY, ARR = "#202020", "#7F7F7F", "#595959"
PBG, PED, ROWED = "#FBFBFD", "#C9C9D2", "#DCDCE3"
NEAR, FAR = "#2E6DA4", "#B04A4A"
SANS, MONO = "Carlito", "DejaVu Sans Mono"

# --------------------------------------------------------- the shipped case --
KS = [0.05, 0.01]          # <half_saturation_constants>, mM
VMAX = [8.0, 0.4]          # <fba_maximum_uptake_flux>, mmol/gDW/h
SUB = ["acetate", "Fe3"]
DECAY = 0.01
VOX = [("fed", {"acetate": 1.00, "Fe3": 0.50}, NEAR),
       ("starved", {"acetate": 0.20, "Fe3": 0.05}, FAR)]


def bound(c):
    return [VMAX[i] * c[SUB[i]] / (KS[i] + c[SUB[i]]) for i in range(2)]


def growth(c):
    raw = net(bound(c))
    return raw, (0.0 if not (raw > 1e-8) else raw)


W, H = 1900, 3350
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
t(40, 48, "Inside the surrogate: one evaluation of the shipped network, on two voxels of one run",
  fs=25, w="bold")
t(40, 88,
  "Every number below is computed from the weights, the two scalings and the training box written into src/surrogateModel.hh, at two compositions the\n"
  "shipped case passes through, so any line can be checked with a calculator. The fed voxel is drawn in blue and the starved one in red throughout.",
  fs=15, c=GREY, it=True, va="top")
for lab, col, x in (("acetate 1.00 mM, Fe3 0.50 mM", NEAR, 40),
                    ("acetate 0.20 mM, Fe3 0.05 mM", FAR, 560)):
    ax.add_patch(FancyBboxPatch((x, 140), 26, 20,
                                boxstyle="round,pad=0,rounding_size=4",
                                fc=col, ec="none", zorder=5))
    t(x + 38, 150, lab, fs=13.5, w="bold")

# ============================================================= START ========
n0 = node(238, 54, "START:  one wet voxel that holds biomass", pill=True)
px, py = panel(186, 120, "the network is asked once per organism per voxel per reaction step, and nothing else on this path is", tie=n0)
t(px + 10, py + 4,
  "No linear program is solved here or anywhere else in the run. Everything this network knows about\n"
  "the cell it learned offline, from a few thousand programs solved before the simulation existed.",
  fs=13, va="top")

# ============================================================= 1 Monod ======
n1 = node(530, 72, "Build the input: the uptake bound\nMonod allows here", "in", num=1)
px, py = panel(340, 380,
               "this is the step readers get wrong: the network is handed a BOUND in mmol/gDW/h, not a concentration in mM", tie=n1)
t(px + 10, py + 2, "substrate   Ks     Vmax", fs=12, c=GREY, f=MONO)
t(px + 400, py + 2, "the fed voxel", fs=12.5, w="bold", ha="center", c=NEAR)
t(px + 760, py + 2, "the starved voxel", fs=12.5, w="bold", ha="center", c=FAR)
for i, s in enumerate(SUB):
    yy = py + 34 + i * 34
    t(px + 10, yy, "%-9s  %-5.2f  %-4.1f" % (s, KS[i], VMAX[i]), fs=12, f=MONO)
    for k, (_, c, col) in enumerate(VOX):
        t(px + 400 + k * 360, yy,
          "%.2f -> %7.4f" % (c[s], bound(c)[i]), fs=12, f=MONO, ha="center", c=col)
t(px + 10, py + 136,
  "v = Vmax C / (Ks + C), in mmol per gram dry weight per hour. The bound is then held under whatever\n"
  "the voxel could actually give up in one step, so a voxel with almost nothing left cannot be asked\n"
  "for a rate it cannot supply. Vmax comes from the metabolic solver's own maximum uptake flux, which\n"
  "is the same field both linear-program paths read, so switching an organism between them does not\n"
  "silently change its limits.",
  fs=12.5, va="top")
t(px + 10, py + 246,
  "The sweep that made this network varied bounds on exchange reactions, so bounds are the units it\n"
  "was fitted in. Handing it 1.00 and 0.50 instead of 7.6190 and 0.3922 would be a category error, and\n"
  "nothing in the code could catch it: both pairs are inside the box and both return a number.",
  fs=12.5, va="top", c=REDEDGE)

# ============================================================= 2 box ========
n2 = node(914, 72, "Hold each input inside the box\nthe fit covered, and count it", "calc", num=2)
px, py = panel(754, 320,
               "the box is written in the file beside the weights, because a network outside it does not fail, it answers", tie=n2)
t(px + 10, py + 2, "input      lowest fitted     highest fitted", fs=12, c=GREY, f=MONO)
for i, s in enumerate(SUB):
    t(px + 10, py + 32 + i * 28, "%-9s  %-15.8g  %-.8g" % (s, LO[i], HI[i]),
      fs=12, f=MONO)
for k, (nm, c, col) in enumerate(VOX):
    t(px + 10, py + 104 + k * 28,
      "%-9s %7.4f  %7.4f   both inside: nothing moved" % (nm, bound(c)[0], bound(c)[1]),
      fs=12, f=MONO, c=col)
t(px + 10, py + 170,
  "One evaluation that moved anything counts once, not once per input: a two-input network counting\n"
  "each value would report two hundred per cent of its evaluations as out of range. The end of the run\n"
  "prints the percentage, and says plainly when it is large enough that the sweep should be widened\n"
  "and the network fitted again.",
  fs=12.5, va="top")

# ============================================================= 3 scale ======
n3 = node(1236, 72, "Put both inputs on the scale\nthe network was trained in", "calc", num=3)
px, py = panel(1108, 256,
               "mapminmax, the same map MATLAB applies: every input on minus one to plus one over the fitted range", tie=n3)
for k, (nm, c, col) in enumerate(VOX):
    b = bound(c); sc = scale(b)
    t(px + 10, py + 12 + k * 56, "%-9s" % nm, fs=13, w="bold", c=col)
    t(px + 140, py + 12 + k * 56,
      "(%7.4f − %.6f) x %.6f − 1  =  %+.6f" % (b[0], xo[0], xg[0], sc[0]),
      fs=12.5, f=MONO, c=col)
    t(px + 140, py + 38 + k * 56,
      "(%7.4f − %.6f) x %.6f − 1  =  %+.6f" % (b[1], xo[1], xg[1], sc[1]),
      fs=12.5, f=MONO, c=col)
t(px + 10, py + 126,
  "The offsets and gains are the two lines that also define the box above: lowest fitted value is the\n"
  "offset, and highest is the offset plus two over the gain. A file that carried the weights without\n"
  "them would be unusable, which is why the loader refuses one.",
  fs=12.5, va="top")

# ============================================================= 4 network ====
n4 = node(1698, 72, "Four layers of ten, then one\nlinear output", "calc", num=4)
px, py = panel(1398, 640,
               "every circle is one unit, shaded by what it actually returns at this voxel: blue positive, red negative", tie=n4)

SIZES = [2, 10, 10, 10, 10, 1]
NW, NH = 370, 300


def drawnet(x0, y0, c, col, label):
    L = layers(bound(c))
    cols = len(SIZES)
    xs = [x0 + i * (NW / float(cols - 1)) for i in range(cols)]
    pos = []
    for i, n in enumerate(SIZES):
        step = NH / 10.0
        top = y0 + (NH - (n - 1) * step) / 2.0
        pos.append([(xs[i], top + j * step) for j in range(n)])
    for i in range(cols - 1):
        for (ax0, ay0) in pos[i]:
            for (ax1, ay1) in pos[i + 1]:
                ax.plot([ax0, ax1], [ay0, ay1], color="#E6E6EC", lw=0.5, zorder=2)
    for i, n in enumerate(SIZES):
        for j in range(n):
            v = float(L[i][j])
            v = max(-1.0, min(1.0, v))
            fc = ("#BDD7EE" if v >= 0 else "#F2B2B2")
            a = 0.18 + 0.82 * abs(v)
            ax.add_patch(FancyBboxPatch((pos[i][j][0] - 9, pos[i][j][1] - 9), 18, 18,
                                        boxstyle="round,pad=0,rounding_size=9",
                                        fc=fc, ec="#8A8A96", lw=0.8, alpha=a,
                                        zorder=4))
            sat = abs(v) > 0.95
            ax.add_patch(FancyBboxPatch((pos[i][j][0] - 9, pos[i][j][1] - 9), 18, 18,
                                        boxstyle="round,pad=0,rounding_size=9",
                                        fc="none",
                                        ec="#2F2F38" if sat else "#8A8A96",
                                        lw=2.2 if sat else 0.8, zorder=5))
    t(x0 + NW / 2, y0 - 22, label, fs=13, w="bold", ha="center", c=col)
    for i, lab in enumerate(["in", "1", "2", "3", "4", "out"]):
        t(xs[i], y0 + NH + 20, lab, fs=11.5, c=GREY, f=MONO, ha="center")
    t(xs[-1] + 16, pos[-1][0][1], "%+.4f" % float(L[-1][0]), fs=12.5, f=MONO,
      w="bold", c=col)
    ns = sum(1 for i in range(1, 5) for v in L[i] if abs(v) > 0.95)
    t(x0 + NW / 2, y0 + NH + 44,
      "%d of the 40 hidden units saturated, %d of them in the last layer"
      % (ns, sum(1 for v in L[4] if abs(v) > 0.95)),
      fs=11.5, c=GREY, it=True, ha="center")
    return L


LA = drawnet(px + 30, py + 50, VOX[0][1], NEAR, "the fed voxel")
LB = drawnet(px + 545, py + 50, VOX[1][1], FAR, "the starved voxel")
t(px + 10, py + 418,
  "Each unit adds up the units before it with its own weights, adds its own bias, and puts the total\n"
  "through the same S-shaped curve, which flattens at minus one and plus one. A unit drawn with a heavy\n"
  "ring has reached one of those flats: it is saturated, and can no longer tell one input from another.\n"
  "Both voxels have some. The difference is where they are. In the starved voxel six of the last layer's\n"
  "ten units are pinned, against three in the fed one, and the output unit, which does no curve at all,\n"
  "lands at %+.4f: exactly the bottom of the range the fit was scaled into. The network is not saying\n"
  "this voxel grows slowly. It is saying the answer is off the bottom of what it was ever shown."
  % float(LB[-1][0]), fs=12.5, va="top")

# ============================================================= 5 back =======
n5 = node(2192, 72, "Undo the output scaling, and floor\nwhat comes out", "out", num=5)
px, py = panel(2032, 320,
               "the reverse of the same map, then one comparison that is doing more work than it looks", tie=n5)
for k, (nm, c, col) in enumerate(VOX):
    raw, mu = growth(c)
    de = (float(layers(bound(c))[-1][0]) - yym) / ygain + yoff
    t(px + 10, py + 14 + k * 62, "%-9s" % nm, fs=13, w="bold", c=col)
    t(px + 140, py + 14 + k * 62,
      "(%+.6f + 1) / %.6f  =  %+.3e" % (float(layers(bound(c))[-1][0]), ygain, de),
      fs=12.5, f=MONO, c=col)
    t(px + 140, py + 40 + k * 62,
      "%s  ->  growth = %s per hour" % (
          "above the floor" if mu > 0 else "not above the floor",
          "%.5f" % mu if mu > 0 else "exactly 0"),
      fs=12.5, f=MONO, c=col)
t(px + 10, py + 150,
  "The floor is written as \"not greater than a hundred-millionth\", which also catches a value that is\n"
  "not a number. The starved voxel is why it is there: its answer de-scales to a small negative rate,\n"
  "which is not a slow cell but an arithmetic leftover, and passing it on would grow biomass backwards.\n"
  "If the file says the output was fitted in logarithms, ten is raised to it first. The shipped network\n"
  "was not, so that line does nothing here. Flux outputs, if the file carries any, are checked and dropped.",
  fs=12.5, va="top")

# ============================================================= 6 use it =====
n6 = node(2531, 72, "Grow, or decay, and take the\nsubstrate the bound asked for", "out", num=6)
px, py = panel(2386, 290,
               "growth and uptake are decided separately on this path, which is what a growth-only network forces", tie=n6)
for k, (nm, c, col) in enumerate(VOX):
    raw, mu = growth(c)
    t(px + 10, py + 14 + k * 58, "%-9s" % nm, fs=13, w="bold", c=col)
    t(px + 140, py + 14 + k * 58,
      ("growth %.5f  ->  biomass increases" % mu) if mu > 0
      else ("growth 0  ->  biomass decays at %.2f" % DECAY),
      fs=12.5, f=MONO, c=col)
    t(px + 140, py + 40 + k * 58,
      "draw on acetate and Fe3  =  %7.4f  %7.4f" % tuple(bound(c)),
      fs=12.5, f=MONO, c=col)
t(px + 10, py + 140,
  "The draw is the bound from step 1, unchanged: the network returned growth and nothing that names a\n"
  "substrate, so the solver keeps its own Monod estimate as the uptake. Summed over every organism in\n"
  "the voxel and then held under what the voxel holds, once, so two organisms cannot each take all of\n"
  "it. This is also why the path can consume but cannot excrete.",
  fs=12.5, va="top")

d7 = dia(2780, 100, "more organisms here,\nmore voxels?", w=450, fs=13.5)
stop = node(2910, 54, "STOP:  back to the solve loop of Figure 1", pill=True)

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

# ===================================================== the learned surface ===
BX, BY, BW_, BH_ = PX, 2700, PW, 600
ax.add_patch(FancyBboxPatch((BX, BY), BW_, BH_,
                            boxstyle="round,pad=0,rounding_size=9",
                            fc=PBG, ec=PED, lw=1.2, zorder=1))
t(BX + 22, BY + 26, "what this network actually learned, drawn across the whole of its box",
  fs=13, c=GREY, w="bold")
ax.plot([n4["cx"] + n4["w"] / 2, BX], [n4["cy"], BY + 120], color="#C0C0CA",
        lw=1.2, ls=(0, (5, 4)), zorder=1)

GX, GY, GW, GH = BX + 96, BY + 62, 700, 268
ax.add_patch(Rectangle((GX, GY), GW, GH, fc="white", ec=EDGE, lw=1.1, zorder=3))
X0, X1, Y1 = 0.0, 10.0, 0.06


def gx(v):
    return GX + (v - X0) / (X1 - X0) * GW


def gy(v):
    return GY + GH - v / Y1 * GH


for v in (0.00, 0.02, 0.04, 0.06):
    ax.plot([GX, GX + GW], [gy(v), gy(v)], color="#E4E4EA", lw=1.0, zorder=3)
    t(GX - 10, gy(v), "%.2f" % v, fs=11.5, f=MONO, c=GREY, ha="right")
for v in (0, 2, 4, 6, 8, 10):
    ax.plot([gx(v), gx(v)], [GY, GY + GH], color="#E4E4EA", lw=1.0, zorder=3)
    t(gx(v), GY + GH + 14, "%g" % v, fs=11.5, f=MONO, c=GREY, ha="center")
N = 200
for nm, c, col in VOX:
    fe = bound(c)[1]
    xs, ys = [], []
    for i in range(N + 1):
        u = X0 + (X1 - X0) * i / float(N)
        b = max(LO[0], min(HI[0], u))
        raw = net([b, fe])
        ys.append(gy(0.0 if not (raw > 1e-8) else raw))
        xs.append(gx(u))
    for i in range(N):
        ax.plot([xs[i], xs[i + 1]], [ys[i], ys[i + 1]], color=col, lw=2.0, zorder=6)
for nm, c, col in VOX:
    b = bound(c); raw, mu = growth(c)
    ax.add_patch(FancyBboxPatch((gx(b[0]) - 7, gy(mu) - 7), 14, 14,
                                boxstyle="round,pad=0,rounding_size=7",
                                fc=col, ec="white", lw=1.4, zorder=8))
    t(gx(b[0]) + (56 if mu > 0 else -46), gy(mu) + (-44 if mu > 0 else -26),
      "%s voxel\n%7.4f -> %s" % (nm, b[0], ("%.5f" % mu) if mu > 0 else "0"),
      fs=11.5, c=col, w="bold", ha="center", va="center")
STEP = 6.8414
ax.plot([gx(STEP), gx(STEP)], [GY, GY + GH], color=GREY, lw=1.2,
        ls=(0, (5, 4)), zorder=5)
t(gx(STEP) - 10, GY + 26, "the step, at 6.84", fs=11.5, c=GREY, it=True, ha="right")
t(GX + GW / 2, GY + GH + 40,
  "the acetate uptake bound handed in, in mmol per gram dry weight per hour",
  fs=12.5, c=GREY, ha="center")
t(GX - 72, GY + GH / 2, "growth, per hour", fs=12.5, c=GREY, ha="center", rot=90)

t(BX + 22, BY + 396,
  "Both curves are the shipped weights, swept across the whole fitted range of the first input, with the second held at each\n"
  "voxel's own value. The two squares are the two voxels above. The shape is the point. Below an acetate bound of about 6.84\n"
  "the network returns the bottom of the range it was fitted into, which de-scales to zero; above it, growth rises almost in a\n"
  "straight line to the edge of the box. The change between the two happens inside five hundredths of a unit, so a voxel can\n"
  "cross it in a single reaction step, and a colony drawing its own acetate down will.\n\n"
  "Whatever the training rows held there, the fit has turned it into a step, and nothing in the file can tell a reader whether\n"
  "that step is the cell's maintenance demand or the edge of the sweep. That is the question to ask of any surrogate before\n"
  "trusting a run that sits near it, and it is why the offline sweep in the yellow band of Figure 3 is worth running.",
  fs=12.5, va="top")

fig.savefig("/home/claude/fig3/m2b.png", dpi=100, facecolor="white")
to_pptx(ax, W, H, "/home/claude/fig3/m2b.pptx", fig=fig)
print("wrote m2b.png and m2b.pptx", W, H)
