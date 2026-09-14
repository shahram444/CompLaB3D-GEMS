#!/usr/bin/env python3
"""
Method 3, companion figure: the search as a flowchart with the formulas drawn.

The flow boxes on the left are deliberately colourless, with only a thin
coloured stripe for the phase. All the colour in the figure belongs to the
formulas on the right, because those are the things being changed. Every step
that does something to a formula shows it.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from draw import Rec, to_pptx
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Rectangle

YEL, BLU, GRN = "#FFE699", "#BDD7EE", "#C6E0B4"      # operator / name / number
PINK = "#F8A5A5"
EDGE, REDEDGE = "#9AA5B1", "#B04A4A"
TXT, GREY, ARR = "#202020", "#7F7F7F", "#595959"
PBG, PED, ROWED = "#FBFBFD", "#C9C9D2", "#DCDCE3"
RING = "#B04A4A"
STRIPE = {"setup": "#8FAADC", "each": "#9CC183", "shelf": "#E8C56A"}
SANS, MONO = "Carlito", "DejaVu Sans Mono"
FILL = {"op": YEL, "var": BLU, "num": GRN}

W, H = 1900, 2880
fig = plt.figure(figsize=(W / 100.0, H / 100.0), dpi=100)
ax = Rec(fig.add_axes([0, 0, 1, 1]))
ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
ax.add_patch(Rectangle((0, 0), W, H, fc="white", ec="none", zorder=0))


def arrow(x0, y0, x1, y1, color=ARR, lw=1.7, z=4, ms=12):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=ms, color=color, lw=lw,
                                 shrinkA=0, shrinkB=0, zorder=z))


def lane(pts, color=ARR, lw=1.7, z=2, ms=12):
    for i in range(len(pts) - 2):
        ax.plot([pts[i][0], pts[i + 1][0]], [pts[i][1], pts[i + 1][1]],
                color=color, lw=lw, solid_capstyle="round", zorder=z)
    arrow(pts[-2][0], pts[-2][1], pts[-1][0], pts[-1][1], color, lw, z, ms)


def t(x, y, s, fs=13, c=TXT, f=SANS, w="normal", ha="left", va="center",
      it=False, z=7, lsp=1.4):
    ax.text(x, y, s, fontsize=fs, color=c, family=f, fontweight=w, ha=ha,
            va=va, style="italic" if it else "normal", zorder=z, linespacing=lsp)


# ------------------------------------------------------------------- trees ---
def N(lab, kind, *kids):
    return (lab, kind, list(kids))


def _place(n, x, y, dx, dy, out):
    if not n[2]:
        out.append((n, x + dx / 2, y)); return x + dx, x + dx / 2
    cx, mids = x, []
    for k in n[2]:
        cx, m = _place(k, cx, y + dy, dx, dy, out); mids.append(m)
    mid = sum(mids) / len(mids); out.append((n, mid, y)); return cx, mid


def _under(sub, n):
    return n is sub or any(_under(k, n) for k in sub[2])


def tree(root, cx, top, bw=58, bh=24, dx=68, dy=38, fs=11, ring=None,
         ringone=None, ringlab=None, z=6):
    out = []
    _place(root, 0, 0, dx, dy, out)
    xs = [p[1] for p in out]
    sh = cx - (min(xs) + max(xs)) / 2
    pos = {id(n): (x + sh, top + y) for n, x, y in out}

    def edges(n):
        for k in n[2]:
            x0, y0 = pos[id(n)]; x1, y1 = pos[id(k)]
            ax.plot([x0, x1], [y0 + bh / 2, y1 - bh / 2], color=EDGE, lw=1.0,
                    zorder=z - 1)
            edges(k)
    edges(root)
    for n, _, _ in out:
        x, y = pos[id(n)]
        ax.add_patch(FancyBboxPatch((x - bw / 2, y - bh / 2), bw, bh,
                                    boxstyle="round,pad=0,rounding_size=4",
                                    fc=FILL[n[1]], ec=EDGE, lw=1.0, zorder=z))
        ax.text(x, y, n[0], ha="center", va="center", fontsize=fs, color=TXT,
                family=MONO, zorder=z + 1)
    tgt = ringone or ring
    if tgt is not None:
        sub = [tgt] if ringone else [n for n, _, _ in out if _under(tgt, n)]
        px = [pos[id(n)][0] for n in sub]; py = [pos[id(n)][1] for n in sub]
        x0, x1 = min(px) - bw / 2 - 7, max(px) + bw / 2 + 7
        y0, y1 = min(py) - bh / 2 - 7, max(py) + bh / 2 + 7
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                    boxstyle="round,pad=0,rounding_size=8",
                                    fc="none", ec=RING, lw=1.6, ls=(0, (4, 3)),
                                    zorder=z + 2))
        if ringlab:
            t((x0 + x1) / 2, y1 + 14, ringlab, fs=11.5, c=RING, it=True,
              ha="center", va="top", z=z + 2)
    return pos


# ------------------------------------------------------------------- nodes ---
FX, FW = 470, 580
PX, PW = 810, 1050
L1, L2, L3 = 140, 92, 44


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
    ax.text(FX + 28 if num is not None else FX + 6, cy, text, ha="center",
            va="center", fontsize=fs, color=TXT, family=SANS,
            fontweight="bold", linespacing=1.3, zorder=7)
    if sub:
        t(FX, cy + h / 2 + 12, sub, fs=12, c=GREY, it=True, ha="center",
          va="top")
    return dict(cx=FX, cy=cy, w=FW, h=h)


def dia(cy, h, text, w=520, fs=14):
    ax.add_patch(Polygon([(FX, cy - h / 2), (FX + w / 2, cy), (FX, cy + h / 2),
                          (FX - w / 2, cy)], closed=True, fc="white", ec=EDGE,
                         lw=1.3, zorder=5))
    ax.text(FX, cy, text, ha="center", va="center", fontsize=fs, color=TXT,
            family=SANS, linespacing=1.3, zorder=7)
    return dict(cx=FX, cy=cy, w=w, h=h)


def panel(y0, h, title, x=PX, w=PW, tie=None):
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


# ------------------------------------------------------------------ header ---
t(40, 48, "Inside the search, drawn: what each step does to a formula", fs=25,
  w="bold")
t(40, 88,
  "The flow is on the left and has no colour on purpose. All the colour belongs to the formulas on the right, because those are what the search changes.",
  fs=15, c=GREY, it=True, va="top")
t(40, 138, "every formula is drawn as boxes:", fs=13.5, w="bold")
kx = 330
for fl, lab in ((YEL, "does one step:  +  -  *  /"),
                (BLU, "a name, read from the voxel"),
                (GRN, "a number, set by fitting")):
    ax.add_patch(FancyBboxPatch((kx, 126), 32, 24,
                                boxstyle="round,pad=0,rounding_size=4",
                                fc=fl, ec=EDGE, lw=1.0, zorder=5))
    t(kx + 42, 138, lab, fs=13)
    kx += 320

# ------------------------------------------------------------------ shapes ---
MONOD = N("/", "op", N("*", "op", N("0.35", "num"), N("a", "var")),
          N("+", "op", N("0.02", "num"), N("a", "var")))
QMARK = N("/", "op", N("*", "op", N("?", "num"), N("a", "var")),
          N("+", "op", N("?", "num"), N("a", "var")))
R1 = N("*", "op", N("a", "var"), N("b", "var"))
R2 = N("+", "op", N("0.7", "num"), N("/", "op", N("b", "var"), N("1.4", "num")))
R3 = N("/", "op", N("a", "var"), N("+", "op", N("a", "var"), N("b", "var")))
R4 = N("-", "op", N("b", "var"), N("0.9", "num"))
R5 = N("*", "op", N("2.1", "num"), N("a", "var"))
PA = N("/", "op", N("b", "var"), N("+", "op", N("0.02", "num"), N("a", "var")))
PB = N("/", "op", N("*", "op", N("0.35", "num"), N("a", "var")), N("b", "var"))
CH = N("/", "op", N("*", "op", N("0.35", "num"), N("a", "var")),
       N("+", "op", N("0.02", "num"), N("a", "var")))
MU = N("*", "op", N("b", "var"), N("+", "op", N("0.02", "num"), N("a", "var")))


# ============================================================ START ==========
n0 = node(230, 54, "START:  the data you bring", pill=True)
px, py = panel(178, 208, tie=n0, title="one row per state of your system: the columns you measured, and the one column you want a formula for")
COLW = 192
for j, (letter, what) in enumerate((
        ("a", "an input you measured\nhere, a concentration"),
        ("b", "a second input\nhere, another concentration"),
        ("y", "the one you want a formula for\nhere, a growth rate"))):
    cx = px + 120 + j * COLW
    t(cx, py + 2, letter, fs=13.5, f=MONO, w="bold", ha="center")
    t(cx, py + 18, what, fs=11.5, c=GREY, it=True, ha="center", va="top")
ax.plot([px + 16, px + 120 + 2 * COLW + 100], [py + 62, py + 62], color=EDGE,
        lw=1.1, zorder=3)
for i, row in enumerate((["0.005", "0.30", "0.0700"], ["0.010", "0.55", "0.1167"],
                         ["0.020", "0.20", "0.1750"], ["...", "...", "..."])):
    for j, c in enumerate(row):
        t(px + 120 + j * COLW, py + 78 + i * 18, c, fs=12, f=MONO, ha="center")
t(px + 690, py + 18,
  "The names are yours. The search only ever\n"
  "sees the numbers, so it cannot be led by\n"
  "what you expect the answer to be.\n\n"
  "Any number of input columns is allowed;\n"
  "two are shown here to keep the pictures small.",
  fs=12.5, va="top")

# ============================================================ 1 ==============
n1 = node(452, 66, "Make a few hundred formulas\nat random", "setup", num=1)
px, py = panel(400, 178, tie=n1, title="a box is drawn from  +  -  *  /  or a name or a number, and whatever hangs under it is drawn the same way")
for k, tr in enumerate((R1, R5, R3, R2, R4)):
    tree(tr, px + 95 + k * 200, py + 18)
t(px, py + 116, "none of them means anything yet", fs=12.5, c=GREY, it=True,
  va="top")

# ============================================================ 2 ==============
n2 = node(672, 66, "Score every one of them", "each", num=2)
px, py = panel(608, 212, tie=n2, title="least squares puts the best possible numbers in, the gap to your data is measured, then a penalty is added per box")
tree(QMARK, px + 140, py + 16)
t(px + 140, py + 118, "the shape, its numbers not chosen yet", fs=11.5, c=GREY,
  it=True, ha="center")
arrow(px + 286, py + 54, px + 368, py + 54, ms=11)
t(px + 327, py + 42, "least squares", fs=11.5, c=GREY, it=True, ha="center",
  va="bottom")
tree(MONOD, px + 510, py + 16)
t(px + 510, py + 118, "the same shape, best numbers in", fs=11.5, c=GREY,
  it=True, ha="center")
arrow(px + 656, py + 54, px + 716, py + 54, ms=11)
t(px + 732, py + 18, "its score", fs=13, w="bold")
ax.add_patch(Rectangle((px + 732, py + 34), 190, 30, fc="#E3A7A7", ec=EDGE,
                       lw=1.0, zorder=5))
ax.add_patch(Rectangle((px + 922, py + 34), 56, 30, fc="#D9D9E0", ec=EDGE,
                       lw=1.0, zorder=5))
t(px + 732, py + 76, "how far it lands from your data", fs=11.5, c=GREY,
  va="top")
t(px + 732, py + 96, "plus a penalty, one per box", fs=11.5, c=GREY, va="top")

# ============================================================ 3 ==============
n3 = node(936, 66, "Put the best of each size\non the shelf", "shelf", num=3)
px, py = panel(858, 256, tie=n3, title="one slot per number of boxes, so a short formula is only ever judged against another short one")
for k, (lab, tr, miss) in enumerate((("3 boxes", R5, "0.075"),
                                     ("5 boxes", R3, "0.056"),
                                     ("7 boxes", MONOD, "0.000"),
                                     ("9 boxes", None, ""))):
    sx = px + 8 + k * 250
    ax.add_patch(FancyBboxPatch((sx, py), 214, 186,
                                boxstyle="round,pad=0,rounding_size=6",
                                fc="white", ec=ROWED, lw=1.2, zorder=3))
    t(sx + 107, py + 20, lab, fs=12, c=GREY, ha="center", z=5)
    if tr is not None:
        tree(tr, sx + 107, py + 50)
        t(sx + 107, py + 164, miss, fs=13, f=MONO, w="bold", ha="center", z=5)
    else:
        t(sx + 107, py + 92, "still empty", fs=12.5, c=GREY, it=True,
          ha="center", z=5)

# ============================================================ round =========
n4 = node(1180, 58, "Start a round with an empty new population", "setup", num=4)
d5 = dia(1292, 100, "is the new population full?", w=500)

# ============================================================ 4 =============
n5 = node(1688, 66, "Make one new formula\nfrom the population", "each", num=5)
px, py = panel(1378, 620, tie=n5, title="one of three moves, chosen at random. The odds are what keep the population both improving and varied")
IW = PW - 44


def chip(x, y, label, body, w=60):
    ax.add_patch(FancyBboxPatch((x, y - 13), w, 26,
                                boxstyle="round,pad=0,rounding_size=13",
                                fc="#4A4A55", ec="none", zorder=6))
    t(x + w / 2, y, label, fs=11.5, c="white", w="bold", ha="center", z=7)
    t(x + w + 14, y, body, fs=12.5, z=7)


def newtag(cx, y):
    ax.add_patch(FancyBboxPatch((cx - 86, y - 13), 172, 26,
                                boxstyle="round,pad=0,rounding_size=6",
                                fc="#E8F0DE", ec="#6E9150", lw=1.2, zorder=6))
    t(cx, y, "the new formula", fs=11.5, c="#3E5A28", w="bold", ha="center", z=7)


# --- what comes in
chip(px, py + 12, "IN", "the population as step 2 left it: every formula in it already has its numbers fitted and its score")

# --- SWAP
S0 = py + 44
ax.add_patch(FancyBboxPatch((px - 10, S0), IW + 20, 236,
                            boxstyle="round,pad=0,rounding_size=7",
                            fc="white", ec=ROWED, lw=1.2, zorder=2))
t(px + 10, S0 + 22, "SWAP", fs=13.5, w="bold", z=4)
t(px + 96, S0 + 22, "6 of every 10", fs=12.5, c=GREY, z=4)
t(px + 250, S0 + 12, "two parents, each the best of three taken at random. Ring one box in\n"
  "each, everything under it included, and let the rings change places",
  fs=12, c=GREY, z=4, va="top")
t(px + 170, S0 + 60, "parent A,  from the population", fs=11.5, w="bold",
  ha="center", z=4)
tree(PA, px + 170, S0 + 84, ring=PA[2][0])
t(px + 170, S0 + 186, "this goes out", fs=11.5, c=RING, it=True, ha="center", z=4)
t(px + 500, S0 + 60, "parent B,  from the population", fs=11.5, w="bold",
  ha="center", z=4)
tree(PB, px + 500, S0 + 84, ring=PB[2][0])
t(px + 500, S0 + 186, "this comes in", fs=11.5, c=RING, it=True, ha="center", z=4)
arrow(px + 630, S0 + 122, px + 700, S0 + 122, ms=11)
t(px + 860, S0 + 60, "the child", fs=11.5, w="bold", ha="center", z=4)
tree(CH, px + 860, S0 + 84, ring=CH[2][0])
newtag(px + 860, S0 + 196)

# --- CHANGE and FRESH
C0 = py + 296
CW_ = 660
ax.add_patch(FancyBboxPatch((px - 10, C0), CW_, 214,
                            boxstyle="round,pad=0,rounding_size=7",
                            fc="white", ec=ROWED, lw=1.2, zorder=2))
t(px + 10, C0 + 22, "CHANGE", fs=13.5, w="bold", z=4)
t(px + 116, C0 + 22, "3 of every 10", fs=12.5, c=GREY, z=4)
t(px + 230, C0 + 22, "one parent from the population, and exactly one box changes",
  fs=12, c=GREY, z=4)
t(px + 140, C0 + 50, "before", fs=11.5, w="bold", ha="center", z=4)
tree(PA, px + 140, C0 + 82, bw=50, bh=21, dx=59, dy=33, fs=10.5, ringone=PA)
arrow(px + 300, C0 + 104, px + 372, C0 + 104, ms=11)
t(px + 500, C0 + 50, "after", fs=11.5, w="bold", ha="center", z=4)
tree(MU, px + 500, C0 + 82, bw=50, bh=21, dx=59, dy=33, fs=10.5, ringone=MU)
t(px + 336, C0 + 124, "only the ringed box\nchanged:  /  became  *", fs=11,
  c=RING, it=True, ha="center", va="top", z=4)
newtag(px + 500, C0 + 190)

ax.add_patch(FancyBboxPatch((px + CW_ - 4, C0), IW - CW_ + 14, 214,
                            boxstyle="round,pad=0,rounding_size=7",
                            fc="white", ec=ROWED, lw=1.2, zorder=2))
t(px + CW_ + 16, C0 + 22, "FRESH", fs=13.5, w="bold", z=4)
t(px + CW_ + 102, C0 + 22, "1 of every 10", fs=12.5, c=GREY, z=4)
t(px + CW_ + 16, C0 + 44, "nothing from the population at all:\na brand new random formula, as in step 1",
  fs=12, c=GREY, va="top", z=4)
tree(R2, px + CW_ + 180, C0 + 96, bw=50, bh=21, dx=59, dy=33, fs=10.5)
newtag(px + CW_ + 180, C0 + 190)

# --- what goes out
chip(px, py + 542, "OUT",
     "whichever of the three ran, exactly one new formula leaves this step: the one tagged above.\n"
     "Step 6 scores it, and the question under step 6 decides whether it earns a slot on the shelf.", w=74)

# ============================================================ 5 =============
n6 = node(2056, 62, "Score it, exactly as in step 2", "each", num=6,
          sub="its shape is new, so the numbers it inherited are refitted before anything is measured")
d7 = dia(2188, 110, "does it beat the shelf's formula\nof the same size?", w=620,
         fs=13.5)

ax.add_patch(FancyBboxPatch((PX + 60, 2158), 520, 62,
                            boxstyle="round,pad=0,rounding_size=8",
                            fc="white", ec=EDGE, lw=1.3, zorder=5))
ax.add_patch(Rectangle((PX + 62, 2160), 8, 58, fc=STRIPE["shelf"], ec="none",
                       zorder=6))
ax.add_patch(FancyBboxPatch((PX + 80, 2174), 30, 30,
                            boxstyle="round,pad=0,rounding_size=15",
                            fc="#4A4A55", ec="none", zorder=7))
t(PX + 95, 2189, "7", fs=13, c="white", w="bold", ha="center", z=8)
t(PX + 345, 2189, "It takes that slot; the old one is gone", fs=15, w="bold",
  ha="center")
t(PX + 325, 2228, "the only place the shelf ever changes", fs=12, c=GREY,
  it=True, ha="center", va="top")
arrow(FX + d7["w"] / 2, 2188, PX + 60, 2188)
t(FX + d7["w"] / 2 + 14, 2170, "yes", fs=12.5, c=GREY)

n8 = node(2328, 58, "Add it to the new population", "each", num=8)
n9 = node(2440, 58, "The new population replaces the old one", "setup", num=9)
d10 = dia(2552, 96, "rounds left?", w=380)
stop = node(2682, 54, "STOP:  write the shelf out, shortest first", pill=True)

px, py = panel(2616, 196, tie=stop, title="the list: one formula per size. Choosing a row off it is your decision, not the tool's")
for k, (lab, tr, miss) in enumerate((("3 boxes", R5, "0.075"),
                                     ("5 boxes", R3, "0.056"),
                                     ("7 boxes", MONOD, "0.000"))):
    sx = px + 30 + k * 330
    t(sx + 100, py + 4, lab, fs=12, c=GREY, ha="center")
    tree(tr, sx + 100, py + 30, bw=52, bh=22, dx=61, dy=33, fs=10.5)
    t(sx + 100, py + 122, miss, fs=12.5, f=MONO, w="bold", ha="center")

# ------------------------------------------------------------------ arrows --
down(n0, n1)
down(n1, n2)
down(n2, n3)
down(n3, n4)
down(n4, d5)
down(d5, n5, "no")
down(n5, n6)
down(n6, d7)
down(d7, n8, "no")
down(n9, d10)
down(d10, stop, "no")

lane([(PX + 320, 2220), (PX + 320, n8["cy"]), (FX + FW / 2, n8["cy"])])

lane([(FX - FW / 2, n8["cy"]), (L1, n8["cy"]), (L1, 1238), (FX - 6, 1238)])
ax.text(L1 + 12, (n8["cy"] + 1238) / 2, "the next formula", rotation=90,
        ha="center", va="center", fontsize=12.5, color=GREY, family=SANS, zorder=6)

lane([(FX - d5["w"] / 2, d5["cy"]), (L2, d5["cy"]), (L2, n9["cy"]),
      (FX - FW / 2, n9["cy"])])
t(FX - d5["w"] / 2 - 12, d5["cy"] - 17, "yes", fs=12.5, c=GREY, ha="right")
ax.text(L2 - 12, (d5["cy"] + n9["cy"]) / 2, "the round ends", rotation=90,
        ha="center", va="center", fontsize=12.5, color=GREY, family=SANS, zorder=6)

lane([(FX - d10["w"] / 2, d10["cy"]), (L3, d10["cy"]), (L3, 1050),
      (FX - 6, 1050)])
t(FX - d10["w"] / 2 - 12, d10["cy"] - 17, "yes", fs=12.5, c=GREY, ha="right")
ax.text(L3 - 12, (d10["cy"] + 1128) / 2, "the next round", rotation=90,
        ha="center", va="center", fontsize=12.5, color=GREY, family=SANS, zorder=6)

fig.savefig("/home/claude/fig3/figB.png", dpi=100, facecolor="white")
to_pptx(ax, W, H, fig=fig, path= "/home/claude/fig3/figB.pptx")
print("wrote fig3e.png", W, H)
