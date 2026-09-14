#!/usr/bin/env python3
"""Method 1, companion: one linear program, drawn, on two voxels of example 09.

Every number is computed here from the four constants the configuration file
carries, so any line can be checked with a calculator, and they agree with the
worked example in Section 8 of the guide.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Rectangle

from draw import Rec, to_pptx

YEL, BLU, GRN, RED = "#FFE699", "#BDD7EE", "#C6E0B4", "#F2B2B2"
PINK, PUR = "#F8A5A5", "#E4D3F5"
EDGE, REDEDGE = "#9AA5B1", "#B04A4A"
TXT, GREY, ARR = "#202020", "#7F7F7F", "#595959"
PBG, PED, ROWED = "#FBFBFD", "#C9C9D2", "#DCDCE3"
FED, DRY = "#2E6DA4", "#B04A4A"
SANS, MONO = "Carlito", "DejaVu Sans Mono"

# ------------------------------------------------------- the shipped case --
KS = 0.05          # <half_saturation_constants>, both substrates, mol/L
VS, VO = 10.0, 4.0  # <fba_maximum_uptake_flux>, mmol/gDW/h
M = 24.6           # <biomass_molar_mass>, gDW/mol
B = 1.0            # biomass in this voxel, mol/L
DT = 0.015         # the reaction step, seconds
VOX = [("fed", 1.00, 1.00, FED), ("drawn down", 0.02, 1.00, DRY)]


def bounds(cs, co):
    return -VS * cs / (KS + cs), -VO * co / (KS + co)


def solve(cs, co):
    ls, lo = bounds(cs, co)
    capS, capO = -ls, -2.0 * lo          # the two ceilings on the growth flux
    g = min(capS, capO)
    return ls, lo, capS, capO, g, ("the donor" if capS < capO else "the acceptor")


def dC(v):
    return v * (B * M) * DT / 3.6e6


W, H = 1900, 3140
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
t(40, 48, "Inside one linear program: the same four reactions, on two voxels of one run",
  fs=25, w="bold")
t(40, 88,
  "Every number below is computed from the four constants in examples/09_fba_glpk, so any line can be checked with a calculator. The model is small on\n"
  "purpose: three metabolites and four reactions, which is few enough that the program can be solved by hand and the solver's answer confirmed.",
  fs=15, c=GREY, it=True, va="top")
for lab, col, x in (("donor 1.00, acceptor 1.00 mol/L", FED, 40),
                    ("donor 0.02, acceptor 1.00 mol/L", DRY, 560)):
    ax.add_patch(FancyBboxPatch((x, 140), 26, 20,
                                boxstyle="round,pad=0,rounding_size=4",
                                fc=col, ec="none", zorder=5))
    t(x + 38, 150, lab, fs=13.5, w="bold")

# ============================================================= START ========
n0 = node(238, 54, "START:  one wet voxel that holds biomass", pill=True)
px, py = panel(186, 120, "the matrix, the capacities and the objective are already built; only two numbers are new", tie=n0)
t(px + 10, py + 4,
  "The program was assembled once at start-up and has been in memory ever since. What arrives here\n"
  "is a pair of concentrations, and what leaves is a growth rate and three exchange fluxes.",
  fs=13, va="top")

# ============================================================= 1 bounds =====
n1 = node(440, 72, "Turn the two concentrations into\ntwo uptake bounds", "in", num=1)
px, py = panel(352, 300,
               "Monod on each substrate, against that substrate's own capacity: equation (3), twice", tie=n1)
t(px + 10, py + 2, "substrate   Ks     capacity", fs=12, c=GREY, f=MONO)
t(px + 410, py + 2, "the fed voxel", fs=12.5, w="bold", ha="center", c=FED)
t(px + 770, py + 2, "the drawn-down voxel", fs=12.5, w="bold", ha="center", c=DRY)
for i, (nm, cap) in enumerate((("donor", VS), ("acceptor", VO))):
    yy = py + 34 + i * 34
    t(px + 10, yy, "%-9s  %-5.2f  %-4.1f" % (nm, KS, cap), fs=12, f=MONO)
    for k, (_, cs, co, col) in enumerate(VOX):
        c = cs if i == 0 else co
        lo = bounds(cs, co)[i]
        t(px + 410 + k * 360, yy, "%.2f -> %9.6f" % (c, lo), fs=12, f=MONO,
          ha="center", c=col)
t(px + 10, py + 116,
  "Negative, because every exchange reaction is written as the metabolite leaving the cell, so taking\n"
  "something up is a flux in the negative direction. The bound is then held above the model's own\n"
  "floor and above what the voxel could actually give up in one step, neither of which bites here.",
  fs=12.5, va="top")
t(px + 10, py + 190,
  "Only the acceptor's bound is unchanged between the two voxels. That single difference is the whole\n"
  "of what the transport solver tells the metabolic model about where it is.",
  fs=12.5, va="top", c=REDEDGE)

# ============================================================= 2 model ======
n2 = node(790, 72, "The program: four reactions and\nthree rows that must balance", "calc", num=2)
px, py = panel(690, 430,
               "the whole model, and then the three balances solved, which is what makes this one checkable by hand", tie=n2)
for i, ln in enumerate([
        "v0   EX_S    S_c ->                 negative: the donor comes IN",
        "v1   EX_O    O_c ->                 negative: the acceptor comes IN",
        "v2   EX_P    P_c ->                 positive: the product goes OUT",
        "v3   BIO     S_c + 0.5 O_c -> P_c   the objective: growth"]):
    t(px + 10, py + 4 + i * 28, ln, fs=12, f=MONO)
ax.plot([px + 10, px + 700], [py + 130, py + 130], color=EDGE, lw=1.1, zorder=4)
for i, ln in enumerate([
        "          v0    v1    v2    v3            row S_c:  -v0 - v3     = 0  ->  v0 = -v3",
        "S_c       -1     0     0    -1            row O_c:  -v1 - 0.5 v3 = 0  ->  v1 = -0.5 v3",
        "O_c        0    -1     0  -0.5            row P_c:  -v2 + v3     = 0  ->  v2 = +v3",
        "P_c        0     0    -1    +1"]):
    t(px + 10, py + 152 + i * 28, ln, fs=12, f=MONO)
t(px + 10, py + 282,
  "Every flux is now written in terms of the growth flux alone, so the program has one unknown left\n"
  "and two bounds on it. A genome-scale model has a few thousand columns and cannot be read this way,\n"
  "but nothing about the procedure changes: the simplex does what the three lines on the right do.",
  fs=12.5, va="top")

# ============================================================= 3 ceilings ===
n3 = node(1280, 72, "Each bound becomes a ceiling on\ngrowth. The lower one wins", "calc", num=3)
px, py = panel(1170, 495,
               "one donor and half an acceptor per turn of the growth reaction, so the acceptor's bound buys twice as much growth", tie=n3)

GX, GY, GW, GH = px + 210, py + 26, 620, 210
ax.add_patch(Rectangle((GX, GY), GW, GH, fc="white", ec=EDGE, lw=1.1, zorder=3))
XMAX = 10.5


def gx(v):
    return GX + v / XMAX * GW


for v in (0, 2, 4, 6, 8, 10):
    ax.plot([gx(v), gx(v)], [GY, GY + GH], color="#E4E4EA", lw=1.0, zorder=3)
    t(gx(v), GY + GH + 14, "%g" % v, fs=11.5, f=MONO, c=GREY, ha="center")
row = 0
for nm, cs, co, col in VOX:
    ls, lo, capS, capO, g, who = solve(cs, co)
    for lab, cap in (("from the donor", capS), ("from the acceptor", capO)):
        yy = GY + 22 + row * 46
        binds = abs(cap - g) < 1e-9
        ax.add_patch(FancyBboxPatch((GX + 2, yy - 12), max(4.0, gx(cap) - GX - 2), 24,
                                    boxstyle="round,pad=0,rounding_size=4",
                                    fc=col if binds else "white", ec=col,
                                    lw=1.4, alpha=1.0 if binds else 1.0, zorder=5))
        t(GX - 12, yy, lab, fs=11.5, c=col, ha="right")
        t(gx(cap) + 10, yy, "%.6f%s" % (cap, "   binds" if binds else ""),
          fs=11.5, f=MONO, c=col, w="bold" if binds else "normal")
        row += 1
    row += 0
t(GX + GW / 2, GY + GH + 40, "the largest growth flux each bound allows, in mmol per gram dry weight per hour",
  fs=12.5, c=GREY, ha="center")
for k, (nm, cs, co, col) in enumerate(VOX):
    t(px + 10, GY + 44 + k * 92, nm, fs=13, w="bold", c=col)
t(px + 10, py + 290,
  "The two ceilings are not the two bounds. A turn of the growth reaction consumes one donor and half\n"
  "an acceptor, so an acceptor bound of 3.809524 permits a growth flux of twice that. It is the\n"
  "stoichiometry, not the size of the two capacities, that decides which substrate limits: 4 is less\n"
  "than 10, yet in the second voxel the donor is the one that binds.\n\n"
  "The limit changes hands at a donor concentration of exactly 0.16 mol per litre, and the panel at\n"
  "the foot of this figure draws the whole of that.",
  fs=12.5, va="top")

# ============================================================= 4 fluxes =====
n4 = node(1730, 72, "The rest of the flux vector\nfollows from the three rows", "calc", num=4)
px, py = panel(1652, 250, "no search is needed once the growth flux is known: the other three are it, scaled", tie=n4)
for k, (nm, cs, co, col) in enumerate(VOX):
    ls, lo, capS, capO, g, who = solve(cs, co)
    yy = py + 16 + k * 84
    t(px + 10, yy, nm, fs=13, w="bold", c=col)
    t(px + 210, yy, "v3 = %+9.6f   growth" % g, fs=12, f=MONO, c=col)
    t(px + 210, yy + 24, "v0 = %+9.6f   donor in      v1 = %+9.6f   acceptor in"
      % (-g, -0.5 * g), fs=12, f=MONO, c=col)
    t(px + 210, yy + 48, "v2 = %+9.6f   product out"
      % g, fs=12, f=MONO, c=col)
t(px + 10, py + 176,
  "In the fed voxel the acceptor is taken up at exactly its bound and the donor is not. In the other\n"
  "the donor is at its bound and the acceptor is well under: a cell has no use for an acceptor it\n"
  "cannot pair with a donor, which is why the uptake falls rather than staying at capacity.",
  fs=12.5, va="top")

# ============================================================= 5 to lattice =
n5 = node(2030, 72, "Fluxes become changes in\nconcentration in this voxel", "out", num=5)
px, py = panel(1952, 290, "the one place the two unit systems meet: per gram of cells per hour, into moles per litre per step", tie=n5)
t(px + 10, py + 4, "dC  =  v  x  ( B x M )  x  dt  /  3.6e6",
  fs=13.5, f=MONO, w="bold")
t(px + 10, py + 32,
  "B the biomass this voxel holds, %.1f mol/L; M the biomass molar mass, %.1f gDW/mol; dt the step, %.3f s"
  % (B, M, DT), fs=12, c=GREY, it=True)
for k, (nm, cs, co, col) in enumerate(VOX):
    g = solve(cs, co)[4]
    yy = py + 74 + k * 72
    t(px + 10, yy, nm, fs=13, w="bold", c=col)
    t(px + 210, yy, "donor    %+.6e      acceptor %+.6e"
      % (dC(-g), dC(-0.5 * g)), fs=12, f=MONO, c=col)
    t(px + 210, yy + 24, "product  %+.6e      all in mol/L"
      % dC(g), fs=12, f=MONO, c=col)
t(px + 10, py + 224,
  "The two checks worth doing every time: the acceptor change is exactly half the donor change,\n"
  "because the reaction says 0.5, and the product change is exactly minus the donor change, because\n"
  "the reaction says one for one. Neither survives a mistake in the unit chain.",
  fs=12.5, va="top")

# ============================================================= 6 biomass ====
n6 = node(2350, 72, "Growth becomes biomass, or,\nif there is none, decay", "out", num=6)
px, py = panel(2272, 250, "the growth flux is a specific growth rate per hour, so this conversion is not the one above", tie=n6)
t(px + 10, py + 4, "dB  =  v3  x  B  x  dt  /  3600", fs=13.5, f=MONO, w="bold")
for k, (nm, cs, co, col) in enumerate(VOX):
    g = solve(cs, co)[4]
    t(px + 10, py + 44 + k * 30, nm, fs=13, w="bold", c=col)
    t(px + 210, py + 44 + k * 30, "%.6f x %.1f x %.3f / 3600  =  %+.6e  mol/L"
      % (g, B, DT, g * B * DT / 3600.0), fs=12, f=MONO, c=col)
t(px + 10, py + 120,
  "A program that came back infeasible, or that returned no growth at all, takes the other branch:\n"
  "the organism loses biomass at its own decay coefficient instead. The two conversions differ by a\n"
  "factor of a thousand and are easy to confuse, which is why growth is per hour and decay per second\n"
  "in the code and each has its own converter.",
  fs=12.5, va="top")

d7 = dia(2610, 100, "more organisms here,\nmore voxels?", w=450, fs=13.5)
stop = node(2740, 54, "STOP:  back to the solve loop of Figure 3", pill=True)

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

# ================================================= where the limit changes ===
BX, BY, BW_, BH_ = PX, 2540, PW, 560
ax.add_patch(FancyBboxPatch((BX, BY), BW_, BH_,
                            boxstyle="round,pad=0,rounding_size=9",
                            fc=PBG, ec=PED, lw=1.2, zorder=1))
t(BX + 22, BY + 26, "where the limit changes hands, drawn across the whole run",
  fs=13, c=GREY, w="bold")
ax.plot([n3["cx"] + n3["w"] / 2, BX], [n3["cy"], BY + 140], color="#C0C0CA",
        lw=1.2, ls=(0, (5, 4)), zorder=1)

CX2, CY2, CW, CH = BX + 96, BY + 62, 700, 250
ax.add_patch(Rectangle((CX2, CY2), CW, CH, fc="white", ec=EDGE, lw=1.1, zorder=3))
X1, Y1 = 1.0, 9.0


def cx_(v):
    return CX2 + v / X1 * CW


def cy_(v):
    return CY2 + CH - v / Y1 * CH


for v in (0, 2, 4, 6, 8):
    ax.plot([CX2, CX2 + CW], [cy_(v), cy_(v)], color="#E4E4EA", lw=1.0, zorder=3)
    t(CX2 - 10, cy_(v), "%g" % v, fs=11.5, f=MONO, c=GREY, ha="right")
for v in (0.0, 0.16, 0.4, 0.6, 0.8, 1.0):
    ax.plot([cx_(v), cx_(v)], [CY2, CY2 + CH], color="#E4E4EA", lw=1.0, zorder=3)
    t(cx_(v), CY2 + CH + 14, "%.2f" % v, fs=11.5, f=MONO, c=GREY, ha="center")
CAPO = 2.0 * VO * 1.0 / (KS + 1.0)
N = 240
xs, ys = [], []
for i in range(N + 1):
    c = X1 * i / float(N)
    xs.append(cx_(c)); ys.append(cy_(min(VS * c / (KS + c), CAPO)))
for i in range(N):
    ax.plot([xs[i], xs[i + 1]], [ys[i], ys[i + 1]], color="#4A4A55", lw=2.2, zorder=6)
ax.plot([cx_(0.16), cx_(0.16)], [CY2, CY2 + CH], color=GREY, lw=1.3,
        ls=(0, (5, 4)), zorder=5)
t(cx_(0.16) - 10, CY2 + 24, "0.16 mol/L", fs=11.5, c=GREY, it=True, ha="right")
t(cx_(0.16) + 12, CY2 + 62, "the donor limits\nleft of here", fs=11.5, c=GREY,
  it=True, ha="left", va="top")
t(cx_(0.16) + 12, CY2 + 24, "the acceptor limits right of here", fs=11.5, c=GREY,
  it=True, ha="left")
for nm, cs, co, col in VOX:
    g = solve(cs, co)[4]
    ax.add_patch(FancyBboxPatch((cx_(cs) - 7, cy_(g) - 7), 14, 14,
                                boxstyle="round,pad=0,rounding_size=7",
                                fc=col, ec="white", lw=1.4, zorder=8))
    t(cx_(cs) + (0 if cs > 0.5 else 18), cy_(g) + (30 if cs > 0.5 else -2),
      "%s\n%.2f -> %.6f" % (nm, cs, g), fs=11.5, c=col, w="bold",
      ha="center" if cs > 0.5 else "left", va="center")
t(CX2 + CW / 2, CY2 + CH + 40, "the donor concentration in this voxel, in moles per litre",
  fs=12.5, c=GREY, ha="center")
t(CX2 - 54, CY2 + CH / 2, "growth flux", fs=12.5, c=GREY, ha="center", rot=90)

t(BX + 22, BY + 380,
  "The acceptor is held at 1.00 mol per litre and the donor swept across everything a run visits. To the right of 0.16 the acceptor\n"
  "limits and growth is flat: adding donor to a voxel that already has enough buys nothing, which is why a concentration field alone\n"
  "cannot tell a modeller whether the medium is rich or merely adequate. To the left the donor limits and growth falls with it,\n"
  "steeply, because the Monod term is still climbing there. The corner is not a numerical artefact: it is the exact concentration at\n"
  "which the two ceilings are equal, and either side of it a different constraint is doing all the work. A run where the corner sits\n"
  "inside the domain is the interesting kind, and nothing in the output says so unless the limiting substrate is printed.",
  fs=12.5, va="top")

fig.savefig("/home/claude/fig3/m1b.png", dpi=100, facecolor="white")
to_pptx(ax, W, H, "/home/claude/fig3/m1b.pptx", fig=fig)
print("wrote m1b.png and m1b.pptx", W, H)
