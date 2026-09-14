#!/usr/bin/env python3
"""
The shared drawing kit for the "path as code" flowcharts.

One of these exists for every rate path in docs/methods. They have to look like
each other or the reader has to learn the notation five times, so everything
that decides how they look lives here and each figure file holds only content.

The palette is the one Figure 1 of every guide already uses.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Rectangle

from draw import Rec, to_pptx

YEL, BLU, GRN, RED = "#FFE699", "#BDD7EE", "#C6E0B4", "#F2B2B2"
PINK, PUR, ORG = "#F8A5A5", "#E4D3F5", "#F8CBAD"
EDGE, REDEDGE = "#9AA5B1", "#B04A4A"
TXT, GREY, ARR = "#202020", "#7F7F7F", "#595959"
PBG, PED, ROWED = "#FBFBFD", "#C9C9D2", "#DCDCE3"
RING = "#B04A4A"
SANS, MONO = "Carlito", "DejaVu Sans Mono"

CX, BW = 560, 700           # the main column
L_LOOP = CX - BW / 2 - 90   # the lane the loop-back arrows run down
PANEL_X, PANEL_W = 1010, 840


class Fig(object):
    """A figure in progress. Everything draws through this."""

    def __init__(self, w, h, title, subtitle):
        self.W, self.H = w, h
        self.fig = plt.figure(figsize=(w / 100.0, h / 100.0), dpi=100)
        self.ax = Rec(self.fig.add_axes([0, 0, 1, 1]))
        self.ax.set_xlim(0, w)
        self.ax.set_ylim(h, 0)
        self.ax.axis("off")
        self.ax.add_patch(Rectangle((0, 0), w, h, fc="white", ec="none", zorder=0))
        self.ax.text(40, 52, title, fontsize=25, fontweight="bold", color=TXT,
                     family=SANS, va="center", zorder=7)
        self.ax.text(40, 100, subtitle, fontsize=15, color=GREY, family=SANS,
                     style="italic", va="top", linespacing=1.45, zorder=7)

    # ---------------------------------------------------------------- atoms --
    def t(self, x, y, s, fs=13, c=TXT, f=SANS, w="normal", ha="left",
          va="center", it=False, z=7, lsp=1.4, rot=0):
        self.ax.text(x, y, s, fontsize=fs, color=c, family=f, fontweight=w,
                     ha=ha, va=va, style="italic" if it else "normal",
                     zorder=z, linespacing=lsp, rotation=rot)

    def note(self, cx, y, s, color=GREY, fs=12.5):
        self.ax.text(cx, y, s, ha="center", va="top", fontsize=fs, color=color,
                     family=SANS, style="italic", linespacing=1.35, zorder=9,
                     bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="none"))

    def arrow(self, x0, y0, x1, y1, color=ARR, lw=1.6, z=2, ms=11):
        self.ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                          mutation_scale=ms, color=color, lw=lw,
                                          shrinkA=0, shrinkB=0, zorder=z))

    def lane(self, pts, color=ARR, lw=1.6, z=2, ms=11):
        for i in range(len(pts) - 2):
            self.ax.plot([pts[i][0], pts[i + 1][0]], [pts[i][1], pts[i + 1][1]],
                         color=color, lw=lw, solid_capstyle="round", zorder=z)
        self.arrow(pts[-2][0], pts[-2][1], pts[-1][0], pts[-1][1], color, lw, z, ms)

    # ---------------------------------------------------------------- nodes --
    def box(self, cy, h, fill, text, sub=None, ec=EDGE, fs=15, w=None, cx=None,
            r=8):
        w = w or BW
        cx = CX if cx is None else cx
        self.ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                         boxstyle="round,pad=0,rounding_size=%g" % r,
                                         fc=fill, ec=ec, lw=1.3, zorder=5))
        self.ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
                     color=TXT, family=SANS, fontweight="bold", linespacing=1.35,
                     zorder=6)
        if sub:
            self.note(cx, cy + h / 2 + 14, sub)
        return dict(cx=cx, cy=cy, w=w, h=h)

    def dia(self, cy, h, text, w=460, fs=14, cx=None):
        cx = CX if cx is None else cx
        self.ax.add_patch(Polygon([(cx, cy - h / 2), (cx + w / 2, cy),
                                   (cx, cy + h / 2), (cx - w / 2, cy)],
                                  closed=True, fc="white", ec=EDGE, lw=1.3,
                                  zorder=5))
        self.ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
                     color=TXT, family=SANS, linespacing=1.3, zorder=6)
        return dict(cx=cx, cy=cy, w=w, h=h)

    def down(self, a, b, label=None, side=1):
        y0, y1 = a["cy"] + a["h"] / 2, b["cy"] - b["h"] / 2
        self.arrow(a["cx"], y0, b["cx"], y1)
        if label:
            self.t(a["cx"] + 13 * side, (y0 + y1) / 2, label, fs=12.5, c=GREY,
                   ha="left" if side > 0 else "right")

    # --------------------------------------------------------------- frames --
    def band(self, y0, y1, colour, name, tag=None, first=None):
        self.ax.add_patch(Rectangle((44, y0), 9, y1 - y0, fc=colour, ec="none",
                                    zorder=2))
        self.t(30, (y0 + y1) / 2, name, fs=15, w="bold", ha="center", rot=90)
        if tag is not None and first is not None:
            self.t(CX - BW / 2, first - 14, tag, fs=12.5, c=GREY, f=MONO,
                   va="bottom")

    def legend(self, rows, y=200, x=PANEL_X, w=PANEL_W, title="what the colours mean"):
        self.t(x, y - 26, title, fs=15, w="bold")
        for fill, lab, ec in rows:
            self.ax.add_patch(FancyBboxPatch((x, y), w, 46,
                                             boxstyle="round,pad=0,rounding_size=6",
                                             fc=fill, ec=ec, lw=1.2, zorder=3))
            self.ax.text(x + w / 2, y + 23, lab, ha="center", va="center",
                         fontsize=14.5, color=TXT, family=SANS, style="italic",
                         zorder=4)
            y += 58
        return y

    def panel(self, x, y0, w, title, italic, rows, anchor=None):
        lh = 46
        h = 64 + (26 if italic else 0) + len(rows) * (lh + 10) + 26
        self.ax.add_patch(FancyBboxPatch((x, y0), w, h,
                                         boxstyle="round,pad=0,rounding_size=8",
                                         fc=PBG, ec=PED, lw=1.3, ls=(0, (6, 5)),
                                         zorder=2))
        yy = y0 + 34
        self.t(x + 26, yy, title, fs=15.5, w="bold")
        yy += 30
        if italic:
            self.t(x + 26, yy, italic, fs=13, c=GREY, it=True, va="top")
            yy += 42
        for r in rows:
            self.ax.add_patch(FancyBboxPatch((x + 26, yy), w - 52, lh,
                                             boxstyle="round,pad=0,rounding_size=5",
                                             fc="white", ec=ROWED, lw=1.0, zorder=3))
            self.ax.text(x + w / 2, yy + lh / 2, r, ha="center", va="center",
                         fontsize=13.5, color=TXT, family=SANS, zorder=4)
            yy += lh + 10
        if anchor:
            self.ax.plot([anchor[0], x], [anchor[1], y0 + h / 2], color=PINK,
                         lw=1.3, ls=(0, (6, 5)), zorder=1)
        return y0 + h

    def handover(self, after, label, note, gap=128):
        y0 = after["cy"] + after["h"] / 2
        y1 = y0 + gap
        self.arrow(CX, y0, CX, y1, lw=2.4, ms=15)
        self.t(CX + 20, y0 + gap * 0.35, label, fs=15.5, w="bold")
        self.t(CX + 20, y0 + gap * 0.6, note, fs=13.5, c=GREY, it=True)
        self.ax.plot([100, self.W - 40], [y0 + gap * 0.47, y0 + gap * 0.47],
                     color="#C9C9D2", lw=1.2, ls=(0, (6, 5)), zorder=1)
        return y1

    def loop(self, frm, to_y, label):
        self.lane([(CX - frm["w"] / 2, frm["cy"]), (L_LOOP, frm["cy"]),
                   (L_LOOP, to_y), (CX - BW / 2, to_y)])
        self.t(CX - frm["w"] / 2 - 14, frm["cy"] - 16, "yes", fs=12.5, c=GREY,
               ha="right")
        self.t(L_LOOP + 11, (frm["cy"] + to_y) / 2, label, fs=12.5, c=GREY,
               ha="center", rot=90)

    # ------------------------------------------------------------------ out --
    def save(self, stem):
        self.fig.savefig(stem + ".png", dpi=100, facecolor="white")
        to_pptx(self.ax, self.W, self.H, stem + ".pptx", fig=self.fig)
        print("wrote %s.png and %s.pptx  (%d x %d)" % (stem, stem, self.W, self.H))
