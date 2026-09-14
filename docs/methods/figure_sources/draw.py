#!/usr/bin/env python3
"""
Record what a figure script draws, then emit it twice: once as a PNG through
matplotlib and once as native PowerPoint shapes.

A figure that ships only as a picture is a figure nobody but its author can
correct. The existing Method guides carry a .pptx of real shapes beside every
guide for exactly that reason, and these two new figures had to join them.

The trick is that nothing in the figure scripts changes. `Rec` wraps the
matplotlib Axes, forwards every call to it, and keeps a copy. `to_pptx` then
walks the copy.

The pixel-to-slide mapping is chosen so that no font size has to be converted:
100 px is one inch, so a matplotlib point at dpi 100 is a PowerPoint point.
"""
import matplotlib.colors as mcolors
from matplotlib.patches import (FancyBboxPatch, FancyArrowPatch, Polygon,
                                Rectangle)

from pptx import Presentation
from pptx.util import Emu, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from lxml import etree

EMU = 9144                      # one pixel, with 100 px to the inch
NS = "http://schemas.openxmlformats.org/drawingml/2006/main"


# ------------------------------------------------------------------ record ---
class Rec(object):
    """Everything the script draws goes to matplotlib and into a list."""

    def __init__(self, ax):
        self.ax = ax
        self.items = []

    def add_patch(self, p):
        self.ax.add_patch(p)
        self.items.append(("patch", p))
        return p

    def text(self, x, y, s, **k):
        h = self.ax.text(x, y, s, **k)
        self.items.append(("text", (x, y, s), dict(k), h))
        return h

    def plot(self, xs, ys, **k):
        self.ax.plot(xs, ys, **k)
        self.items.append(("plot", (list(xs), list(ys)), dict(k)))

    def __getattr__(self, name):
        return getattr(self.ax, name)


# ------------------------------------------------------------------- pptx ----
def _hex(c):
    if c is None:
        return None
    if isinstance(c, (tuple, list)):
        if len(c) == 4 and c[3] == 0:
            return None
        r, g, b = c[0], c[1], c[2]
    else:
        if c in ("none", "None"):
            return None
        r, g, b = mcolors.to_rgb(c)
    return RGBColor(int(round(r * 255)), int(round(g * 255)), int(round(b * 255)))


def _dashed(style):
    if style in (None, "-", "solid"):
        return False
    if isinstance(style, tuple):
        return True
    return style in ("--", ":", "-.", "dashed", "dotted", "dashdot")


def _line(shape, color, lw, dash=False, head=False):
    ln = shape.line
    if color is None:
        ln.fill.background()
        return
    ln.color.rgb = color
    ln.width = Pt(max(0.5, lw * 0.75))
    el = ln._get_or_add_ln()
    if dash:
        d = etree.SubElement(el, "{%s}prstDash" % NS)
        d.set("val", "dash")
    if head:
        e = etree.SubElement(el, "{%s}tailEnd" % NS)
        e.set("type", "triangle")
        e.set("w", "med")
        e.set("len", "med")


def _fill(shape, color):
    if color is None:
        shape.fill.background()
    else:
        shape.fill.solid()
        shape.fill.fore_color.rgb = color


def _box(slide, kind, x, y, w, h, adj=None):
    s = slide.shapes.add_shape(kind, Emu(int(x * EMU)), Emu(int(y * EMU)),
                               Emu(int(max(1, w) * EMU)), Emu(int(max(1, h) * EMU)))
    s.shadow.inherit = False
    if adj is not None and len(s.adjustments):
        s.adjustments[0] = adj
    return s


def _rounding(p):
    """The corner radius a FancyBboxPatch was built with, in pixels."""
    bs = p.get_boxstyle()
    return float(getattr(bs, "rounding_size", 0) or 0)


def _connector(slide, x0, y0, x1, y1, color, lw, dash, head):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Emu(int(x0 * EMU)),
                                   Emu(int(y0 * EMU)), Emu(int(x1 * EMU)),
                                   Emu(int(y1 * EMU)))
    _line(c, color, lw, dash, head)
    return c


def _textbox(slide, x, y, s, k, box):
    """`box` is the text's real rectangle, measured from the rendered figure.

    Sizing the shape to the text rather than to a generous guess is what makes
    the placement survive PowerPoint, LibreOffice and Keynote alike: a box that
    is exactly as wide as its contents cannot be aligned wrongly."""
    fs = float(k.get("fontsize", 12))
    ha = k.get("ha", k.get("horizontalalignment", "left"))
    rot = float(k.get("rotation", 0) or 0)
    lines = str(s).split("\n")
    align = {"center": PP_ALIGN.CENTER, "right": PP_ALIGN.RIGHT}.get(
        ha, PP_ALIGN.LEFT)

    if rot:
        w = max(len(ln) for ln in lines) * fs * 0.62 + 10
        h = len(lines) * fs * 1.45 + 6
        bx, by = x - w / 2, y - h / 2
    else:
        x0, y0, x1, y1 = box
        bx, by = x0 - 6, y0 - 4
        w, h = (x1 - x0) + 12, (y1 - y0) + 8

    tb = slide.shapes.add_textbox(Emu(int(bx * EMU)), Emu(int(by * EMU)),
                                  Emu(int(max(4, w) * EMU)),
                                  Emu(int(max(4, h) * EMU)))
    tf = tb.text_frame
    tf.word_wrap = False
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0

    fam = str(k.get("family", "Calibri"))
    if "Carlito" in fam:
        fam = "Calibri"
    if "Mono" in fam:
        fam = "Consolas"
    bold = k.get("fontweight", "normal") in ("bold", "heavy", 700)
    ital = k.get("style", "normal") == "italic"
    col = _hex(k.get("color", "#000000"))

    for i, line in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run()
        r.text = line
        f = r.font
        f.size = Pt(fs)
        f.bold = bold
        f.italic = ital
        f.name = fam
        if col is not None:
            f.color.rgb = col

    if rot:
        tb.rotation = (-rot) % 360
    return tb


def to_pptx(rec, W, H, path, fig=None):
    return to_deck([(rec, W, H, fig)], path)


def to_deck(figures, path):
    """One deck, one slide per figure. A presentation has a single slide size,
    so the tallest figure sets it and the others sit at the top of theirs."""
    W = max(f[1] for f in figures)
    H = max(f[2] for f in figures)
    prs = Presentation()
    prs.slide_width = Emu(int(W * EMU))
    prs.slide_height = Emu(int(H * EMU))
    for rec, _w, _h, fig in figures:
        _one(prs, rec, H, fig)
    prs.save(path)
    return path


def _one(prs, rec, H, fig):
    slide = prs.slides.add_slide(prs.slide_layouts[6])   # blank

    if fig is not None:
        fig.canvas.draw()
    ren = fig.canvas.get_renderer() if fig is not None else None
    H = fig.get_size_inches()[1] * fig.dpi if fig is not None else H

    def z(item):
        if item[0] == "patch":
            return item[1].get_zorder()
        return item[2].get("zorder", 3)

    for item in sorted(rec.items, key=z):
        if item[0] == "patch":
            p = item[1]
            fc = _hex(p.get_facecolor())
            ec = _hex(p.get_edgecolor())
            lw = float(p.get_linewidth() or 1.0)
            dash = _dashed(p.get_linestyle())
            if isinstance(p, FancyArrowPatch):
                (x0, y0), (x1, y1) = p._posA_posB
                _connector(slide, x0, y0, x1, y1, ec or fc, lw, dash, True)
            elif isinstance(p, Polygon):
                pts = p.get_xy()
                xs = [q[0] for q in pts]
                ys = [q[1] for q in pts]
                s = _box(slide, MSO_SHAPE.DIAMOND, min(xs), min(ys),
                         max(xs) - min(xs), max(ys) - min(ys))
                _fill(s, fc)
                _line(s, ec, lw, dash)
            elif isinstance(p, FancyBboxPatch):
                r = _rounding(p)
                w, h = p.get_width(), p.get_height()
                if r <= 1:
                    s = _box(slide, MSO_SHAPE.RECTANGLE, p.get_x(), p.get_y(), w, h)
                else:
                    adj = min(0.5, r / max(1.0, min(w, h)))
                    s = _box(slide, MSO_SHAPE.ROUNDED_RECTANGLE, p.get_x(),
                             p.get_y(), w, h, adj)
                _fill(s, fc)
                _line(s, ec, lw, dash)
            elif isinstance(p, Rectangle):
                s = _box(slide, MSO_SHAPE.RECTANGLE, p.get_x(), p.get_y(),
                         p.get_width(), p.get_height())
                _fill(s, fc)
                _line(s, ec, lw, dash)
        elif item[0] == "plot":
            (xs, ys), k = item[1], item[2]
            col = _hex(k.get("color", "#000000"))
            lw = float(k.get("lw", k.get("linewidth", 1.0)))
            dash = _dashed(k.get("ls", k.get("linestyle")))
            for i in range(len(xs) - 1):
                _connector(slide, xs[i], ys[i], xs[i + 1], ys[i + 1], col, lw,
                           dash, False)
        else:
            (x, y, s), k, handle = item[1], item[2], item[3]
            bb = handle.get_window_extent(renderer=ren)
            _textbox(slide, x, y, s, k,
                     (bb.x0, H - bb.y1, bb.x1, H - bb.y0))
