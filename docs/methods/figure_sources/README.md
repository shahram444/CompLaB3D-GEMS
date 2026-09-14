# Rebuilding Figures 4 and 5 of the Method 3 guide

Figures 4 and 5 of `Method_3_Symbolic_guide.docx` are generated, not drawn by
hand. This directory holds what generates them, so a correction can be made once,
in the place the picture actually comes from, instead of being made twice and
drifting.

```
draw.py     records what a figure script draws and emits it as PowerPoint shapes
figA.py     Figure 4, the symbolic path as code
figB.py     Figure 5, inside the search
insert.py   puts both into the guide, with captions, and renumbers what follows
```

## Running them

```bash
python3 figA.py          # writes figA.png and figA.pptx
python3 figB.py          # writes figB.png and figB.pptx
```

The two-slide deck that ships beside the guide as
`Method_3_Symbolic_guide-flowcharts.pptx` comes from handing both recordings to
`to_deck`:

```python
import runpy
from draw import to_deck
recs = [(g["ax"], g["W"], g["H"], g["fig"])
        for g in (runpy.run_path("figA.py"), runpy.run_path("figB.py"))]
to_deck(recs, "Method_3_Symbolic_guide-flowcharts.pptx")
```

Needs `matplotlib`, `python-pptx` and `lxml`. Carlito and DejaVu Sans Mono are
the fonts the PNGs are drawn with; the slides ask for Calibri and Consolas,
which are what the other Method decks use.

## Putting them back into the guide

The guide is hand-edited and is never rebuilt from a script. `insert.py` unzips
it, adds the two pictures with their captions and one paragraph each, renumbers
the figures that now come after them, and zips it back. Its paths point at a
working copy; read the top of the file before running it.

If a figure changes but its size does not, replacing `word/media/image10.png` or
`image11.png` inside the docx is enough and `insert.py` is not needed at all.

## Why the PowerPoint is worth the machinery

A figure that ships only as a picture is a figure nobody but its author can
correct, and every other Method guide already carries a deck of real shapes
beside it. `draw.py` keeps that promise without asking the figure scripts to be
written twice: `Rec` wraps the matplotlib axes, passes every call through to it
and keeps a copy, and `to_pptx` walks the copy afterwards. Rounded boxes become
rounded rectangles, the decision diamonds become diamonds, arrows become
connectors with a triangular head, and every label becomes a text box.

Two details are worth knowing before changing anything here.

**The slide size is not arbitrary.** It is the figure's pixel size divided by a
hundred, in inches, so that a hundred pixels is one inch. At matplotlib's dpi of
100 that makes a point in the drawing a point on the slide, which means no font
size has to be converted and none of the labels have to be nudged back into
place after the fact.

**Text boxes are sized to their text.** The width comes from measuring the
rendered label rather than from a generous guess. A box exactly as wide as its
contents cannot be aligned wrongly, and the first version of this, which used a
fixed wide box and trusted the alignment setting, came out centred in
LibreOffice and left-aligned in PowerPoint.

## If you edit the .pptx instead

That is fine for a one-off tweak, but the next run of `figA.py` will not know
about it. Either fold the change back into the script, or accept that the deck
and the script have parted company and say so where it matters.
