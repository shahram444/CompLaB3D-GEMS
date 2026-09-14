import runpy, matplotlib
matplotlib.use("Agg")
from draw import to_deck
figs = []
for stem in ("m1a", "m1c", "m1b"):
    g = runpy.run_path(stem + ".py")
    if "F" in g:
        F = g["F"]; figs.append((F.ax, F.W, F.H, F.fig))
    else:
        figs.append((g["ax"], g["W"], g["H"], g["fig"]))
p = to_deck(figs, "/home/claude/deploy/docs/methods/Method_1_FBA_guide-flowcharts.pptx")
print("wrote", p)
