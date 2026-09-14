#!/usr/bin/env python3
"""Method 2 guide: the same pass Methods 3, 4 and 5 had.

Out: the end-to-end band chart, the two-phases plate, the offline pipeline
plate, the sample-split plate, the run-time pipeline plate and the two-forms
comparison. Each is now said better by one of the two new flowcharts, and the
"how to read Figure N" sections that went with them. Sections 8, 10 and 11 go,
as 8 to 11 did in Method 3.

In: the path as code, and one evaluation of the shipped network drawn on two
real voxels. Section 7 gains the five equations the guide never wrote down and
the order the code applies all thirteen in.
"""
import os
import re
import shutil
import sys

sys.path.insert(0, "/home/claude/fig3")
from edit_m4 import (P_HEAD, P_GAP, esc, run, eq, where, body, mono, block)  # noqa

ROOT = ("/tmp/claude-0/-home-claude/907c0fba-0d6e-5185-95b4-57fae18fda84"
        "/scratchpad/m2/u")
DOC = os.path.join(ROOT, "word", "document.xml")
RELS = os.path.join(ROOT, "word", "_rels", "document.xml.rels")
MEDIA = os.path.join(ROOT, "word", "media")


P_EQ2 = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
         '<w:spacing w:before="40" w:after="140"/><w:jc w:val="center"/></w:pPr>%s</w:p>')


def run2(text, kind="n", sz=None):
    """run(), plus 'p' for a superscript."""
    if kind != "p":
        return run(text, kind, sz)
    pr = '<w:rPr><w:i/><w:iCs/>%s<w:vertAlign w:val="superscript"/></w:rPr>' % (
        '<w:sz w:val="%d"/><w:szCs w:val="%d"/>' % (sz, sz) if sz else "")
    return '<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r>' % (pr, esc(text))


def eq2(parts, label):
    return P_EQ2 % ("".join(run2(t, k, 24) for t, k in parts)
                    + run("        (%s)" % label, "n", 24))



P_MONO_END = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
              '<w:spacing w:line="252" w:lineRule="auto"/></w:pPr><w:r><w:rPr>'
              '<w:rFonts w:ascii="Consolas" w:eastAsia="Consolas" '
              'w:hAnsi="Consolas" w:cs="Consolas"/><w:sz w:val="17"/>'
              '<w:szCs w:val="17"/></w:rPr><w:t xml:space="preserve">%s</w:t>'
              '</w:r></w:p>')
P_GAP_END = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
             '<w:spacing w:after="60"/></w:pPr></w:p>')


def monoblock(lines):
    """A run of monospace lines that stays together but may be broken after."""
    return mono(lines[:-1]) + P_MONO_END % lines[-1]



def para_start(doc, needle):
    """Offset of the <w:p> that contains `needle`."""
    i = doc.index(needle)
    return max(doc.rfind("<w:p ", 0, i), doc.rfind("<w:p>", 0, i))


def cut(doc, first, last):
    """Remove everything from the paragraph holding `first` up to `last`.

    Deleting a section paragraph by paragraph is wrong twice over when the
    section's content is a table: the table shell survives with every cell
    emptied, and an empty <w:tc> is invalid OOXML that Word refuses to open.
    Worse, cell paragraphs are often byte-identical to each other, so
    replace(p, "", 1) removes the FIRST match in the document rather than the
    one intended, and takes an unrelated table apart.
    """
    return doc[:para_start(doc, first)] + doc[para_start(doc, last):]


CAP_A = ("Figure 2. The surrogate path as code, from a sweep of linear programs "
         "to a change in one voxel. The side panels open the two steps that are "
         "otherwise taken on trust: what one row of the sweep contains, and what "
         "the run-time file actually gives back.")

TXT_A = ("The three bands are three different time scales. The yellow band runs "
         "once, on a workstation, before any simulation exists, and it is the "
         "only place in this whole method that a linear program is solved. The "
         "blue band runs once when the solver starts. The green band runs in "
         "every wet voxel that holds biomass, at every reaction step, and it is "
         "the box drawn in red in Figure 1. Each band carries the name of the "
         "source file it lives in, so a reader who wants the code can go "
         "straight to it. The dashed line across the page is the handover, and "
         "the .srg file is the only thing that crosses it: nothing in the run "
         "can reach back to the metabolic model, which is the whole point and "
         "also the whole cost. Two boxes are drawn in red because they are the "
         "two places where a mistake does not announce itself. A file that "
         "arrives without a complete training box is refused at start-up rather "
         "than loaded, because a network asked about conditions it was never "
         "shown does not fail, it answers; and the clamp in the green band is "
         "counted once per evaluation rather than once per input, so the "
         "percentage printed at the end of the run means what a reader would "
         "assume it means.")

CAP_B = ("Figure 6. One evaluation of the shipped network, on two voxels of the "
         "same run. Every number is computed from the weights, the two scalings "
         "and the training box written into surrogateModel.hh, at two "
         "compositions the shipped case passes through.")

TXT_B = ("The two voxels differ in one thing only: the fed one holds the acetate "
         "and iron the case starts with, and the starved one holds what is left "
         "after a colony has drawn them down. That difference reaches the "
         "network as a change of about one and a quarter in the acetate uptake "
         "bound, and the network turns it into the difference between growing "
         "and not growing at all. The two pictures of the units are the reason. "
         "Every unit puts its total through the same S-shaped curve, which "
         "flattens at minus one and plus one, and a unit that has reached one of "
         "those flats can no longer tell one input from another. In the starved "
         "voxel six of the last layer's ten units are pinned there, and the "
         "output unit, which applies no curve at all, lands on exactly the "
         "bottom of the range the fit was scaled into. De-scaled, that is a "
         "small negative number rather than a small growth rate, which is why "
         "the floor of equation (10) is not housekeeping: without it the "
         "solver would be handed a rate that grows biomass backwards. The panel "
         "at the foot of the figure sweeps the same weights across the whole "
         "fitted range of the first input and shows what the network learned as "
         "a shape rather than as a number, which is the one check worth making "
         "on any surrogate before trusting a run whose voxels sit near the step.")

EQNS = (
    P_HEAD % "7.3  The five equations the sections above do not write down"
    + body("Equations (1) to (6) carry a value from a pair of concentrations to "
           "a number in physical units. Five more turn that number into "
           "something the transport solver can use, and they are the ones a "
           "reader is most likely to assume rather than check.")
    + body("A file may record that the fit was done on the logarithm of the "
           "first output, which is worth doing when growth spans decades. When "
           "it does, the first output and only the first is raised back:")
    + eq2([("o", "i"), ("0", "s"), (" ← 10", "i"), ("o0", "p"),
           ("        applied only when the file carries logoutput 1", "i")], "9")
    + body("Then the floor, which is the shortest line in the path and does the "
           "most work:")
    + eq2([("μ = o", "i"), ("0", "s"), ("   if   o", "i"), ("0", "s"),
           (" > 10", "i"), ("−8", "p"), (" ,      and   μ = 0   otherwise", "i")], "10")
    + where([("the comparison is written the way round that also catches a value "
              "which is not a number, since any comparison against one is false. "
              "Growth is floored; a flux is not, because a flux may legitimately "
              "be zero or negative and clamping one would corrupt it.", "n")])
    + body("If the thermodynamic gate of Method 5 is switched on and this "
           "organism carries a reaction block, every rate the network returned "
           "and the growth with them are scaled by the factor that method "
           "computes:")
    + eq([("v", "i"), ("k", "s"), (" ← F", "i"), ("T", "s"), (" v", "i"),
          ("k", "s"), (" ,      μ ← F", "i"), ("T", "s"), (" μ", "i")], "11")
    + where([("F", "i"), ("T", "s"),
             (" is one when the gate is off or the organism has no block, so an "
              "existing case is bit-identical. The factor is applied to the "
              "answer rather than to the inputs, because the network was fitted "
              "on a metabolic model that knows nothing about the local energy "
              "balance: what it predicts is the rate the organism could run if "
              "energy allowed, and the factor is the fraction of that which "
              "energy does allow.", "n")])
    + body("The growth rate then becomes a biomass increment, or, when it is "
           "zero, a loss at the organism's own decay coefficient:")
    + eq([("ΔB = μ B Δt      if μ > 0 ,      ΔB = − b B Δt      "
           "if μ = 0", "i")], "12")
    + where([("b", "i"), (" the decay coefficient. The test is on the growth "
             "rate, not on the biomass, so a voxel whose network answer fell to "
             "the floor loses biomass at the same rate whether it is starving "
             "slightly or completely. ", "n")])
    + body("And the substrate draw, which on the growth-only path is the bound "
           "of equation (1) rather than anything the network said:")
    + eq([("ΔC", "i"), ("k", "s"), (" = − Δt Σ", "i"), ("m", "s"),
          (" v", "i"), ("m,k", "s"), (" B", "i"), ("m", "s"), (" M", "i"),
          (" ,      held under what the voxel holds", "i")], "13")
    + where([("the sum runs over every organism present in the voxel and the "
              "limit is applied once, to the total. Applying it inside the loop "
              "instead lets two organisms each take the whole supply, and the "
              "lattice goes negative on the first step. ", "n"),
             ("M", "i"), (" converts biomass to the units the lattice carries.", "n")])
)

ORDER = (
    P_HEAD % "7.4  The order in which the code applies them"
    + body("Read down. The shape of this list is the method: everything "
           "expensive is above the last group, and the last group is the only "
           "part that runs more than once.", keep=True)
    + P_GAP
    + monoblock([
        "before the run                 choose the substrates and the sweep range",
        "  generateTrainingData.py   7  place the points, evenly in every decade",
        "                            1  one linear program per point, same bounds",
        "                               an infeasible solve is kept, as zeros",
        "",
        "  trainSurrogate.py         2  the input scaling, a pair per input",
        "                            6  the output scaling, a pair per output",
        "                          3,4  the architecture the fit will use",
        "                            8  minimise the error, from five restarts",
        "                               write the .srg: weights, scalings, box",
        "",
        "once, at start-up              read it; check each layer against its shape",
        "  complab3d_surrogate.hh       refuse it if the training box is incomplete",
        "                               bind each input to a substrate by name",
    ])
    + P_GAP_END
    + monoblock([
        "per organism, per voxel,",
        "per reaction step           1  the uptake bound Monod allows here",
        "  complab3d_processors_        hold each bound inside the box, and count it",
        "  surrogate.hh              2  rescale both inputs onto minus one to plus one",
        "                          3,4  four hidden layers of ten",
        "                            5  the output layer: weights, a bias, no curve",
        "                            6  back into physical units",
        "                            9  ten to it, if the fit was in logarithms",
        "                           10  the floor: not above 1e-8 becomes exactly zero",
        "                           11  the thermodynamic factor, if the gate is on",
        "                           12  grow the biomass, or decay it",
        "                           13  the substrate draw, summed then floored",
    ])
    + P_GAP_END
    + body("Two matrix products of ten by ten, two of ten by two and one of one "
           "by ten, forty hyperbolic tangents, and a dozen lines of arithmetic "
           "around them. That is the entire run-time cost of this path, and it "
           "is why a case that would not finish with a linear program in every "
           "voxel finishes with this one. Everything that made the cost "
           "worthwhile is in the first two groups, and none of it is repeated.")
)


def main():
    doc = open(DOC, encoding="utf-8").read()
    rels = open(RELS, encoding="utf-8").read()

    for src, name, rid in (("m2a.png", "imageA.png", "rId90"),
                           ("m2b.png", "imageB.png", "rId91")):
        shutil.copy("/home/claude/fig3/" + src, os.path.join(MEDIA, name))
        rels = rels.replace("</Relationships>",
                            '<Relationship Id="%s" Type="http://schemas.openxml'
                            'formats.org/officeDocument/2006/relationships/image"'
                            ' Target="media/%s"/></Relationships>' % (rid, name))

    paras = re.findall(r'<w:p [^>]*>.*?</w:p>|<w:p/>', doc, re.S)

    def txt(p):
        return re.sub(r'<[^>]+>', '', p).strip()

    def find(prefix):
        for i, p in enumerate(paras):
            if txt(p).startswith(prefix):
                return i
        raise SystemExit("not found: " + prefix)

    # anchors, taken before anything is deleted
    fig3_at = paras[find("Figure 3. The two phases") - 2]

    # --- out: the six plates that the two new flowcharts now carry ----------
    for cap in ("Figure 2. How this method works",
                "Figure 3. The two phases",
                "Figure 6. The offline pipeline",
                "Figure 7. How the samples are divided",
                "Figure 8. The run-time pipeline"):
        i = find(cap)
        doc = doc.replace(paras[i - 1], "", 1)
        doc = doc.replace(paras[i], "", 1)

    # --- out: the three "how to read" sections ------------------------------
    for head, stop in (("5.2 How to read Figures 3 and 4", "6. Building the network"),
                       ("6.1 How to read Figure 6", "6.2 Where the sample points go"),
                       ("7.1 How to read Figure 8", "7.2 The range check")):
        doc = cut(doc, head, stop)

    # --- out: section 8, and sections 10 and 11 -----------------------------
    doc = cut(doc, "10. Running it", "12. References")

    # --- renumber what is left ---------------------------------------------
    for old, new in (("6.2 Where the sample points go", "6.1 Where the sample points go"),
                     ("6.3 What the fit minimises", "6.2 What the fit minimises"),
                     ("6.4 Which number goes in the paper", "6.3 Which number goes in the paper"),
                     ("7.2 The range check", "7.1 The range check"),
                     ("7.3 The ordering", "7.2 The ordering"),
                     ("12. References", "@@10@@ References")):
        doc = doc.replace(old, new, 1)
    doc = doc.replace("The whole of Section 7.2 is about",
                      "The whole of Section 7.1 is about")

    # the one assumption keeps its equation, but it now comes after the five
    # this pass adds, so it takes the last number rather than the first free one
    doc = doc.replace('<w:t xml:space="preserve">        (9)</w:t>',
                      '<w:t xml:space="preserve">        (14)</w:t>', 1)
    doc = doc.replace("equation (9)", "equation (14)")
    doc = doc.replace("Equation (9)", "Equation (14)")

    # figure numbers: 4 and 5 move up to 3 and 4, and the two new ones take
    # 2 and 5. Renumber BEFORE inserting, or the new captions get renumbered.
    doc = doc.replace("Figure 4. The single-output form", "@@F3@@ The single-output form")
    doc = doc.replace("Figure 5. The multi-output form", "@@F4@@ The multi-output form")
    doc = doc.replace("Figure 9. The two forms compared", "@@F5@@ The two forms compared")
    doc = doc.replace("Figure 4 has one", "Figure 3 has one")
    doc = doc.replace("Figure 5 has", "Figure 4 has")
    doc = re.sub(r"@@F(\d)@@ ", lambda m: "Figure %s. " % m.group(1), doc)
    doc = re.sub(r"@@(\d+)@@ ", lambda m: "%s. " % m.group(1), doc)

    # --- in: the two flowcharts, and the equations --------------------------
    h = next(p for p in re.findall(r'<w:p [^>]*>.*?</w:p>', doc, re.S)
             if txt(p).startswith("8. The one assumption"))
    doc = doc.replace(h, EQNS + ORDER + h, 1)

    doc = doc.replace(fig3_at, fig3_at + block("rId90", 621, 2900, CAP_A, TXT_A,
                                               cx=4450000), 1)
    sec9_at = next(p for p in re.findall(r'<w:p [^>]*>.*?</w:p>', doc, re.S)
                   if txt(p).startswith("9.2 Step 1: from concentrations"))
    doc = doc.replace(sec9_at, block("rId91", 622, 3350, CAP_B, TXT_B,
                                     cx=4000000) + sec9_at, 1)

    # Safety net. Any cell left without a paragraph gets an empty one, because
    # a <w:tc> with no block-level child is invalid OOXML: Word declines to open
    # the file and LibreOffice repairs it silently, which is how such a file
    # gets shipped without anybody noticing.
    def fill(m):
        c = m.group(0)
        if "<w:p" in c:
            return c
        return c.replace("</w:tc>", "<w:p/></w:tc>")

    doc, n = re.subn(r'<w:tc>.*?</w:tc>', fill, doc, flags=re.S)
    bad = [c for c in re.findall(r'<w:tc>.*?</w:tc>', doc, re.S) if "<w:p" not in c]
    assert not bad, "%d table cells still have no paragraph" % len(bad)

    open(DOC, "w", encoding="utf-8").write(doc)
    open(RELS, "w", encoding="utf-8").write(rels)
    print("figures:", sorted(set(re.findall(r'Figure (\d)\. ', doc))))
    print("sections:", sorted(set(re.findall(r'>(\d{1,2})\.  ?[A-Z]', doc)), key=int))
    print("7.x:", sorted(set(re.findall(r'>(7\.\d)  ?', doc))))
    print("9.x:", sorted(set(re.findall(r'>(9\.\d) [A-Z]', doc))))


if __name__ == "__main__":
    main()
