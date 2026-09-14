#!/usr/bin/env python3
"""Method 4 guide: the same pass Method 3 had.

Out: the end-to-end band chart, the two-halves-of-a-round plate, the six-steps
plate, the two-readouts plate and the offline-row plate. Each is now said better
by one of the two new flowcharts, and a guide that shows the same thing twice
teaches the reader to skim.

In: the path as code, and what one evaluation does to the graph. Section 5 gains
the equations that produced the file, which it never had, and the order the code
applies all of them in. Sections 7 to 9 go, as 8 to 11 did in Method 3.
"""
import os
import re
import shutil

ROOT = "/tmp/claude-0/-home-claude/907c0fba-0d6e-5185-95b4-57fae18fda84/scratchpad/m4/u"
DOC = os.path.join(ROOT, "word", "document.xml")
RELS = os.path.join(ROOT, "word", "_rels", "document.xml.rels")
MEDIA = os.path.join(ROOT, "word", "media")
EMU_W = 5334000

P_BODY = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>%s'
          '<w:spacing w:after="140" w:line="288" w:lineRule="auto"/></w:pPr>%s</w:p>')
P_WHERE = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
           '<w:spacing w:after="140" w:line="288" w:lineRule="auto"/>'
           '<w:ind w:left="260"/></w:pPr>%s</w:p>')
P_EQ = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
        '<w:spacing w:before="40" w:after="140"/><w:jc w:val="center"/></w:pPr>%s</w:p>')
P_HEAD = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
          '<w:spacing w:before="160" w:after="80"/></w:pPr>'
          '<w:r><w:rPr><w:b/><w:bCs/></w:rPr>'
          '<w:t xml:space="preserve">%s</w:t></w:r></w:p>')
P_MONO = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
          '<w:keepNext/><w:keepLines/><w:spacing w:line="252" w:lineRule="auto"/>'
          '</w:pPr><w:r><w:rPr><w:rFonts w:ascii="Consolas" w:eastAsia="Consolas" '
          'w:hAnsi="Consolas" w:cs="Consolas"/><w:sz w:val="17"/>'
          '<w:szCs w:val="17"/></w:rPr><w:t xml:space="preserve">%s</w:t></w:r></w:p>')
P_GAP = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
         '<w:keepNext/><w:spacing w:after="60"/></w:pPr></w:p>')
PIC = (
    '<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
    '<w:spacing w:before="140" w:after="40"/><w:jc w:val="center"/></w:pPr>'
    '<w:r><w:rPr><w:noProof/></w:rPr><w:drawing>'
    '<wp:inline distT="0" distB="0" distL="0" distR="0">'
    '<wp:extent cx="%(cx)d" cy="%(cy)d"/><wp:effectExtent l="0" t="0" r="0" b="0"/>'
    '<wp:docPr id="%(id)d" name="Picture %(id)d"/><wp:cNvGraphicFramePr>'
    '<a:graphicFrameLocks xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" noChangeAspect="1"/>'
    '</wp:cNvGraphicFramePr>'
    '<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
    '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
    '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
    '<pic:nvPicPr><pic:cNvPr id="0" name=""/><pic:cNvPicPr/></pic:nvPicPr>'
    '<pic:blipFill><a:blip r:embed="%(rid)s"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
    '<pic:spPr bwMode="auto"><a:xfrm><a:off x="0" y="0"/><a:ext cx="%(cx)d" cy="%(cy)d"/></a:xfrm>'
    '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
    '</pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>')
CAP = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
       '<w:spacing w:after="120"/><w:jc w:val="center"/></w:pPr>'
       '<w:r><w:rPr><w:i/><w:iCs/><w:sz w:val="18"/><w:szCs w:val="18"/></w:rPr>'
       '<w:t xml:space="preserve">%s</w:t></w:r></w:p>')


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def run(text, kind="n", sz=None):
    rpr = []
    if kind in ("i", "s"):
        rpr.append("<w:i/><w:iCs/>")
    if kind == "b":
        rpr.append("<w:b/><w:bCs/>")
    if sz:
        rpr.append('<w:sz w:val="%d"/><w:szCs w:val="%d"/>' % (sz, sz))
    if kind == "s":
        rpr.append('<w:vertAlign w:val="subscript"/>')
    pr = "<w:rPr>%s</w:rPr>" % "".join(rpr) if rpr else ""
    return '<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r>' % (pr, esc(text))


def eq(parts, label):
    return P_EQ % ("".join(run(t, k, 24) for t, k in parts)
                   + run("        (%s)" % label, "n", 24))


def where(parts):
    return P_WHERE % (run("where ", "b") + "".join(run(t, k) for t, k in parts))


def body(text, keep=False):
    return P_BODY % ("<w:keepNext/><w:keepLines/>" if keep else "", run(text))


def mono(lines):
    return "".join(P_MONO % ln for ln in lines)


def block(rid, pid, px_h, cap, txt, cx=EMU_W):
    cy = int(round(cx * px_h / 1900.0))
    return (PIC % dict(cx=cx, cy=cy, id=pid, rid=rid) + CAP % esc(cap)
            + body(txt))



def para_start(doc, needle):
    """Offset of the <w:p> that contains `needle`."""
    i = doc.index(needle)
    return max(doc.rfind("<w:p ", 0, i), doc.rfind("<w:p>", 0, i))


def cut(doc, first, last):
    """Remove everything from the paragraph holding `first` up to `last`.

    Deleting a section paragraph by paragraph is wrong twice over when the
    section's content is a table: the table shell survives with every cell
    emptied, and an empty <w:tc> is invalid OOXML that Word refuses to open
    while LibreOffice repairs it silently. Cell paragraphs are also often
    byte-identical, so replace(p, "", 1) can remove the first match in the
    document rather than the intended one and take another table apart.
    """
    return doc[:para_start(doc, first)] + doc[para_start(doc, last):]


# ------------------------------------------------------------- new text -----
CAP_A = ("Figure 2. The graph network path as code, from a table of states to a "
         "change in one voxel. The side panels open the two steps that are "
         "otherwise taken on trust: one epoch of the fit, and one evaluation in "
         "one voxel.")

TXT_A = ("The three bands are the three time scales this method lives on. The "
         "yellow band happens once, on a workstation, and none of it is inside "
         "CompLaB3D. The blue band happens once when the solver starts. The "
         "green band happens in every wet voxel that holds biomass, at every "
         "reaction step, and is the box drawn in red in Figure 1. Each band "
         "carries the name of the source file it lives in. The dashed line "
         "across the page is the handover, and the .gnn file is the only thing "
         "that crosses it. Two boxes are drawn in red because they are the two "
         "places where a mistake does not announce itself: a species name the "
         "run does not recognise stops the run at start-up rather than binding "
         "the network to the wrong lattice, and a concentration arriving from "
         "outside the box the network was fitted over is pulled back to the "
         "edge and counted.")

CAP_B = ("Figure 4. Inside the graph network: what one evaluation does, drawn on "
         "the shipped anaerobic methane oxidation case. Four species, one "
         "reaction, six numbers per node, two rounds.")

TXT_B = ("The flow down the left is drawn without colour so that the only "
         "coloured objects are the network itself. A species node is blue, a "
         "reaction node is yellow, and a number the fit chose is green. The "
         "numbers on the edges are the stoichiometry and they are not fitted: "
         "they are read from the file, and they decide which messages exist at "
         "all. That is the difference between this path and a dense network "
         "told about stoichiometry through a penalty term in its loss. Steps 3 "
         "and 4 together are one round, and the loop back from the question "
         "under step 4 is what the rounds line in the file controls. Step 5 is "
         "the step the whole method turns on: one number read off each reaction "
         "node rather than one off each species node, which is what makes the "
         "ratios between the species rates exact instead of approximately "
         "right.")

S_OFFLINE = (
    P_HEAD % "5.1  The offline equations: what the fit minimises"
    + body("The five subsections after this one are what the solver evaluates. "
           "This one is what produced the file it evaluates. None of it runs "
           "during a simulation, and it is here because the numbers in a .gnn "
           "file cannot be read honestly without knowing what was minimised to "
           "arrive at them.")
    + body("Every concentration column is put on one scale first, so that a "
           "species measured in millimoles and one measured in moles carry the "
           "same weight in what follows:")
    + eq([("x", "i"), ("i", "s"), (" = ( C", "i"), ("i", "s"), (" − lo", "i"),
          ("i", "s"), (" ) · 2 / ( hi", "i"), ("i", "s"), (" − lo", "i"),
          ("i", "s"), (" ) − 1", "i")], "S1")
    + where([("lo", "i"), ("i", "s"), (" and ", "n"), ("hi", "i"), ("i", "s"),
             (" the smallest and largest value of species i anywhere in the "
              "training table. Those two numbers are written into the file as "
              "the trainmin and trainmax lines, and they are the box the run "
              "enforces later in equation (1).", "n")])
    + body("The rates are then turned into extents, one per reaction, by the "
           "pseudo-inverse of the stoichiometric matrix:")
    + eq([("Ξ = Y S", "i"), ("+", "s"), ("      and      T = Ξ / max | Ξ |", "i")],
         "S2")
    + where([("Y", "i"), (" the table of species rates; ", "n"), ("S", "i"),
             ("+", "s"),
             (" the pseudo-inverse of the stoichiometry; the division is a pure "
              "scale with no offset, because an offset does not commute with "
              "the sum in equation (8) and would destroy the exactness that "
              "readout exists for. If more than five per cent of the training "
              "rates do not lie in the column space of the stoichiometry the "
              "fit says so, because a rate vector that is not a set of reaction "
              "rates cannot be represented at all.", "n")])
    + body("What is minimised is then one line, with no weighting, no "
           "regularisation and no penalty term for the chemistry:")
    + eq([("L = mean over samples and outputs of ( network − T )", "i"),
          ("2", "s")], "S3")
    + body("and the gradient is taken by finite difference, on two dozen weights "
           "drawn at random from each array every epoch:")
    + eq([("g", "i"), ("k", "s"), (" = ( L( w", "i"), ("k", "s"),
          (" + h ) − L( w", "i"), ("k", "s"), (" − h ) ) / 2h ,"
                                               "        h = 10", "i"),
          ("−5", "i")], "S4")
    + body("Adam then moves every weight, and five epochs without improvement "
           "halve the step and rewind to the best weights seen. There is no "
           "held-out set: the best training loss is what is kept. That is worth "
           "saying plainly, because it means the fit cannot tell you whether it "
           "generalises. The honest check is the one example 18 performs "
           "afterwards, on three hundred samples drawn fresh from the law the "
           "training table came from.")
)

ORDER = (
    P_HEAD % "5.7  The order in which the code applies them"
    + body("Read down. The three groups are three different time scales, and the "
           "gap between them is the cost argument for this whole path.", keep=True)
    + P_GAP
    + mono([
        "once, offline               S1   put every concentration on one scale",
        "  tools/method_4_graphnet/train_graphnet.py   S2   project the rates onto one extent per reaction",
        "                            S3   the loss: plain mean squared error",
        "                            S4   the gradient, by finite difference",
        "",
        "once, at start-up                parse the file and check every array's shape",
        "  complab3d_graphnet.hh          bind every species name to its lattice",
        "                                 store the training box and the unit scale",
        "",
        "per voxel, per reaction step",
        "  complab3d_processors_       1  clamp each input to the training box",
        "    graphnet.hh               2  scale it, and encode it onto its species node",
        "                              3  species speak to reactions",
        "                              4  reactions speak back, and again",
        "                              7  read one extent off each reaction node",
        "                              8  build the species rates from the stoichiometry",
        "                              9  growth, off the average species node",
        "                             10  divide by 3600, gate, multiply by biomass",
        "                             11  hand the increments to transport, floored",
    ])
    + P_GAP
    + body("Nothing in the first group runs during a simulation, and nothing in "
           "the second runs more than once. Only the third is inside the loop, "
           "and for the shipped four-species network it is eighty-four "
           "hyperbolic tangents and a few hundred multiplications. That is the "
           "trade this path makes against a linear program solved per voxel: the "
           "expensive work is the work that never runs again.")
)


def main():
    doc = open(DOC, encoding="utf-8").read()
    rels = open(RELS, encoding="utf-8").read()

    for src, name, rid in (("m4a.png", "imageA.png", "rId90"),
                           ("m4b.png", "imageB.png", "rId91")):
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

    # --- the five plates the flowcharts replace --------------------------------
    for cap in ("Figure 2. How this method works",
                "Figure 4. The two halves of one round",
                "Figure 5. The six steps of one evaluation",
                "Figure 6. The two readouts drawn side by side",
                "Figure 7. The offline row and the file"):
        i = find(cap)
        doc = doc.replace(paras[i - 1], "", 1)
        doc = doc.replace(paras[i], "", 1)

    # --- sections 7, 8 and 9; the references become section 7 ------------------
    doc = cut(doc, "7.  Fitting it offline", "10.  References")
    doc = doc.replace("10.  References", "7.  References")

    # --- two sentences pointed into what has gone ------------------------------
    doc = doc.replace("which is why Section 7 recommends it",
                      "which is why it is the default and the other mode exists "
                      "only to read files written before it did")
    doc = doc.replace("and it is the reason Section 8 exists",
                      "and it is the reason the species names are the binding "
                      "and the training box is enforced rather than advised")

    # --- section 5 gains the offline equations and the running order -----------
    for old, new in ((5, 6), (4, 5), (3, 4), (2, 3), (1, 2)):
        doc = doc.replace("5.%d  " % old, "@@5.%d@@  " % new)
    doc = re.sub(r"@@5\.(\d)@@  ", lambda m: "5.%s  " % m.group(1), doc)

    h = next(p for p in re.findall(r'<w:p [^>]*>.*?</w:p>', doc, re.S)
             if txt(p).startswith("5.2  The clamp"))
    doc = doc.replace(h, S_OFFLINE + h, 1)
    h6 = next(p for p in re.findall(r'<w:p [^>]*>.*?</w:p>', doc, re.S)
              if txt(p).startswith("6.  A worked example"))
    doc = doc.replace(h6, ORDER + h6, 1)

    # --- the two new figures ---------------------------------------------------
    a_at = next(p for p in re.findall(r'<w:p [^>]*>.*?</w:p>', doc, re.S)
                if txt(p).startswith("The whole run-time cost")
                or txt(p).startswith("One time step of CompLaB3D"))
    doc = doc.replace(a_at, a_at + block("rId90", 601, 2630, CAP_A, TXT_A,
                                         cx=4300000), 1)
    b_at = next(p for p in re.findall(r'<w:p [^>]*>.*?</w:p>', doc, re.S)
                if txt(p).startswith("4.5  What one voxel actually does"))
    doc = doc.replace(b_at, b_at + block("rId91", 602, 2630, CAP_B, TXT_B,
                                         cx=4300000), 1)

    # Safety net: a <w:tc> with no block-level child is invalid OOXML.
    def fill(m):
        c = m.group(0)
        return c if "<w:p" in c else c.replace("</w:tc>", "<w:p/></w:tc>")

    doc = re.sub(r'<w:tc>.*?</w:tc>', fill, doc, flags=re.S)
    bad = [c for c in re.findall(r'<w:tc>.*?</w:tc>', doc, re.S) if "<w:p" not in c]
    assert not bad, "%d table cells still have no paragraph" % len(bad)

    open(DOC, "w", encoding="utf-8").write(doc)
    open(RELS, "w", encoding="utf-8").write(rels)
    print("figures:", sorted(set(re.findall(r'Figure (\d)\. ', doc))))
    print("sections:", sorted(set(re.findall(r'>(\d{1,2})\.  [A-Z]', doc)), key=int))
    print("5.x:", sorted(set(re.findall(r'>(5\.\d)  ', doc))))


if __name__ == "__main__":
    main()
