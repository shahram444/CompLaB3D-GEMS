#!/usr/bin/env python3
"""Method 1 guide: the same pass Methods 2, 3, 4 and 5 had.

Out: nine plates. The end-to-end band chart, the one-voxel plate, the two
pipeline plates for GLPK and COBRApy, the pipeline-with-numbers plate, the
chain and competition pipeline plates, and the two end-to-end charts for the
refinements. Three new flowcharts carry all of them, and a guide that shows the
same thing twice teaches the reader to skim. Section 10 goes, as the pointer
sections did in Methods 2, 3, 4 and 5.

In: the path as code; the two refinements as code; and one linear program drawn
on two voxels of example 09, where the limiting substrate changes hands. A new
Section 8 writes down the four equations the guide used but never stated, and
lists all eleven in the order the code applies them.
"""
import os
import re
import shutil
import sys

sys.path.insert(0, "/home/claude/fig3")
from edit_m4 import (P_HEAD, P_GAP, P_MONO, esc, run, eq, where, body, mono,
                     block)  # noqa

ROOT = ("/tmp/claude-0/-home-claude/907c0fba-0d6e-5185-95b4-57fae18fda84"
        "/scratchpad/m1/u")
DOC = os.path.join(ROOT, "word", "document.xml")
RELS = os.path.join(ROOT, "word", "_rels", "document.xml.rels")
MEDIA = os.path.join(ROOT, "word", "media")

P_MONO_END = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
              '<w:spacing w:line="252" w:lineRule="auto"/></w:pPr><w:r><w:rPr>'
              '<w:rFonts w:ascii="Consolas" w:eastAsia="Consolas" '
              'w:hAnsi="Consolas" w:cs="Consolas"/><w:sz w:val="17"/>'
              '<w:szCs w:val="17"/></w:rPr><w:t xml:space="preserve">%s</w:t>'
              '</w:r></w:p>')
P_GAP_END = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
             '<w:spacing w:after="60"/></w:pPr></w:p>')


def monoblock(lines):
    """Monospace lines that stay together but may be broken after.

    mono() interpolates straight into the run, so a line carrying an XML tag
    name -- and this block names two of them -- has to be escaped first.
    """
    ln = [esc(x) for x in lines]
    return mono(ln[:-1]) + P_MONO_END % ln[-1]


# ------------------------------------------------------------- captions -----
CAP_A = ("Figure 4. The flux balance path as code, from a metabolic model on "
         "disk to a change in one voxel. The side panels open the two steps "
         "that are otherwise taken on trust: what the program actually says, "
         "and why the answer is sometimes solved again.")

TXT_A = ("The three bands are three different time scales, and the shape of "
         "them is the reason this method is affordable at all. The yellow band "
         "happens once, before any simulation exists, and none of it is code. "
         "The blue band happens once when the solver starts, and it is where "
         "the expensive part is paid for: the matrix, the fixed bounds and the "
         "objective are the same in every voxel of every step, so they are "
         "turned into one problem object that survives the whole run. The green "
         "band happens in every wet voxel that holds biomass, at every reaction "
         "step, and it is the box drawn in red in Figure 3. What it does to "
         "that problem object is small: it moves one number per substrate, "
         "solves, and reads two things back. The purple box is the only step "
         "where the two back ends differ, and they differ in cost rather than "
         "in answer, which is what makes a disagreement between them a defect "
         "in one of the couplings rather than a property of either solver. The "
         "two boxes drawn in red are the two places where a mistake would "
         "otherwise be silent: a substrate the model has no column for would be "
         "consumed by nothing for the whole run, and two organisms in one voxel "
         "would each be granted the whole supply unless their draws were added "
         "up afterwards and checked.")

CAP_C = ("Figure 5. The two refinements of the solve step, as code. Neither "
         "replaces the method of Figure 4: the first replaces the single "
         "program inside it with a chain, and the second wraps whichever solve "
         "is configured and runs it once per carbon source.")

TXT_C = ("They answer two different complaints, and they can be used together. "
         "The chain answers the one about reproducibility: a single program "
         "pins the growth rate but often leaves the individual exchange fluxes "
         "free to trade against one another, and the solver returns whichever "
         "member of that set its pivoting rule reached. Each stage of the chain "
         "pins one more quantity until nothing is free. The competition answers "
         "the one about behaviour: an organism offered several substrates eats "
         "the best one and switches when it runs out, which no single program "
         "and no rate law will do on its own. The two boxes drawn in red are "
         "the places where each one goes quietly wrong. In the chain, the "
         "growth rate must be read out of the last stage's flux vector rather "
         "than the first stage's optimum, because with any retain fraction "
         "below one those two are different numbers and reporting the first "
         "while handing back the second's fluxes stops mass balancing. In the "
         "competition, shutting an option means forbidding its substrate's "
         "uptake and never its release, because the compound one option "
         "excretes is what the next option will live on.")

CAP_B = ("Figure 8. One linear program, drawn, on two voxels of the same run. "
         "Every number is computed from the four constants in "
         "examples/09_fba_glpk, and the model is small enough that the "
         "program can be solved by hand and the solver's answer confirmed.")

TXT_B = ("The two voxels differ in the donor concentration and in nothing else, "
         "and that one difference moves which constraint is doing the work. In "
         "the fed voxel the acceptor limits, and the reason is stoichiometry "
         "rather than size: a turn of the growth reaction consumes one donor "
         "and half an acceptor, so an acceptor bound of 3.809524 buys twice as "
         "much growth as its own size suggests, and the two ceilings come out "
         "at 9.523810 and 7.619048. In the drawn-down voxel the donor's own "
         "ceiling has fallen to 2.857143 and it is the donor that binds; the "
         "acceptor uptake falls with it rather than staying at its capacity, "
         "because a cell has no use for an acceptor it cannot pair with a "
         "donor. The panel at the foot of the figure sweeps the donor across "
         "the whole range a run visits and shows where the two ceilings cross: "
         "at 0.16 moles per litre exactly, with the acceptor limiting above it "
         "and the donor below. A run whose domain straddles that concentration "
         "is a run in which the limiting substrate changes from one place to "
         "another, and nothing in a concentration field says so.")

# ------------------------------------------------------------ equations -----
EQNS = (
    P_HEAD % "8.1  The four equations the sections above do not write down"
    + body("Equations (1) to (7) describe the program and the competition "
           "around it. Four more carry the answer back to the transport solver, "
           "and they are the ones a reader is most likely to assume rather than "
           "check.")
    + body("The chain of Section 6 is written there as a list of programs. As "
           "an equation it is one inequality per stage, imposed on the stages "
           "that come after it:")
    + eq([("v", "i"), ("j", "s"), (" ≥ a", "i"), ("j", "s"), (" f", "i"),
          ("j", "s"), (" ,      for every stage j already solved", "i")], "8")
    + where([("f", "i"), ("j", "s"),
             (" what stage j reached when it was the objective, and ", "n"),
             ("a", "i"), ("j", "s"),
             (" the fraction of it later stages must preserve. Because every "
              "objective in the chain is a single reaction column, this is a "
              "bound on a column rather than a constraint on a weighted sum, so "
              "no row is added, the matrix never changes, and the warm start "
              "survives. That is the whole reason n stages cost about 1.6 "
              "single solves rather than n of them.", "n")])
    + body("The fluxes come back in millimoles per gram of dry cells per hour, "
           "and the lattice carries moles per litre. One line reconciles them, "
           "over every organism present in the voxel:")
    + eq([("ΔC", "i"), ("k", "s"), (" = Δt Σ", "i"), ("m", "s"), (" v", "i"),
          ("m,k", "s"), (" B", "i"), ("m", "s"), (" M", "i"),
          (" / 3.6 × 10", "i"), ("6", "s")], "9")
    + where([("B", "i"), ("m", "s"),
             (" the biomass of organism m in this voxel and ", "n"),
             ("M", "i"),
             (" its biomass molar mass, so their product is grams of dry cells "
              "per litre of water, which is what turns a per-gram flux into a "
              "per-litre rate. The sum runs over every organism, and the total "
              "is checked against what the voxel holds before it is applied: if "
              "it overruns, the organisms that wanted that substrate are given "
              "one pooled supply and solved again, at most three times, and "
              "then the draw is clamped rather than refused.", "n")])
    + body("If the thermodynamic gate of Method 5 is switched on and this "
           "organism carries a reaction block, the solved fluxes and the growth "
           "are scaled before equation (9) is applied for the last time:")
    + eq([("v", "i"), ("k", "s"), (" ← F", "i"), ("T", "s"), (" v", "i"),
          ("k", "s"), (" ,      μ ← F", "i"), ("T", "s"), (" μ", "i")], "10")
    + where([("the factor multiplies the solved answer and not the bounds that "
              "produced it. A program asked for a smaller bound would "
              "redistribute its whole flux pattern and could return a different "
              "set of by-products, which is a modelling change; scaling the "
              "solution keeps the metabolic answer and asks only what fraction "
              "of it the local energy balance permits. It is also why the "
              "budget check above cannot be broken by the gate: scaling can "
              "only make the draw smaller.", "n")])
    + body("Last, the growth flux becomes a biomass increment, or, when the "
           "program returned none, a loss:")
    + eq([("ΔB = μ B Δt / 3600      if μ ≠ 0 ,      ΔB = − b B Δt      "
           "if μ = 0", "i")], "11")
    + where([("b", "i"), (" the decay coefficient. The two conversions differ "
             "by a factor of a thousand and are easy to confuse: the growth "
             "flux is per hour, because that is the unit a metabolic model "
             "works in, and the decay coefficient is per second, because that "
             "is the unit the transport solver works in. They have separate "
             "converters in the code for exactly that reason. An infeasible "
             "program takes the decay branch too, and is counted so that the "
             "run can say how often it happened.", "n")])
)

ORDER = (
    P_HEAD % "8.2  The order in which the code applies them"
    + body("Read down. The shape of this list is the argument for the method: "
           "everything that costs anything to build is above the last group, "
           "and the last group is the only part that runs more than once.",
           keep=True)
    + P_GAP
    + monoblock([
        "before the run                   write or obtain the metabolic model",
        "                                 name each exchange, and set its capacity",
        "",
        "once, at start-up           1,2  the matrix and the fixed bounds, once",
        "  complab3d_metabolic.hh      4  the objective column, once",
        "                                 the column of every named exchange, by name",
    ])
    + P_GAP_END
    + monoblock([
        "per organism, per voxel,",
        "per reaction step             3  the uptake bound Monod allows here",
        "  complab3d_processors_          held above the model's floor and the supply",
        "  fba.hh                    5,6  the weights, if <cybernetic> is on",
        "                                 one solve per source, the others shut",
        "                              4  maximise growth",
        "                              8  or a chain, each stage pinning one more",
        "                              7  blend the options, if <cybernetic> is on",
        "                              9  what every organism here would draw",
        "                                 if the total overruns, share and re-solve",
        "                             10  the thermodynamic factor, if the gate is on",
        "                              9  the draws again, from the gated fluxes",
        "                             11  growth into biomass, or decay",
    ])
    + P_GAP_END
    + body("One matrix, built once. One number per substrate, moved per voxel. "
           "Two numbers read back. Everything else in this document is about "
           "what happens when that is not quite enough: when the answer is not "
           "unique, when the organism has a choice of substrate, when two "
           "organisms share a voxel, or when the energy does not add up.")
)


def main():
    doc = open(DOC, encoding="utf-8").read()
    rels = open(RELS, encoding="utf-8").read()

    for src, name, rid in (("m1a.png", "imageA.png", "rId90"),
                           ("m1c.png", "imageC.png", "rId92"),
                           ("m1b.png", "imageB.png", "rId91")):
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

    # anchors, taken before anything is deleted or renumbered
    figA_at = paras[find("Figure 6. Pipeline a") - 2]
    figC_at = paras[find("Figure 9. Pipeline c") - 2]

    # --- out: the nine plates the three new flowcharts now carry ------------
    for cap in ("Figure 2. How this method works",
                "Figure 4. One voxel, end to end",
                "Figure 6. Pipeline a",
                "Figure 7. The same pipeline with numbers",
                "Figure 8. Pipeline b",
                "Figure 9. Pipeline c",
                "Figure 11. Multi-step flux balance analysis",
                "Figure 12. Pipeline d",
                "Figure 14. Cybernetic substrate switching"):
        i = find(cap)
        doc = doc.replace(paras[i - 1], "", 1)
        doc = doc.replace(paras[i], "", 1)

    # the paragraph that walks the deleted numbers plate goes with it
    i = find("Figure 7 is small enough to check by hand")
    doc = doc.replace(paras[i], "", 1)

    # --- out: the pointer section -------------------------------------------
    # Cut this one by offset rather than paragraph by paragraph. Its content is
    # a table, so a paragraph-list deletion would leave the empty table shell
    # behind, and worse: its cells hold paragraphs whose XML is byte-identical
    # to cells of the troubleshooting table above, so replace(p, "", 1) removes
    # the earlier one and takes that table apart instead.
    def para_start(needle):
        i = doc.index(needle)
        j = doc.rfind("<w:p ", 0, i)
        k = doc.rfind("<w:p>", 0, i)
        return max(j, k)

    doc = doc[:para_start("10. Where to go next")] + doc[para_start("11. References"):]

    # --- renumber: a new Section 8 opens up, and everything after moves one --
    for old, new in (("9. What to check when it has run", "@@10@@ What to check when it has run"),
                     ("8. A worked example", "@@9@@ A worked example")):
        doc = doc.replace(old, new, 1)
    for a, b in (("8.1", "9.1"), ("8.2", "9.2"), ("8.3", "9.3"), ("8.4", "9.4"),
                 ("8.5", "9.5"), ("8.6", "9.6"), ("8.7", "9.7")):
        doc = doc.replace(">" + a + " ", ">" + b + " ")
    doc = doc.replace("Section 8.5 is where", "Section 9.5 is where")
    doc = doc.replace("the two bounds from Section 8.3", "the two bounds from Section 9.3")

    # --- renumber the figures, before anything new is inserted --------------
    doc = doc.replace("Figure 3. The program, block by block", "@@F2@@ The program, block by block")
    doc = doc.replace("Figure 5. The solve loop", "@@F3@@ The solve loop")
    doc = doc.replace("Figure 10. The chain on the network", "@@F6@@ The chain on the network")
    doc = doc.replace("Figure 13. The competition on a network", "@@F7@@ The competition on a network")
    doc = doc.replace("the small network of Figure 10", "the small network of Figure 6")
    doc = doc.replace("Figure 13 is the whole mechanism", "Figure 7 is the whole mechanism")
    doc = re.sub(r"@@F(\d)@@ ", lambda m: "Figure %s. " % m.group(1), doc)
    doc = re.sub(r"@@(\d+)@@ ", lambda m: "%s. " % m.group(1), doc)

    # --- in: the equations, then the three flowcharts ------------------------
    h = next(p for p in re.findall(r'<w:p [^>]*>.*?</w:p>', doc, re.S)
             if txt(p).startswith("9. A worked example"))
    doc = doc.replace(h, P_HEAD % "8.  Every equation, in the order the code applies them"
                      + EQNS + ORDER + h, 1)

    doc = doc.replace(figA_at, figA_at + block("rId90", 641, 2760, CAP_A, TXT_A,
                                               cx=4450000), 1)
    doc = doc.replace(figC_at, figC_at + block("rId92", 642, 2380, CAP_C, TXT_C,
                                               cx=4600000), 1)
    step1 = next(p for p in re.findall(r'<w:p [^>]*>.*?</w:p>', doc, re.S)
                 if txt(p).startswith("9.3 Step 1: turn the concentrations"))
    doc = doc.replace(step1, block("rId91", 643, 3140, CAP_B, TXT_B,
                                   cx=4100000) + step1, 1)

    # Every run of monospace lines in this document is one block that has to be
    # read together, and the ones that came with the guide carry no keepNext, so
    # a matrix or a flux vector splits across a page break. Give each line
    # keepNext except the last of its run, which lets the block move as a unit
    # while still allowing a break after it.
    paras2 = re.findall(r'<w:p [^>]*>.*?</w:p>|<w:p/>', doc, re.S)
    ismono = [("Consolas" in p) for p in paras2]
    fixed = 0
    for i, p in enumerate(paras2):
        if not ismono[i] or (i + 1 < len(ismono) and not ismono[i + 1]):
            continue
        if "<w:keepNext/>" in p:
            continue
        if "<w:pPr>" in p:
            q = p.replace("<w:pPr>", "<w:pPr><w:keepNext/><w:keepLines/>", 1)
        else:
            q = p.replace(">", "><w:pPr><w:keepNext/><w:keepLines/></w:pPr>", 1)
        doc = doc.replace(p, q, 1)
        fixed += 1
    print("held %d monospace lines to the line below them" % fixed)

    # A table that butts straight up against the paragraph before it renders
    # with that paragraph pulled into its first cell, which is how the GLPK and
    # COBRApy comparison has been coming out. One empty paragraph in front of
    # each table fixes it, and costs four points of space.
    spacer = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
              '<w:spacing w:after="60"/></w:pPr></w:p>')
    doc, n = re.subn(r'</w:p><w:tbl>', '</w:p>' + spacer + '<w:tbl>', doc)
    print("spaced %d tables off the paragraph above them" % n)

    open(DOC, "w", encoding="utf-8").write(doc)
    open(RELS, "w", encoding="utf-8").write(rels)
    print("figures:", sorted(set(re.findall(r'Figure (\d+)\. ', doc)), key=int))
    print("sections:", sorted(set(re.findall(r'>(\d{1,2})\.  ?[A-Z]', doc)), key=int))
    print("8.x:", sorted(set(re.findall(r'>(8\.\d)  ?[A-Z]', doc))))
    print("9.x:", sorted(set(re.findall(r'>(9\.\d) [A-Z]', doc))))
    print("eq labels:", re.findall(r'\s{4,}\((\d+)\)', re.sub(r'<[^>]+>', '', doc)))


if __name__ == "__main__":
    main()
