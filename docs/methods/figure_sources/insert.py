#!/usr/bin/env python3
"""Put the two flowcharts into Method_3_Symbolic_guide.docx, in place.

The guide is hand-edited, so it is never rebuilt from a script: it is unzipped,
its document.xml is changed, and it is zipped back. This adds two pictures with
their captions and one paragraph of explanation each, and renumbers the four
figures that now come after them. Nothing else in the document is touched.
"""
import os
import re
import shutil

ROOT = "/tmp/claude-0/-home-claude/907c0fba-0d6e-5185-95b4-57fae18fda84/scratchpad/m3b/u"
DOC = os.path.join(ROOT, "word", "document.xml")
RELS = os.path.join(ROOT, "word", "_rels", "document.xml.rels")
MEDIA = os.path.join(ROOT, "word", "media")

EMU_W = 5334000                              # what every other figure uses

PIC = (
    '<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
    '<w:spacing w:before="140" w:after="40"/><w:jc w:val="center"/></w:pPr>'
    '<w:r><w:rPr><w:noProof/></w:rPr><w:drawing>'
    '<wp:inline distT="0" distB="0" distL="0" distR="0">'
    '<wp:extent cx="%(cx)d" cy="%(cy)d"/>'
    '<wp:effectExtent l="0" t="0" r="0" b="0"/>'
    '<wp:docPr id="%(id)d" name="Picture %(id)d"/><wp:cNvGraphicFramePr>'
    '<a:graphicFrameLocks xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" noChangeAspect="1"/>'
    '</wp:cNvGraphicFramePr>'
    '<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
    '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
    '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
    '<pic:nvPicPr><pic:cNvPr id="0" name=""/><pic:cNvPicPr/></pic:nvPicPr>'
    '<pic:blipFill><a:blip r:embed="%(rid)s"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
    '<pic:spPr bwMode="auto"><a:xfrm><a:off x="0" y="0"/>'
    '<a:ext cx="%(cx)d" cy="%(cy)d"/></a:xfrm>'
    '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
    '</pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>'
)

CAP = (
    '<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
    '<w:spacing w:after="120"/><w:jc w:val="center"/></w:pPr>'
    '<w:r><w:rPr><w:i/><w:iCs/><w:sz w:val="18"/><w:szCs w:val="18"/></w:rPr>'
    '<w:t xml:space="preserve">%s</w:t></w:r></w:p>'
)

BODY = (
    '<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
    '<w:spacing w:after="140" w:line="288" w:lineRule="auto"/></w:pPr>'
    '<w:r><w:t xml:space="preserve">%s</w:t></w:r></w:p>'
)


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


CAP4 = ("Figure 4. The symbolic path as code, from a table of numbers to a "
        "change in one voxel. The side panels open the two steps that are "
        "otherwise taken on trust: one generation of the search, and one walk "
        "of one expression tree in one voxel.")

TXT4 = ("The three bands are the three time scales this method lives on. The "
        "yellow band happens once, on a workstation, and none of it is inside "
        "CompLaB3D. The blue band happens once when the solver starts. The "
        "green band happens in every wet voxel at every reaction step, and is "
        "the box drawn in red in Figure 1. Each band carries the name of the "
        "source file it lives in, so a box in the chart can be traced to the "
        "code that performs it. The dashed line across the page is the "
        "handover, and the .sym file is the only thing that crosses it. Two "
        "boxes are drawn in red because they are the two places where a "
        "mistake does not announce itself: a variable name the run does not "
        "recognise stops the run at start-up rather than becoming a silent "
        "zero, and a concentration arriving from outside the range the fit was "
        "shown is pulled back to the edge and counted.")

CAP5 = ("Figure 5. Inside the search: what each step does to a formula. Step 5 "
        "is drawn out in full because it is the step that carries the method, "
        "and the one whose description elsewhere has to be taken on trust.")

TXT5 = ("The flow down the left is drawn without colour so that the only "
        "coloured objects in the figure are the formulas, which are what the "
        "search changes. A formula is held as a picture of boxes: a yellow box "
        "performs one arithmetic step, a blue box holds a variable name the "
        "run fills in per voxel, and a green box holds a number the fitting is "
        "free to move. The count of those boxes is the length that Figure 6 "
        "indexes its list by. Steps 2 to 5 of the chart are one round, and the "
        "three lanes down the left edge are the three loops: one adds the next "
        "formula to the population, one closes the round, one begins the next. "
        "Step 5 states what enters it and what leaves it, because that is "
        "where the description usually goes vague: in, the population as step "
        "2 left it, every member already fitted and scored; out, exactly one "
        "new formula, whichever of the three moves was drawn. The shelf keeps "
        "one slot per number of boxes, which is what stops a short formula "
        "from being discarded for being less accurate than a long one.")


def main():
    doc = open(DOC, encoding="utf-8").read()
    rels = open(RELS, encoding="utf-8").read()

    # --- the two new images, and relationships for them ----------------------
    plan = [("figA.png", "image10.png", "rId80", 2320, 501),
            ("figB.png", "image11.png", "rId81", 2880, 502)]
    add = []
    for src, name, rid, px_h, _id in plan:
        shutil.copy("/home/claude/fig3/" + src, os.path.join(MEDIA, name))
        add.append('<Relationship Id="%s" Type="http://schemas.openxmlformats.org'
                   '/officeDocument/2006/relationships/image" Target="media/%s"/>'
                   % (rid, name))
    rels = rels.replace("</Relationships>", "".join(add) + "</Relationships>")

    def block(rid, pid, px_h, cap, txt, cx=EMU_W):
        cy = int(round(cx * px_h / 1900.0))
        return (PIC % dict(cx=cx, cy=cy, id=pid, rid=rid)
                + CAP % esc(cap) + BODY % esc(txt))

    # --- the four figures that moved down ------------------------------------
    for old in (7, 6, 5, 4):
        doc = doc.replace("<w:t>Figure %d. " % old,
                          "<w:t>Figure %d. " % (old + 2))
        doc = doc.replace('<w:t xml:space="preserve">Figure %d. ' % old,
                          '<w:t xml:space="preserve">Figure %d. ' % (old + 2))

    paras = re.findall(r'<w:p [^>]*>.*?</w:p>|<w:p/>', doc, re.S)

    def find(prefix):
        for i, p in enumerate(paras):
            t = re.sub(r'<[^>]+>', '', p)
            if t.strip().startswith(prefix):
                return i
        raise SystemExit("not found: " + prefix)

    # Figure 4 closes section 3; Figure 5 opens section 4, before its list.
    after_sec3 = paras[find("That last point is the same reasoning")]
    before_list = paras[find("A symbolic regression run does not return one")]

    doc = doc.replace(after_sec3,
                      after_sec3 + block("rId80", 501, 2320, CAP4, TXT4), 1)
    doc = doc.replace(before_list,
                      before_list + block("rId81", 502, 2880, CAP5, TXT5,
                                          cx=4765833), 1)

    open(DOC, "w", encoding="utf-8").write(doc)
    open(RELS, "w", encoding="utf-8").write(rels)
    print("inserted; captions now:",
          sorted(set(re.findall(r'Figure (\d)\. ', doc))))


if __name__ == "__main__":
    main()
