#!/usr/bin/env python3
"""Method 5 guide: the same pass Methods 3 and 4 had."""
import os
import re
import shutil
import sys

sys.path.insert(0, "/home/claude/fig3")
from edit_m4 import (P_HEAD, P_GAP, esc, run, eq, where, body, mono, block)  # noqa

ROOT = "/tmp/claude-0/-home-claude/907c0fba-0d6e-5185-95b4-57fae18fda84/scratchpad/m5/u"
DOC = os.path.join(ROOT, "word", "document.xml")
RELS = os.path.join(ROOT, "word", "_rels", "document.xml.rels")
MEDIA = os.path.join(ROOT, "word", "media")

CAP_A = ("Figure 2. The thermodynamic gate as code, from a file of energies to a "
         "rate that is allowed to run. The gate makes no rate of its own, which "
         "is why the green band sits inside the box marked in red in Figure 1 "
         "rather than beside it.")

TXT_A = ("The three bands are three different time scales, and the first of them "
         "is optional. Nothing has to be fitted for this path: the .thm file is "
         "written by hand from the literature, at the temperature of the case, "
         "and the sweep in the yellow band only reads it back. That sweep is "
         "worth running anyway, because a gate that is open everywhere and a "
         "gate that is shut everywhere are both invisible in the output "
         "afterwards. The blue band happens once when the solver starts, and it "
         "is where a species name that is not in the run, a second block on an "
         "organism that already has one, or a floor of zero under the activities "
         "stops the run rather than becoming a silent wrong answer. The green "
         "band happens per gated organism per voxel per reaction step, and it "
         "sits inside whichever rate path that organism already uses: the gate "
         "does not care whether the number it is about to multiply came from "
         "compiled kinetics, a linear program, a surrogate, a symbolic law or a "
         "graph network.")

CAP_B = ("Figure 4. The four equations on two voxels of one aggregate, with the "
         "factor drawn. Every number is computed from the constants in the "
         "shipped file at the two compositions the offline sweep runs between.")

TXT_B = ("The two voxels differ in one thing only: deep inside the aggregate the "
         "sulfide the reaction has already made has built up, so the same "
         "reaction is close to its own equilibrium. That moves the reaction "
         "quotient by about eleven units of its logarithm, which moves the free "
         "energy by twenty-five kilojoules per mole, which is enough to carry it "
         "across the threshold. The curve in the fourth panel is the whole of the "
         "method: left of the threshold the reaction can pay for its ATP and the "
         "factor climbs to one, right of it nothing runs, and the climb between "
         "them is not a step. Its width is chi R T ln 2, here 1.6 kilojoules per "
         "mole, which is why a reaction front computed this way has a thickness "
         "rather than an edge, and why that thickness is set by the chemistry and "
         "the temperature rather than by a parameter anyone chose.")

YIELD = (
    P_HEAD % "5.6  The optional growth yield, and what it does not do"
    + body("A block may also carry a yield line. When it does, the same free "
           "energy of equation (2) is turned into a biomass yield by the "
           "correlation of Heijnen and van Dijken, which prices the dissipation "
           "a cell cannot avoid:")
    + eq([("dG", "i"), ("dis", "s"), (" = 200 + 18 | 6 − NoC |", "i"),
          ("1.8", "i"), (" + exp [ ( ( 3.8 − γ )", "i"), ("2", "s"),
          (" )", "i"), ("0.16", "i"), (" ( 3.6 + 0.4 NoC ) ]", "i")], "7")
    + eq([("Y = dG / − ( dG", "i"), ("ana", "s"), (" + dG", "i"), ("dis", "s"),
          (" )", "i")], "8")
    + where([("NoC", "i"),
             (" the number of carbon atoms in the carbon source and ", "n"),
             ("γ", "i"), (" its degree of reduction; ", "n"), ("dG", "i"),
             ("ana", "s"),
             (" the energy of building biomass from it, which is negative; ", "n"),
             ("dG", "i"), ("dis", "s"),
             (" a positive dissipation, so the denominator is negative and Y "
              "comes out positive, in C-moles of biomass per mole of reaction. "
              "A dGdis line in the file overrides the correlation entirely. The "
              "absolute value in the first term matters: without it any carbon "
              "source with more than six carbons raises a negative number to a "
              "fractional power, and the yield came back as not-a-number.", "n")])
    + body("One thing about this has to be said plainly, because it is the kind of "
           "thing a reader would otherwise assume. The function that evaluates "
           "equation (8) is not called by any rate path. With a yield line "
           "present the number appears in the start-up log and in the offline "
           "sweep's table, and it is exercised by the tests, but the biomass a "
           "run produces still comes from whatever yield that organism's own rate "
           "path uses. In the shipped example that is the constant in the "
           "compiled kinetics file, not equation (8). Read it as a reporting "
           "quantity until that changes.")
)

ORDER = (
    P_HEAD % "5.7  The order in which the code applies them"
    + body("Read down. The gate is cheap, and the reason is in the shape of this "
           "list: nothing above the last group runs more than once.", keep=True)
    + P_GAP
    + mono([
        "before the run                   write the .thm file, from the literature",
        "  offline/thermo_curve.py        optional: sweep the compositions the run visits",
        "",
        "once, at start-up                parse the file; scale the energies to kJ",
        "  complab3d_thermo.hh            bind every species of the reaction by name",
        "                                 map each organism to its block, or to none",
        "",
        "per gated organism, per voxel,",
        "per reaction step             1  the reaction quotient, as a sum of logarithms",
        "  inside each rate path's     2  the free energy available here",
        "  own processor               3  what is left once the ATP is paid for",
        "                              4  the factor",
        "                            5,6  multiply every substrate rate, and the growth",
    ])
    + P_GAP
    + body("Four lines of arithmetic and one exponential, over as many species as "
           "the reaction has. That is the whole run-time cost, and it is why the "
           "gate can sit in front of a linear program without being noticed in "
           "the timings. Equations (7) and (8) are not in this list because "
           "nothing in the loop calls them.")
)


def main():
    doc = open(DOC, encoding="utf-8").read()
    rels = open(RELS, encoding="utf-8").read()

    for src, name, rid in (("m5a.png", "imageA.png", "rId90"),
                           ("m5b.png", "imageB.png", "rId91")):
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

    fig2_at = paras[find("Figure 2. How the gate works") - 2]
    fig5_at = paras[find("Figure 5. The four equations") - 2]

    for cap in ("Figure 2. How the gate works",
                "Figure 3. Every biotic rate path returns",
                "Figure 5. The four equations",
                "Figure 7. The file shipped with",
                "Figure 8. The energy-based growth yield"):
        i = find(cap)
        doc = doc.replace(paras[i - 1], "", 1)
        doc = doc.replace(paras[i], "", 1)

    s7, s9 = find("7.  The energetics file"), find("9.  References")
    for p in paras[s7:s9]:
        doc = doc.replace(p, "", 1)
    doc = doc.replace("9.  References", "7.  References")

    doc = doc.replace("and Section 8.3 explains why the difference is the point",
                      "and evaluating the gate at the bulk composition instead is "
                      "a different question with a different answer")
    doc = doc.replace("See Section 7.1.", "See Section 6.5.")

    # figure numbers: 4 and 6 move up to 3 and 5
    doc = doc.replace("Figure 4. Left, a slice", "@@3@@ Left, a slice")
    doc = doc.replace("Figure 6. The same aggregate", "@@5@@ The same aggregate")
    doc = doc.replace("The middle panel of Figure 4 is", "The middle panel of Figure 3 is")
    doc = re.sub(r"@@(\d)@@ ", lambda m: "Figure %s. " % m.group(1), doc)

    h = next(p for p in re.findall(r'<w:p [^>]*>.*?</w:p>', doc, re.S)
             if txt(p).startswith("6.  A worked example"))
    doc = doc.replace(h, YIELD + ORDER + h, 1)

    doc = doc.replace(fig2_at, fig2_at + block("rId90", 611, 2480, CAP_A, TXT_A,
                                               cx=4450000), 1)
    doc = doc.replace(fig5_at, fig5_at + block("rId91", 612, 2360, CAP_B, TXT_B,
                                               cx=4700000), 1)

    open(DOC, "w", encoding="utf-8").write(doc)
    open(RELS, "w", encoding="utf-8").write(rels)
    print("figures:", sorted(set(re.findall(r'Figure (\d)\. ', doc))))
    print("sections:", sorted(set(re.findall(r'>(\d{1,2})\.  [A-Z]', doc)), key=int))
    print("5.x:", sorted(set(re.findall(r'>(5\.\d)  ', doc))))


if __name__ == "__main__":
    main()
