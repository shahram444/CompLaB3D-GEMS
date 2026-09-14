#!/usr/bin/env python3
"""Second pass on Method_3_Symbolic_guide.docx.

Three figures go: the end-to-end band chart, the complexity-versus-error plot
and the expression-tree plate. Each of them is now said better by one of the two
flowcharts, and a guide that shows the same thing twice teaches the reader to
skim. Their captions go with them and the figures that follow are renumbered.

Section 6 grows. It used to hold only the five equations the solver evaluates,
which left every number that arrives in a .sym file unexplained: what was
minimised to produce it, and what the solver works out for itself before the
first step. Those are added as (S1) to (S6) and (D1), and the section ends with
the order in which all of them actually run.
"""
import os
import re

ROOT = "/tmp/claude-0/-home-claude/907c0fba-0d6e-5185-95b4-57fae18fda84/scratchpad/m3c/u"
DOC = os.path.join(ROOT, "word", "document.xml")

P_BODY = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
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
          '<w:keepNext/><w:keepLines/>'
          '<w:spacing w:line="252" w:lineRule="auto"/></w:pPr>'
          '<w:r><w:rPr><w:rFonts w:ascii="Consolas" w:eastAsia="Consolas" '
          'w:hAnsi="Consolas" w:cs="Consolas"/><w:sz w:val="17"/>'
          '<w:szCs w:val="17"/></w:rPr><w:t xml:space="preserve">%s</w:t></w:r></w:p>')
P_GAP = ('<w:p w:rsidR="00000000" w:rsidRDefault="00000000"><w:pPr>'
         '<w:spacing w:after="60"/></w:pPr></w:p>')


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def run(text, kind="n", sz=None):
    """kind: n plain, i italic, s italic subscript, b bold."""
    rpr = []
    if kind in ("i", "s"):
        rpr.append("<w:i/><w:iCs/>")
    if kind == "b":
        rpr.append("<w:b/><w:bCs/>")
    if sz:
        rpr.append("<w:sz w:val=\"%d\"/><w:szCs w:val=\"%d\"/>" % (sz, sz))
    if kind == "s":
        rpr.append('<w:vertAlign w:val="subscript"/>')
    pr = "<w:rPr>%s</w:rPr>" % "".join(rpr) if rpr else ""
    return '<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r>' % (pr, esc(text))


def eq(parts, label):
    """parts is a list of (text, kind); label is the equation number."""
    body = "".join(run(t, k, 24) for t, k in parts)
    return P_EQ % (body + run("        (%s)" % label, "n", 24))


def where(parts):
    return P_WHERE % (run("where ", "b") + "".join(run(t, k) for t, k in parts))


def body(text):
    return P_BODY % run(text)


def mono(lines):
    return "".join(P_MONO % ln for ln in lines)


# ---------------------------------------------------------------- section 6 --
S61 = (
    P_HEAD % "6.1  The offline equations: what the search minimises"
    + body("The four subsections after this one are what the solver evaluates. "
           "This one is what produced the file it evaluates, and it belongs here "
           "because a number in a .sym file cannot be read honestly without "
           "knowing what was minimised to arrive at it. None of it runs during a "
           "simulation.")
    + body("Every sample is weighted before anything is squared. A rate law is "
           "used across orders of magnitude, and plain least squares on the "
           "absolute error only cares about the large values, so the fit comes "
           "out excellent at the top of the range and can be wrong by several "
           "hundred per cent at the bottom, which is exactly where a substrate is "
           "running out and the interesting chemistry is:")
    + eq([("r", "i"), ("i", "s"), (" = ( E( x", "i"), ("i", "s"),
          (" , c ) − y", "i"), ("i", "s"), (" ) · w", "i"), ("i", "s")], "S1")
    + where([("E", "i"), (" the expression being tried; ", "n"), ("c", "i"),
             (" its constants; ", "n"), ("x", "i"), ("i", "s"),
             (" the inputs of sample i; ", "n"), ("y", "i"), ("i", "s"),
             (" the rate that sample is to reproduce; ", "n"), ("w", "i"),
             ("i", "s"),
             (" the weight, which is 1 under absolute loss and 1 / max( |y", "n"),
             ("i", "s"),
             ("| , f ) under relative loss, the default. The floor f is one per "
              "cent of the root-mean-square rate, and without it a sample whose "
              "true rate is essentially zero would carry unbounded weight and "
              "drag the whole fit onto a corner of the domain that carries no "
              "information.", "n")])
    + body("The numbers inside one shape are then fitted by Levenberg−Marquardt, "
           "on the residuals of (S1) with a small ridge attached:")
    + eq([("S( c ) = Σ", "i"), ("i", "s"), (" r", "i"), ("i", "s"),
          ("² + ρ² Σ", "i"), ("j", "s"), (" ( c", "i"), ("j", "s"),
          (" / c̄ )²", "i")], "S2")
    + where([("ρ", "i"),
             (" ten to the minus seven times the span of the data, and ", "n"),
             ("c̄", "i"),
             (" the larger of that span and one. The second sum is a ridge and it "
              "exists for one purpose: when a shape carries more constants than "
              "the data can pin down, the best fit is not a point but a valley, "
              "and without the ridge two runs land in different places in it and "
              "the same law comes back spelled two different ways. It is dropped "
              "again before the fit is scored, so it is never counted as error.", "n")])
    + body("What is scored is the first sum alone, as a mean over the samples:")
    + eq([("M = ( 1 / N ) Σ", "i"), ("i", "s"), (" r", "i"), ("i", "s"),
          ("²", "i")], "S3")
    + body("and what the search ranks on adds a price for every box the shape uses:")
    + eq([("F = M / M", "i"), ("0", "s"), (" + p · n ,        M", "i"),
          ("0", "s"), (" = ( 1 / N ) Σ", "i"), ("i", "s"), (" ( y", "i"),
          ("i", "s"), (" w", "i"), ("i", "s"), (" )² ,        p = 4×10", "i"),
          ("−3", "i")], "S4")
    + where([("n", "i"), (" the number of boxes in the shape and ", "n"),
             ("M", "i"), ("0", "s"),
             (" the error of predicting zero everywhere, which is what makes F a "
              "pure number. Without the second term the search buys accuracy with "
              "boxes until nothing on the list can be read, and the whole point "
              "of this path is that the answer is something a person can argue "
              "with.", "n")])
    + body("The two percentages printed beside each row of the list are the "
           "typical and the worst error over the samples, judged the same way the "
           "fit was:")
    + eq([("typical = 100 √( mean d² ) ,        worst = 100 max d ,"
           "        d", "i"), ("i", "s"), (" = | E( x", "i"), ("i", "s"),
          (" , c ) − y", "i"), ("i", "s"), (" | · w", "i"), ("i", "s")], "S5")
    + body("and the payoff column, which is the one to read down, is how much the "
           "error fell for each extra box:")
    + eq([("g = ( ln M", "i"), ("prev", "s"), (" − ln M ) / ( n − n", "i"),
          ("prev", "s"), (" )", "i")], "S6")
    + body("The automatic pick is the largest g. It is a rule of thumb and not an "
           "answer, which is why the whole list is printed and why --pick exists: "
           "on the shipped training table the rule lands on the crude bilinear "
           "formula and skips the one that recovers both half-saturation "
           "constants. Choosing a row is a scientific act, not a numerical one.")
)

S62 = (
    P_HEAD % "6.2  What the solver works out for itself at start-up"
    + body("If the file states a reaction, a yield and a biomass name, the "
           "substrate lines are not read from the file at all. They are derived "
           "once, before the first step, and printed in the log marked as derived:")
    + eq([("coefficient( i ) = ( ν", "i"), ("i", "s"),
          (" / | ν", "i"), ("key", "s"), (" | ) / Y", "i")], "D1")
    + where([("ν", "i"), ("i", "s"),
             (" the stoichiometric coefficient of species i in the reaction line, "
              "negative for consumed and positive for produced; ", "n"),
             ("ν", "i"), ("key", "s"),
             (" the coefficient of the species the yield is quoted against; ", "n"),
             ("Y", "i"),
             (" the yield, in moles of biomass per mole of that species.", "n")])
    + body("For the shipped example, acetate at −1 and oxygen at −2 with a "
           "yield of 0.4 give −2.5 and −5.0 exactly. The ratio between the two "
           "substrate lines is then arithmetic on the reaction rather than two "
           "numbers typed twice, and it cannot disagree with the chemistry. "
           "Mistype one digit in a hand-written pair and the run still finishes, "
           "the fields still look smooth, and the model quietly creates or "
           "destroys oxygen every step for the rest of the simulation. A species "
           "may have a rate line or a place in the reaction, never both; two "
           "sources of truth for one number is the thing this removes.")
)

ORDER = (
    P_HEAD % "6.7  The order in which the code applies them"
    + P_BODY.replace('<w:pPr>', '<w:pPr><w:keepNext/><w:keepLines/>')
      % run("Read down. The three groups are three different time scales, and "
            "the gap between them is the cost argument for this whole path.")
    + P_GAP.replace('<w:pPr>', '<w:pPr><w:keepNext/>')
    + mono([
        "once, offline                 S1   weight every sample, so the small rates count",
        "  tools/method_3_symbolic/fit_symbolic.py       S2   fit the numbers inside one shape",
        "                              S3   score that shape against the data",
        "                              S4   rank it, with a price for every box",
        "                              S5   report typical and worst error for the list",
        "                              S6   the payoff column, and the automatic pick",
        "",
        "once, at start-up             D1   derive one substrate coefficient per species",
        "  complab3d_symbolic.hh            parse each rate line into a tree",
        "                                   bind every variable name to a lattice",
        "",
        "per voxel, per reaction step   1   clamp each input to its range, and count it",
        "  complab3d_integration.hh     2   walk each tree, in file order",
        "                               3   put the result in the solver's time unit",
        "                               4   apply the substrate increments, floored",
        "                               5   apply the biomass increment",
    ])
    + P_GAP
    + body("Nothing in the first group runs during a simulation, and nothing in "
           "the second runs more than once. Only the third is inside the loop, "
           "and it is five lines of arithmetic over a tree of a few dozen boxes. "
           "That is the trade this path makes: the expensive equations are the "
           "ones that never run again.")
)


def main():
    doc = open(DOC, encoding="utf-8").read()
    paras = re.findall(r'<w:p [^>]*>.*?</w:p>|<w:p/>', doc, re.S)

    def text_of(p):
        return re.sub(r'<[^>]+>', '', p).strip()

    def find(prefix, start=0):
        for i in range(start, len(paras)):
            if text_of(paras[i]).startswith(prefix):
                return i
        raise SystemExit("not found: " + prefix)

    # --- the three figures that are now said better elsewhere ----------------
    drop = []
    for cap in ("Figure 2. How this method works",
                "Figure 6. The list the search returns",
                "Figure 7. One line of a .sym file"):
        i = find(cap)
        drop += [paras[i - 1], paras[i]]          # the picture, then its caption
    for d in drop:
        doc = doc.replace(d, "", 1)

    # --- one sentence pointed at a figure that has gone ----------------------
    doc = doc.replace(
        "The count of those boxes is the length that Figure 6 indexes its list by.",
        "The count of those boxes is what the list calls length.")

    # --- renumber what is left: 3,4,5,8,9 become 2,3,4,5,6 -------------------
    for old, new in ((3, 2), (4, 3), (5, 4), (8, 5), (9, 6)):
        doc = doc.replace("Figure %d. " % old, "@@%d@@ " % new)
    doc = re.sub(r"@@(\d)@@ ", lambda m: "Figure %s. " % m.group(1), doc)

    # --- section 6 gains three subsections; the four it had shift down -------
    for old, new in ((4, 6), (3, 5), (2, 4), (1, 3)):
        doc = doc.replace("6.%d  " % old, "@@6.%d@@  " % new)
    doc = re.sub(r"@@6\.(\d)@@  ", lambda m: "6.%s  " % m.group(1), doc)
    doc = doc.replace("The unit scale of Section 6.3 is",
                      "The unit scale of Section 6.5 is")

    head63 = paras_after(doc, "6.3  What is handed in, and the clamp")
    doc = doc.replace(head63, S61 + S62 + head63, 1)

    # section 5 lost its opening plate, so it needs a sentence of its own
    h5 = paras_after(doc, "5.  The expression, as the solver holds it")
    doc = doc.replace(h5, h5 + body(
        "The picture of a formula that Figure 4 introduces is exactly what the "
        "solver holds: one node per box, built once at start-up and walked per "
        "voxel thereafter. This section is the language those boxes are written "
        "in, and the single rule that governs the order they are read in."), 1)

    tail = paras_after(doc, "7.  A worked example")
    doc = doc.replace(tail, ORDER + tail, 1)

    open(DOC, "w", encoding="utf-8").write(doc)
    print("figures now:", sorted(set(re.findall(r'Figure (\d)\. ', doc))))
    print("6.x now:", sorted(set(re.findall(r'>(6\.\d)  ', doc))))


def paras_after(doc, prefix):
    for p in re.findall(r'<w:p [^>]*>.*?</w:p>|<w:p/>', doc, re.S):
        if re.sub(r'<[^>]+>', '', p).strip().startswith(prefix):
            return p
    raise SystemExit("not found: " + prefix)


if __name__ == "__main__":
    main()
