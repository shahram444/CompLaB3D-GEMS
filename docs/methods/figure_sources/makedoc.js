// Build the new Method 3 flowchart document.
//
// Three sections. The first keeps the guide's page setup for running text; the
// two that follow pull the margins in to 0.6 inch so each flowchart gets as
// much of the sheet as it can. Both are tall figures and height, not width, is
// what limits them, so every tenth of an inch of margin is a tenth of an inch
// of legibility.
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, ImageRun,
  AlignmentType, PageOrientation,
} = require("docx");

const BODY = "Cambria";
const TW = 12240, TH = 15840;                 // US Letter, in twentieths of a point
const TEXT_MARGIN = { top: 1440, bottom: 1440, left: 1620, right: 1620 };
const PLATE_MARGIN = { top: 864, bottom: 864, left: 864, right: 864 };

const page = (margin) => ({
  page: {
    size: { width: TW, height: TH, orientation: PageOrientation.PORTRAIT },
    margin,
  },
});

const p = (text, o = {}) =>
  new Paragraph({
    spacing: { after: o.after === undefined ? 160 : o.after, line: 288 },
    children: [new TextRun({ text, font: BODY, size: o.size || 22,
                             bold: !!o.bold, italics: !!o.it })],
  });

const title = (text) =>
  new Paragraph({
    spacing: { after: 220 },
    children: [new TextRun({ text, font: BODY, size: 32, bold: true })],
  });

const heading = (text) =>
  new Paragraph({
    spacing: { before: 280, after: 120 },
    children: [new TextRun({ text, font: BODY, size: 24, bold: true })],
  });

const picture = (file, w, h) =>
  new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 0, after: 0 },
    children: [new ImageRun({ type: "png", data: fs.readFileSync(file),
                              transformation: { width: w, height: h } })],
  });

const caption = (text) =>
  new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 140, after: 0, line: 264 },
    children: [new TextRun({ text, font: BODY, size: 18, italics: true })],
  });

const doc = new Document({
  sections: [
    // ------------------------------------------------------------ the text --
    {
      properties: page(TEXT_MARGIN),
      children: [
        title("Rate Path 3, Symbolic Rate Laws: the two flowcharts"),

        p("This document holds the two flowcharts that go with Method 3. The " +
          "first answers where the symbolic path sits in the code and what " +
          "happens at each point along it. The second opens up the one step the " +
          "first can only name, the search itself, and shows what it does to a " +
          "formula. They are meant to be read in that order."),

        p("Both use the same visual language. A formula is drawn as a picture of " +
          "boxes: a yellow box does one arithmetic step, a blue box holds a name " +
          "the run fills in per voxel, and a green box holds a number the fitting " +
          "is free to change. Count the boxes in a formula and you have its " +
          "length, which is the only thing that word means in this method and the " +
          "quantity the list of candidates is indexed by."),

        heading("How to read Figure 1"),

        p("Read it straight down. The yellow band runs once, before any simulation " +
          "exists, and none of it is inside CompLaB3D. The blue band runs once when " +
          "the solver starts. The green band runs in every wet voxel at every " +
          "reaction step, and is the box drawn in red in Figure 1 of the guide. " +
          "Each band names the source file it lives in, so a box can be traced to " +
          "the code that performs it."),

        p("The dashed line across the page is the boundary between the two phases, " +
          "and the .sym file is the only thing that crosses it. That is what lets " +
          "the fitting code, the metabolic model and the search stay outside the " +
          "solver entirely. The two red boxes are the guard rails: a name the run " +
          "does not know stops the run before the first step rather than becoming a " +
          "silent zero, and a concentration arriving from outside the range the fit " +
          "was shown is pulled back to the edge and counted, because a formula used " +
          "outside its range does not fail, it returns a confident number."),

        heading("How to read Figure 2"),

        p("The flow down the left has no colour, so that all the colour in the " +
          "figure belongs to the formulas on the right, which are the things being " +
          "changed. The steps are numbered because the cross-references are real: " +
          "step 6 scores a formula exactly as step 2 did, and FRESH draws one " +
          "exactly as step 1 did. Each picture is tied to its step by a dashed " +
          "leader, and the three lanes down the far left are the three loops: one " +
          "adds the next formula to the population, one closes the round, one " +
          "starts the next round."),

        p("Step 5 is the step that carries the method, so it states what comes in " +
          "and what goes out. In: the population as step 2 left it, every formula " +
          "already fitted and scored. Out: exactly one new formula, whichever of " +
          "the three moves ran, each of them tagging its own result. The shelf " +
          "holds one slot per number of boxes, which is what stops a short formula " +
          "a reader can argue with from being discarded for being less accurate " +
          "than a long one nobody can."),

        heading("Editing them"),

        p("Both flowcharts ship as editable PowerPoint shapes in " +
          "Method_3_Symbolic_flowcharts-figures.pptx, one slide each, in the same " +
          "way the other Method guides carry their figures. Nothing in either " +
          "picture is a bitmap: every box, diamond, arrow and label is a shape that " +
          "can be selected, retyped, recoloured or moved. Change a figure there and " +
          "export it again rather than editing the image."),

        p("The slides are 19 by 28.8 inches, which is unusual but deliberate. At a " +
          "hundred pixels to the inch a point on the slide is a point in the " +
          "drawing, so every label keeps the size it was designed at and nothing " +
          "has to be rescaled by eye.", { after: 0 }),
      ],
    },

    // --------------------------------------------------------- figure one ---
    {
      properties: page(PLATE_MARGIN),
      children: [
        picture("/home/claude/fig3/figA.png", 710, 867),
        caption(
          "Figure 1. The symbolic path as code: from a table of numbers to a change " +
          "in one voxel. The side panels open up the two steps that are otherwise " +
          "easy to take on trust: one generation of the search, and one walk of one " +
          "expression tree in one voxel."),
      ],
    },

    // --------------------------------------------------------- figure two ---
    {
      properties: page(PLATE_MARGIN),
      children: [
        picture("/home/claude/fig3/figB.png", 563, 853),
        caption(
          "Figure 2. Inside the search: what each step does to a formula. Step 5 is " +
          "drawn out in full because it is the step that carries the method, and the " +
          "one a reader is most likely to have to take on trust anywhere else."),
      ],
    },
  ],
});

Packer.toBuffer(doc).then((b) => {
  fs.writeFileSync("/home/claude/fig3/Method_3_Symbolic_flowcharts.docx", b);
  console.log("wrote Method_3_Symbolic_flowcharts.docx");
});
