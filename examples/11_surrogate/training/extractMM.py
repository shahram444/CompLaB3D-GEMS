#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extractMM.py  --  metabolic model converter for CompLaB3D
=========================================================

WHAT THIS SCRIPT DOES (for readers who do not program)
------------------------------------------------------
A genome-scale metabolic model is a big table of numbers describing every
chemical reaction a microbe can run.  The community distributes those models as
SBML (.xml), MATLAB (.mat) or JSON (.json) files.  The CompLaB3D solver cannot
read those formats directly: it wants one small, flat XML file holding just the
numbers the flux-balance solver (GLPK) needs.

This script is the translator.  You run it ONCE per microbe, before the
simulation, and it writes the flat XML that you then name in CompLaB.xml.

    python3 extractMM.py iAF987.xml                  ->  iAF987_mm.xml
    python3 extractMM.py iAF987.mat -o geobacter.xml ->  geobacter.xml
    python3 extractMM.py --help                      ->  usage information

THE OUTPUT FORMAT (do not change it -- the C++ loader depends on it)
--------------------------------------------------------------------
    <Metabolic_Model>
        <name>    ... short text label
        <nmet>    ... number of metabolites  (= number of rows of S)
        <nrxn>    ... number of reactions    (= number of columns of S)
        <objLoc>  ... which reaction is the biomass/objective reaction (counting from 0)
        <S>       ... the stoichiometric matrix, FLATTENED ROW-MAJOR, i.e.
                      metabolite-major: all of metabolite 0's coefficients (one
                      per reaction), then all of metabolite 1's, and so on.
                      The C++ side reads it back as S[irow*nrxn + icol]
                      (see load_metabolic_models() in src/complab_functions.hh).
        <b>       ... right-hand side of S*v = b; all zeros at steady state
        <c>       ... objective coefficients, one per reaction
        <lb>,<ub> ... lower / upper flux bounds, one per reaction

    Every numeric list is written as space-separated text with ONE leading and
    ONE trailing space -- the Palabos XMLreader expects that padding.

WHAT CHANGED RELATIVE TO THE 2D CompLaB VERSION
-----------------------------------------------
Five real defects in the 2D script were fixed.  Every fix is marked in the code
with a "# [FIX-3D]" comment so you can find it:

  1. objLoc was only set when some reaction had objective coefficient EXACTLY
     1.0.  Otherwise the variable never existed and the script died with a
     confusing NameError near the end.  Now the objective is detected properly
     (largest nonzero coefficient, falling back to model.objective) with a clear
     error message if there genuinely is no objective.
  2. The output file name was derived from the input file name, so an SBML input
     called "model.xml" SILENTLY OVERWROTE ITSELF with its own converted output
     -- destroying the original model.  The script now refuses to overwrite its
     input, writes "<name>_mm.xml" by default, and accepts -o/--output.
  3. Numbers were written with "%g", which truncates to 6 significant digits.
     That is lossy for stoichiometry and for flux bounds.  Full precision now.
  4. Command-line handling is done with argparse, so there is a real --help.
  5. A sanity report and a round-trip check of S are printed at the end so you
     can confirm the conversion actually worked.

REQUIREMENTS
------------
    pip install cobra numpy
"""

import os
import sys
import argparse

import numpy as np
from xml.dom import minidom


# ---------------------------------------------------------------------------
# NUMBER FORMATTING
#   [FIX-3D] The 2D script used "%g", which keeps only 6 significant digits:
#   a bound of -1000.0000001 became -1000, and a stoichiometric coefficient of
#   0.0060205 became 0.0060205 but 1.234567891e-3 became 0.00123457.  For flux
#   balance analysis that is a silent change of the model.
#   "%.17g" is the shortest format that is GUARANTEED to reproduce an IEEE-754
#   double exactly when read back, so the XML is a lossless copy of the model.
#   TRADEOFF: the file gets roughly 2-3x larger than with "%g" (a genome-scale
#   model such as iAF987 goes from a few MB to roughly 10-20 MB).  That is a
#   one-off cost paid on disk at conversion time; the solver reads the file once
#   at startup, so it does not affect run speed.  If disk space is genuinely
#   tight you may lower this to "%.12g" -- but never back to "%g".
# ---------------------------------------------------------------------------
NUMFMT = "%.17g"

# [FIX-3D] Infinite bounds must NOT be written as "inf".  The C++ loader reads these
# lists through Palabos' XMLreader, which is std::stringstream extraction, and
# libstdc++'s num_get<double> does not accept "inf" -- the stream fails at the first
# one, the list comes back short, and load_metabolic_models3D() aborts the run with a
# length error that says nothing about the real cause.  Map non-finite bounds onto a
# large finite sentinel instead; complab3d_glpkcpp.hh treats anything at or beyond
# this magnitude as a free bound.
BIG = 1.0e30

def _finite(x):
    """Clamp +/-inf and NaN onto the finite sentinel the C++ side understands."""
    v = float(x)
    if v != v:            # NaN
        return 0.0
    if v > BIG:
        return BIG
    if v < -BIG:
        return -BIG
    return v


def fmt_list(values):
    """Turn a list of numbers into the space-padded text the C++ reader wants."""
    return " " + " ".join(NUMFMT % _finite(x) for x in values) + " "


# ---------------------------------------------------------------------------
# COMMAND LINE
#   [FIX-3D] The 2D script read sys.argv[1] directly: no --help, no way to
#   choose the output name, and an unhelpful IndexError when run with no
#   arguments.  argparse gives all three for free.
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="extractMM.py",
        description=(
            "Convert an SBML (.xml), MATLAB (.mat) or JSON (.json) genome-scale "
            "metabolic model into the flat XML that CompLaB3D reads."),
        epilog=(
            "Examples:\n"
            "  python3 extractMM.py iAF987.xml\n"
            "      -> writes iAF987_mm.xml\n"
            "  python3 extractMM.py iAF987.mat -o models/geobacter.xml\n"
            "      -> writes models/geobacter.xml\n"
            "  python3 extractMM.py iAF987.xml --objective BIOMASS_Gm_GS15_core_79p20M\n"
            "      -> forces which reaction is the growth (objective) reaction\n"),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input",
                   help="metabolic model file to convert (.xml / .mat / .json)")
    p.add_argument("-o", "--output", default=None,
                   help="output XML file (default: <input basename>_mm.xml). "
                        "The script refuses to overwrite its own input.")
    p.add_argument("-n", "--name", default=None,
                   help="text label written into <name> (default: input basename)")
    p.add_argument("--objective", default=None,
                   help="reaction ID to use as the objective (biomass) reaction. "
                        "Only needed when auto-detection fails or picks the wrong one.")
    p.add_argument("-f", "--force", action="store_true",
                   help="overwrite the output file if it already exists")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# MODEL LOADING
# ---------------------------------------------------------------------------
def load_model(fin):
    """Read the model with cobrapy, dispatching on the file extension."""
    import cobra  # imported here so that --help works without cobra installed

    ext = os.path.splitext(fin)[1].lower()
    if ext == ".xml" or ext == ".sbml":
        return cobra.io.read_sbml_model(fin)
    elif ext == ".mat":
        return cobra.io.load_matlab_model(fin)
    elif ext == ".json":
        return cobra.io.load_json_model(fin)
    else:
        sys.exit("ERROR: '%s' has extension '%s'. Only .xml (SBML), .mat and "
                 ".json files are supported." % (fin, ext if ext else "(none)"))


# ---------------------------------------------------------------------------
# OBJECTIVE DETECTION
#   [FIX-3D] The 2D script did:
#
#       objcoef = model.reactions[j].objective_coefficient
#       if objcoef == 1:
#           objloc = j
#
#   so objloc existed ONLY if some reaction had a coefficient of exactly 1.0.
#   Models whose biomass reaction is weighted 0.5, or -1 (minimisation), or
#   whose objective lives only in model.objective and never propagates to the
#   per-reaction coefficients, left objloc undefined -- and the script then died
#   with "NameError: name 'objloc' is not defined" ~30 lines later, long after
#   the expensive model parsing, with no hint of what was wrong.
#
#   The replacement tries, in order:
#     (a) an explicit --objective REACTION_ID given by the user,
#     (b) the reaction with the largest |objective_coefficient| (this covers the
#         old exact-1.0 case as a special case, and also 0.5, -1, etc.),
#     (c) cobra's linear_reaction_coefficients(model), which reads the optlang
#         objective expression directly (some SBML files carry the objective
#         only there),
#   and otherwise stops with a clear, actionable message.
# ---------------------------------------------------------------------------
def find_objective(model, c, forced_id=None):
    """Return (objloc, reaction_id, how_it_was_found)."""
    rxn_ids = [r.id for r in model.reactions]

    # (a) user override -------------------------------------------------------
    if forced_id is not None:
        if forced_id not in rxn_ids:
            sys.exit("ERROR: --objective '%s' is not a reaction in this model.\n"
                     "       Reactions whose IDs contain 'BIOMASS' or 'biomass':\n"
                     "       %s" % (forced_id,
                                    [r for r in rxn_ids if "iomass" in r or "IOMASS" in r]
                                    or "(none found)"))
        j = rxn_ids.index(forced_id)
        return j, forced_id, "given on the command line with --objective"

    # (b) largest nonzero objective coefficient -------------------------------
    nonzero = [j for j, v in enumerate(c) if v != 0.0]
    if nonzero:
        j = max(nonzero, key=lambda k: abs(c[k]))
        how = "largest nonzero objective coefficient (c = %s)" % (NUMFMT % c[j])
        if len(nonzero) > 1:
            how += "; WARNING: %d reactions have a nonzero coefficient %s -- " \
                   "check this is the growth reaction, or force it with --objective" \
                   % (len(nonzero), [rxn_ids[k] for k in nonzero])
        return j, rxn_ids[j], how

    # (c) fall back on model.objective ---------------------------------------
    try:
        from cobra.util.solver import linear_reaction_coefficients
        lin = linear_reaction_coefficients(model)          # {Reaction: coefficient}
        if lin:
            rxn = max(lin, key=lambda r: abs(lin[r]))
            j = rxn_ids.index(rxn.id)
            # keep c consistent with what we detected, so the XML <c> vector
            # actually expresses the objective the solver should maximise
            c[j] = float(lin[rxn])
            return j, rxn.id, "model.objective expression (coefficient %s)" % (NUMFMT % c[j])
    except Exception as exc:                                # noqa: BLE001
        print("  note: could not inspect model.objective (%s)" % exc)

    # (d) genuinely no objective ---------------------------------------------
    candidates = [r for r in rxn_ids if "iomass" in r or "IOMASS" in r or r.lower().startswith("bio")]
    sys.exit(
        "ERROR: this model has NO objective reaction.\n"
        "       Every reaction has objective_coefficient == 0 and model.objective\n"
        "       is empty, so CompLaB3D would have nothing to maximise.\n"
        "       Fix it by naming the growth reaction explicitly, e.g.\n"
        "           python3 extractMM.py <input> --objective <REACTION_ID>\n"
        "       Reactions that look like biomass reactions in this model:\n"
        "           %s" % (candidates if candidates else "(none found -- check the model)"))


# ---------------------------------------------------------------------------
# OUTPUT PATH
#   [FIX-3D] The 2D script built its output name as
#
#       name           = fin.split('.', 1)[0]
#       save_path_file = name + ".xml"
#
#   Two bugs in two lines:
#     * split('.', 1) cuts at the FIRST dot, so a path like "../models/iAF987.xml"
#       became ".." and the model was written to "...xml" in the wrong place;
#     * for the normal SBML case the input IS "<name>.xml", so the output path
#       equalled the input path and the script SILENTLY OVERWROTE THE INPUT
#       MODEL with its own converted output.  The original SBML was destroyed,
#       and re-running the script then failed because the converted file is not
#       a valid SBML model.
#
#   Now: split off only the LAST extension, of only the file name; default to
#   "<name>_mm.xml"; and hard-refuse any output path that resolves to the input.
# ---------------------------------------------------------------------------
def resolve_output(fin, requested, force):
    base = os.path.splitext(os.path.basename(fin))[0]
    if requested is None:
        fout = os.path.join(os.path.dirname(os.path.abspath(fin)), base + "_mm.xml")
    else:
        fout = os.path.abspath(requested)

    if os.path.abspath(fin) == fout:
        sys.exit("ERROR: refusing to overwrite the input model '%s' with its own\n"
                 "       converted output. Choose a different name with -o/--output." % fin)
    if os.path.exists(fout) and not force:
        sys.exit("ERROR: output file '%s' already exists.\n"
                 "       Use -f/--force to overwrite it, or -o to pick another name." % fout)
    return base, fout


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)
    fin = args.input
    if not os.path.isfile(fin):
        sys.exit("ERROR: input file '%s' does not exist." % fin)

    base, fout = resolve_output(fin, args.output, args.force)
    name = args.name if args.name is not None else base

    print("reading metabolic model : %s" % os.path.abspath(fin))
    model = load_model(fin)

    # --- the stoichiometric matrix: rows = metabolites, columns = reactions ---
    import cobra.util.array
    S = cobra.util.array.create_stoichiometric_matrix(model)
    S = np.asarray(S, dtype=float)

    # --- per-reaction vectors -------------------------------------------------
    lb, ub, c, b = [], [], [], []
    for rxn in model.reactions:
        lb.append(float(rxn.lower_bound))
        ub.append(float(rxn.upper_bound))
        c.append(float(rxn.objective_coefficient))
    # [FIX-3D] b is the right-hand side of the mass-balance constraint S*v = b, so
    # it has ONE ENTRY PER METABOLITE (one per ROW of S). The 2D script appended a
    # zero inside the loop over REACTIONS, giving it nrxn entries instead of nmet.
    # That went unnoticed only because genome-scale models usually have more
    # reactions than metabolites and every entry is zero; the C++ side reads
    # b[0..nrow-1] with nrow = nmet (see glpk_call() in src/complab_glpkcpp.hh),
    # so any model with nmet > nrxn read past the end of the vector.
    b = [0.0] * S.shape[0]          # steady state: S*v = 0, one entry per METABOLITE

    nmet, nrxn = S.shape
    if len(lb) != nrxn:
        sys.exit("ERROR: internal inconsistency -- S has %d columns but the model "
                 "has %d reactions." % (nrxn, len(lb)))

    objloc, objid, objhow = find_objective(model, c, args.objective)

    # --- flatten S row-major (metabolite-major), exactly as the C++ expects ---
    S_flat = S.flatten()            # numpy default order='C' == row-major
    str_S = fmt_list(S_flat)

    # -----------------------------------------------------------------------
    # BUILD THE XML
    # -----------------------------------------------------------------------
    root = minidom.Document()
    model_elm = root.createElement("Metabolic_Model")
    root.appendChild(model_elm)

    def add(tag, text):
        child = root.createElement(tag)
        child.appendChild(root.createTextNode(text))
        model_elm.appendChild(child)

    add("name",   name)
    add("nmet",   str(nmet))
    add("nrxn",   str(nrxn))
    add("objLoc", str(objloc))
    add("S",      str_S)
    add("b",      fmt_list(b))
    add("c",      fmt_list(c))
    add("lb",     fmt_list(lb))
    add("ub",     fmt_list(ub))

    with open(fout, "w") as f:
        f.write(root.toprettyxml(indent="\t"))

    # -----------------------------------------------------------------------
    # ROUND-TRIP CHECK OF S
    #   [FIX-3D] The 2D script never verified the flattening. A wrong order
    #   (metabolite-major vs reaction-major) or a lossy number format produces a
    #   file that loads fine and gives silently wrong fluxes. Here we parse the
    #   text we just wrote back into numbers, reshape it the way the C++ loader
    #   does -- S[irow*nrxn + icol], see load_metabolic_models() in
    #   src/complab_functions.hh -- and require it to equal the original matrix
    #   EXACTLY (which "%.17g" guarantees for IEEE-754 doubles).
    # -----------------------------------------------------------------------
    back = np.array([float(t) for t in str_S.split()], dtype=float)
    ok = (back.size == nmet * nrxn)
    if ok:
        back = back.reshape((nmet, nrxn))       # row-major, same as C++
        ok = np.array_equal(back, S)
    roundtrip = "PASS" if ok else "FAIL"

    # -----------------------------------------------------------------------
    # SANITY REPORT
    #   [FIX-3D] New: the 2D script printed nothing at all, so a conversion that
    #   produced garbage looked exactly like one that worked. Read these numbers
    #   and compare them against what you expect from the published model.
    # -----------------------------------------------------------------------
    nnz = int(np.count_nonzero(S))
    finite_lb = [v for v in lb if np.isfinite(v)]
    finite_ub = [v for v in ub if np.isfinite(v)]
    n_rev = sum(1 for v in lb if v < 0.0)
    n_fixed0 = sum(1 for l, u in zip(lb, ub) if l == 0.0 and u == 0.0)

    print("")
    print("=================== CONVERSION REPORT ===================")
    print("  input                : %s" % os.path.abspath(fin))
    print("  output               : %s" % fout)
    print("  <name>               : %s" % name)
    print("  metabolites  (nmet)  : %d" % nmet)
    print("  reactions    (nrxn)  : %d" % nrxn)
    print("  nonzeros in S        : %d  (density %.4f %%)"
          % (nnz, 100.0 * nnz / float(nmet * nrxn) if nmet * nrxn else 0.0))
    print("  |S| range (nonzero)  : %.6g .. %.6g"
          % ((np.min(np.abs(S[S != 0.0])), np.max(np.abs(S))) if nnz else (0.0, 0.0)))
    print("  objective reaction   : '%s'  (index %d, counting from 0)" % (objid, objloc))
    print("      detected via     : %s" % objhow)
    print("  lower bounds lb      : %.6g .. %.6g   (%d reversible, i.e. lb < 0)"
          % (min(finite_lb) if finite_lb else 0.0,
             max(finite_lb) if finite_lb else 0.0, n_rev))
    print("  upper bounds ub      : %.6g .. %.6g" % (min(finite_ub) if finite_ub else 0.0,
                                                     max(finite_ub) if finite_ub else 0.0))
    print("  blocked reactions    : %d  (lb == ub == 0)" % n_fixed0)
    print("  b vector             : %d entries, all zero (steady state)" % len(b))
    print("  S round-trip check   : %s" % roundtrip)
    print("  file size            : %.2f MB" % (os.path.getsize(fout) / (1024.0 * 1024.0)))
    print("=========================================================")
    print("")
    # [FIX-3D] CompLaB.xml needs <exchange_reaction_indices>, one POSITIONAL index
    # into model.reactions per substrate.  The converted file contains no reaction
    # IDs, so without this table the user has no way to recover those indices -- and
    # an off-by-one here is silent: the flux of the wrong reaction is attributed to
    # the substrate.  This is the highest-risk manual step in the whole workflow, so
    # print the candidates.
    ex_idx = []
    for j, rxn in enumerate(model.reactions):
        rid = str(rxn.id)
        up = rid.upper()
        if up.startswith("EX_") or up.startswith("R_EX_") or up.startswith("DM_") \
           or up.startswith("SINK_") or up.startswith("EXCHANGE"):
            ex_idx.append((j, rid))
    print("")
    print("EXCHANGE REACTIONS  (for <exchange_reaction_indices> in CompLaB.xml)")
    print("---------------------------------------------------------------")
    if ex_idx:
        for j, rid in ex_idx:
            print("  index %-6d  %s" % (j, rid))
        print("")
        print("  %d exchange reaction(s) found. In CompLaB.xml give ONE index per substrate," % len(ex_idx))
        print("  in substrate order, using -1 for a substrate this organism does not exchange.")
    else:
        print("  None matched the usual EX_ / DM_ / SINK_ naming. List them yourself with:")
        print("    python3 -c \"import cobra;m=cobra.io.read_sbml_model('%s');"
              "print([(j,r.id) for j,r in enumerate(m.reactions) if r.boundary])\"" % args.input)
    print("")
    print("Next step: put the OUTPUT file name into CompLaB.xml, in this microbe's")
    print("<model_filename> entry (without the .xml), and fill in")
    print("<exchange_reaction_indices> from the table above.")

    if not ok:
        sys.exit("ERROR: the S matrix did NOT survive the flatten/reshape round trip.\n"
                 "       The written file is unusable. Please report this.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
