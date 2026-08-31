#!/usr/bin/env python3
"""
Check every example against the rules the solver actually enforces.

CompLB3D validates its input at start-up and stops with a message when
something does not add up: a vector of the wrong length, a solver that needs a
tag that is not there, a reaction type whose switch is off. Those checks are
good, but you only see them after a build and a run. This script applies the
same rules up front, so a broken example is caught in a second rather than
after a compile.

What it checks

  1. every XML is well formed, and no comment contains a double hyphen
     (which XML forbids and which this codebase has tripped over repeatedly)
  2. every tag name appears in the set the SOURCE actually reads, extracted
     mechanically from CompLB3D rather than typed out here, so a tag that was
     renamed upstream shows up as an unknown tag rather than as silence
  3. the vector-length and presence rules that terminate a run
  4. geometry files: the right number of voxels, and the material numbers the
     XML refers to actually present
  5. metabolic model files: nmet, nrxn, and the lengths of S, b, c, lb, ub
  6. the kinetics headers COMPILE, at -Wall -Wextra, and their own array
     guards are consistent with the number of substrates in the XML

What it does NOT check: this is not Palabos. It does not build lattices or
run anything. It is a fast structural check, not a substitute for the run.

    python3 validateExamples.py [path to CompLB3D]

Exits non-zero if anything fails.

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
"""

import os
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.abspath(__file__))

# tags that are numbered rather than fixed
DYNAMIC = re.compile(r'^(substrate|microbe|phase|species)\d+$')

# Keep this in step with rxntype::parse() in src/complab3d_metabolic.hh.  If the
# solver accepts a spelling and this does not, a perfectly good case fails the
# check; if this accepts one the solver does not, a broken case passes it.
VALID_REACTION = {
    "none", "no", "0", "kinetics", "kns", "1", "2", "3", "4", "5", "6", "7",
    "8", "9", "10", "11",
    "glpk", "fba", "fba_glpk", "surrogate", "srg", "ann", "cobrapy", "cpy",
    "fba_cobrapy", "glpk_and_kinetics", "glpk+kinetics",
    "surrogate_and_kinetics", "surrogate+kinetics",
    "cobrapy_and_kinetics", "cobrapy+kinetics",
    "symbolic", "sym", "expression",
    "symbolic_and_kinetics", "symbolic+kinetics",
    "graphnet", "gnn", "graph_network",
    "graphnet_and_kinetics", "graphnet+kinetics",
}
USES_SYM = {"symbolic", "sym", "expression", "8",
            "symbolic_and_kinetics", "symbolic+kinetics", "9"}
USES_GNN = {"graphnet", "gnn", "graph_network", "10",
            "graphnet_and_kinetics", "graphnet+kinetics", "11"}
VALID_SOLVER = {"fd", "finite difference", "finite_difference",
                "ca", "cellular automata", "cellular_automata",
                "lbm", "lb", "lattice boltzmann"}
VALID_BC = {"dirichlet", "neumann"}

USES_GLPK = {"glpk", "fba", "fba_glpk", "2", "glpk_and_kinetics", "glpk+kinetics", "5"}
USES_CPY = {"cobrapy", "cpy", "fba_cobrapy", "4", "cobrapy_and_kinetics",
            "cobrapy+kinetics", "7"}
USES_SRG = {"surrogate", "srg", "ann", "3", "surrogate_and_kinetics",
            "surrogate+kinetics", "6"}
USES_KNS = {"kinetics", "kns", "1", "glpk_and_kinetics", "glpk+kinetics", "5",
            "surrogate_and_kinetics", "surrogate+kinetics", "6",
            "cobrapy_and_kinetics", "cobrapy+kinetics", "7",
            "symbolic_and_kinetics", "symbolic+kinetics", "9",
            "graphnet_and_kinetics", "graphnet+kinetics", "11"}
USES_FBA = USES_GLPK | USES_CPY


def known_tags(src_dir):
    """Every XML tag the C++ actually reads, extracted from the source.

    Typed by hand this list rots. Extracted, it tracks the code."""
    tags = set()
    if not os.path.isdir(src_dir):
        return None
    for base, _, files in os.walk(src_dir):
        for fn in files:
            if not fn.endswith((".hh", ".cpp")):
                continue
            text = open(os.path.join(base, fn), errors="ignore").read()
            tags |= set(re.findall(r'\[\s*"([a-zA-Z_0-9]+)"\s*\]', text))
            # keys built in tables rather than written inline
            tags |= set(re.findall(r'"(enable_[a-z_]+|fix_[a-z_]+)"', text))
    return tags


class Case(object):
    def __init__(self, folder):
        self.folder = folder
        self.name = os.path.basename(folder)
        self.errors = []
        self.warnings = []

    def err(self, msg):
        self.errors.append(msg)

    def warn(self, msg):
        self.warnings.append(msg)


def text_of(node, path, default=None):
    el = node.find(path)
    if el is None or el.text is None:
        return default
    return el.text.strip()


def nums(node, path):
    t = text_of(node, path)
    if t is None:
        return None
    return t.split()


def truthy(s):
    return s is not None and s.strip().lower() in ("true", "yes", "1", "on")


# ===========================================================================
def check_xml_syntax(c, path):
    raw = open(path, "rb").read()
    # XML forbids "--" inside a comment. Python's parser reports it, but with a
    # message that does not say which comment, so find it first.
    for m in re.finditer(rb"<!--(.*?)-->", raw, re.S):
        if b"--" in m.group(1):
            line = raw[:m.start()].count(b"\n") + 1
            c.err("comment starting at line %d contains a double hyphen, which "
                  "XML forbids" % line)
    try:
        return ET.parse(path).getroot()
    except ET.ParseError as e:
        c.err("not well formed: %s" % e)
        return None


def check_tag_names(c, root, allowed):
    if allowed is None:
        return
    for el in root.iter():
        tag = el.tag
        if tag == "parameters" or DYNAMIC.match(tag):
            continue
        if tag not in allowed:
            c.err("unknown tag <%s>: nothing in the CompLB3D source reads it "
                  "(typo, or the tag was renamed upstream)" % tag)


def check_geometry(c, root):
    fn = text_of(root, "LB_numerics/domain/filename")
    nx = int(text_of(root, "LB_numerics/domain/nx"))
    ny = int(text_of(root, "LB_numerics/domain/ny"))
    nz = int(text_of(root, "LB_numerics/domain/nz"))
    path = os.path.join(c.folder, text_of(root, "path/input_path", "input"), fn)
    if not os.path.exists(path):
        c.err("geometry %s not found" % path)
        return set()
    vals = [int(v) for v in open(path).read().split()]
    if len(vals) != nx * ny * nz:
        c.err("geometry has %d values but nx*ny*nz = %d" % (len(vals), nx * ny * nz))
    present = set(vals)

    mats = root.find("LB_numerics/domain/material_numbers")
    if mats is not None:
        for el in mats:
            for v in (el.text or "").split():
                if int(v) not in present:
                    if el.tag in ("solid", "bounce_back", "pore"):
                        # Declaring a code that happens not to occur is fine:
                        # channel.dat has no solid voxels, and <solid>0</solid>
                        # still has to be there to say what 0 would mean.
                        continue
                    c.err("material_numbers/%s refers to %s, which does not "
                          "occur in the geometry, so that microbe would be "
                          "seeded nowhere" % (el.tag, v))
    return present


def check_model_file(c, folder, root, name, nsub, exch):
    path = os.path.join(folder, text_of(root, "path/input_path", "input"), name + ".xml")
    if not os.path.exists(path):
        c.err("model file %s not found" % path)
        return
    try:
        m = ET.parse(path).getroot()
    except ET.ParseError as e:
        c.err("model file %s is not well formed: %s" % (name, e))
        return
    if m.tag != "Metabolic_Model":
        c.err("model file %s has root <%s>; the loader looks for "
              "<Metabolic_Model>" % (name, m.tag))
        return
    nmet = int(text_of(m, "nmet"))
    nrxn = int(text_of(m, "nrxn"))
    if nmet <= 0 or nrxn <= 0:
        c.err("model %s: nmet and nrxn must both be positive" % name)
        return
    for tag, want in (("S", nmet * nrxn), ("b", nmet), ("c", nrxn),
                      ("lb", nrxn), ("ub", nrxn)):
        got = len((text_of(m, tag) or "").split())
        if got != want:
            c.err("model %s: <%s> has %d entries, expected %d"
                  % (name, tag, got, want))
    obj = text_of(m, "objLoc")
    if obj is not None and not (0 <= int(obj) < nrxn):
        c.err("model %s: objLoc %s is outside 0..%d" % (name, obj, nrxn - 1))
    for k, j in enumerate(exch or []):
        if int(j) >= nrxn:
            c.err("exchange_reaction_indices entry %s for substrate %d is "
                  "outside 0..%d of model %s" % (j, k, nrxn - 1, name))


def check_case(c, allowed):
    xml = os.path.join(c.folder, "CompLaB.xml")
    if not os.path.exists(xml):
        c.err("no CompLaB.xml")
        return
    root = check_xml_syntax(c, xml)
    if root is None:
        return
    check_tag_names(c, root, allowed)
    present = check_geometry(c, root)

    # ---- chemistry ------------------------------------------------------
    nsub = int(text_of(root, "chemistry/number_of_substrates"))
    chem = root.find("chemistry")
    declared = [el.tag for el in chem if el.tag.startswith("substrate")]
    if len(declared) != nsub:
        c.err("number_of_substrates is %d but %d <substrateN> blocks are "
              "present" % (nsub, len(declared)))
    for i in range(nsub):
        if "substrate%d" % i not in declared:
            c.err("substrate%d is missing; the blocks must be numbered 0..%d "
                  "with no gaps" % (i, nsub - 1))
    immobile = set()
    for i in range(nsub):
        s = chem.find("substrate%d" % i)
        if s is None:
            continue
        if truthy(text_of(s, "immobile")):
            immobile.add(i)
        for side in ("left", "right"):
            bt = text_of(s, "%s_boundary_type" % side)
            if bt is not None and bt.lower() not in VALID_BC:
                c.err("substrate%d %s_boundary_type '%s' is neither Dirichlet "
                      "nor Neumann" % (i, side, bt))
        if text_of(s, "substrate_diffusion_coefficients/in_pore") is None:
            c.err("substrate%d has no diffusion coefficient in_pore" % i)
    for tag in ("fix_concentration", "fix_lower_bounds"):
        v = nums(chem, tag)
        if v is not None and len(v) != nsub:
            c.err("<%s> has %d entries but there are %d substrates"
                  % (tag, len(v), nsub))

    # ---- microbiology ---------------------------------------------------
    biotic = truthy(text_of(root, "simulation_mode/biotic_mode"))
    micro = root.find("microbiology")
    nmic = int(text_of(root, "microbiology/number_of_microbes", "0"))
    if not biotic and nmic > 0:
        c.warn("biotic_mode is false but number_of_microbes is %d; the microbes "
               "will be ignored" % nmic)

    mats = root.find("LB_numerics/domain/material_numbers")
    seeded = {}
    if mats is not None:
        for el in mats:
            if el.tag.startswith("microbe"):
                seeded[el.tag] = (el.text or "").split()

    sw = {
        "glpk": truthy(text_of(root, "simulation_mode/enable_fba_glpk")),
        "cobrapy": truthy(text_of(root, "simulation_mode/enable_fba_cobrapy")),
        "surrogate": truthy(text_of(root, "simulation_mode/enable_surrogate")),
    }
    kinetics_on = truthy(text_of(root, "simulation_mode/enable_kinetics"))

    ca_count = 0
    visc_count = 0
    used = {"glpk": 0, "cobrapy": 0, "surrogate": 0, "kinetics": 0}

    if biotic:
        for i in range(nmic):
            m = micro.find("microbe%d" % i) if micro is not None else None
            if m is None:
                c.err("microbe%d is missing" % i)
                continue
            tag = "microbe%d" % i
            is_biofilm = tag in seeded

            solver = (text_of(m, "solver_type") or "").lower()
            if solver not in VALID_SOLVER:
                c.err("microbe%d solver_type '%s' is not FD, CA or LBM" % (i, solver))
            if solver.startswith("ca") or "cellular" in solver:
                ca_count += 1
            is_fd = solver in ("fd", "finite difference", "finite_difference")

            rxn = (text_of(m, "reaction_type") or "").lower()
            if rxn not in VALID_REACTION:
                c.err("microbe%d reaction_type '%s' is not a spelling the "
                      "parser accepts" % (i, rxn))
            # A symbolic or graphnet organism is useless without a file to read.
            # Unlike the metabolic solvers these need no <simulation_mode> switch:
            # naming the file is what turns them on, so that is what to check.
            for word, uses, block, key in (("symbolic", USES_SYM, "symbolic", "expressions_file"),
                                           ("graphnet", USES_GNN, "graphnet", "network_file")):
                if rxn in uses:
                    blk = root.find(block)
                    named = blk is not None and blk.find(key) is not None
                    if not named:
                        c.err("microbe%d asks for reaction_type %s but <%s><%s> "
                              "is not set, so it has no rate law to read"
                              % (i, word, block, key))
                    elif (text_of(blk, "enabled") or "").strip().lower() not in ("true", "1", "yes"):
                        c.err("microbe%d asks for reaction_type %s but <%s><enabled> "
                              "is not true" % (i, word, block))

            for key, s in (("glpk", USES_GLPK), ("cobrapy", USES_CPY),
                           ("surrogate", USES_SRG), ("kinetics", USES_KNS)):
                if rxn in s:
                    used[key] += 1
                    if key != "kinetics" and not sw[key]:
                        c.err("microbe%d asks for reaction_type %s but "
                              "<enable_fba_%s> / <enable_%s> is false"
                              % (i, rxn, key, key))
                    if key == "kinetics" and not kinetics_on:
                        c.err("microbe%d asks for a kinetics reaction_type but "
                              "<enable_kinetics> is false" % i)

            if is_fd and not is_biofilm:
                c.err("microbe%d uses solver_type FD but has no <microbe%d> "
                      "entry under material_numbers. FD is rejected for "
                      "planktonic microbes." % (i, i))
            if is_fd and text_of(m, "biomass_diffusion_coefficients/in_pore") is None:
                c.err("microbe%d uses FD, which requires "
                      "<biomass_diffusion_coefficients>" % i)

            if text_of(m, "viscosity_ratio_in_biofilm") is not None:
                visc_count += 1
                if not is_biofilm:
                    c.err("microbe%d has <viscosity_ratio_in_biofilm> but no "
                          "<microbe%d> entry under material_numbers. The parser "
                          "counts these and compares against the number of "
                          "seeded microbes." % (i, i))
            elif is_biofilm:
                c.err("microbe%d is seeded under material_numbers but has no "
                      "<viscosity_ratio_in_biofilm>. The counts must match." % i)

            b0 = nums(m, "initial_densities")
            if is_biofilm and b0 is not None and len(b0) != len(seeded[tag]):
                c.err("microbe%d has %d initial_densities but %d material "
                      "numbers; they must correspond one to one"
                      % (i, len(b0), len(seeded[tag])))

            for t in ("half_saturation_constants", "maximum_uptake_flux",
                      "fba_maximum_uptake_flux", "substrate_lower_bounds",
                      "substrate_upper_bounds", "exchange_reaction_indices"):
                v = nums(m, t)
                if v is not None and len(v) != nsub:
                    c.err("microbe%d <%s> has %d entries but there are %d "
                          "substrates" % (i, t, len(v), nsub))

            if rxn in USES_FBA:
                mf = text_of(m, "model_filename")
                if mf is None:
                    c.err("microbe%d uses a flux balance reaction_type but has "
                          "no <model_filename>" % i)
                else:
                    check_model_file(c, c.folder, root, mf, nsub,
                                     nums(m, "exchange_reaction_indices"))
                if nums(m, "exchange_reaction_indices") is None:
                    c.err("microbe%d uses flux balance analysis but has no "
                          "<exchange_reaction_indices>" % i)
                od = (text_of(m, "objective_direction") or "maximize").lower()
                if od not in ("maximize", "max", "minimize", "min"):
                    c.err("microbe%d objective_direction '%s' is neither "
                          "maximize nor minimize" % (i, od))

        if ca_count > 0 and text_of(root, "microbiology/thrd_biofilm_fraction") is None:
            c.err("a microbe uses the cellular automaton, which requires "
                  "<thrd_biofilm_fraction>")
        cam = text_of(root, "microbiology/CA_method")
        if ca_count > 0 and cam is None:
            c.warn("a microbe uses the cellular automaton but <CA_method> is "
                   "not set")

    # a switch on with nothing using it is wasted work, not an error
    for key in ("glpk", "cobrapy", "surrogate"):
        if sw[key] and used[key] == 0:
            c.warn("<enable_%s> is true but no microbe uses it" % key)

    # ---- equilibrium ----------------------------------------------------
    eq = root.find("equilibrium")
    if eq is not None and truthy(text_of(eq, "enabled")):
        comps = (text_of(eq, "components") or "").split()
        if not comps:
            c.err("equilibrium is enabled but <components> is empty")
        names = [text_of(chem.find("substrate%d" % i), "name_of_substrates")
                 for i in range(nsub)]
        for comp in comps:
            if comp not in names:
                c.err("equilibrium component '%s' does not match any "
                      "<name_of_substrates>" % comp)
        st = eq.find("stoichiometry")
        if st is None:
            c.err("equilibrium is enabled but there is no <stoichiometry>")
        else:
            for el in st:
                row = (el.text or "").split()
                if len(row) != len(comps):
                    c.err("stoichiometry/%s has %d numbers but there are %d "
                          "components" % (el.tag, len(row), len(comps)))

    # ---- precipitation and dissolution ----------------------------------
    pr = root.find("precipitation")
    if pr is not None and truthy(text_of(pr, "enabled")):
        ss = text_of(pr, "solid_substrate")
        if ss is None:
            c.err("precipitation is enabled but <solid_substrate> is not set")
        else:
            if not (0 <= int(ss) < nsub):
                c.err("precipitation solid_substrate %s is outside 0..%d"
                      % (ss, nsub - 1))
            elif int(ss) not in immobile:
                c.err("precipitation points at substrate%s, which is not "
                      "<immobile>true</immobile>. A mineral that advects will "
                      "not accumulate where it forms." % ss)
        if text_of(pr, "max_precipRho") is None:
            c.err("precipitation is enabled but <max_precipRho> is not set")

    ds = root.find("dissolution")
    if ds is not None and truthy(text_of(ds, "enabled")):
        phases = [el for el in ds if el.tag.startswith("phase")]
        if not phases:
            c.err("dissolution is enabled but no <phaseN> block is declared, "
                  "so nothing can dissolve")
        for ph in phases:
            sub = text_of(ph, "substrate")
            if sub is None:
                c.err("%s has no <substrate>" % ph.tag)
            elif not (0 <= int(sub) < nsub):
                c.err("%s substrate %s is outside 0..%d" % (ph.tag, sub, nsub - 1))
            elif int(sub) not in immobile:
                c.err("%s points at substrate%s, which is not "
                      "<immobile>true</immobile>" % (ph.tag, sub))
            fd = text_of(ph, "full_density")
            if fd is None or float(fd) <= 0:
                c.err("%s needs a positive <full_density>" % ph.tag)
            init = text_of(ph, "initial_fill")
            if init is not None and fd is not None and float(init) > float(fd):
                c.err("%s initial_fill %s exceeds full_density %s"
                      % (ph.tag, init, fd))
            mn = text_of(ph, "material_number")
            if mn is not None and int(mn) not in present:
                c.warn("%s has material_number %s, which does not occur in the "
                       "geometry. That is correct for a precipitate that grows "
                       "into existing solid, and wrong for an original grain."
                       % (ph.tag, mn))
        rf = text_of(ds, "reopen_fraction")
        if rf is not None and not (0 < float(rf) < 1):
            c.err("reopen_fraction %s must lie strictly between 0 and 1" % rf)

    return nsub


def check_kinetics_headers(c, nsub):
    """Compile the headers, and check their own guards against the XML.

    A header written for 3 substrates and dropped into a 2-substrate case reads
    past the end of the array. The shipped uranium network guards with
    C.size() < 14, so this catches exactly the mistake of forgetting to replace
    it."""
    stub = """#ifndef PLBSTUB
#define PLBSTUB
namespace plb { typedef long long plint; }
#endif
"""
    for fn in ("defineKinetics.hh", "defineAbioticKinetics.hh"):
        path = os.path.join(c.folder, fn)
        if not os.path.exists(path):
            c.err("%s is missing; the shipped uranium network would be compiled "
                  "in instead, and it indexes C[0] to C[94]" % fn)
            continue
        text = open(path).read()
        for m in re.finditer(r'C\.size\(\)\s*<\s*(\d+)', text):
            need = int(m.group(1))
            if need > nsub:
                c.err("%s guards on C.size() < %d but this case has only %d "
                      "substrates; the file was written for different chemistry"
                      % (fn, need, nsub))

        tmp = tempfile.mkdtemp(prefix="kincheck_")
        with open(os.path.join(tmp, "plbstub.h"), "w") as f:
            f.write(stub)
        src = os.path.join(tmp, "t.cpp")
        with open(src, "w") as f:
            f.write('#include "plbstub.h"\n#include <vector>\n#include "%s"\nint main(){return 0;}\n'
                    % path)
        r = subprocess.run([os.environ.get("CXX", "g++"), "-fsyntax-only",
                            "-std=c++11", "-Wall", "-Wextra", src],
                           capture_output=True, text=True)
        if r.returncode != 0:
            first = next((l for l in r.stderr.splitlines() if "error:" in l), "")
            c.err("%s does not compile: %s" % (fn, first.strip()))
        elif r.stderr.strip():
            first = r.stderr.strip().splitlines()[0]
            c.warn("%s compiles with a warning: %s" % (fn, first))


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(ROOT), "CompLB3D")
    allowed = known_tags(os.path.join(src, "src")) if os.path.isdir(src) else None
    if allowed:
        allowed |= known_tags(src) or set()
        print("Tag whitelist: %d names read from %s" % (len(allowed), src))
    else:
        print("CompLB3D source not found at %s; skipping the tag-name check.\n"
              "Pass the path as the first argument to enable it." % src)

    folders = sorted(d for d in os.listdir(ROOT)
                     if re.match(r'^\d\d_', d) and os.path.isdir(os.path.join(ROOT, d)))
    if not folders:
        sys.exit("No example folders found. Run makeExamples.py first.")

    bad = 0
    total_warn = 0
    for f in folders:
        c = Case(os.path.join(ROOT, f))
        try:
            nsub = check_case(c, allowed)
            if nsub:
                check_kinetics_headers(c, nsub)
        except Exception as e:                       # noqa: BLE001
            c.err("validator crashed on this case: %r" % e)
        status = "FAIL" if c.errors else ("ok  " if not c.warnings else "ok* ")
        print("  [%s] %s" % (status, c.name))
        for e in c.errors:
            print("      ERROR   %s" % e)
        for w in c.warnings:
            print("      note    %s" % w)
        bad += 1 if c.errors else 0
        total_warn += len(c.warnings)

    print("\n%d case(s) checked, %d failed, %d note(s)."
          % (len(folders), bad, total_warn))
    print("Structural check only. It does not build lattices and does not "
          "replace running the cases.")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
