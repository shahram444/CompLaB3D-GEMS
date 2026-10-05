#!/usr/bin/env python3
"""Make the v1.3.2 Damkohler banner work for a case that starts empty.

WHY THIS PATCH EXISTS
---------------------
v1.3.2 prints, next to the Peclet number, the Damkohler number the run actually
achieved:

    Da_r = k_r L^2 / D          k_r = r(Cref) / Cref

It builds Cref out of <initial_concentration>. Every case in the ABC sweep sets
every <initial_concentration> to 0.0 on purpose: nothing is pre-loaded, A enters
at x = 0 and B enters at x = nx-1 and the two travel toward each other. So
r(Cref) = 0 for every species, and the banner prints

    Da_r: not reported. No substrate is consumed at the initial composition.

for all 270 runs. The number the whole sweep is organised around would never be
reported by the solver, and the only record of it would be the one make_campaign.py
wrote into params.json -- which is the arithmetic being checked, not a check of it.

WHAT IT CHANGES
---------------
One thing: where the reference composition comes from. A species that starts at
zero falls back to its own Dirichlet feed value, whichever end it enters from:

    Cref[i] = C0[i]                              if C0[i] > 0
            = max(left BC, right BC) over the ends that are Dirichlet, otherwise

and the heading says which of the two it used. A case that does pre-load its
species is untouched, byte for byte: the fallback only fires where C0 is zero.

This is the right reference for this sweep on physical grounds as well as
practical ones. The feed value is the concentration the reaction sees where it
actually runs, once the plumes have arrived. Zero is the concentration it sees
before anything has happened, which is not a state worth naming a number after.

USAGE
-----
    python3 apply_damkohler_feed_reference.py /path/to/src/complab.cpp

Idempotent: run twice and the second run reports that it is already applied and
changes nothing. Writes complab.cpp.orig beside the file the first time.
"""
import os
import shutil
import sys

MARK = "[ABC] REFERENCE COMPOSITION FOR THE DAMKOHLER BANNER"

OLD_REF = """        std::vector<double> c0ref(vec_c0.begin(), vec_c0.end());
        std::vector<double> rate0((size_t) num_of_substrates, 0.0);"""

NEW_REF = """        /* """ + MARK + """
         *
         * Stock v1.3.2 read <initial_concentration> here. A case where nothing is
         * pre-loaded -- every species entering through a Dirichlet plane and starting
         * from an empty domain -- then has r(Cref) = 0 for every species and reports no
         * Damkohler number at all, which is exactly the case a Damkohler sweep is.
         *
         * So a species that starts at zero falls back to its own feed value, taken from
         * whichever end holds it. A species that is pre-loaded is unaffected: the
         * fallback fires only where C0 is zero, so an existing case prints the same
         * line it printed before.
         *
         * max() rather than the inlet end specifically, because CompLaB has no notion of
         * which end is the inlet: <Peclet> sets the flow direction and either x-normal
         * plane may be Dirichlet. A counter-current case feeds one species from each end
         * and both are the feed value of the species that uses them. */
        std::vector<double> c0ref((size_t) num_of_substrates, 0.0);
        bool cref_from_feed = false;
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            double cref = (double) vec_c0[iS];
            if (!(cref > 0)) {
                double feed = 0.0;
                if (vec_left_btype[iS]  == 0 && (double) vec_left_bcondition[iS]  > feed)
                    feed = (double) vec_left_bcondition[iS];
                if (vec_right_btype[iS] == 0 && (double) vec_right_bcondition[iS] > feed)
                    feed = (double) vec_right_bcondition[iS];
                if (feed > 0) { cref = feed; cref_from_feed = true; }
            }
            c0ref[(size_t) iS] = cref;
        }
        std::vector<double> rate0((size_t) num_of_substrates, 0.0);"""

OLD_C = """            const double c = (double) vec_c0[iS];
            if (!(r > 0) || !(c > 0)) continue;          /* produced, inert, or absent at t = 0    */"""

NEW_C = """            const double c = c0ref[(size_t) iS];     /* [ABC] the reference built above       */
            if (!(r > 0) || !(c > 0)) continue;          /* produced, inert, or absent everywhere */"""

OLD_HEAD_ASCII = '                pcout << "│ DAMKOHLER  Da_r = k_r L^2 / D,  k_r = r(C0)/C0,  L = "'
NEW_HEAD_ASCII = ('                pcout << "│ DAMKOHLER  Da_r = k_r L^2 / D,  k_r = r(Cref)/Cref"\n'
                  '                      << (cref_from_feed ? ",  Cref = feed" : ",  Cref = C0")\n'
                  '                      << ",  L = "')

OLD_NONE = """            pcout << "│ Damkohler: not reported. No substrate is consumed at the initial\\n";
            pcout << "│   composition, so there is no rate to build k_r from.\\n";"""

NEW_NONE = """            pcout << "│ Damkohler: not reported. No substrate is consumed at the reference\\n";
            pcout << "│   composition, so there is no rate to build k_r from. If every\\n";
            pcout << "│   <initial_concentration> is zero, check that the species that should\\n";
            pcout << "│   be consumed really do enter through a Dirichlet plane.\\n";"""


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[-4].strip())
        print("usage: apply_damkohler_feed_reference.py <path to src/complab.cpp>")
        return 2
    path = argv[1]
    if not os.path.isfile(path):
        print("not a file: %s" % path)
        return 1
    src = open(path, encoding="utf8").read()

    if MARK in src:
        print("already applied: %s" % path)
        return 0

    missing = [name for name, frag in
               (("reference vector", OLD_REF),
                ("per-species reference", OLD_C),
                ("banner heading", OLD_HEAD_ASCII),
                ("empty-reference message", OLD_NONE))
               if frag not in src]
    if missing:
        print("this does not look like CompLaB3D v1.3.2 complab.cpp.")
        print("could not find: %s" % ", ".join(missing))
        print("nothing was changed.")
        return 1

    out = src.replace(OLD_REF, NEW_REF)
    out = out.replace(OLD_C, NEW_C)
    out = out.replace(OLD_HEAD_ASCII, NEW_HEAD_ASCII)
    out = out.replace(OLD_NONE, NEW_NONE)

    backup = path + ".orig"
    if not os.path.exists(backup):
        shutil.copy2(path, backup)
        print("kept the original as %s" % backup)
    open(path, "w", encoding="utf8").write(out)
    print("patched %s" % path)
    print("  Cref now falls back to the Dirichlet feed where C0 is zero")
    print("  the banner says which of the two it used")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
