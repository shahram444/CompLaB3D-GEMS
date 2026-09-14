#!/usr/bin/env python3
"""
check_wizard.py -- the guided helper, driven by a script instead of a person.

Part of the CompLaB program.  GNU Affero General Public License v3 or later.

make_rate_law.py exists so that somebody who has never seen a command line option can
still produce a rate law file.  That makes it the one tool here whose failure mode is
silent: nobody runs it in a hurry, so a break would sit unnoticed until the next new user
met it.  This drives it with a canned set of answers, the way a person would type them,
and checks that what comes out the far end is a file the solver accepts.

Both kinds are covered, because they take different questions: an organism with a yield,
and a reaction in open water with none.

    python3 check_wizard.py
"""

import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WIZ = os.path.join(HERE, "..", "tools", "make_rate_law.py")

fails = 0


def ok(cond, what):
    global fails
    print("  %-58s %s" % (what, "ok" if cond else "FAIL"))
    if not cond:
        fails += 1


def run(answers, cwd):
    """Feed the helper one answer per line, exactly as a person would type them."""
    p = subprocess.Popen([sys.executable, WIZ], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         cwd=cwd, universal_newlines=True)
    out, _ = p.communicate("\n".join(answers) + "\n")
    return p.returncode, out


def main():
    print("the guided helper")
    tmp = tempfile.mkdtemp()

    # --- a table with an organism -------------------------------------------------------
    bio = os.path.join(tmp, "bio.csv")
    with open(bio, "w") as f:
        f.write("acetate,o2,growth\n")
        for i in range(1, 121):
            a = 5e-3 * i / 120.0
            o = 2e-3 * (121 - i) / 120.0
            f.write("%.8g,%.8g,%.8g\n" % (a, o, 0.35 * a / (0.05 + a) * o / (0.01 + o)))

    out_bio = os.path.join(tmp, "bio.sym")
    rc, log = run(["bio.csv",        # which file
                   "", "",           # target and inputs: take the defaults
                   "yes", "Bug",     # an organism, called Bug
                   "-1", "-2",       # one acetate, two oxygen, both used up
                   "no",             # nothing else takes part
                   "acetate", "0.4", # the yield, against acetate
                   "per_hour",
                   "quick", "no", "1",
                   out_bio, "yes"], tmp)
    ok(rc == 0, "an organism's law: the helper finishes")
    ok(os.path.isfile(out_bio), "an organism's law: the file is written")
    text = open(out_bio).read() if os.path.isfile(out_bio) else ""
    ok("reaction acetate -1 o2 -2" in text, "the reaction it was told is in the file")
    ok("yield    acetate 0.4" in text, "the yield is in the file")
    ok("biomass  Bug" in text, "the organism is named in the file")
    ok("rate   growth =" in text, "the fitted line is called growth")
    ok("rate   acetate" not in text, "no substrate line was written by hand")
    ok("--reaction" in log and "fit_symbolic.py" in log,
       "it prints the command that would repeat this without questions")

    # --- a table with no organism -------------------------------------------------------
    ab = os.path.join(tmp, "ab.csv")
    with open(ab, "w") as f:
        f.write("Fe,HS,extent\n")
        for i in range(1, 121):
            fe = 1e-3 * i / 120.0
            hs = 1e-3 * (121 - i) / 120.0
            f.write("%.8g,%.8g,%.8g\n" % (fe, hs, 160.0 * fe * hs))

    out_ab = os.path.join(tmp, "ab.sym")
    rc, log = run(["ab.csv",
                   "", "",
                   "no",             # no organism
                   "-1", "-1",       # Fe and HS both used up
                   "yes", "FeS", "1",  # and FeS is made
                   "no",
                   "per_second",
                   "quick", "no", "1",
                   out_ab, "yes"], tmp)
    ok(rc == 0, "a law with no organism: the helper finishes")
    ok(os.path.isfile(out_ab), "a law with no organism: the file is written")
    text = open(out_ab).read() if os.path.isfile(out_ab) else ""
    ok("reaction Fe -1 HS -1 FeS 1" in text, "the product it was told about is in the reaction")
    ok("rate   extent =" in text, "the fitted line is called extent")
    ok("yield" not in text and "biomass" not in text,
       "no yield and no biomass, which is what makes it abiotic")

    # --- it must not write anything when told not to ------------------------------------
    out_no = os.path.join(tmp, "never.sym")
    rc, log = run(["bio.csv", "", "", "yes", "Bug", "-1", "-2", "no", "acetate", "0.4",
                   "per_hour", "quick", "no", "1", out_no, "no"], tmp)
    ok(rc == 0 and not os.path.isfile(out_no),
       "answering no at the end writes nothing at all")

    # --- and it must refuse a table it cannot read --------------------------------------
    bad = os.path.join(tmp, "bad.csv")
    with open(bad, "w") as f:
        f.write("a,b\n1,2\nnot,numbers\n")
    rc, log = run(["bad.csv"], tmp)
    ok(rc != 0, "a table that is not numbers is refused, with a reason")

    print("\n%s" % ("FAILED" if fails else "all checks passed"))
    return fails


if __name__ == "__main__":
    sys.exit(main())
