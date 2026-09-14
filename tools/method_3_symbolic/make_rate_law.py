#!/usr/bin/env python3
"""
make_rate_law.py -- find a rate law from your own measurements, by answering questions.

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
Meile Lab, University of Georgia.

WHAT THIS IS FOR

    fit_symbolic.py does the work, but it has fourteen options and you have to know which
    ones your case needs before you can start.  That is a fine tool for the second time
    and a poor one for the first.

    This asks instead.  Point it at a table of your own measurements and it reads the
    columns, shows you what it found, and asks a short list of questions in ordinary
    words.  At the end it writes a rate law file the solver can read, and prints the
    single command that would produce the same file again without the questions, so the
    result stays reproducible and you never have to answer them twice.

    Nothing is written until the last step, and it says what it is about to write first.

    python3 make_rate_law.py                     asks for everything, including the file
    python3 make_rate_law.py mydata.csv          starts from that table

WHAT IT NEEDS FROM YOU

    A table with one row per measurement.  One column is the thing you want a formula
    for, and the rest are what it depends on:

        acetate      o2           growth
        1.577e-05    3.147e-06    3.457e-08
        9.599e-05    1.611e-05    1.072e-06

    Where that table came from is your business: a laboratory experiment, a sweep of a
    slower model, or output from flux balance analysis.  This does not care.
"""

import os as _os

# Same reasoning as fit_symbolic.py: the thread pools have to be pinned before numpy
# loads, or the same seed does not give the same answer twice.  See the long note at the
# top of that file.
COMPLAB_THREADS = _os.environ.get("COMPLAB_THREADS", "1")
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    _os.environ[_var] = COMPLAB_THREADS

import subprocess
import sys

HERE = _os.path.dirname(_os.path.abspath(__file__))

try:
    import numpy as np
except ImportError:
    sys.exit("This needs numpy.  Install it with:  pip install numpy\n"
             "scipy is strongly recommended too:   pip install scipy")


# ------------------------------------------------------------------------------------------------
#  asking
#
#  Every question shows its answer in brackets and takes that answer when you press Enter,
#  so the whole thing can be run by holding Enter down to see what it would do.
# ------------------------------------------------------------------------------------------------
def ask(question, default=None, note=None):
    print("")
    if note:
        for line in note.split("\n"):
            print("      " + line)
    prompt = "  %s " % question
    if default is not None:
        prompt += "[%s] " % default
    while True:
        try:
            answer = input(prompt).strip()
        except (EOFError, KeyboardInterrupt):
            print("\n  stopped.  Nothing was written.")
            sys.exit(1)
        if answer:
            return answer
        if default is not None:
            return str(default)
        print("      an answer is needed here.")


def ask_yes_no(question, default=True, note=None):
    d = "yes" if default else "no"
    while True:
        a = ask(question + " (yes/no)", d, note).lower()
        if a in ("y", "yes"):
            return True
        if a in ("n", "no"):
            return False
        print("      please answer yes or no.")


def ask_number(question, default=None, note=None, positive=False):
    while True:
        a = ask(question, default, note)
        try:
            v = float(a)
        except ValueError:
            print("      that is not a number.")
            continue
        if positive and not v > 0:
            print("      that has to be greater than zero.")
            continue
        return v


def ask_choice(question, choices, default, note=None):
    """choices is a list of (key, description)."""
    if note:
        for line in note.split("\n"):
            print("      " + line)
    for k, d in choices:
        print("      %-10s %s" % (k, d))
    keys = [k for k, _ in choices]
    while True:
        a = ask(question, default).lower()
        if a in keys:
            return a
        print("      please answer one of: %s" % ", ".join(keys))


def heading(n, total, title):
    print("")
    print("  " + "-" * 68)
    print("  Step %d of %d.  %s" % (n, total, title))
    print("  " + "-" * 68)


# ------------------------------------------------------------------------------------------------
#  reading the table
# ------------------------------------------------------------------------------------------------
def read_table(path):
    """Return (column names, values as an array).  Accepts commas, semicolons or tabs."""
    with open(path) as f:
        first = f.readline()
    sep = max((",", ";", "\t"), key=first.count)
    header = [h.strip() for h in first.strip().split(sep)]
    if not header or any(h == "" for h in header):
        sys.exit("The first line of %s does not look like column names." % path)
    try:
        data = np.loadtxt(path, delimiter=sep, skiprows=1, ndmin=2)
    except ValueError as e:
        sys.exit("Could not read the numbers in %s.\n  %s\n"
                 "Every line after the first has to be numbers separated by '%s'."
                 % (path, e, sep))
    if data.shape[1] != len(header):
        sys.exit("%s has %d column names but %d columns of numbers."
                 % (path, len(header), data.shape[1]))
    return header, data


def looks_like_a_rate(name):
    """A guess at which column is the thing to be predicted, used only to pick a default."""
    n = name.lower()
    return any(w in n for w in ("growth", "rate", "extent", "flux", "mu", "velocity", "speed"))


def describe_columns(header, data):
    print("")
    print("      %-14s %-13s %-13s %s" % ("column", "smallest", "largest", "what it might be"))
    for i, h in enumerate(header):
        col = data[:, i]
        kind = "the thing to predict" if looks_like_a_rate(h) else "something it depends on"
        print("      %-14s %-13.4g %-13.4g %s" % (h, col.min(), col.max(), kind))
    print("")
    print("      %d rows." % data.shape[0])


# ------------------------------------------------------------------------------------------------
def main():
    print("")
    print("  ====================================================================")
    print("    Finding a rate law from your measurements")
    print("  ====================================================================")
    print("")
    print("    This asks you some questions and writes a rate law file that")
    print("    CompLaB3D can read.  Press Enter to accept anything in brackets.")
    print("    Nothing is written until the last step.")

    # ---------------------------------------------------------------------------- 1. the table
    heading(1, 7, "Your measurements")
    path = sys.argv[1] if len(sys.argv) > 1 else ask(
        "Which file holds your measurements?",
        "training/growth_samples.csv",
        "A table with one row per measurement and a name on top of each column.\n"
        "Commas, semicolons or tabs between the columns, all of them fine.")
    if not _os.path.isfile(path):
        sys.exit("  There is no file called %s here." % path)
    header, data = read_table(path)
    describe_columns(header, data)

    # ------------------------------------------------------------------- 2. what to find a law for
    heading(2, 7, "What do you want a formula for?")
    guess = next((h for h in header if looks_like_a_rate(h)), header[-1])
    target = ask("Which column is the thing you want a formula for?", guess,
                 "This is the one you cannot predict and want the search to explain.\n"
                 "For an organism it is usually the growth rate.  For a reaction with no\n"
                 "organism it is the rate of the reaction itself.")
    if target not in header:
        sys.exit("  There is no column called '%s'.  The file has: %s" % (target, ", ".join(header)))

    rest = [h for h in header if h != target]
    inputs = ask("Which columns does it depend on?", ",".join(rest),
                 "Separate them with commas.  Leave out anything that is a label, a\n"
                 "repeat, or a quantity the solver will not know at run time.")
    inputs = [s.strip() for s in inputs.replace(" ", ",").split(",") if s.strip()]
    for nm in inputs:
        if nm not in header:
            sys.exit("  There is no column called '%s'." % nm)
    if not inputs:
        sys.exit("  A formula needs at least one thing to depend on.")

    # ------------------------------------------------------------------------- 3. organism or not
    heading(3, 7, "Is an organism doing this?")
    biotic = ask_yes_no("Is this reaction carried out by an organism?", True,
                        "Yes, if something living eats the reactants and grows.  The rate\n"
                        "then depends on how much of that organism is present.\n"
                        "No, if it happens in the water on its own: a mineral forming or\n"
                        "dissolving, a redox reaction, sorption.")

    organism, yield_species, yield_value = "", "", 0.0
    if biotic:
        organism = ask("What is the organism called?", "Bug",
                       "Use the same name as <name_of_microbes> in your CompLaB.xml, or the\n"
                       "name you intend to give it there.  They have to match exactly.")

    # -------------------------------------------------------------------------- 4. the reaction
    heading(4, 7, "The reaction")
    print("      Now the chemistry, which is not fitted.  You already know it, so it is")
    print("      typed rather than guessed, and the solver works the rest out from it.")
    print("")
    print("      For each substance, say how many of it take part in one turn of the")
    print("      reaction.  Use a minus sign for something used up and a plus sign for")
    print("      something made.  Example, for  1 acetate + 2 oxygen -> biomass:")
    print("")
    print("            acetate   -1")
    print("            oxygen    -2")
    print("")

    stoich = []
    for nm in inputs:
        v = ask_number("How many %s per turn of the reaction?" % nm, "-1")
        if v == 0:
            print("      0 means it takes no part, so it is left out of the reaction.")
            continue
        stoich.append((nm, v))

    while ask_yes_no("Is anything else made or used that is not a column in your table?", False,
                     "A product, usually.  Nothing was measured about it, so it is not in\n"
                     "the table, but the solver still has to know it appears."):
        nm = ask("What is it called?", None,
                 "The same name as in <name_of_substrates> in your CompLaB.xml.")
        v = ask_number("How many %s per turn?" % nm, "1")
        if v != 0:
            stoich.append((nm, v))

    if not stoich:
        sys.exit("  A reaction with nothing in it cannot be written.")

    if biotic:
        names = [n for n, _ in stoich]
        yield_species = ask("Which substance is the yield measured against?",
                            names[0],
                            "The yield says how much organism is built per unit of one\n"
                            "substance eaten.  Say which substance that is.")
        if yield_species not in names:
            sys.exit("  '%s' is not in the reaction." % yield_species)
        yield_value = ask_number(
            "How much %s is made per %s used?" % (organism, yield_species),
            default="0.4",
            note="In the same units as the biomass in your CompLaB.xml.  0.4 means that\n"
                 "four tenths of what is eaten becomes organism and the rest is respired.",
            positive=True)

    # ----------------------------------------------------------------------------- 5. time units
    heading(5, 7, "Time")
    units = ask_choice("Are the rates in your table per second or per hour?",
                       [("per_second", "the rate columns are per second"),
                        ("per_hour", "the rate columns are per hour")],
                       "per_hour",
                       "Get this wrong and everything is out by a factor of 3600, and the\n"
                       "run will still finish and still look sensible.  Check the table.")

    # ---------------------------------------------------------------------------- 6. how long
    heading(6, 7, "How hard should it look?")
    effort = ask_choice("How thorough a search?",
                        [("quick",    "about a minute.  Enough to see whether it works at all"),
                         ("normal",   "a few minutes.  A reasonable answer"),
                         ("thorough", "ten minutes or more.  Use this for anything you will publish")],
                        "normal",
                        "The search tries a great many candidate formulas and keeps the ones\n"
                        "that match your table.  Trying more of them finds better formulas\n"
                        "and takes longer.  Nothing else changes.")
    pop, gens = {"quick": (150, 10), "normal": (400, 30), "thorough": (600, 60)}[effort]

    verify = ask_yes_no("Run it twice and check the two agree?", True,
                        "This doubles the time.  It is the only way to know the answer is\n"
                        "repeatable on your machine, and it is worth it before a number\n"
                        "goes into a paper.")

    seed = int(ask_number("Any whole number, to make the run repeatable", "1"))

    # ---------------------------------------------------------------------------- 7. where to write
    heading(7, 7, "Where to put the answer")
    # The default goes into input/ beside the case's other input files, and carries
    # "discovered" in its name, so that a first run cannot land on top of a rate law the case
    # already ships.  offline.sh makes the same choice for the same reason: a fitted law is
    # something to compare against the shipped one, not to silently replace it with.
    default_out = _os.path.join("input", "growth_discovered.sym" if biotic
                                         else "reaction_discovered.sym")
    while True:
        out = ask("What should the rate law file be called?", default_out,
                  "A plain text file you can read and keep.  Point CompLaB.xml at it with\n"
                  + ("<expressions_file>." if biotic else "<abiotic_file>."))
        d = _os.path.dirname(_os.path.abspath(out))
        if not _os.path.isdir(d):
            print("      there is no folder called %s." % _os.path.dirname(out))
            continue
        if _os.path.exists(out):
            # A rate law is a result somebody may have spent an afternoon choosing off a list.
            # Overwriting one without asking is the kind of small rudeness that loses work.
            print("")
            print("      %s already exists." % out)
            if not ask_yes_no("Overwrite it?", False,
                              "Answer no to give the new one a different name and keep both,\n"
                              "which is what you want if you are comparing two laws."):
                default_out = out
                continue
        break

    # ------------------------------------------------------------------- build the command and show it
    rate_name = "growth" if biotic else "extent"
    reaction = "  ".join("%s %g" % (n, v) for n, v in stoich)

    cmd = [sys.executable, _os.path.join(HERE, "fit_symbolic.py"),
           "--data", path, "--target", target, "--inputs", ",".join(inputs),
           "--reaction", reaction, "--units", units, "--rate-name", rate_name,
           "--pop", str(pop), "--gens", str(gens), "--depth", "6", "--seed", str(seed),
           "--out", out]
    if biotic:
        cmd += ["--yield", "%s %g" % (yield_species, yield_value), "--biomass", organism]
    if verify:
        cmd += ["--verify"]

    print("")
    print("  " + "=" * 68)
    print("  Here is what I have.  Nothing has been written yet.")
    print("  " + "=" * 68)
    print("")
    print("      measurements     %s, %d rows" % (path, data.shape[0]))
    print("      finding          a formula for %s" % target)
    print("      from             %s" % ", ".join(inputs))
    print("      kind             %s" % ("carried out by %s" % organism if biotic
                                         else "no organism, happens in the water"))
    print("      reaction         %s" % reaction)
    if biotic:
        print("      yield            %g %s per %s" % (yield_value, organism, yield_species))
    print("      rates are        %s" % units.replace("_", " "))
    print("      search           %s (%d candidates, %d rounds)" % (effort, pop, gens))
    print("      writes           %s" % out)
    print("")
    print("      To do exactly this again without the questions:")
    print("")
    quoted = " ".join(('"%s"' % c) if " " in c else c for c in cmd[1:])
    print("        python %s" % quoted)
    print("")

    if not ask_yes_no("Run the search now?", True):
        print("\n  Nothing was written.  The command above will do it when you are ready.")
        return 0

    print("")
    print("  " + "=" * 68)
    print("  Searching.  What follows comes from fit_symbolic.py.")
    print("  " + "=" * 68)
    rc = subprocess.call(cmd)
    if rc != 0:
        print("\n  The search stopped with an error.  Nothing usable was written.")
        return rc

    # ------------------------------------------------------------------------------ what to do next
    print("")
    print("  " + "=" * 68)
    print("  Done.  %s is written." % _os.path.abspath(out))
    print("  " + "=" * 68)
    print("")
    print("      READ THE LIST ABOVE BEFORE YOU TRUST THE FILE.  The search returns one")
    print("      formula for each length, and it picked one for you by a rule of thumb")
    print("      that is often not the one you want.")
    print("")
    print("      Read down the list and watch the error.  It falls a lot at first, then")
    print("      barely moves.  The last formula where it still fell by something worth")
    print("      having is the one to take.  Longer ones are describing the noise in your")
    print("      measurements rather than the chemistry, and they will be wrong at")
    print("      concentrations you never measured.")
    print("")
    print("      To take a different one, add --pick and the length from the list:")
    print("")
    print("        python %s --pick 13" % quoted)
    print("")
    print("      Then point CompLaB.xml at it:")
    print("")
    if biotic:
        print("        <symbolic>")
        print("            <enabled>true</enabled>")
        print("            <expressions_file>%s</expressions_file>" % out)
        print("        </symbolic>")
        print("")
        print("      and give the organism  <reaction_type>symbolic</reaction_type>.")
    else:
        print("        <symbolic>")
        print("            <enabled>true</enabled>")
        print("            <abiotic_file>%s</abiotic_file>" % out)
        print("        </symbolic>")
    print("")
    print("      The substrate lines are NOT in the file and do not need to be.  The")
    print("      solver works them out from the reaction you gave and prints them in the")
    print("      log, so their proportions cannot disagree with the chemistry.")
    print("")
    return 0


if __name__ == "__main__":
    sys.exit(main())
