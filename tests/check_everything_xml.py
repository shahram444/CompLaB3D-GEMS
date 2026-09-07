#!/usr/bin/env python3
"""Which XML tags the solver reads but config/CompLaB.everything.xml does not list.

That file claims to carry one line for every ability the solver has. A claim like
that decays the first time someone adds a tag and forgets to index it, and the
decay is invisible: a missing line looks exactly like a feature nobody uses.

So the claim is checked rather than trusted. This walks every source file for the
XML paths the solver actually opens -- doc["a"]["b"] and getOpt(doc, "a", "b") --
and prints, space separated, any that the everything file does not mention.
Silence means the index is complete.

Called by tests/check_repo.sh. Run it on its own while editing the parser:

    python3 tests/check_everything_xml.py

Two kinds of name are deliberately not required. Blocks built from a loop index
(substrate0, microbe1, phase2, species3, input4) are represented in the file by
one worked instance, not by all sixty-four. And the flat metabolic-model format
(<Metabolic_Model><S>...) is a data file the exporter writes, not a setting a
user edits.
"""
import os
import re
import sys

EVERYTHING = os.path.join("config", "CompLaB.everything.xml")
SRC = "src"

# Values and container names that appear in the same syntactic position as a tag
# but are not tags a user sets.
NOT_A_TAG = {
    "parameters", "true", "false",
    "Metabolic_Model", "S", "b", "c", "lb", "ub", "nmet", "nrxn", "objLoc",
}
INDEXED = re.compile(r"(substrate|microbe|phase|species|input)\d+$")


def listed_tags(path):
    text = open(path, encoding="utf8").read()
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    return set(re.findall(r"<([A-Za-z_][A-Za-z0-9_]*)>", text))


def tags_the_solver_reads(root):
    found = set()
    for dirpath, _, files in os.walk(root):
        for f in files:
            if not f.endswith((".hh", ".cpp")):
                continue
            s = open(os.path.join(dirpath, f), encoding="utf8", errors="replace").read()
            # doc["parameters"]["chemistry"]["number_of_substrates"]
            #
            # Every ["name"] in the file, not only the ones preceded by a
            # closing bracket. Requiring the bracket looks tighter and is
            # wrong: re.findall does not overlap, so in a chain the closing
            # ] of one key is consumed by its own match and the NEXT key is
            # never seen. That silently reduced this check to "the first key
            # of every chain", which is the one that is never the new tag.
            found |= set(re.findall(r'\[\s*"([A-Za-z_][A-Za-z0-9_]*)"\s*\]', s))
            # getOpt(doc, "parameters", "model_source", 0, tmp)
            for m in re.finditer(r"getOpt\s*\(\s*\w+\s*,(.*?)\)\s*;", s, re.S):
                found |= set(re.findall(r'"([A-Za-z_][A-Za-z0-9_]*)"', m.group(1)))
    return found


def main():
    if not os.path.exists(EVERYTHING):
        print("MISSING-FILE")
        return 1
    listed = listed_tags(EVERYTHING)
    read = tags_the_solver_reads(SRC)
    missing = sorted(
        t for t in read
        if t not in listed and t not in NOT_A_TAG and not INDEXED.match(t)
    )
    print(" ".join(missing))
    return 0


if __name__ == "__main__":
    sys.exit(main())
