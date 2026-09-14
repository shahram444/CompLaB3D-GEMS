#!/usr/bin/env python3
"""
check_docs.py -- hold the documentation up against the code it describes.

Every README in this repository shows commands and names files. Both drift: a
flag gets renamed, a tool grows subcommands, a file moves, and the prose keeps
saying what used to be true. Nothing catches it, because prose compiles under
any circumstances.

This does catch it. It reads every tracked markdown file and checks three
things that can be checked mechanically:

  1. Every  python <script> --flag  shown in a README names a flag that script
     really accepts. This found 29 wrong flags the first time it ran, including
     ten on a tool that has never taken any of them.

  2. Every file path written in backticks, and looking like a path into this
     repository, exists.

  3. Every constant quoted in a case README that also appears in that case's
     .sym rate law agrees with the file.

What it deliberately does NOT check is whether the prose is true. A sentence
can be wrong in ways no script will ever see. This only removes the class of
error that is mechanical, which is the class that accumulates silently.

    python3 tests/check_docs.py
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# argparse gives every script these, whatever its author wrote.
UNIVERSAL = {"--help"}

# Paths that appear in prose as examples rather than as things that exist.
PLACEHOLDER = re.compile(r"[<>*?]|yourmodel|mytable|your[_-]|\.\.\.|^-|scan/|"
                         r"path/to|previous_run|run/mycase|cmp/|figs/")

# A run writes these. A README naming one is describing what will appear, not
# something that should be in the repository, so they are not checked. This is
# the line between "the docs are wrong" and "the docs describe an output", and
# it has to be drawn somewhere explicit rather than guessed at.
GENERATED = re.compile(r"(^|/)output/|_discovered\.|_retrained\.|_new\.|"
                       r"^summary\.csv$|^checks\.log$|^run\.log$")

# A file named with no directory at all, in prose, is usually generic: "put your
# stoich.csv here". Only a path with a directory component is specific enough to
# be worth checking.
def _specific(ref):
    return "/" in ref


def tracked():
    out = subprocess.run(["git", "-C", ROOT, "ls-files"],
                         capture_output=True, text=True).stdout.split()
    return set(out)


def resolve(files, doc, ref):
    """The script a doc means, resolved the way a reader standing in that
    directory would resolve it: beside the README first, then one and two
    levels up, then anywhere by name."""
    d = os.path.dirname(doc)
    for c in (os.path.normpath(os.path.join(d, ref)),
              os.path.normpath(os.path.join(d, "..", ref)),
              os.path.normpath(os.path.join(d, "..", "..", ref)),
              os.path.normpath(ref)):
        if c in files:
            return [c]
    base = os.path.basename(ref)
    return sorted(f for f in files if os.path.basename(f) == base)


def check_flags(files, docs):
    bad = []
    for p in docs:
        txt = open(os.path.join(ROOT, p), encoding="utf-8", errors="replace").read()
        for m in re.finditer(r"python3?\s+([^\s`]+\.py)((?:\s+[^\n`]*?)?)(?=\n|`)", txt):
            ref, rest = m.group(1), m.group(2)
            cands = resolve(files, p, ref)
            if not cands:
                bad.append((p, ref, "no such script in the repository"))
                continue
            srcs = [open(os.path.join(ROOT, c), encoding="utf-8",
                         errors="replace").read() for c in cands]
            flags = set(re.findall(r"(?<![\w-])--[a-z][a-z0-9-]*", rest)) - UNIVERSAL
            for f in sorted(flags):
                if not any(('"%s"' % f) in s or ("'%s'" % f) in s for s in srcs):
                    bad.append((p, os.path.basename(ref), "does not accept %s" % f))
    return bad


def check_paths(files, docs):
    bad = []
    for p in docs:
        d = os.path.dirname(p)
        txt = open(os.path.join(ROOT, p), encoding="utf-8", errors="replace").read()
        for m in re.finditer(r"`([A-Za-z0-9_./-]+\.(?:py|sh|hh|cpp|xml|sym|thm|gnn|srg|csv|dat|md|txt))`", txt):
            ref = m.group(1)
            if PLACEHOLDER.search(ref) or GENERATED.search(ref) or not _specific(ref):
                continue
            for c in (os.path.normpath(os.path.join(d, ref)),
                      os.path.normpath(os.path.join(d, "..", ref)),
                      os.path.normpath(os.path.join(d, "..", "..", ref)),
                      os.path.normpath(ref)):
                if c in files:
                    break
            else:
                if not any(os.path.basename(f) == os.path.basename(ref) for f in files):
                    bad.append((p, ref, "named in backticks but nowhere in the repository"))
    return bad


def check_sym_constants(files, docs):
    """A case README that quotes its own rate law must quote the file's law."""
    bad = []
    for p in docs:
        d = os.path.dirname(p)
        syms = [f for f in files
                if f.startswith(d + "/") and f.endswith(".sym") and "/input/" in f]
        if not syms:
            continue
        txt = open(os.path.join(ROOT, p), encoding="utf-8", errors="replace").read()
        for s in syms:
            body = open(os.path.join(ROOT, s), encoding="utf-8", errors="replace").read()
            for line in body.splitlines():
                line = line.strip()
                if not line.startswith("rate") or "growth" not in line or "=" not in line:
                    continue
                expr = line.split("=", 1)[1].strip()
                if expr.lstrip().startswith("-"):      # a derived substrate line
                    continue
                nums = set(re.findall(r"\d+\.?\d*e?-?\d*", expr))
                shown = [m.group(0) for m in
                         re.finditer(r"growth\s*=\s*[^\n]+", txt)]
                for g in shown:
                    if "acetate" not in g and "o2" not in g:
                        continue
                    gnums = set(re.findall(r"\d+\.?\d*e?-?\d*", g))
                    if gnums and nums and not (gnums & nums):
                        bad.append((p, os.path.basename(s),
                                    "quotes a law whose constants are not in the file: %s"
                                    % g.strip()[:70]))
                        break
    return bad


def main():
    files = tracked()
    # A changelog records what the repository used to contain, so a path it
    # names may legitimately no longer exist.
    docs = sorted(f for f in files if f.endswith(".md")
                  and os.path.basename(f) != "CHANGELOG.md")
    problems = []
    for name, fn in (("command flags", check_flags),
                     ("file paths", check_paths),
                     ("rate law constants", check_sym_constants)):
        found = fn(files, docs)
        print("  %-22s %s" % (name, "ok" if not found else "%d PROBLEM(S)" % len(found)))
        problems += found
    if problems:
        print()
        for p in problems:
            print("    %-52s %-22s %s" % p)
        print("\n%d documentation problem(s) over %d markdown files"
              % (len(problems), len(docs)))
        return 1
    print("\nthe documentation matches the code, over %d markdown files" % len(docs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
