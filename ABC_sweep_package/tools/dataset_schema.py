#!/usr/bin/env python3
"""The layout of the campaign database file, written down once.

Three collectors write this file: collect_complab_output.py for a campaign our
own builder produced, collect_foreign_complab.py for any other folder of
CompLaB runs, and collect_to_h5.py on the CompLaB side. Before this module
existed each of them decided the layout for itself, and they drifted: one
stored the input file as raw text, another parsed it into tables, and a reader
had to cope with both. The layout is now stated here, in one place, as data.

    python3 dataset_schema.py check   dataset.h5
    python3 dataset_schema.py upgrade dataset.h5 --out dataset_v1.h5

`check` reports every entry the document asks for and says whether the file has
it. `upgrade` writes a new file that does, without re-running the campaign.
Collectors call `finalise()` so that a file is right the moment it is written.

WHAT IS AND IS NOT INVENTED

An entry that can be derived from what the file already holds is derived. An
entry that cannot is left absent and the absence is recorded in writing, in
ancillary/absent. Nothing is ever filled with zeros to make the layout look
complete, because "this field was absent" and "this field was measured as zero"
would then be the same numbers, and the second is a result while the first is
a gap.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

try:
    import h5py
except ImportError:                                       # pragma: no cover
    sys.exit("h5py is required: pip install h5py")


SCHEMA = "complab-dataset"
SCHEMA_VERSION = "1.0"


# ===========================================================================
#  the document, as data
#
#  Each row is (name, required, meaning). `required` is True for an entry every
#  campaign must have, False for one that exists only when the campaign used
#  the thing it describes, such as a mineral or an evolving pore space.
# ===========================================================================

ROOT_LABELS = [
    ("species",           True,  "the field names, in order"),
    ("species_role",      True,  "dissolved, mineral or microbe"),
    ("reactions",         True,  "the name of each reaction channel"),
    ("param_names",       True,  "what each condition is"),
    ("param_units",       True,  "the unit of each one"),
    ("shape",             True,  "the grid, always three numbers"),
    ("spacing",           True,  "dx, dy, dz"),
    ("spacing_unit",      True,  "m, mm or um"),
    ("mode",              True,  "steady or transient"),
    ("n_times",           True,  "how many snapshots per run"),
    ("n_samples",         True,  "how many runs"),
    ("n_geometries",      True,  "how many rocks"),
    ("structure_evolves", True,  "did the pore space change"),
]

GEOM = [
    ("gid",      True,  "identity of each rock"),
    ("material", True,  "0 solid, 1 wall, 2 pore"),
    ("gdf",      True,  "distance from the inlet, through pore space"),
    ("edt",      True,  "distance to the nearest solid, in a straight line"),
    ("biofilm",  False, "where biofilm was at t=0. Only when the campaign "
                        "carried a microbe"),
    ("mineral",  False, "the solid phase chemical. 3D campaigns only"),
]

SAMPLES = [
    ("geom_index", True,  "which rock each run used"),
    ("run_id",     True,  "identity of each run"),
    ("params",     True,  "the conditions of each run"),
    ("t_norm",     True,  "snapshot times, 0 to 1"),
    ("t_seconds",  True,  "the same times in seconds"),
    ("conc",       True,  "every field"),
    ("conc_scale", True,  "one divisor per field"),
    ("velocity",   True,  "the flow field"),
    ("rate",       False, "how fast each reaction ran. Only when the runs "
                          "were told to write it"),
    ("rate_scale", False, "one divisor per rate channel"),
    ("material",   False, "only when structure_evolves"),
    ("biofilm",    False, "only when structure_evolves"),
    ("mineral",    False, "only when structure_evolves"),
    ("gdf",        False, "only when structure_evolves"),
]

GRID = [
    ("boundary_conditions", True, "per field, inlet and outlet"),
]

INPUTS = [
    ("xml",      True, "every setting in the input file"),
    ("kinetics", True, "every rate constant"),
    ("order",    True, "the chemical order the kinetics file expects"),
    ("files",    True, "which file each run got, or where it was looked for"),
]

# The conditions the document lists, in its order. A campaign records the ones
# it has. A Damkohler number belongs to a reaction, so a campaign that ran only
# a purely chemical reaction has da_abio and no da_bio, and one that ran only a
# microbial reaction has the opposite. The generic name `da` does not say which
# it is, so it is made specific from the da_column label wherever that is the
# only Damkohler number in the file.
CONDITIONS = [
    ("pe",          "Peclet number, flow against diffusion"),
    ("da_bio",      "Damkohler number for the microbial reaction"),
    ("da_abio",     "Damkohler number for the chemical reaction"),
    ("ks_ac_norm",  "half saturation for the donor, over its inlet "
                    "concentration"),
    ("ks_a_norm",   "half saturation for the acceptor, over its inlet "
                    "concentration"),
    ("y_norm",      "yield, times the donor inlet concentration, over the "
                    "initial biomass"),
]

ANCILLARY = [
    ("run_note",     True, "comments about a run"),
    ("run_name",     True, "the folder it came from"),
    ("dataset_note", True, "comments about the campaign"),
]


def _vs():
    return h5py.string_dtype("utf-8")


# The four input tables, with exactly the columns the document names.
def _dtypes():
    vs = _vs()
    return {
        "xml":      np.dtype([("name", vs), ("value", vs), ("number", "f8"),
                              ("unit", vs), ("note", vs)]),
        "kinetics": np.dtype([("name", vs), ("value", "f8"), ("unit", vs),
                              ("note", vs), ("file", vs)]),
        "order":    np.dtype([("name", vs), ("position", "i4"), ("file", vs)]),
        "files":    np.dtype([("what", vs), ("found", vs), ("where", vs)]),
    }


# ===========================================================================
#  reading what a file already has
# ===========================================================================

def _s(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.ndarray):
        if v.dtype.kind == "S":
            return [x.decode("utf-8", "replace") for x in v.ravel().tolist()]
        return v.tolist()
    return v


def _num_or_nan(text):
    try:
        return float(str(text).strip())
    except Exception:
        return float("nan")


def _xml_rows(text):
    """(path, value, comment) for every leaf setting in one input file.

    ElementTree drops comments and the comments in these files are where a
    person says why a number is that number, so the comment that sits on the
    line before a setting is recovered from the raw text and attached to it.
    """
    import re
    rows, comments = [], {}
    for m in re.finditer(r"<!--(.*?)-->\s*<([A-Za-z_][\w.\-]*)>", text, re.S):
        comments[m.group(2)] = " ".join(m.group(1).split())
    try:
        import xml.etree.ElementTree as ET
        root = ET.fromstring(text)
    except Exception:
        return rows

    def walk(node, prefix):
        kids = list(node)
        if not kids:
            val = (node.text or "").strip()
            rows.append(("%s/%s" % (prefix, node.tag) if prefix else node.tag,
                         val, comments.get(node.tag, "")))
            return
        for k in kids:
            walk(k, "%s/%s" % (prefix, node.tag) if prefix else node.tag)

    for k in list(root):
        walk(k, "")
    return rows


def _kinetics_constants(text):
    """(name, default, note, env_name) for every rate constant in a header.

    Two ways these headers state a constant, and both have to be caught:

        double k = 1.2e-3;                       compiled in, same every run
        q.k_abio = envd("PRT_KABIO", 1.0);       read from the environment

    The second is the one that matters, because it is how a sweep varies the
    rate constant from run to run: the number in the header is only the
    fallback, and the value a run actually used is in its own environment
    record. env_name is empty for the first form and the variable name for the
    second, which is what lets the table resolve it per run.
    """
    import re
    out, seen = [], set()
    env = re.compile(
        r"([A-Za-z_]\w*)\s*=\s*envd\s*\(\s*\"([A-Za-z_]\w*)\"\s*,"
        r"\s*([-+0-9.eE]+)\s*\)\s*;\s*(?://(.*))?")
    for m in env.finditer(text or ""):
        try:
            v = float(m.group(3))
        except ValueError:
            continue
        if m.group(1) in seen:
            continue
        seen.add(m.group(1))
        out.append((m.group(1), v, " ".join((m.group(4) or "").split()),
                    m.group(2)))
    # A running total declared `static double x = 0.0;` inside the kinetics
    # function is a counter the solver accumulates into, not a constant the
    # chemistry was given. Keeping it would put a value of zero in the rate
    # constant table and read as a rate constant that was set to zero.
    ACCUMULATOR = re.compile(r"^(iter_|n_)|_(total|count|sum|calls)$")
    plain = re.compile(
        r"^\s*(?:static\s+)?(?:const\s+)?(?:double|float)\s+"
        r"([A-Za-z_]\w*)\s*=\s*([-+0-9.eE]+)\s*;\s*(?://(.*))?$", re.M)
    for m in plain.finditer(text or ""):
        try:
            v = float(m.group(2))
        except ValueError:
            continue
        if m.group(1) in seen or ACCUMULATOR.search(m.group(1)):
            continue
        seen.add(m.group(1))
        out.append((m.group(1), v, " ".join((m.group(3) or "").split()), ""))
    return out


def _env_values(text):
    """{VARIABLE: number} from a run's environment record."""
    import re
    out = {}
    for m in re.finditer(r"^\s*(?:export\s+)?([A-Z][A-Z0-9_]*)\s*=\s*"
                         r"([-+0-9.eE]+)\s*$", text or "", re.M):
        try:
            out[m.group(1)] = float(m.group(2))
        except ValueError:
            pass
    return out


def _kinetics_order(text):
    """(name, position) for the chemical order a kinetics header expects.

    Written in these headers as a comment mapping each slot of the
    concentration vector onto a chemical, for example  C[0]=A  C[1]=B.
    """
    import re
    out = []
    for m in re.finditer(r"\b[A-Za-z]+\s*\[\s*(\d+)\s*\]\s*=\s*([A-Za-z]\w*)",
                         text or ""):
        name, pos = m.group(2), int(m.group(1))
        if (name, pos) not in out:
            out.append((name, pos))
    return out


# ===========================================================================
#  checking a file against the document
# ===========================================================================

RESULTS_LAYOUT = [
    ("design",      "the campaign: the pore spaces, the run table, the scales"),
    ("models",      "one group per fitted model: its settings, its split, its "
                    "fitting history and its weights"),
    ("results",     "the error of every model on train, val and test, pooled "
                    "and broken down, with the held out volumes"),
    ("comparisons", "the geometry ablation and the speed comparison"),
    ("notes",       "the written notes that travel with the file"),
]


def file_kind(path):
    """'dataset', 'results' or 'unknown'.

    Two different files come out of this project and they are not two versions
    of one layout: a campaign dataset holds simulations, a results file holds
    what the fitted models did with them. Checking one against the other's
    layout reports everything as missing, which is alarming and wrong, so the
    kind is settled before anything is checked.
    """
    with h5py.File(path, "r") as h:
        if str(_s(h.attrs.get("schema", ""))) == "prt-results":
            return "results"
        if "samples/conc" in h and "geom/material" in h:
            return "dataset"
        if "results" in h and "models" in h:
            return "results"
    return "unknown"


def validate(path):
    """Return (rows, problems). One row per documented entry.

    Only the campaign dataset is described here. A results file is reported
    on its own terms rather than against a layout that was never meant for it.
    """
    if file_kind(path) == "results":
        rows, problems = [], []
        with h5py.File(path, "r") as h:
            for name, meaning in RESULTS_LAYOUT:
                present = name in h
                rows.append(("results file", name, "yes" if present else "no",
                             "must", meaning, ""))
                if not present:
                    problems.append("%s is missing" % name)
        return rows, problems

    rows, problems = [], []
    with h5py.File(path, "r") as h:
        a = h.attrs
        evolves = bool(a.get("structure_evolves", False))
        has_rate = "samples/rate" in h

        def row(group, name, required, present, meaning, note=""):
            rows.append((group, name, "yes" if present else "no",
                         "must" if required else "when used", meaning, note))
            if required and not present:
                problems.append("%s/%s is missing" % (group or "root", name))

        for n, req, mean in ROOT_LABELS:
            row("root label", n, req, n in a, mean)

        for n, req, mean in GEOM:
            need = req or (n == "biofilm" and _has_microbe(h)) \
                or (n == "mineral" and _has_mineral(h))
            row("geom", n, need, ("geom/%s" % n) in h, mean)

        for n, req, mean in SAMPLES:
            need = req
            if n in ("rate", "rate_scale"):
                need = has_rate
            if n in ("material", "biofilm", "mineral", "gdf"):
                need = evolves
            row("samples", n, need, ("samples/%s" % n) in h, mean)

        for n, req, mean in GRID:
            row("grid", n, req, ("grid/%s" % n) in h, mean)

        dt = _dtypes()
        for n, req, mean in INPUTS:
            p = "inputs/%s" % n
            ok = p in h
            note = ""
            if ok:
                got = h[p].dtype
                want = dt[n]
                if got.names is None:
                    ok, note = False, "stored as raw text, not as a table"
                elif list(got.names) != list(want.names):
                    note = "columns are %s" % (", ".join(got.names),)
            row("inputs", n, req, ok, mean, note)

        for n, req, mean in ANCILLARY:
            p = "ancillary/%s" % n
            present = p in h or n in h.get("ancillary", h).attrs
            row("ancillary", n, req, present, mean)

        # velocity has a time axis in the document
        if "samples/velocity" in h:
            nd = h["samples/velocity"].ndim
            ok = nd == 6
            rows.append(("shape", "velocity", "yes" if ok else "no", "must",
                         "(S,Tv,3,nx,ny,nz)",
                         "" if ok else "is %dD, missing the Tv axis" % nd))
            if not ok:
                problems.append("samples/velocity has no Tv axis")
    return rows, problems


def _has_microbe(h):
    roles = _s(h.attrs.get("species_role", np.array([], "S1")))
    return any("microbe" in str(r) for r in (roles or []))


def _has_mineral(h):
    roles = _s(h.attrs.get("species_role", np.array([], "S1")))
    return any("mineral" in str(r) for r in (roles or []))


def report(path, out=sys.stdout):
    kind = file_kind(path)
    rows, problems = validate(path)
    if kind == "results":
        print("%s\n%s\n" % (path, "=" * len(path)), file=out)
        print("This is a results file, not a campaign dataset. It has its own "
              "layout,\nand the campaign layout does not apply to it.\n",
              file=out)
        for g, n, has, _w, mean, _nt in rows:
            print("   %-14s %-4s %s" % (n, has, mean), file=out)
        print("", file=out)
        print("This file matches the results layout." if not problems
              else "missing: %s" % ", ".join(problems), file=out)
        return problems
    w = [max(len(str(r[i])) for r in rows + [("GROUP", "ENTRY", "HAS", "WHEN",
                                              "MEANING", "NOTE")])
         for i in range(5)]
    print("%s\n%s\n" % (path, "=" * len(path)), file=out)
    print("%-*s  %-*s  %-*s  %-*s  %s" % (w[0], "GROUP", w[1], "ENTRY",
                                          3, "HAS", 9, "WHEN", "MEANING"),
          file=out)
    print("-" * (w[0] + w[1] + 40), file=out)
    for g, n, has, when, mean, note in rows:
        print("%-*s  %-*s  %-3s  %-9s  %s%s"
              % (w[0], g, w[1], n, has, when, mean,
                 ("   <- " + note) if note else ""), file=out)
    print("", file=out)
    if problems:
        print("%d entr%s the document asks for and this file does not have:"
              % (len(problems), "y" if len(problems) == 1 else "ies"), file=out)
        for p in problems:
            print("   %s" % p, file=out)
        print("\nRun  dataset_schema.py upgrade %s --out fixed.h5" % path,
              file=out)
    else:
        print("This file matches the document.", file=out)
    return problems


# ===========================================================================
#  bringing a file up to the document
# ===========================================================================

def finalise(h, xml_texts=None, kin_texts=None, kin_names=None,
             file_records=None, reactions=None, boundaries=None,
             run_notes=None, dataset_note=None, verbose=True):
    """Fill in everything the document asks for that can be derived.

    Collectors call this just before closing the file. Everything it writes is
    either derived from what is already there or passed in by the caller.
    Anything it cannot supply is recorded in ancillary/absent rather than
    invented.
    """
    if str(_s(h.attrs.get("schema", ""))) == "prt-results":
        if verbose:
            print("   schema: this is a results file, which has its own "
                  "layout. Nothing to do.")
        return [], []

    said = []
    absent = []
    vs = _vs()
    dt = _dtypes()

    def note(msg):
        said.append(msg)
        if verbose:
            print("   schema: %s" % msg)

    h.attrs["schema"] = SCHEMA
    h.attrs["schema_version"] = SCHEMA_VERSION

    S = int(h.attrs.get("n_samples", 0)) or (
        len(h["samples/run_id"]) if "samples/run_id" in h else 0)

    # ---- root labels ------------------------------------------------------
    if "reactions" not in h.attrs:
        species = [_dec(x) for x in np.atleast_1d(h.attrs["species"])] \
            if "species" in h.attrs else []
        n_rate = None
        if "samples/rate" in h and h["samples/rate"].ndim > 2:
            n_rate = int(h["samples/rate"].shape[2])

        if reactions:
            names = list(reactions)
        elif n_rate is not None and species and n_rate == len(species):
            # The label has to describe the array it sits beside. This
            # campaign's rate array has one channel per CHEMICAL, not one per
            # reaction: with a single reaction the solver wrote each
            # chemical's own net rate of change rather than the reaction
            # extent. Naming the channels after the reaction would say there
            # are three reactions, which there are not.
            names = ["net rate of %s" % x for x in species]
            note("rate carries one channel per chemical, not per reaction, "
                 "so the channels are named for the chemicals")
        elif "reaction" in h.attrs:
            names = [_dec(h.attrs["reaction"])]
        else:
            names = []

        if names:
            h.attrs["reactions"] = np.array([n.encode() for n in names])
            note("reactions label written, %d channel(s)" % len(names))
            if "reaction" in h.attrs:
                h.attrs["chemistry"] = h.attrs["reaction"]
        else:
            absent.append("reactions: the campaign did not say what its "
                          "reaction channels are called")

    if "mode" not in h.attrs and "n_times" in h.attrs:
        h.attrs["mode"] = "steady" if int(h.attrs["n_times"]) == 1 \
            else "transient"
        note("mode derived from n_times")

    if "structure_evolves" not in h.attrs:
        h.attrs["structure_evolves"] = False
        note("structure_evolves written as false")

    # ---- which conditions this campaign actually has ---------------------
    if "param_names" in h.attrs:
        names = [_dec(x) for x in np.atleast_1d(h.attrs["param_names"])]
        da = [i for i, n in enumerate(names) if n in ("da", "Da")]
        if len(da) == 1 and "da_column" in h.attrs:
            # The file already records which Damkohler number this is. Saying
            # it in the name as well means a reader never has to find the
            # other label to know what the column means.
            specific = _dec(h.attrs["da_column"])
            if specific in [c for c, _ in CONDITIONS]:
                names[da[0]] = specific
                h.attrs["param_names"] = np.array(
                    [n.encode() for n in names])
                note("condition %d renamed from da to %s, from the da_column "
                     "label" % (da[0], specific))
        known = [c for c, _ in CONDITIONS]
        missing = [(c, w) for c, w in CONDITIONS if c not in names]
        if missing and all(n in known for n in names):
            for c, w in missing:
                absent.append("condition %s (%s): this campaign did not vary "
                              "it, and it is left out rather than stored as "
                              "zero, because a zero would read as a value "
                              "that was set rather than one that does not "
                              "apply" % (c, w))

    if "param_units" not in h.attrs and "param_names" in h.attrs:
        n = len(_s(h.attrs["param_names"]))
        h.attrs["param_units"] = np.array([b"1"] * n)
        note("param_units written as dimensionless")

    # ---- scales, as datasets ---------------------------------------------
    # The document asks for these as arrays beside the fields they rescale,
    # not as labels on the group, so that a reader finds them the same way it
    # finds everything else.
    if "samples" in h:
        for key, src in (("conc_scale", "conc_scale"),
                         ("rate_scale", "rate_scale")):
            if key in h["samples"]:
                continue
            if src in h["samples"].attrs:
                v = np.asarray(h["samples"].attrs[src], np.float64)
                h["samples"].create_dataset(key, data=v)
                h["samples"][key].attrs["how_to_read_this"] = (
                    b"Multiply the stored field by this to get it back in its "
                    b"own unit. One entry per channel.")
                note("samples/%s written as a dataset" % key)

    # ---- grid/boundary_conditions ----------------------------------------
    if "grid/boundary_conditions" not in h:
        bc = boundaries
        if bc is None and "boundaries" in h.attrs:
            bc = str(_s(h.attrs["boundaries"]))
        if bc:
            g = h.require_group("grid")
            species = [str(x) for x in (_s(h.attrs.get("species", [])) or [])]
            BC = np.dtype([("field", vs), ("inlet", vs), ("outlet", vs)])
            arr = np.empty(max(len(species), 1), BC)
            arr[...] = ("", "", "")
            # One description covers every field unless the campaign recorded
            # them separately, which the ABC sweep did not.
            for i, nm in enumerate(species):
                arr[i] = (nm, str(bc), str(bc))
            d = g.create_dataset("boundary_conditions", data=arr)
            d.attrs["how_to_read_this"] = (
                b"One row per field, with what was imposed at the inlet and "
                b"at the outlet. Where the campaign recorded one description "
                b"for the whole domain rather than one per field, that same "
                b"description appears on every row.")
            note("grid/boundary_conditions written from the boundaries label")
        else:
            absent.append("grid/boundary_conditions: the campaign did not "
                          "record what was imposed at the faces")

    # ---- velocity gets its time axis --------------------------------------
    if "samples/velocity" in h and h["samples/velocity"].ndim == 5:
        old = h["samples/velocity"]
        data = old[...]
        at = dict(old.attrs)
        del h["samples/velocity"]
        d = h["samples"].create_dataset(
            "velocity", data=data[:, None], dtype=data.dtype,
            compression="gzip", compression_opts=4)
        for k, v in at.items():
            d.attrs[k] = v
        d.attrs["how_to_read_this"] = (
            b"(run, flow snapshot, component, grid). A rock that never "
            b"changes shape needs one flow field for the whole run, so the "
            b"second axis has length 1. It becomes as long as the snapshot "
            b"axis only when structure_evolves is true.")
        note("samples/velocity reshaped to carry its Tv axis")

    # ---- the four input tables -------------------------------------------
    g = h.require_group("inputs")

    if xml_texts is not None and _needs_table(h, "inputs/xml", dt["xml"]):
        _write_xml_table(h, g, xml_texts, dt["xml"], note)
    elif "inputs/xml" in h and _needs_table(h, "inputs/xml", dt["xml"]):
        raw = [_dec(x) for x in h["inputs/xml"][:]]
        _write_xml_table(h, g, raw, dt["xml"], note, rename_raw=True)
    elif "inputs/xml" not in h:
        absent.append("inputs/xml: no input file was found for any run")

    if kin_texts is not None and _needs_table(h, "inputs/kinetics",
                                              dt["kinetics"]):
        _write_kin_tables(h, g, kin_texts, kin_names, dt, S, note,
                          env_texts=_env_texts(h))
    elif "inputs/kinetics" in h and _needs_table(h, "inputs/kinetics",
                                                 dt["kinetics"]):
        raw = [_dec(x) for x in h["inputs/kinetics"][:]]
        nm = ([_dec(x) for x in h["inputs/kinetics_name"][:]]
              if "inputs/kinetics_name" in h else None)
        _write_kin_tables(h, g, raw, nm, dt, S, note, rename_raw=True,
                          env_texts=_env_texts(h))
    elif "inputs/kinetics" not in h:
        absent.append("inputs/kinetics: no kinetics file was found")

    if "inputs/files" not in h:
        FILES = dt["files"]
        if file_records:
            arr = np.array(file_records, FILES)
        else:
            what = []
            for label, key in (("the input file", "inputs/xml_raw"),
                               ("the kinetics file", "inputs/kinetics_raw")):
                what.append((label, "yes" if key in h else "no",
                             "stored in this file" if key in h
                             else "not recorded by this collector"))
            arr = np.empty((max(S, 1), len(what)), FILES)
            for k in range(max(S, 1)):
                for i, row in enumerate(what):
                    arr[k, i] = row
        d = g.create_dataset("files", data=arr)
        d.attrs["how_to_read_this"] = (
            b"CompLaB reads its input file from the folder it is launched in "
            b"and never copies it next to the results, so output can arrive "
            b"with nothing describing it. A no here is that, not a fault, and "
            b"the where column says exactly where it was looked for.")
        note("inputs/files written")

    # ---- ancillary --------------------------------------------------------
    a = h.require_group("ancillary")
    if "run_note" not in a:
        a.create_dataset("run_note",
                         data=np.array(list(run_notes) if run_notes
                                       else [""] * max(S, 1), dtype=object),
                         dtype=vs)
        a["run_note"].attrs["how_to_read_this"] = (
            b"Free text about one run. Empty where nobody wrote any.")
        note("ancillary/run_note written")
    if "dataset_note" not in a:
        txt = dataset_note if dataset_note is not None else \
            str(_s(a.attrs.get("dataset_note", "")))
        a.create_dataset("dataset_note", data=np.array([txt], dtype=object),
                         dtype=vs)
        note("ancillary/dataset_note written")

    # ---- what could not be supplied ---------------------------------------
    if absent:
        if "absent" in a:
            del a["absent"]
        a.create_dataset("absent", data=np.array(absent, dtype=object),
                         dtype=vs)
        a["absent"].attrs["how_to_read_this"] = (
            b"Entries the layout describes that this campaign genuinely does "
            b"not have, each with the reason. Nothing here was filled with "
            b"zeros, because a measured zero and a missing value must stay "
            b"distinguishable.")
        if verbose:
            for x in absent:
                print("   schema: absent, %s" % x)
    return said, absent


def _dec(x):
    """One element of a variable length string array, as text.

    h5py hands these back as `object` arrays of bytes rather than as a fixed
    width "S" dtype, so str() on one would give the b'...' repr rather than
    the string it holds.
    """
    if isinstance(x, bytes):
        return x.decode("utf-8", "replace")
    return str(x)


def _label(note, env_name):
    """The note column: the comment, and where the value came from."""
    if env_name and note:
        return "%s. Set per run by %s" % (note, env_name)
    if env_name:
        return "set per run by %s" % env_name
    return note


def _env_texts(h):
    """Each run's environment record, where the collector stored one.

    This is where the value a run actually used lives, for every constant the
    kinetics header reads through getenv rather than hard coding.
    """
    for key in ("inputs/env", "inputs/environment"):
        if key in h:
            return [_dec(x) for x in h[key][:]]
    return []


def _needs_table(h, path, want):
    if path not in h:
        return True
    got = h[path].dtype
    return got.names is None or list(got.names) != list(want.names)


def _write_xml_table(h, g, texts, SETTING, note, rename_raw=False):
    """Part 5 of the document: every setting, with its value as a number."""
    try:
        from collect_foreign_complab import describe_setting
    except Exception:
        def describe_setting(path):
            return path.rsplit("/", 1)[-1].replace("_", " "), "", ""

    S = len(texts)
    order, seen = [], {}
    parsed = []
    for t in texts:
        rows = _xml_rows(t or "")
        parsed.append(rows)
        for p, _, _ in rows:
            if p not in seen:
                seen[p] = len(order)
                order.append(p)
    P = max(len(order), 1)
    arr = np.empty((max(S, 1), P), SETTING)
    arr[...] = ("", "", np.nan, "", "")
    for i, p in enumerate(order):
        nm, unit, _mean = describe_setting(p)
        for k in range(max(S, 1)):
            arr[k, i] = (nm, "", np.nan, unit, "")
    for k, rows in enumerate(parsed):
        for p, v, comment in rows:
            j = seen[p]
            nm, unit, mean = describe_setting(p)
            n = comment or mean
            if len(n) > 60:
                n = n[:57].rstrip(" ,.;") + "..."
            arr[k, j] = (nm, v, _num_or_nan(v), unit, n)

    if rename_raw and "xml" in g and "xml_raw" not in g:
        g.move("xml", "xml_raw")
        g["xml_raw"].attrs["how_to_read_this"] = (
            b"The input file of each run, exactly as it was written. The "
            b"table in inputs/xml is this parsed into columns.")
    elif "xml" in g:
        del g["xml"]
    d = g.create_dataset("xml", data=arr)
    d.attrs["how_to_read_this"] = (
        b"One row per setting, one table per run. The number column is a real "
        b"number wherever the setting is one, so it can be sorted, differenced "
        b"between two runs, or fed straight into anything else without "
        b"parsing the value column. Not-a-number means the setting is not a "
        b"number, or that run had no input file.")
    note("inputs/xml written as a table, %d settings x %d runs" % (len(order), S))


def _write_kin_tables(h, g, texts, names, dt, S, note, rename_raw=False,
                      env_texts=None):
    """Part 5: every rate constant, and the chemical order it expects.

    The header is compiled into the executable, so it is the same for every
    run. The VALUES are not: a sweep varies the rate constant by setting an
    environment variable per run, and the number in the header is only the
    fallback used when nothing set it. So the table carries one row per run,
    and where a constant is environment backed its value is read from that
    run's own environment record. A constant that is not, or a run with no
    record, keeps the header's own number.
    """
    CONSTANT, ORDER = dt["kinetics"], dt["order"]
    names = [_dec(x) for x in (names if names else
                               ["kinetics" for _ in texts])]
    envs = [_env_values(t) for t in (env_texts or [])]

    keys, kseen = [], {}
    per_file = []
    for t, fam in zip(texts, names):
        cs = _kinetics_constants(t)
        per_file.append((fam, cs))
        for nm, _v, nt, ev in cs:
            if (fam, nm) not in kseen:
                kseen[(fam, nm)] = len(keys)
                keys.append((nm, fam, nt, ev))
    K = max(len(keys), 1)
    arr = np.empty((max(S, 1), K), CONSTANT)
    arr[...] = ("", np.nan, "", "", "")
    n_from_env = 0
    for i, (nm, fam, nt, ev) in enumerate(keys):
        label = _label(nt, ev)
        for k in range(max(S, 1)):
            arr[k, i] = (nm, np.nan, "", label[:60], fam)
    for fam, cs in per_file:
        for nm, v, nt, ev in cs:
            i = kseen[(fam, nm)]
            label = _label(nt, ev)
            for k in range(max(S, 1)):
                val = v
                if ev and k < len(envs) and ev in envs[k]:
                    val = envs[k][ev]
                    n_from_env += 1
                arr[k, i] = (nm, val, "", label[:60], fam)
    if n_from_env:
        note("%d rate constant values read from the per run environment, "
             "not from the header default" % n_from_env)

    if rename_raw and "kinetics" in g and "kinetics_raw" not in g:
        g.move("kinetics", "kinetics_raw")
    elif "kinetics" in g:
        del g["kinetics"]
    d = g.create_dataset("kinetics", data=arr)
    d.attrs["how_to_read_this"] = (
        b"The rate constants, with the note written beside each one in the "
        b"source. These live in the .hh files and never in the input file, "
        b"which is why a result cannot be reproduced without them.")
    note("inputs/kinetics written as a table, %d constants" % len(keys))

    rows = []
    for t, fam in zip(texts, names):
        fam = _dec(fam)
        for nm, pos in _kinetics_order(t):
            if (fam, nm) not in [(f, n) for n, _p, f in rows]:
                rows.append((nm, pos, fam))
    if "order" in g:
        del g["order"]
    orr = np.empty(max(len(rows), 1), ORDER)
    orr[...] = ("", -1, "")
    for i, r in enumerate(rows):
        orr[i] = r
    d = g.create_dataset("order", data=orr)
    d.attrs["how_to_read_this"] = (
        b"The slot each chemical occupies in the concentration vector the "
        b"kinetics file was written against. Comparing this with the species "
        b"label is the one check the solver will not do for you: get it wrong "
        b"and every rate is computed on the wrong chemical, with nothing "
        b"anywhere complaining.")
    note("inputs/order written, %d entries" % len(rows))


# ===========================================================================
#  upgrading a file that already exists
# ===========================================================================

def upgrade(src, dst, verbose=True):
    """Copy a file and bring the copy up to the document."""
    if os.path.abspath(src) == os.path.abspath(dst):
        raise SystemExit("write the upgraded file somewhere else, so the "
                         "original is still there if anything goes wrong")
    if verbose:
        print("copying %s -> %s" % (src, dst))
    with h5py.File(src, "r") as a, h5py.File(dst, "w") as b:
        for k, v in a.attrs.items():
            b.attrs[k] = v
        for name in a:
            a.copy(name, b, name=name)
        if verbose:
            print("bringing it up to the layout:")
        said, absent = finalise(b, verbose=verbose)
    if verbose:
        print("\n%d change%s, %d entr%s genuinely absent"
              % (len(said), "" if len(said) == 1 else "s",
                 len(absent), "y" if len(absent) == 1 else "ies"))
    return dst


def fix(path, verbose=True):
    """Run the layout pass on a file in place.

    Everything finalise writes is small: labels, four tables and a couple of
    arrays. The one exception is reshaping velocity, which rewrites a large
    array, and HDF5 does not give back the space the old one occupied. So this
    refuses to do that one and sends you to upgrade instead, which writes a
    fresh file and leaves no hole in it.
    """
    with h5py.File(path, "r") as h:
        if "samples/velocity" in h and h["samples/velocity"].ndim == 5:
            raise SystemExit(
                "this file still needs its velocity reshaped, which would "
                "leave a hole in it.\n"
                "Use:  dataset_schema.py upgrade %s --out new.h5" % path)
    with h5py.File(path, "r+") as h:
        said, absent = finalise(h, verbose=verbose)
    if verbose:
        print("\n%d change%s" % (len(said), "" if len(said) == 1 else "s"))
    return said


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check", help="report this file against the document")
    c.add_argument("file")
    u = sub.add_parser("upgrade", help="write a copy that matches it")
    u.add_argument("file")
    u.add_argument("--out", required=True)
    f = sub.add_parser("fix", help="run the layout pass in place, for a file "
                                   "that needs no large array rewritten")
    f.add_argument("file")
    a = ap.parse_args(argv)
    if a.cmd == "check":
        return 1 if report(a.file) else 0
    if a.cmd == "fix":
        fix(a.file)
        print()
        return 1 if report(a.file) else 0
    upgrade(a.file, a.out)
    print()
    return 1 if report(a.out) else 0


if __name__ == "__main__":
    sys.exit(main())
