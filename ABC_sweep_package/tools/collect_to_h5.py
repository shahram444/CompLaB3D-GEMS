#!/usr/bin/env python3
"""Turn a finished ABC campaign into the .h5 PRT-DeepONet-3D reads.

    python3 tools/collect_to_h5.py --campaign campaign --geometries geometries \\
                                   --out dataset

Writes dataset/dataset.h5, dataset/failures.csv and dataset/report.md.

YOU MAY NOT NEED THIS FILE
--------------------------
3D/tools/collect_complab_output.py in the PRT-DeepONet-GUI already does this job,
and the campaign is built to satisfy it exactly: every params.json carries all
thirteen keys it requires, the geometry directories are in the layout it reads,
and the species names match the .vti base names CompLaB writes. If you have that
tool, use it:

    python3 3D/tools/collect_complab_output.py \\
        --campaign campaign --geometries geometries --out dataset

This file exists so the package stands on its own, and so the .h5 can be built on
a machine that has only NumPy and h5py. It writes the SAME layout, byte for byte
in structure, and tools/verify_dataset.py will accept either.

THE LAYOUT
----------
    /geom/gid            (G,)                 int32
    /geom/material       (G,nx,ny,nz)         uint8    0 solid, 1 wall, 2 pore
    /geom/gdf            (G,nx,ny,nz)         float32  geodesic distance from inlet
    /geom/edt            (G,nx,ny,nz)         float32  distance to nearest solid
    /samples/geom_index  (S,)                 int32    row into /geom
    /samples/run_id      (S,)                 int32
    /samples/params      (S,2)                float32  [Pe, Da]
    /samples/t_norm      (S,T)                float32  0 .. 1 of this run's length
    /samples/conc        (S,T,C,nx,ny,nz)     float16  RAW, in mol/L
    /samples/velocity    (S,3,nx,ny,nz)       float16  the steady flow field
    /samples.attrs/conc_scale   (C,)          float32  per-species max

    root attrs: species, param_names, param_layout, da_column, mode, shape,
                n_samples, n_geometries

conc is stored RAW and conc_scale is recorded beside it, because
dataset_reader.py divides by conc_scale when it hands a target to the model. A
file that stored the already-divided field would be normalised twice.

Storing raw in float16 is only safe because the feed concentration is 1.0: five
decades of the field then sit inside float16's normal range. See the note in
make_campaign.py.

THE SPLIT IS BY GEOMETRY
------------------------
It is recorded per sample as /samples/split so the training script does not have
to re-derive it, and it is derived from the POROSITY, not from a random draw. A
network tested on a pore space it trained on is being scored on memory rather
than transfer, and two packings at the same porosity are close enough to each
other that splitting between them would be almost as bad.

WHAT IT REFUSES
---------------
A run is dropped, with its reason, into failures.csv when it never ran, exited
non-zero, is missing a species, holds a non-finite value, is identically zero, or
has more than 5 per cent of its voxels negative. That last one is the important
one and it is the reason tools/check_campaign.py clamp exists: a campaign whose
high-Damkohler half went negative produces a dataset quietly smaller than the
campaign, and nothing else says so.

The VTI decoding below is adapted from 3D/tools/collect_complab_output.py, which
has been through several real campaigns; the two conventions for an inline
base64 payload are not hypothetical.
"""
import argparse
import base64
import csv
import json
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import re
import struct
import sys
import zlib
from collections import OrderedDict

import numpy as np

PORE = 2

# ============================================================== VTI reading ===
_HDR = {"UInt32": ("<I", 4), "UInt64": ("<Q", 8)}


def _decode(raw, dtype, header_type, compressed):
    fmt, hs = _HDR[header_type]
    if not compressed:
        n = struct.unpack(fmt, raw[:hs])[0]
        return np.frombuffer(raw[hs:hs + n], dtype=dtype)
    nblocks = struct.unpack(fmt, raw[:hs])[0]
    off = hs * 3
    sizes = [struct.unpack(fmt, raw[off + i * hs: off + (i + 1) * hs])[0]
             for i in range(nblocks)]
    off += nblocks * hs
    buf = b""
    for s in sizes:
        buf += zlib.decompress(raw[off:off + s])
        off += s
    return np.frombuffer(buf, dtype=dtype)


def _decode_inline(b64, dtype, header_type, compressed):
    """An inline <DataArray format="binary"> payload.

    Two conventions exist and a reader has to survive both: one base64 stream
    over header and data together, or the header encoded separately and the two
    strings concatenated, which is what the VTK library itself writes. Under the
    second, decoding the concatenation as one stream returns the header alone and
    then silently yields an empty array, so the RESULT is checked here rather
    than only the exception.
    """
    b64 = re.sub(r"\s", "", b64)
    fmt, hs = _HDR[header_type]
    try:
        out = _decode(base64.b64decode(b64), dtype, header_type, compressed)
        if out.size:
            return out
    except Exception:
        pass
    if not compressed:
        hchars = 4 * ((hs + 2) // 3)
        return _decode(base64.b64decode(b64[:hchars])
                       + base64.b64decode(b64[hchars:]), dtype, header_type, compressed)
    probe = base64.b64decode(b64[:24] + "==")
    nblocks = struct.unpack(fmt, probe[:hs])[0]
    hbytes = (3 + nblocks) * hs
    hchars = 4 * ((hbytes + 2) // 3)
    return _decode(base64.b64decode(b64[:hchars])
                   + base64.b64decode(b64[hchars:]), dtype, header_type, compressed)


def _read_vti_manual(path):
    with open(path, "rb") as f:
        blob = f.read()
    head_end = blob.find(b"<AppendedData")
    header = (blob if head_end < 0 else blob[:head_end]).decode("utf-8", "replace")

    m = re.search(r'WholeExtent="([-\d\s]+)"', header)
    ext = [int(v) for v in m.group(1).split()]
    dims = (ext[1] - ext[0] + 1, ext[3] - ext[2] + 1, ext[5] - ext[4] + 1)
    n = int(np.prod(dims))

    compressed = "vtkZLibDataCompressor" in header
    ht = re.search(r'header_type="(\w+)"', header)
    header_type = ht.group(1) if ht else "UInt32"

    appended = None
    if head_end >= 0:
        enc = re.search(r'<AppendedData\s+encoding="(\w+)"',
                        blob[head_end:head_end + 200].decode("utf-8", "replace")).group(1)
        start = blob.find(b"_", head_end) + 1
        stop = blob.find(b"</AppendedData>", start)
        appended = (enc, blob[start:stop])

    np_of = {"Float64": np.float64, "Float32": np.float32,
             "Int32": np.int32, "UInt8": np.uint8, "Int64": np.int64}
    out = {}
    for mt in re.finditer(
            r'<DataArray\b([^>]*)>(.*?)</DataArray>|<DataArray\b([^>]*)/>', header, re.S):
        attrs = mt.group(1) or mt.group(3)
        body = mt.group(2) or ""
        gn = re.search(r'Name="([^"]+)"', attrs)
        gt = re.search(r'type="(\w+)"', attrs)
        gf = re.search(r'format="(\w+)"', attrs)
        if not (gn and gt and gf):
            continue
        name = gn.group(1)
        dtype = np_of.get(gt.group(1), np.float64)
        fmt = gf.group(1)
        ncomp = int((re.search(r'NumberOfComponents="(\d+)"', attrs)
                     or re.match(r"(1)", "1")).group(1))
        if fmt == "ascii":
            arr = np.fromstring(body.strip(), sep=" ", dtype=np.float64).astype(dtype)
        elif fmt == "binary":
            arr = _decode_inline(body, dtype, header_type, compressed)
        elif fmt == "appended":
            off = int(re.search(r'offset="(\d+)"', attrs).group(1))
            enc, data = appended
            raw = data[off:] if enc == "raw" \
                else base64.b64decode(re.sub(rb"\s", b"", data))[off:]
            arr = _decode(raw, dtype, header_type, compressed)
        else:
            continue
        if arr.size == n and ncomp == 1:
            out[name] = arr.reshape(dims[::-1]).transpose(2, 1, 0).astype(np.float32)
        elif arr.size == 3 * n:
            out[name] = arr.reshape(dims[::-1] + (3,)).transpose(3, 2, 1, 0).astype(np.float32)
    return out, dims


def _read_vti_pyvista(path):
    import pyvista as pv
    m = pv.read(path)
    ext = m.extent
    dims = (ext[1] - ext[0] + 1, ext[3] - ext[2] + 1, ext[5] - ext[4] + 1)
    out = {}
    n = int(np.prod(dims))
    for name in m.array_names:
        a = np.asarray(m[name])
        if a.size == n:
            out[name] = a.reshape(dims[::-1]).transpose(2, 1, 0)
        elif a.size == 3 * n:
            out[name] = a.reshape(dims[::-1] + (3,)).transpose(3, 2, 1, 0)
    return out, dims


def read_vti(path):
    """Every readable array in a .vti, keyed by name.

    pyvista first because it is faster; our own reader second. An EMPTY result
    from pyvista counts as a failure and falls through, because pyvista drops any
    array whose size does not match the extent rather than saying so.
    """
    out = None
    try:
        out, dims = _read_vti_pyvista(path)
    except Exception:
        out = None
    if not out:
        out, dims = _read_vti_manual(path)
    if not out:
        raise ValueError("no readable data array in %s" % os.path.basename(path))
    return out, dims


# ================================================================ one run ====
SNAP_RE_FMT = r"^%s_(\d{7})\.vti$"


def snapshots(out_dir, name):
    """[(iteration, path)] for every <name>_<7 digits>.vti, in iteration order.

    The ITERATION is kept, not just the ordering, so t_norm is the run's own
    fraction of its own length rather than a position in a list.
    """
    rx = re.compile(SNAP_RE_FMT % re.escape(name))
    hits = []
    if not os.path.isdir(out_dir):
        return hits
    for f in os.listdir(out_dir):
        m = rx.match(f)
        if m:
            hits.append((int(m.group(1)), os.path.join(out_dir, f)))
    return sorted(hits)


def pick_scalar(arrs, species, shape):
    """CompLaB names it 'Density'. The fallbacks cost nothing and turn a total
    failure into a successful read if a different Palabos version names it after
    the species instead."""
    for key in ("Density", species, species.lower(), species.upper(),
                "concentration", "Concentration"):
        if key in arrs and arrs[key].shape == shape:
            return arrs[key]
    cands = [v for v in arrs.values() if v.shape == shape]
    return cands[0] if len(cands) == 1 else None


def load_run(rdir, species, shape, max_conc, drop_first):
    """(conc (T,C,nx,ny,nz) float32, t_norm (T,), velocity (3,...) or None)."""
    od = os.path.join(rdir, "output")
    if not os.path.isdir(od):
        raise ValueError("no output directory")

    per = OrderedDict()
    for nm in species:
        s = snapshots(od, nm)
        if not s:
            raise ValueError("no .vti for species '%s'" % nm)
        per[nm] = s

    ntimes = min(len(v) for v in per.values())
    if ntimes < 2:
        raise ValueError("only %d snapshot(s)" % ntimes)

    iters = [it for it, _ in per[species[0]][:ntimes]]

    # THE NEAR-DUPLICATE LAST FRAME.
    #
    # The solver writes inside its loop whenever iT % save_VTK_interval == 0, so
    # with ade_max_iT = V*20 + 1 that is 21 evenly spaced frames ending at 20V.
    # Then PHASE 7 writes ONE MORE at iT = ade_max_iT, unconditionally, which is
    # a single timestep after the one before it. Keeping it costs a sample that
    # is a copy of its neighbour and makes the time axis non-uniform at the end,
    # where t_norm would read 0.9999 and then 1.0.
    #
    # Dropped when the final gap is under half the typical one, which catches
    # exactly that case and leaves a genuinely uneven series alone.
    if ntimes >= 3:
        gaps = [iters[i + 1] - iters[i] for i in range(len(iters) - 1)]
        typical = sorted(gaps)[len(gaps) // 2]
        if gaps[-1] * 2 < typical:
            ntimes -= 1
            iters = iters[:ntimes]

    span = max(iters[-1], 1)
    idx = list(range(ntimes))
    if drop_first and iters[0] == 0:
        # The frame at iT = 0 is identically zero in every case: nothing starts
        # in the domain. Keeping it spends a twenty-first of the dataset on a
        # sample whose answer is "all zero" for every geometry and every
        # parameter, which the network learns in one step and then has to keep
        # being told.
        idx = idx[1:]
    t_norm = np.array([iters[i] for i in idx], np.float32) / span

    conc = np.zeros((len(idx), len(species)) + shape, np.float32)
    for ci, nm in enumerate(species):
        for ti, k in enumerate(idx):
            path = per[nm][k][1]
            arrs, _ = read_vti(path)
            a = pick_scalar(arrs, nm, shape)
            if a is None:
                raise ValueError("no usable scalar in %s" % os.path.basename(path))
            if a.shape != shape:
                raise ValueError("shape %s != %s in %s"
                                 % (a.shape, shape, os.path.basename(path)))
            conc[ti, ci] = a

    if not np.all(np.isfinite(conc)):
        raise ValueError("non-finite values in the concentration fields")
    if float(np.abs(conc).max()) == 0.0:
        raise ValueError("every concentration field is identically zero")
    if float(conc.max()) > max_conc:
        raise ValueError("concentration %.3g exceeds the bound %.3g (blow-up)"
                         % (float(conc.max()), max_conc))
    frac_neg = float((conc < -1e-12).mean())
    if frac_neg > 0.05:
        raise ValueError("%.1f%% of voxels negative; check the [CLAMP] line"
                         % (100 * frac_neg))

    vel = None
    for base in ("nsLattice", "velocity", "vel", "flow"):
        s = snapshots(od, base)
        if not s:
            continue
        arrs, _ = read_vti(s[-1][1])            # the flow is steady; take the last
        for nm in ("velocity", "Velocity", "u", "U"):
            if nm in arrs and arrs[nm].shape == (3,) + shape:
                vel = arrs[nm]
                break
        if vel is not None:
            break

    return conc, t_norm, vel


# ================================================================== collect ===
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--campaign", default="campaign")
    ap.add_argument("--geometries", default="geometries")
    ap.add_argument("--out", default="dataset")
    ap.add_argument("--max-conc", type=float, default=1e3,
                    help="reject a run whose peak exceeds this (blow-up guard)")
    ap.add_argument("--compress", default="gzip", choices=["gzip", "lzf", "none"])
    ap.add_argument("--keep-first-frame", action="store_true",
                    help="keep the all-zero snapshot at iteration 0")
    ap.add_argument("--no-velocity", action="store_true")
    args = ap.parse_args(argv)

    import h5py

    os.makedirs(args.out, exist_ok=True)
    rd = os.path.join(args.campaign, "runs")
    if not os.path.isdir(rd):
        sys.exit("no runs/ under %s" % args.campaign)
    cases = sorted(d for d in os.listdir(rd) if d.startswith("run_"))
    if not cases:
        sys.exit("no run_* directories under %s" % rd)
    print("collecting %d cases from %s" % (len(cases), args.campaign))

    # ---- triage -----------------------------------------------------------
    good, failures = [], []
    for d in cases:
        cdir = os.path.join(rd, d)
        try:
            p = json.load(open(os.path.join(cdir, "params.json")))
        except Exception as e:
            failures.append(dict(run=d, gid=-1, stage="params", reason=str(e)))
            continue
        sp = os.path.join(cdir, "status.json")
        if not os.path.isfile(sp):
            failures.append(dict(run=d, gid=p["gid"], stage="run",
                                 reason="never ran (no status.json)"))
            continue
        try:
            st = json.load(open(sp))
        except Exception as e:
            failures.append(dict(run=d, gid=p["gid"], stage="run", reason=str(e)))
            continue
        if st.get("state") != "ok":
            failures.append(dict(run=d, gid=p["gid"], stage="run",
                                 reason=st.get("reason") or "marked failed",
                                 exit_code=st.get("exit_code"),
                                 wall_s=st.get("wall_s")))
            continue
        good.append((d, cdir, p, st))

    if not good:
        sys.exit("no successful runs: %d failed. Run "
                 "`make_campaign.py status` to see why." % len(failures))

    shape = (good[0][2]["nx"], good[0][2]["ny"], good[0][2]["nz"])
    species = list(good[0][2]["species"])
    C = len(species)

    # ---- geometries -------------------------------------------------------
    gids = sorted({p["gid"] for _, _, p, _ in good})
    gindex = {g: i for i, g in enumerate(gids)}
    G = len(gids)
    mat = np.zeros((G,) + shape, np.uint8)
    gdf = np.zeros((G,) + shape, np.float32)
    edt = np.zeros((G,) + shape, np.float32)
    for g in gids:
        f = os.path.join(args.geometries, "geom_%04d" % g, "geom_%04d.npz" % g)
        if not os.path.isfile(f):
            sys.exit("missing %s -- run tools/make_geometries.py into --geometries" % f)
        z = np.load(f)
        i = gindex[g]
        mat[i] = z["material"]
        gdf[i] = np.nan_to_num(z["gdf"])
        edt[i] = z["edt"]

    # ---- read the fields, and write each one as it is read ----------------
    #
    # STREAMED, NOT ACCUMULATED. Holding every run's field until the end costs
    # (runs) x (frames) x (species) x (voxels) x 4 bytes plus the velocity: about
    # 2.2 GB at 270 runs, all of it live while h5py allocates its own buffers on
    # top. That is a MemoryError on a login node after an hour of decoding, with
    # nothing written. Each run is written the moment it is read instead, and only
    # the per-run metadata, which is a few hundred bytes each, is kept in memory.
    comp = None if args.compress == "none" else args.compress
    h5p = os.path.join(args.out, "dataset.h5")
    N = len(good)                     # an upper bound; trimmed at the end
    meta, T, scale = [], None, None
    dc = dv = None
    n_vel = 0

    with h5py.File(h5p, "w") as h:
        gg = h.create_group("geom")
        gg.create_dataset("gid", data=np.array(gids, np.int32))
        gg.create_dataset("material", data=mat, compression=comp)
        gg.create_dataset("gdf", data=gdf, compression=comp)
        gg.create_dataset("edt", data=edt, compression=comp)
        sg = h.create_group("samples")

        for d, cdir, p, st in good:
            try:
                conc, tn, vel = load_run(cdir, species, shape, args.max_conc,
                                         not args.keep_first_frame)
            except Exception as e:
                failures.append(dict(run=d, gid=p["gid"], stage="data", reason=str(e),
                                     wall_s=st.get("wall_s")))
                continue
            if T is None:
                T = conc.shape[0]
                scale = np.full(C, 1e-30, np.float64)
                dc = sg.create_dataset("conc", shape=(N, T, C) + shape,
                                       maxshape=(None, T, C) + shape,
                                       dtype=np.float16, compression=comp,
                                       chunks=(1, 1, 1) + shape)
                if not args.no_velocity:
                    dv = sg.create_dataset("velocity", shape=(N, 3) + shape,
                                           maxshape=(None, 3) + shape,
                                           dtype=np.float16, compression=comp,
                                           chunks=(1, 3) + shape)
            if conc.shape[0] != T:
                failures.append(dict(run=d, gid=p["gid"], stage="data",
                                     reason="%d snapshots, expected %d. Every case "
                                            "should hold the same number; check "
                                            "ade_max_iT and save_VTK_interval."
                                            % (conc.shape[0], T)))
                continue

            i = len(meta)
            dc[i] = conc.astype(np.float16)
            if dv is not None and vel is not None:
                dv[i] = vel.astype(np.float16)
                n_vel += 1
            for ci in range(C):
                m = float(np.abs(conc[:, ci]).max())
                if m > scale[ci]:
                    scale[ci] = m
            meta.append((d, p, tn))
            print("  %s  gid=%-3d Pe=%-6g Da=%-6g  max=%.4g  neg=%.3f%%"
                  % (d, p["gid"], p["pe"], p["da_abio"], float(conc.max()),
                     100.0 * float((conc < -1e-12).mean())))
            del conc, vel

        if not meta:
            sys.exit("every run failed validation; see %s/failures.csv" % args.out)

        S = len(meta)
        if S != N:
            dc.resize(S, axis=0)
            if dv is not None:
                dv.resize(S, axis=0)

        sg.create_dataset("geom_index",
                          data=np.array([gindex[m[1]["gid"]] for m in meta], np.int32))
        sg.create_dataset("run_id",
                          data=np.array([m[1]["run_id"] for m in meta], np.int32))
        # [Pe, Da]. Nothing in this sweep is biotic, so the single Da column is
        # the abiotic one and da_column says so.
        sg.create_dataset("params", data=np.array(
            [[float(m[1]["pe"]), float(m[1]["da_abio"])] for m in meta], np.float32))
        sg.create_dataset("t_norm",
                          data=np.array([m[2] for m in meta], np.float32))
        sg.create_dataset("split", data=np.array(
            [m[1]["split"].encode() for m in meta], dtype="S8"))
        sg.create_dataset("porosity", data=np.array(
            [float(m[1]["porosity"]) for m in meta], np.float32))
        sg.attrs["conc_scale"] = scale.astype(np.float32)

        if dv is not None and n_vel < S:
            print()
            print("!" * 72)
            print("only %d of %d runs had a readable flow field. Without it the" % (n_vel, S))
            print("--flow-proxy and --dim-free switches cannot be used: they replace")
            print("the geometry with the flow, so they need it present. Look for")
            print("nsLattice_*.vti in each run\'s output/.")
            print("!" * 72)

        h.attrs["species"] = np.array(species, dtype="S8")
        h.attrs["param_names"] = np.array(["pe", "da"], dtype="S16")
        h.attrs["param_layout"] = b"pe_da"
        h.attrs["da_column"] = b"da_abio"
        h.attrs["mode"] = "transient"
        h.attrs["shape"] = np.array(shape, np.int32)
        h.attrs["n_samples"] = S
        h.attrs["n_geometries"] = G
        h.attrs["reaction"] = b"A + B -> C,  R = k[A][B]"
        h.attrs["boundaries"] = str(meta[0][1].get("boundaries", "unknown")).encode()

        # ---- the documented layout, in one place -------------------------
        # Three collectors write this file and each used to decide the layout
        # for itself, which is how they drifted apart. dataset_schema holds
        # the layout now, and fills in anything derivable that is still
        # missing. An entry a campaign genuinely does not have is recorded in
        # ancillary/absent with its reason rather than filled with zeros.
        try:
            import dataset_schema
            dataset_schema.finalise(h)
        except ImportError:
            print("   note: dataset_schema.py not beside this script, so the "
                  "file was written without the final layout pass. Run "
                  "dataset_schema.py upgrade on it to bring it up.")


    # ---- the two reports --------------------------------------------------
    fcsv = os.path.join(args.out, "failures.csv")
    if failures:
        keys = sorted({k for f in failures for k in f})
        with open(fcsv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(failures)
    elif os.path.exists(fcsv):
        os.remove(fcsv)

    splits = {}
    for _, p, _ in meta:
        splits.setdefault(p["split"], set()).add(p["gid"])
    with open(os.path.join(args.out, "report.md"), "w") as f:
        f.write("# ABC dataset\n\n")
        f.write("- %d samples from %d runs over %d geometries\n" % (S, S, G))
        f.write("- %d snapshots each, species %s\n" % (T, ", ".join(species)))
        f.write("- shape %s, conc_scale %s\n" % (shape, scale.round(6).tolist()))
        f.write("- %d run(s) dropped%s\n"
                % (len(failures), " (see failures.csv)" if failures else ""))
        f.write("\n## split, by geometry\n\n")
        for k in sorted(splits):
            f.write("- %s: %d geometries, gids %s\n"
                    % (k, len(splits[k]), sorted(splits[k])))
        if failures:
            f.write("\n## dropped\n\n")
            for x in failures:
                f.write("- %s (gid %s): %s\n" % (x["run"], x.get("gid"), x["reason"]))

    print()
    print("wrote %s" % h5p)
    print("  %d samples, %d geometries, %d snapshots, species %s"
          % (S, G, T, ", ".join(species)))
    print("  conc_scale %s" % scale.round(6).tolist())
    for k in sorted(splits):
        print("  %-9s %d geometries" % (k, len(splits[k])))
    if failures:
        print("  %d run(s) dropped, see %s" % (len(failures), fcsv))
    print()
    print("next:  python3 tools/verify_dataset.py --dataset %s" % h5p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
