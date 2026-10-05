#!/usr/bin/env python3
"""A window for looking at the two kinds of .h5 this project produces.

    python3 prt_h5_viewer.py                      then File, Open
    python3 prt_h5_viewer.py some_file.h5         open it straight away

THE TWO KINDS, AND HOW THEY ARE TOLD APART

  campaign dataset      what CompLaB3D produced, collected into one file.
                        Recognised by /samples/conc beside /geom/material.
                        Holds the pore spaces, the concentration, reaction rate
                        and velocity volumes of every run, and the complete
                        input record of each one.

  results               what the fitted models do, written by make_results_h5.py.
                        Recognised by the root attribute schema="prt-results".
                        Holds the accuracy on all three splits, the predicted
                        and simulated volumes of the held out set, the design of
                        the campaign, the settings every model was fitted with,
                        the weights and the notes.

Every tab is offered for both kinds. A tab whose data a particular file does not
carry says so in place of its picture rather than disappearing, so the window
looks the same whatever was opened. Anything else still opens too: the Browse tab
is a plain HDF5 tree and works on any file, so an unrecognised one is
inspectable rather than refused.

NOTHING IS LOADED UNTIL IT IS ASKED FOR. A 586 MB dataset opens instantly
because only the attributes and the small index arrays are read; a volume is
read one snapshot at a time, 32 kB at a time, when a control changes.

REQUIREMENTS
    numpy, h5py, matplotlib, and tkinter, which ships with most CPython builds.
    scikit-image is optional and only used for the grain surface in the 3D view.
"""
from __future__ import annotations

import io
import json
import os
import sys
import traceback

import numpy as np

try:
    import h5py
except ImportError:                                   # pragma: no cover
    sys.exit("h5py is required: pip install h5py")

import matplotlib
matplotlib.use("TkAgg") if "pytest" not in sys.modules else None
import matplotlib.pyplot as plt                                      # noqa: E402

try:
    from skimage import measure as _measure
except Exception:                                     # pragma: no cover
    _measure = None


# ===========================================================================
#  the file
# ===========================================================================

KIND_DATASET = "campaign dataset"
KIND_RESULTS = "results"
KIND_UNKNOWN = "unrecognised"


def _s(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.ndarray) and v.dtype.kind == "S":
        return [x.decode("utf-8", "replace") for x in v.ravel().tolist()]
    if isinstance(v, np.ndarray):
        return v.tolist()
    return v


class Doc:
    """One open file, with the handful of accessors the window needs.

    Everything here is lazy. The constructor reads attributes and the small
    index arrays and nothing else, so opening a 586 MB dataset is instant.
    """

    def __init__(self, path):
        self.path = path
        self.f = h5py.File(path, "r")
        self.kind = self._detect()
        self._cache = {}

        a = self.f.attrs
        self.species = [str(x) for x in (_s(a["species"]) if "species" in a
                                         else ["A", "B", "C"])]
        self.param_names = [str(x) for x in (_s(a["param_names"])
                                             if "param_names" in a else [])]
        self.shape = tuple(int(v) for v in a["shape"]) if "shape" in a else \
            (tuple(int(v) for v in a["grid"]) if "grid" in a else ())
        self.dim = int(a["dim"]) if "dim" in a else len(self.shape)

        self.conc_scale = None
        if "samples" in self.f and "conc_scale" in self.f["samples"].attrs:
            self.conc_scale = np.asarray(self.f["samples"].attrs["conc_scale"],
                                         np.float64)
        elif "design/scales/conc_scale" in self.f:
            self.conc_scale = self.f["design/scales/conc_scale"][:].astype(float)

        if self.kind == KIND_DATASET:
            self.gid = self.f["geom/gid"][:]
            self.geom_index = self.f["samples/geom_index"][:]
            self.params = self.f["samples/params"][:]
            self.run_id = self.f["samples/run_id"][:]
            self.t_norm = self.f["samples/t_norm"][:]
            self.n_runs = len(self.run_id)
            self.n_times = self.t_norm.shape[1]
        else:
            self.gid = (self.f["design/geometry/gid"][:]
                        if "design/geometry/gid" in self.f else np.array([]))
            self.n_runs = int(self.f.attrs.get("n_runs_in_dataset", 0))
            self.n_times = 0

    # ------------------------------------------------------------------ meta
    def _detect(self):
        a = self.f.attrs
        if str(_s(a.get("schema", ""))) == "prt-results":
            return KIND_RESULTS
        if "samples/conc" in self.f and "geom/material" in self.f:
            return KIND_DATASET
        return KIND_UNKNOWN

    def close(self):
        try:
            self.f.close()
        except Exception:
            pass

    def header(self):
        a = self.f.attrs
        L = ["FILE     %s" % os.path.basename(self.path),
             "SIZE     %.1f MB" % (os.path.getsize(self.path) / 1e6),
             "KIND     %s" % self.kind, ""]
        if self.kind == KIND_DATASET:
            L += ["This is what CompLaB3D produced, collected into one file.",
                  "It holds the pore spaces, and for every simulation the",
                  "concentration, reaction rate and velocity volumes at every",
                  "stored time, plus the complete input record of each run.", ""]
        elif self.kind == KIND_RESULTS:
            L += ["This is what the fitted models do, written by",
                  "make_results_h5.py. It holds the accuracy of every model on",
                  "the training, validation and test sets, the design of the",
                  "campaign, the settings each model was fitted with, the",
                  "fitting history, the weights and the notes.", ""]
        for k in sorted(a):
            if k.startswith("dataset_sha256") and k.endswith("bytes"):
                continue
            v = _s(a[k])
            if isinstance(v, str) and len(v) > 300:
                v = v[:300] + " ..."
            L.append("%-24s %s" % (k, v))
        return "\n".join(L)

    # -------------------------------------------------------------- geometry
    def geom_rows(self):
        """(gid, porosity, tortuosity) per pore space, where known."""
        if self.kind == KIND_DATASET:
            key = "dsgeom"
            if key not in self._cache:
                por = np.full(len(self.gid), np.nan)
                tor = np.full(len(self.gid), np.nan)
                for s in range(len(self.geom_index)):
                    g = int(self.geom_index[s])
                    if np.isnan(por[g]):
                        try:
                            pj = json.loads(_s(self.f["inputs/params_json"][s]))
                        except Exception:
                            pj = {}
                        por[g] = float(pj.get("porosity", np.nan))
                        tor[g] = float(pj.get("tortuosity",
                                              pj.get("tau_geom", np.nan)))
                self._cache[key] = (self.gid, por, tor)
            return self._cache[key]
        g = self.f["design/geometry"]
        return (g["gid"][:], g["porosity"][:], g["tortuosity"][:])

    def has_geometry_volumes(self):
        return ("geom/material" in self.f
                or "design/geometry/material" in self.f)

    def geometry_volume(self, row, what):
        """`what` is material, geodesic or euclidean."""
        if self.kind == KIND_DATASET:
            key = {"material": "geom/material", "geodesic": "geom/gdf",
                   "euclidean": "geom/edt"}[what]
        else:
            key = {"material": "design/geometry/material",
                   "geodesic": "design/geometry/geodesic_distance",
                   "euclidean": "design/geometry/euclidean_distance"}[what]
        if key not in self.f:
            return None
        return np.asarray(self.f[key][row], np.float32)

    # ------------------------------------------------------- dataset volumes
    def run_label(self, s):
        p = " ".join("%s %g" % (n, v)
                     for n, v in zip(self.param_names, self.params[s]))
        return "run %-5d  gid %-3d  %s" % (self.run_id[s],
                                           self.gid[self.geom_index[s]], p)

    def conc(self, s, t, c):
        v = np.asarray(self.f["samples/conc"][s, t, c], np.float32)
        if self.conc_scale is not None:
            v = v * float(self.conc_scale[c])
        return v

    def rate(self, s, t, c):
        if "samples/rate" not in self.f:
            return None
        v = np.asarray(self.f["samples/rate"][s, t, c], np.float32)
        rs = self.f["samples"].attrs.get("rate_scale")
        return v * float(rs[c]) if rs is not None else v

    def speed(self, s):
        """Flow speed of one run, from either velocity layout.

        The documented layout carries a flow snapshot axis,
        (run, flow snapshot, component, grid), because a pore space that seals
        or reopens has its flow re-solved. Files written before that axis
        existed are (run, component, grid). Both are read here, so an older
        campaign file opens unchanged.
        """
        d = self.f["samples/velocity"]
        v = np.asarray(d[s], np.float32)
        if d.ndim == 6:                       # (S, Tv, 3, nx, ny, nz)
            v = v[-1]                         # the flow the last snapshot saw
        return np.sqrt((v ** 2).sum(0))

    def input_choices(self):
        """What inputs/ actually holds in this file, in a sensible order."""
        if "inputs" not in self.f:
            return []
        order = ["xml", "kinetics", "order", "files", "xml_raw",
                 "kinetics_raw", "params_json", "env", "time_step"]
        got = list(self.f["inputs"])
        return [k for k in order if k in got] + \
               [k for k in got if k not in order]

    def run_text(self, s, which):
        """One entry of inputs/, rendered for reading.

        Four of these are tables with named columns rather than text, which is
        the whole point of storing them that way, but a table printed with
        str() comes out as one unbroken line of braces. So a compound dtype is
        laid out in aligned columns here, and only a genuine string is printed
        as it stands.
        """
        key = "inputs/%s" % which
        if key not in self.f:
            return "(not in this file)"
        o = self.f[key]
        if o.dtype.names is None:
            v = o[s] if o.ndim and o.shape[0] > s else o[()]
            return _s(v) if not isinstance(v, np.ndarray) else str(_s(v))
        rows = o[s] if o.ndim == 2 else o[:]
        return _table_text(rows, o.dtype.names, title=key,
                           note=_s(o.attrs.get("how_to_read_this", "")))

    # ------------------------------------------------------- results accessors
    def models(self):
        return sorted(self.f["results"]) if "results" in self.f else []

    def splits(self, tag):
        order = ["train", "val", "test"]
        got = list(self.f["results/%s" % tag]) if "results" in self.f else []
        return [s for s in order if s in got] + [s for s in got if s not in order]

    def pooled(self, tag, split):
        p = "results/%s/%s/pooled" % (tag, split)
        if p not in self.f:
            return None
        g = self.f[p]
        return {k: float(g[k][()]) for k in g if g[k].shape == ()}

    def table(self, path):
        if path not in self.f:
            return None, None
        g = self.f[path]
        cols = [str(c) for c in _s(g.attrs["columns"])]
        return cols, {c: g[c][:] for c in cols}

    def has_fields(self, tag, split):
        return "results/%s/%s/fields" % (tag, split) in self.f

    def field_rows(self, tag, split):
        g = self.f["results/%s/%s/fields" % (tag, split)]
        return g["run_id"][:], g["gid"][:], g["t_norm"][:]

    def field(self, tag, split, row):
        g = self.f["results/%s/%s/fields" % (tag, split)]
        fac = float(g.attrs.get("to_feed_factor", 1.0))
        t = np.asarray(g["truth"][row], np.float32) * fac
        p = np.asarray(g["pred"][row], np.float32) * fac
        return t, p

    # ---------------------------------------------------- settings accessors
    def model_args(self, tag):
        """The complete argument record a model was fitted with."""
        k = "models/%s/args" % tag
        if k not in self.f:
            return {}
        try:
            return json.loads(_s(self.f[k][()]))
        except Exception:
            return {}

    def model_checkpoint(self, tag):
        k = "models/%s/training/checkpoint" % tag
        if k not in self.f:
            return {}, ""
        g = self.f[k]
        out = {n: (float(g[n][()]) if g[n].dtype.kind == "f"
                   else int(g[n][()])) for n in g if g[n].shape == ()}
        return out, _s(g.attrs.get("stopped_because", ""))

    def dataset_settings(self, s=0):
        """(name, value) for one run of a campaign dataset.

        These two arrays are variable length strings, which h5py hands back as
        `object` arrays of bytes rather than as a fixed width "S" dtype, so
        each element is decoded here rather than by the general helper.
        """
        if "inputs/settings_names" not in self.f:
            return []

        def dec(x):
            return x.decode("utf-8", "replace") if isinstance(x, bytes) else str(x)

        names = [dec(x) for x in self.f["inputs/settings_names"][:]]
        vals = [dec(x) for x in self.f["inputs/settings_values"][s]]
        return list(zip(names, vals))

    def design_runs(self):
        """The run table of a results file, as (columns, {column: array})."""
        return self.table("design/runs")

    # -------------------------------------------------- results-side volumes
    def field_pairs(self):
        """Every (model, split) in a results file that carries volumes."""
        return [(t, sp) for t in self.models() for sp in self.splits(t)
                if self.has_fields(t, sp)]

    def material_for_gid(self, gid_value):
        if "design/geometry/gid" not in self.f:
            return None
        gids = self.f["design/geometry/gid"][:]
        w = np.where(gids == gid_value)[0]
        if not len(w) or "design/geometry/material" not in self.f:
            return None
        return np.asarray(self.f["design/geometry/material"][int(w[0])])

    def notes(self):
        if "notes" not in self.f:
            return {}
        return {k: _s(self.f["notes"][k][()]) for k in self.f["notes"]}

    def figures(self):
        if "figures" not in self.f:
            return []
        return sorted(self.f["figures"])


# ===========================================================================
#  what every setting means, in words
#
#  The Settings tab prints one row per setting: its name, the value this file
#  was made with, and a sentence saying what it does. The sentences below are
#  written for someone who works on reactive transport rather than on code, so
#  they say what the number changes about the science, not what it does to the
#  program.
# ===========================================================================

NOT_USED = "Optional feature. Switched off for this campaign."

PRT_SETTINGS = [
    # (name, short group, explanation)
    ("species", "what is being predicted",
     "Which chemical this model predicts. There is one separate model per "
     "chemical, so A, B and C each have their own fitted weights."),
    ("distance", "what is being predicted",
     "Which distance field the model is shown. gdf is the geodesic distance, "
     "the shortest path a molecule can actually swim through the pore space. "
     "edt is the straight line distance, which ignores grains in the way. "
     "none withholds the distance altogether. Comparing the three is the "
     "geometry ablation."),
    ("reaction", "what is being predicted",
     "The file describing the reaction, so the model knows which chemicals "
     "are fed in and which are produced."),
    ("data", "what is being predicted",
     "The campaign file the model was fitted to."),
    ("out", "what is being predicted",
     "Where the fitted weights and the training log were written."),

    ("split_file", "how the data was divided",
     "The explicit list of which pore spaces are used for training, for "
     "validation and for testing. Using a fixed list rather than a random "
     "draw is what makes the division repeatable, and it is applied by pore "
     "space, so no geometry appears on two sides."),
    ("val_frac", "how the data was divided",
     "Fraction of pore spaces that would be held back for validation if no "
     "split file were given. The split file above was given, so this was not "
     "used."),
    ("test_frac", "how the data was divided",
     "The same for the test set, and likewise unused here."),
    ("two_way_split", "how the data was divided",
     "If on, there would be no separate validation set and stopping would be "
     "decided on the test set. Off, which is what keeps the test score "
     "honest."),
    ("target_scale", "how the data was divided",
     "What the chemical values are divided by before fitting, so that all "
     "three chemicals are the same size to the optimiser. Set to train, "
     "meaning the scale is measured on the training pore spaces only, so "
     "nothing from the held out sets can leak into the fit."),

    ("n_points", "how the fitting was done",
     "How many pore voxels are drawn from each snapshot at each step. "
     "Neighbouring voxels of a smooth field carry almost the same "
     "information, so sampling a subset costs a fraction of the time and "
     "loses very little."),
    ("batch_size", "how the fitting was done",
     "How many snapshots are handled together before the weights are "
     "adjusted once."),
    ("epochs", "how the fitting was done",
     "The largest number of passes through the training set that was "
     "allowed. The run normally stops well before this, see patience."),
    ("lr", "how the fitting was done",
     "Learning rate: how large a step the optimiser takes each time it "
     "adjusts the weights. Too large and the fit never settles, too small "
     "and it takes far longer than it needs to."),
    ("patience", "how the fitting was done",
     "Stop once the validation error has failed to improve for this many "
     "passes in a row. This is what actually ended the run, and the weights "
     "kept are those of the best pass, not the last one."),
    ("seed", "how the fitting was done",
     "The random seed. Running the same command with the same seed "
     "reproduces the same model."),
    ("workers", "how the fitting was done",
     "How many processes read data at the same time. This changes how long "
     "the run takes and nothing about the answer."),

    ("with_velocity", "optional inputs, all off here", NOT_USED +
     " It would hand the model the flow field as an extra input."),
    ("flow_proxy", "optional inputs, all off here", NOT_USED +
     " It would replace the flow field with a cheaper stand in."),
    ("flow_mode", "optional inputs, all off here", NOT_USED +
     " It selects which stand in flow_proxy would use."),
    ("u_floor", "optional inputs, all off here", NOT_USED +
     " It is the smallest velocity the flow stand in would admit."),
    ("velocity_informed", "optional inputs, all off here", NOT_USED +
     " It would let the flow field shape the trunk network."),
    ("geom_features", "optional inputs, all off here", NOT_USED +
     " It would add summary numbers of the pore space as extra inputs."),
    ("keep_geometry_channel", "optional inputs, all off here", NOT_USED +
     " It would keep the raw grain map as an input channel alongside the "
     "distance field."),
    ("distance_convention", "optional inputs, all off here",
     "Which sign and origin convention the distance field uses. Left at the "
     "project default."),
    ("with_time", "optional inputs, all off here", NOT_USED +
     " It would override how time is handed to the model."),
    ("dim_free", "optional inputs, all off here", NOT_USED +
     " It would let one set of weights serve both 2D and 3D."),

    ("transfer_2d", "starting from another model, all off here", NOT_USED +
     " It would start the 3D fit from a model already fitted in 2D."),
    ("transfer_2d_frac", "starting from another model, all off here", NOT_USED +
     " It is how much of the 2D model transfer_2d would carry over."),
    ("init_from", "starting from another model, all off here", NOT_USED +
     " It would start from an existing checkpoint instead of from scratch."),
    ("freeze_trunk", "starting from another model, all off here", NOT_USED +
     " It would hold the trunk network fixed and fit only the branches."),
]

COMPLAB_SETTINGS = [
    ("Peclet", "the two numbers that define a case",
     "How much of the transport is carried by the flow against how much by "
     "diffusion. Below one, diffusion wins and the chemicals spread as much "
     "across the flow as along it."),
    ("characteristic_length", "the two numbers that define a case",
     "The length used in the Peclet and Damkohler numbers, in voxels. Every "
     "case in a campaign uses the same one so the numbers are comparable."),
    ("abiotic_rate_scale", "the two numbers that define a case",
     "A plain multiplier on the reaction rate. Left at 1 so the rate "
     "constant is set in one place only. Setting both would multiply the "
     "rate twice."),

    ("nx", "the block", "Grid size across the block, in voxels."),
    ("ny", "the block", "Grid size across the block, in voxels."),
    ("nz", "the block", "Grid size along the flow, in voxels."),
    ("dx", "the block", "The size of one voxel."),
    ("unit", "the block", "The unit the voxel size is given in."),
    ("filename", "the block", "The geometry file this case was run on."),
    ("solid", "the block",
     "The material number that marks a grain voxel in the geometry file."),
    ("pore", "the block",
     "The material number that marks an open pore voxel."),
    ("bounce_back", "the block",
     "The wall condition. Fluid that meets a grain is sent back the way it "
     "came, which is what makes the grain surface a no slip wall."),

    ("tau", "the flow",
     "The lattice relaxation time, which is what sets the fluid viscosity in "
     "the flow solver."),
    ("delta_P", "the flow",
     "The pressure difference applied across the block, which is what drives "
     "the flow."),
    ("ns_max_iT1", "the flow",
     "The largest number of flow iterations allowed in the first stage."),
    ("ns_max_iT2", "the flow",
     "The same for the second stage."),
    ("ns_converge_iT1", "the flow",
     "How still the flow has to become before the first stage is called "
     "converged."),
    ("ns_converge_iT2", "the flow",
     "The same for the second stage, a looser tolerance because the field is "
     "already close."),
    ("ns_update_interval", "the flow",
     "How often the flow is recomputed while the chemistry runs. Set very "
     "large here because the grains do not move, so the flow is solved once "
     "and reused."),

    ("number_of_substrates", "the chemistry",
     "How many chemicals are carried through the block."),
    ("name_of_substrates", "the chemistry",
     "What those chemicals are called."),
    ("in_pore", "the chemistry",
     "The diffusion coefficient of a chemical in the open pore space."),
    ("in_biofilm", "the chemistry",
     "The diffusion coefficient inside biofilm, where there is any."),
    ("initial_concentration", "the chemistry",
     "What every chemical starts at, everywhere in the block, at time zero."),
    ("left_boundary_type", "the chemistry",
     "What is imposed at the inlet face. Held fixed, or left free."),
    ("left_boundary_condition", "the chemistry",
     "The value held at the inlet face, when it is held."),
    ("right_boundary_type", "the chemistry",
     "The same for the outlet face."),
    ("right_boundary_condition", "the chemistry",
     "The value held at the outlet face, when it is held."),
    ("outer_faces", "the chemistry",
     "What happens at the four side faces, the ones not along the flow."),
    ("enable_abiotic_kinetics", "the chemistry",
     "Switches on the reaction that needs no cells. This is the path the ABC "
     "sweep used."),
    ("enable_kinetics", "the chemistry",
     "Switches on the reaction that needs cells. Off for the ABC sweep, and "
     "the path an AOM run would use instead."),
    ("biotic_mode", "the chemistry",
     "Whether biomass is carried and allowed to change."),
    ("components", "the chemistry",
     "Any extra chemical components beyond the substrates."),

    ("ade_max_iT", "time",
     "The total number of transport steps the case runs for. This is what "
     "sets how far in time the simulation reaches."),
    ("ade_update_interval", "time",
     "How often the transport solver takes a step, relative to the flow."),
    ("ade_converge_iT", "time",
     "Stop early if the chemistry stops changing. Zero means run the full "
     "number of steps regardless."),
    ("interval", "time",
     "How many steps between one stored snapshot and the next."),
    ("save_VTK_interval", "time",
     "How often a snapshot is written to disk."),
    ("save_CHK_interval", "time",
     "How often a restart file is written. Zero means never."),

    ("tolerance", "housekeeping", "The general convergence tolerance."),
    ("enabled", "housekeeping", "Whether this block of settings is active."),
    ("track_performance", "housekeeping",
     "Whether timing information is recorded."),
    ("enable_validation_diagnostics", "housekeeping",
     "Whether extra self checks are run during the case."),
    ("read_ADE_file", "housekeeping",
     "Whether to restart the chemistry from a saved file."),
    ("read_NS_file", "housekeeping",
     "Whether to restart the flow from a saved file."),
    ("ns_rerun_iT0", "housekeeping",
     "Step at which the flow would be solved again."),
    ("input_path", "housekeeping", "Where the case reads its inputs from."),
    ("output_path", "housekeeping", "Where the case writes its outputs."),
    ("src_path", "housekeeping", "Where the source tree sits."),
    ("summary_csv", "housekeeping", "The one line summary written per case."),
    ("mask_filename", "housekeeping", "Name of the mask lattice file."),
    ("bio_filename", "housekeeping", "Name of the biomass lattice file."),
    ("ns_filename", "housekeeping", "Name of the flow lattice file."),
    ("subs_filename", "housekeeping", "Name of the substrate lattice file."),
    ("species0", "housekeeping",
     "The starting value of the first chemical, when it is set separately."),
]


# ===========================================================================
#  the layout, as the write-up draws it
#
#  Each row is (indent, path in the file, the label the write-up puts beside
#  it). The Structure tab walks this and fills in what the open file actually
#  has, so the picture a reader has in their head and the file in front of
#  them are checked against each other rather than taken on trust.
# ===========================================================================

DATASET_TREE = [
    ("head", "", "(labels on the file)", "not data, just facts about it"),
    ("attr", "species",          "species        (C)",  "the field names, in order"),
    ("attr", "species_role",     "species_role   (C)",  "dissolved, mineral or microbe"),
    ("attr", "reactions",        "reactions      (nr)", "the name of each reaction channel"),
    ("attr", "param_names",      "param_names    (P)",  "what each condition is"),
    ("attr", "param_units",      "param_units    (P)",  "the unit of each one"),
    ("attr", "shape",            "shape",               "the grid, always three numbers"),
    ("attr", "spacing",          "spacing",             "dx, dy, dz"),
    ("attr", "spacing_unit",     "spacing_unit",        "m, mm or um"),
    ("attr", "mode",             "mode",                "steady or transient"),
    ("attr", "n_times",          "n_times        (T)",  "how many snapshots per run"),
    ("attr", "n_samples",        "n_samples      (S)",  "how many runs"),
    ("attr", "n_geometries",     "n_geometries   (G)",  "how many rocks"),
    ("attr", "structure_evolves", "structure_evolves",  "did the pore space change"),

    ("head", "", "grid/", ""),
    ("set",  "grid/boundary_conditions", "boundary_conditions",
     "per field, inlet and outlet"),

    ("head", "", "geom/", "the ROCKS, as they started"),
    ("set",  "geom/gid",      "gid      (G)",            "identity of each rock"),
    ("set",  "geom/material", "material (G,nx,ny,nz)",   "0 solid, 1 wall, 2 pore"),
    ("opt",  "geom/biofilm",  "biofilm  (G,nx,ny,nz)",   "where biofilm was at t=0"),
    ("opt",  "geom/mineral",  "mineral  (G,nx,ny,nz)",   "the solid phase chemical, 3D only"),
    ("set",  "geom/gdf",      "gdf      (G,nx,ny,nz)",   "distance from the inlet"),
    ("set",  "geom/edt",      "edt      (G,nx,ny,nz)",   "distance to the nearest solid"),

    ("head", "", "samples/", "the RUNS, one entry each"),
    ("set",  "samples/geom_index", "geom_index (S)",     "which rock each run used"),
    ("set",  "samples/run_id",     "run_id     (S)",     "identity of each run"),
    ("set",  "samples/params",     "params     (S,P)",   "the conditions of each run"),
    ("set",  "samples/t_norm",     "t_norm     (S,T)",   "snapshot times, 0 to 1"),
    ("set",  "samples/t_seconds",  "t_seconds  (S,T)",   "the same times in seconds"),
    ("set",  "samples/conc",       "conc  (S,T,C,nx,ny,nz)",  "every field"),
    ("set",  "samples/conc_scale", "conc_scale (C)",     "one divisor per field"),
    ("opt",  "samples/rate",       "rate  (S,T,nr,nx,ny,nz)", "how fast each reaction ran"),
    ("opt",  "samples/rate_scale", "rate_scale (nr)",    "one divisor per rate channel"),
    ("set",  "samples/velocity",   "velocity (S,Tv,3,nx,ny,nz)", "the flow field"),
    ("opt",  "samples/material",   "material (S,T,...)", "only when structure_evolves"),
    ("opt",  "samples/biofilm",    "biofilm  (S,T,...)", "only when structure_evolves"),
    ("opt",  "samples/mineral",    "mineral  (S,T,...)", "only when structure_evolves"),
    ("opt",  "samples/gdf",        "gdf      (S,T,...)", "only when structure_evolves"),

    ("head", "", "inputs/", "what each run was TOLD to do"),
    ("set",  "inputs/xml",      "xml       TABLE", "every setting in the input file"),
    ("set",  "inputs/kinetics", "kinetics  TABLE", "every rate constant"),
    ("set",  "inputs/order",    "order     TABLE", "the chemical order it expects"),
    ("set",  "inputs/files",    "files     TABLE", "which file each run got"),

    ("head", "", "ancillary/", "free text"),
    ("set",  "ancillary/run_note",     "run_note     (S)", "comments about a run"),
    ("set",  "ancillary/run_name",     "run_name     (S)", "the folder it came from"),
    ("set",  "ancillary/dataset_note", "dataset_note",     "comments about the campaign"),
]


def structure_text(doc, width=116):
    """The tree, with what this file actually has filled in beside it."""
    h = doc.f
    L = ["%s" % os.path.basename(doc.path),
         "=" * len(os.path.basename(doc.path)), "",
         "The layout on the left is what the write-up draws. The column on "
         "the right is",
         "what this file actually holds. An entry marked `when used` is one a "
         "campaign",
         "only has if it used the thing it describes, so a blank there is not "
         "a fault.",
         ""]
    L.append("%-44s %-19s %s" % ("LAYOUT", "IN HERE", "WHAT IT IS"))
    L.append("-" * width)
    missing = []
    for kind, path, drawn, meaning in DATASET_TREE:
        if kind == "head":
            L.append("")
            L.append("+-- %-40s %-19s %s" % (drawn, "", meaning))
            continue
        if kind == "attr":
            present = path in h.attrs
            got = ""
            if present:
                v = _s(h.attrs[path])
                if isinstance(v, list):
                    got = "(%d)" % len(v)
                else:
                    got = str(v)
                    if len(got) > 19:
                        got = got[:18] + "\u2026"
        else:
            present = path in h
            got = str(h[path].shape).replace(" ", "") if present else ""
            if present and h[path].dtype.names is not None:
                got = "table"
        if not present:
            got = "--" if kind == "opt" else "ABSENT"
            if kind == "set":
                missing.append(path)
        L.append("|     %-38s %-19s %s" % (drawn, got, meaning))
    # The labels say how many of each thing there are. The arrays say it
    # again, in their own shapes. Where the two disagree, one of them is
    # wrong, and which it is cannot be decided from here, so it is reported
    # rather than quietly corrected.
    def count(name):
        v = _s(h.attrs[name]) if name in h.attrs else None
        return len(v) if isinstance(v, list) else None

    mismatch = []
    pairs = [("species", "samples/conc", 2, "chemicals"),
             ("reactions", "samples/rate", 2, "reaction channels"),
             ("param_names", "samples/params", 1, "conditions")]
    for label, path, axis, what in pairs:
        n = count(label)
        if n is None or path not in h or h[path].ndim <= axis:
            continue
        got = int(h[path].shape[axis])
        if got != n:
            mismatch.append(
                "%s says %d %s, but %s carries %d on axis %d"
                % (label, n, what, path, got, axis))
    L.append("")
    if mismatch:
        L.append("%d label%s that does not agree with the array it describes:"
                 % (len(mismatch), "" if len(mismatch) == 1 else "s"))
        for m in mismatch:
            L.append("   %s" % m)
        L.append("")
    if missing:
        L.append("%d entr%s the layout asks for and this file does not have:"
                 % (len(missing), "y" if len(missing) == 1 else "ies"))
        for m in missing:
            L.append("   %s" % m)
        L.append("")
        L.append("dataset_schema.py upgrade <file> --out fixed.h5  writes a "
                 "copy that does.")
    elif not mismatch:
        L.append("This file matches the layout.")
    return "\n".join(L)


def _table_text(rows, names, title="", note="", width=118):
    """A compound dataset as aligned columns, one row per line.

    Column widths are measured from the content rather than fixed, so a table
    of short flags does not get the same spacing as one of file paths, and a
    value too wide for its share is cut with an ellipsis rather than allowed
    to shove every column after it out of line.
    """
    def cell(v):
        if isinstance(v, bytes):
            return v.decode("utf-8", "replace")
        if isinstance(v, (float, np.floating)):
            return "" if not np.isfinite(v) else ("%g" % v)
        return str(v)

    body = [[cell(r[n]) for n in names] for r in np.atleast_1d(rows)]
    body = [r for r in body if any(x.strip() for x in r)]
    if not body:
        return "%s\n\n(every row of this table is empty for this run)" % title
    w = [max(len(names[i]), max(len(r[i]) for r in body))
         for i in range(len(names))]
    # share the line out fairly when the natural widths do not fit
    while sum(w) + 2 * len(w) > width:
        j = w.index(max(w))
        if w[j] <= 12:
            break
        w[j] -= 1
    def line(cells):
        out = []
        for i, c in enumerate(cells):
            c = c if len(c) <= w[i] else c[:w[i] - 1] + "\u2026"
            out.append("%-*s" % (w[i], c))
        return "  ".join(out).rstrip()

    L = []
    if title:
        L += [title, "=" * len(title), ""]
    if note:
        L += _wrap(note, width) + [""]
    L += [line([n.upper() for n in names]), "-" * min(sum(w) + 2 * len(w), width)]
    L += [line(r) for r in body]
    L += ["", "%d row%s" % (len(body), "" if len(body) == 1 else "s")]
    return "\n".join(L)


def _wrap(text, width):
    """Wrap on spaces, and hard break anything with none.

    A file path has no spaces in it, so wrapping on whitespace alone lets it
    run past its column and shove the next one sideways. Breaking the long
    token keeps the three columns square whatever is in them.
    """
    out, line = [], ""
    for w in str(text).split():
        while len(w) > width:
            if line:
                out.append(line); line = ""
            out.append(w[:width]); w = w[width:]
        if not w:
            continue
        if len(line) + len(w) + 1 > width:
            out.append(line); line = w
        else:
            line = (line + " " + w).strip()
    if line:
        out.append(line)
    return out or [""]


def settings_table(pairs, glossary, title, lead, name_w=28, val_w=26,
                   exp_w=62):
    """Three aligned columns: the setting, the value, and what it means.

    Settings are printed in the order of the glossary and grouped by its
    second field, because the order a program happens to declare its options
    in is not the order a reader wants to meet them. Anything in the file
    that the glossary does not cover is printed at the end rather than
    dropped, so the table is always complete.
    """
    have = dict(pairs)
    known = [g[0] for g in glossary]
    rule = "-" * (name_w + val_w + exp_w + 4)
    L = [title, "=" * len(title), ""]
    L += _wrap(lead, name_w + val_w + exp_w + 4) + [""]
    L += ["%-*s  %-*s  %s" % (name_w, "SETTING", val_w, "VALUE",
                              "WHAT IT MEANS"), rule]
    group = None
    for name, grp, expl in glossary:
        if name not in have:
            continue
        if grp != group:
            group = grp
            L += ["", "[ %s ]" % grp.upper(), ""]
        val = str(have[name])
        if val in ("None", "none", "null", ""):
            val = "not set"
        vl = _wrap(val, val_w)
        el = _wrap(expl, exp_w)
        for i in range(max(len(vl), len(el))):
            L.append("%-*s  %-*s  %s" % (
                name_w, name if i == 0 else "",
                val_w, vl[i] if i < len(vl) else "",
                el[i] if i < len(el) else ""))
        L.append("")
    extra = [k for k in have if k not in known]
    if extra:
        L += ["", "[ EVERYTHING ELSE IN THE FILE ]", ""]
        for k in sorted(extra):
            L.append("%-*s  %s" % (name_w, k, str(have[k])[:val_w + exp_w]))
    return "\n".join(L)


# ===========================================================================
#  drawing. Every function here takes a Figure and nothing from tkinter, so
#  each one can be rendered and checked without a window.
# ===========================================================================

GRAIN = "#d9d9d9"


def _mask_to_nan(vol, material):
    """Grain is absent, not zero. Blanking it keeps it out of every colour
    limit, which is the difference between a readable picture and one whose
    range is set by voxels that hold no value at all."""
    if material is None:
        return vol
    out = np.array(vol, np.float32)
    out[material < 2] = np.nan
    return out


MATERIAL_CMAP = matplotlib.colors.ListedColormap(
    ["#3b3b3b", "#9a9a9a", "#f2e8c9"])          # grain, grain surface, pore


def draw_ortho(fig, vol, material, title, cmap="viridis", unit="",
               vmin=None, vmax=None, idx=None, discrete=False):
    """Three orthogonal slices through the middle, and the distribution.

    Two rows rather than one. Three square panels across a wide canvas leave
    most of the height empty, and the fourth panel earns its place: it says how
    the values are distributed, which a slice through the middle cannot.
    """
    fig.clear()
    v = _mask_to_nan(vol, material)
    if np.all(~np.isfinite(v)):
        fig.text(.5, .5, "nothing to show", ha="center"); return
    vmin = np.nanmin(v) if vmin is None else vmin
    vmax = np.nanmax(v) if vmax is None else vmax
    if vmax <= vmin:
        vmax = vmin + 1e-9
    nx, ny, nz = v.shape
    i, j, k = idx or (nx // 2, ny // 2, nz // 2)
    panes = [(v[i, :, :], "x = %d" % i, "y", "z"),
             (v[:, j, :], "y = %d" % j, "x", "z"),
             (v[:, :, k], "z = %d" % k, "x", "y")]
    axes = fig.subplots(2, 2).ravel()
    cm = MATERIAL_CMAP if discrete else cmap
    if discrete:
        vmin, vmax = -0.5, 2.5
    im = None
    for ax, (sl, ttl, xl, yl) in zip(axes, panes):
        ax.set_facecolor(GRAIN)
        im = ax.imshow(sl.T, origin="lower", cmap=cm, vmin=vmin, vmax=vmax,
                       interpolation="nearest")
        ax.set_title(ttl, fontsize=9)
        ax.set_xlabel(xl, fontsize=8); ax.set_ylabel(yl, fontsize=8)
        ax.tick_params(labelsize=7)
    cb = fig.colorbar(im, ax=axes[:3].tolist(), fraction=.045, pad=.02,
                      shrink=.9)
    if discrete:
        cb.set_ticks([0, 1, 2])
        cb.set_ticklabels(["grain", "grain surface", "pore"])
    cb.set_label("" if discrete else unit, fontsize=8)
    cb.ax.tick_params(labelsize=7)

    ax = axes[3]
    good = v[np.isfinite(v)]
    if discrete:
        n = [int(np.sum(np.asarray(vol) == c)) for c in (0, 1, 2)]
        ax.bar(["grain", "grain\nsurface", "pore"], n,
               color=["#3b3b3b", "#9a9a9a", "#d9c99a"])
        tot = max(sum(n), 1)
        for x, c in enumerate(n):
            ax.text(x, c, " %.1f%%" % (100.0 * c / tot), ha="center",
                    va="bottom", fontsize=7.5)
        ax.set_ylabel("voxels", fontsize=8)
        ax.margins(y=.16)
    else:
        ax.hist(good, bins=60, color="#5a5a5a")
        ax.set_yscale("log")
        ax.set_xlabel(unit, fontsize=8)
        ax.set_ylabel("pore voxels", fontsize=8)
    ax.set_title("voxel counts over the whole grid" if discrete
                 else "over the whole volume", fontsize=9)
    ax.tick_params(labelsize=7)
    fig.suptitle(title, fontsize=10)


def _cut_mask(shape, cut):
    """Which voxels to keep, so the inside of the volume can be seen.

    Drawing every pore voxel of a 32 cubed domain puts about nine thousand
    opaque markers in front of each other: the outside hides the inside and the
    picture reads as noise. Removing a corner, or half, leaves the structure
    visible and shows what is in the middle, which is where the chemistry is.
    """
    nx, ny, nz = shape
    X, Y, Z = np.ogrid[:nx, :ny, :nz]
    if cut == "corner":
        return ~((X >= nx // 2) & (Y >= ny // 2) & (Z >= nz // 2))
    if cut == "half":
        return Y < ny // 2
    if cut == "quarter":
        return ~((X >= nx // 2) & (Y >= ny // 2))
    return np.ones(shape, bool)


def draw_3d(fig, vol, material, title, cmap="viridis", unit="",
            vmin=None, vmax=None, grain=True, cut="corner", top=None):
    """The pore voxels in three dimensions, coloured by value.

    A scatter rather than a volume render, because it needs nothing beyond
    matplotlib and shows the pore network honestly. The colour limits are taken
    from the WHOLE volume, not from the part left after the cut, so changing
    the cut never changes what a colour means. The grain surface is drawn
    faintly behind it when scikit-image is present.
    """
    fig.clear()
    ax = fig.add_subplot(111, projection="3d")
    v = _mask_to_nan(vol, material)
    if not np.any(np.isfinite(v)):
        fig.text(.5, .5, "nothing to show", ha="center"); return
    vmin = float(np.nanmin(v)) if vmin is None else vmin
    vmax = float(np.nanmax(v)) if vmax is None else vmax
    if vmax <= vmin:
        vmax = vmin + 1e-9

    keep = np.isfinite(v) & _cut_mask(v.shape, cut)
    if top:                       # only the highest values, for a hot spot
        thr = np.nanpercentile(v, 100.0 - float(top))
        keep &= (v >= thr)
    pts = np.argwhere(keep)
    if pts.size == 0:
        fig.text(.5, .5, "the cut removed everything", ha="center"); return
    val = v[pts[:, 0], pts[:, 1], pts[:, 2]]

    if grain and material is not None and _measure is not None:
        try:
            m = (np.asarray(material) >= 1).astype(np.float32)
            m = np.where(_cut_mask(m.shape, cut), m, 0.0)
            verts, faces, _, _ = _measure.marching_cubes(m, 0.5)
            ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2],
                            color=GRAIN, alpha=.18, linewidth=0, shade=True)
        except Exception:
            pass

    p = ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=val, cmap=cmap,
                   vmin=vmin, vmax=vmax, s=9, marker="s", depthshade=True,
                   linewidths=0, alpha=.95)
    cb = fig.colorbar(p, ax=ax, fraction=.03, pad=.02)
    cb.set_label(unit, fontsize=8); cb.ax.tick_params(labelsize=7)
    ax.set_xlabel("x", fontsize=8); ax.set_ylabel("y", fontsize=8)
    ax.set_zlabel("z", fontsize=8)
    ax.tick_params(labelsize=7)
    for a, n in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), v.shape):
        a(0, n)
    try:
        ax.set_box_aspect(v.shape)
    except Exception:
        pass
    ax.view_init(elev=22, azim=-125)
    sub = {"corner": "a corner removed", "half": "cut in half",
           "quarter": "a quarter removed", "none": "whole"}.get(cut, "")
    if top:
        sub += ", only the highest %g%% of values" % top
    fig.suptitle("%s        %s        colour limits from the whole volume"
                 % (title, sub), fontsize=10)


def draw_compare(fig, truth, pred, material, title, cmap="viridis",
                 unit="fraction of the feed", k=None):
    """Simulated, predicted and the difference, on two cuts.

    The first two share one colour scale, so a difference between them is a
    difference in the field and not in the scale. The third gets a symmetric
    scale about zero, so the colour says the SIGN of the error at a glance.
    Two cuts rather than one, because a single plane can flatter a prediction.
    """
    fig.clear()
    t = _mask_to_nan(truth, material); p = _mask_to_nan(pred, material)
    d = p - t
    nx, ny, nz = t.shape
    k = nz // 2 if k is None else k
    j = ny // 2
    vmin = float(np.nanmin([np.nanmin(t), np.nanmin(p)]))
    vmax = float(np.nanmax([np.nanmax(t), np.nanmax(p)]))
    if vmax <= vmin:
        vmax = vmin + 1e-9
    lim = float(np.nanmax(np.abs(d))) or 1e-9
    axes = fig.subplots(2, 3)
    rows = [((t[:, :, k], p[:, :, k], d[:, :, k]), "z = %d" % k, "x", "y"),
            ((t[:, j, :], p[:, j, :], d[:, j, :]), "y = %d" % j, "x", "z")]
    for r, (sls, cut, xl, yl) in enumerate(rows):
        for c, (sl, ttl, cm, a, b) in enumerate(
                [(sls[0], "simulated", cmap, vmin, vmax),
                 (sls[1], "predicted", cmap, vmin, vmax),
                 (sls[2], "predicted minus simulated", "RdBu_r", -lim, lim)]):
            ax = axes[r, c]
            ax.set_facecolor(GRAIN)
            im = ax.imshow(sl.T, origin="lower", cmap=cm, vmin=a, vmax=b,
                           interpolation="nearest")
            ax.set_title("%s,  %s" % (ttl, cut) if r == 0 else cut, fontsize=8.5)
            ax.set_xlabel(xl, fontsize=8); ax.set_ylabel(yl, fontsize=8)
            ax.tick_params(labelsize=7)
            cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.02)
            cb.ax.tick_params(labelsize=6)
            if r == 0 and c == 2:
                cb.set_label("red high, blue low", fontsize=7)
    ok = np.isfinite(t) & np.isfinite(p)
    rmse = float(np.sqrt(np.mean((p[ok] - t[ok]) ** 2))) if ok.any() else np.nan
    bias = float(np.mean(p[ok] - t[ok])) if ok.any() else np.nan
    fig.suptitle("%s        over the whole volume: rmse %.5f, bias %+.5f  (%s)"
                 % (title, rmse, bias, unit), fontsize=10)


def draw_timeseries(fig, t, mean, mx, mn, title, unit):
    fig.clear()
    ax = fig.add_subplot(111)
    ax.fill_between(t, mn, mx, color="#cfcfcf", label="range over pore voxels")
    ax.plot(t, mean, "-o", color="#1f1f1f", lw=1.5, ms=3, label="mean")
    ax.set_xlabel("normalised time, t / t_end", fontsize=9)
    ax.set_ylabel(unit, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=.25, lw=.6)


# ===========================================================================
#  the window
# ===========================================================================

def _run_gui(initial=None, on_ready=None):
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg, NavigationToolbar2Tk)
    from matplotlib.figure import Figure

    CMAPS = ["viridis", "turbo", "plasma", "inferno", "magma", "cividis",
             "coolwarm", "Spectral_r", "gray"]

    class Panel(ttk.Frame):
        """A figure with a toolbar, plus a row of controls above it."""

        def __init__(self, master):
            super().__init__(master)
            self.controls = ttk.Frame(self)
            self.controls.pack(side="top", fill="x", padx=6, pady=4)
            self.fig = Figure(figsize=(9, 5.2), dpi=100,
                              constrained_layout=True)
            self.canvas = FigureCanvasTkAgg(self.fig, master=self)
            self.canvas.get_tk_widget().pack(side="top", fill="both",
                                             expand=True)
            tb = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
            tb.update()
            tb.pack(side="bottom", fill="x")

        def redraw(self):
            self.canvas.draw_idle()

    class TextPanel(ttk.Frame):
        def __init__(self, master, mono=True):
            super().__init__(master)
            self.txt = tk.Text(self, wrap="none", font=
                               ("Consolas" if mono else "Segoe UI", 9),
                               background="#ffffff", relief="flat")
            ys = ttk.Scrollbar(self, orient="vertical", command=self.txt.yview)
            xs = ttk.Scrollbar(self, orient="horizontal",
                               command=self.txt.xview)
            self.txt.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
            self.txt.grid(row=0, column=0, sticky="nsew")
            ys.grid(row=0, column=1, sticky="ns")
            xs.grid(row=1, column=0, sticky="ew")
            self.rowconfigure(0, weight=1); self.columnconfigure(0, weight=1)

        def set(self, s):
            self.txt.configure(state="normal")
            self.txt.delete("1.0", "end")
            self.txt.insert("1.0", s)
            self.txt.configure(state="disabled")

    # ----------------------------------------------------------------- app
    class App(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title("PRT campaign file viewer")
            self.geometry("1280x860")
            self.doc = None

            bar = ttk.Frame(self)
            bar.pack(side="top", fill="x", padx=8, pady=6)
            ttk.Button(bar, text="Open a .h5 file",
                       command=self.open_dialog).pack(side="left")
            self.lbl = ttk.Label(bar, text="no file open",
                                 font=("Segoe UI", 10, "bold"))
            self.lbl.pack(side="left", padx=12)
            self.kindlbl = ttk.Label(bar, text="", foreground="#555555")
            self.kindlbl.pack(side="left")

            self.nb = ttk.Notebook(self)
            self.nb.pack(side="top", fill="both", expand=True, padx=8, pady=4)

            self.status = ttk.Label(self, text="ready", anchor="w",
                                    relief="sunken")
            self.status.pack(side="bottom", fill="x")

        # -------------------------------------------------------- plumbing
        def say(self, s):
            self.status.configure(text=s)
            self.update_idletasks()

        def open_dialog(self):
            p = filedialog.askopenfilename(
                title="Open a .h5 file",
                filetypes=[("HDF5", "*.h5 *.hdf5"), ("all files", "*.*")])
            if p:
                self.open(p)

        def open(self, path):
            try:
                if self.doc:
                    self.doc.close()
                self.say("opening %s" % os.path.basename(path))
                self.doc = Doc(path)
            except Exception as e:
                messagebox.showerror("Could not open the file",
                                     "%s\n\n%s" % (path, e))
                return
            self.lbl.configure(text=os.path.basename(path))
            self.kindlbl.configure(text="   %s   %.1f MB"
                                   % (self.doc.kind,
                                      os.path.getsize(path) / 1e6))
            self.build_tabs()
            self.say("%s, %s" % (os.path.basename(path), self.doc.kind))

        def build_tabs(self):
            for t in self.nb.tabs():
                self.nb.forget(t)
            d = self.doc
            self.tab_overview()
            if d.has_geometry_volumes():
                self.tab_pore_spaces()
            if d.kind == KIND_DATASET:
                self.tab_simulations()
                self.tab_inputs()
                self.tab_settings()
                self.tab_structure()
            elif d.kind == KIND_RESULTS:
                self.tab_simulations_results()
                self.tab_inputs_results()
                self.tab_predictions()
                self.tab_settings()
                if d.figures():
                    self.tab_figures()
            self.tab_notes()
            self.tab_browse()

        # -------------------------------------------------------- overview
        def tab_overview(self):
            f = ttk.Frame(self.nb); self.nb.add(f, text="Overview")
            tp = TextPanel(f); tp.pack(fill="both", expand=True)
            d = self.doc
            out = [d.header()]
            if d.kind == KIND_DATASET:
                gid, por, tor = d.geom_rows()
                out += ["", "PORE SPACES  %d" % len(gid),
                        "  porosity   %.3f to %.3f" % (np.nanmin(por),
                                                       np.nanmax(por)),
                        "  tortuosity %.2f to %.2f" % (np.nanmin(tor),
                                                       np.nanmax(tor)),
                        "", "SIMULATIONS  %d, each with %d stored times"
                        % (d.n_runs, d.n_times)]
                for j, n in enumerate(d.param_names):
                    lv = sorted(set(np.round(d.params[:, j], 6).tolist()))
                    out.append("  %-4s levels %s" % (n, lv))
            elif "results" in d.f:
                out += ["", "MODELS  %s" % ", ".join(d.models()),
                        "",
                        "The numbers each model scored are in the file and can",
                        "be read with read_results_h5.py, or in the Browse tab",
                        "under results. This window shows the fields and the",
                        "pore spaces rather than the scores."]
            tp.set("\n".join(out))

        # ------------------------------------------------------ pore spaces
        def tab_pore_spaces(self):
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Pore spaces")
            p = Panel(f); p.pack(fill="both", expand=True)
            gid, por, tor = d.geom_rows()

            labels = ["gid %d   porosity %.3f   tortuosity %.2f"
                      % (g, por[i], tor[i]) for i, g in enumerate(gid)]
            which = tk.StringVar(value=labels[0] if labels else "")
            what = tk.StringVar(value="material")
            view = tk.StringVar(value="slices")
            cmap = tk.StringVar(value="viridis")
            cut = tk.StringVar(value="corner")

            def draw(*_):
                try:
                    row = labels.index(which.get())
                except ValueError:
                    return
                self.say("reading pore space %d" % gid[row])
                mat = d.geometry_volume(row, "material")
                w = what.get()
                vol = mat if w == "material" else d.geometry_volume(row, w)
                if vol is None:
                    p.fig.clear()
                    p.fig.text(.5, .5, "%s is not in this file" % w,
                               ha="center")
                    p.redraw(); return
                unit = ("0 grain, 1 grain surface, 2 pore" if w == "material"
                        else "distance, voxels")
                ttl = "pore space %d, %s" % (gid[row], w)
                m = None if w == "material" else mat
                if view.get() == "slices":
                    draw_ortho(p.fig, vol, m, ttl, cmap.get(), unit,
                               discrete=(w == "material"))
                else:
                    draw_3d(p.fig, vol, m if m is not None else
                            (mat if w != "material" else None), ttl,
                            cmap.get(), unit, grain=(w != "material"),
                            cut=cut.get())
                p.redraw(); self.say("ready")

            c = p.controls
            ttk.Label(c, text="pore space").pack(side="left")
            cb = ttk.Combobox(c, textvariable=which, values=labels, width=42,
                              state="readonly")
            cb.pack(side="left", padx=6); cb.bind("<<ComboboxSelected>>", draw)
            ttk.Label(c, text="show").pack(side="left", padx=(12, 0))
            cb2 = ttk.Combobox(c, textvariable=what, state="readonly", width=12,
                               values=["material", "geodesic", "euclidean"])
            cb2.pack(side="left", padx=6); cb2.bind("<<ComboboxSelected>>", draw)
            for v, t in (("slices", "slices"), ("3d", "3D")):
                ttk.Radiobutton(c, text=t, value=v, variable=view,
                                command=draw).pack(side="left", padx=4)
            ttk.Label(c, text="3D cut").pack(side="left", padx=(10, 0))
            cb4 = ttk.Combobox(c, textvariable=cut, width=8, state="readonly",
                               values=["corner", "quarter", "half", "none"])
            cb4.pack(side="left", padx=4); cb4.bind("<<ComboboxSelected>>", draw)
            ttk.Label(c, text="colours").pack(side="left", padx=(12, 0))
            cb3 = ttk.Combobox(c, textvariable=cmap, values=CMAPS, width=11,
                               state="readonly")
            cb3.pack(side="left", padx=6); cb3.bind("<<ComboboxSelected>>", draw)
            draw()

        # ------------------------------------------------------ simulations
        def tab_simulations(self):
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Simulations")
            p = Panel(f); p.pack(fill="both", expand=True)

            labels = [d.run_label(s) for s in range(d.n_runs)]
            which = tk.StringVar(value=labels[0])
            field = tk.StringVar(value="concentration")
            sp = tk.StringVar(value=d.species[0])
            view = tk.StringVar(value="slices")
            cmap = tk.StringVar(value="viridis")
            cut = tk.StringVar(value="corner")
            tix = tk.IntVar(value=d.n_times - 1)

            def draw(*_):
                try:
                    s = labels.index(which.get())
                except ValueError:
                    return
                t = int(tix.get())
                c = d.species.index(sp.get())
                mat = np.asarray(d.f["geom/material"][int(d.geom_index[s])])
                self.say("reading run %d" % d.run_id[s])
                if field.get() == "concentration":
                    vol = d.conc(s, t, c); unit = "c / c_in"
                elif field.get() == "reaction rate":
                    vol = d.rate(s, t, c); unit = "rate"
                elif field.get() == "velocity magnitude":
                    vol = d.speed(s); unit = "speed, lattice units"
                else:
                    vol = None; unit = ""
                if field.get() == "time series":
                    pore = mat == 2
                    mean, mx, mn = [], [], []
                    for k in range(d.n_times):
                        v = d.conc(s, k, c)[pore]
                        mean.append(v.mean()); mx.append(v.max())
                        mn.append(v.min())
                    draw_timeseries(p.fig, d.t_norm[s], np.array(mean),
                                    np.array(mx), np.array(mn),
                                    "%s   species %s" % (labels[s], sp.get()),
                                    "c / c_in")
                    p.redraw(); self.say("ready"); return
                if vol is None:
                    p.fig.clear()
                    p.fig.text(.5, .5, "not in this file", ha="center")
                    p.redraw(); return
                ttl = "%s\n%s %s   t/t_end %.2f" % (
                    labels[s], field.get(), sp.get()
                    if field.get() != "velocity magnitude" else "",
                    d.t_norm[s, t])
                if view.get() == "slices":
                    draw_ortho(p.fig, vol, mat, ttl, cmap.get(), unit)
                else:
                    draw_3d(p.fig, vol, mat, ttl, cmap.get(), unit,
                            cut=cut.get())
                p.redraw(); self.say("ready")

            c = p.controls
            ttk.Label(c, text="run").pack(side="left")
            cb = ttk.Combobox(c, textvariable=which, values=labels, width=38,
                              state="readonly")
            cb.pack(side="left", padx=6); cb.bind("<<ComboboxSelected>>", draw)
            cb2 = ttk.Combobox(c, textvariable=field, state="readonly", width=18,
                               values=["concentration", "reaction rate",
                                       "velocity magnitude", "time series"])
            cb2.pack(side="left", padx=6); cb2.bind("<<ComboboxSelected>>", draw)
            cb3 = ttk.Combobox(c, textvariable=sp, values=d.species, width=4,
                               state="readonly")
            cb3.pack(side="left", padx=4); cb3.bind("<<ComboboxSelected>>", draw)
            ttk.Label(c, text="time").pack(side="left", padx=(12, 2))
            sc = ttk.Scale(c, from_=0, to=d.n_times - 1, orient="horizontal",
                           length=180,
                           command=lambda v: (tix.set(int(float(v))), draw()))
            sc.set(d.n_times - 1); sc.pack(side="left")
            for v, t in (("slices", "slices"), ("3d", "3D")):
                ttk.Radiobutton(c, text=t, value=v, variable=view,
                                command=draw).pack(side="left", padx=4)
            cb5 = ttk.Combobox(c, textvariable=cut, width=8, state="readonly",
                               values=["corner", "quarter", "half", "none"])
            cb5.pack(side="left", padx=4); cb5.bind("<<ComboboxSelected>>", draw)
            cb4 = ttk.Combobox(c, textvariable=cmap, values=CMAPS, width=11,
                               state="readonly")
            cb4.pack(side="left", padx=6); cb4.bind("<<ComboboxSelected>>", draw)
            draw()

        # ----------------------------------------------------------- inputs
        def tab_inputs(self):
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Run inputs")
            top = ttk.Frame(f); top.pack(side="top", fill="x", padx=6, pady=4)
            tp = TextPanel(f); tp.pack(fill="both", expand=True)
            labels = [d.run_label(s) for s in range(d.n_runs)]
            which = tk.StringVar(value=labels[0])
            choices = d.input_choices() or ["xml"]
            what = tk.StringVar(value=choices[0])

            def show(*_):
                try:
                    s = labels.index(which.get())
                except ValueError:
                    s = 0
                tp.set(d.run_text(s, what.get()))

            ttk.Label(top, text="run").pack(side="left")
            cb = ttk.Combobox(top, textvariable=which, values=labels, width=38,
                              state="readonly")
            cb.pack(side="left", padx=6); cb.bind("<<ComboboxSelected>>", show)
            cb2 = ttk.Combobox(top, textvariable=what, state="readonly",
                               width=14, values=choices)
            cb2.pack(side="left", padx=6); cb2.bind("<<ComboboxSelected>>", show)
            ttk.Label(top, text="   the complete record of what this run was "
                                "told to do").pack(side="left")
            show()

        # ------------------------------------------------------ predictions
        def tab_predictions(self):
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Predictions")
            p = Panel(f); p.pack(fill="both", expand=True)
            pairs = [(t, s) for t in d.models() for s in d.splits(t)
                     if d.has_fields(t, s)]
            if not pairs:
                return
            names = ["%s, %s" % x for x in pairs]
            which = tk.StringVar(value=names[0])
            view = tk.StringVar(value="compare")
            cmap = tk.StringVar(value="viridis")
            cut = tk.StringVar(value="corner")
            rowv = tk.IntVar(value=0)
            state = {}

            def reload_rows(*_):
                t, s = pairs[names.index(which.get())]
                rid, gid, tn = d.field_rows(t, s)
                state["t"], state["s"] = t, s
                state["rid"], state["gid"], state["tn"] = rid, gid, tn
                state["labels"] = ["run %d  gid %d  t/t_end %.2f"
                                   % (rid[i], gid[i], tn[i])
                                   for i in range(len(rid))]
                sc.configure(to=len(rid) - 1)
                # Start at the first fully developed snapshot. Row 0 is t = 0,
                # where every field is still essentially empty, which is a poor
                # thing to meet first.
                last = int(np.where(tn == tn.max())[0][0])
                rowv.set(last); sc.set(last)
                draw()

            def material_for(gid_value):
                if "design/geometry/gid" not in d.f:
                    return None
                gids = d.f["design/geometry/gid"][:]
                w = np.where(gids == gid_value)[0]
                if not len(w) or "design/geometry/material" not in d.f:
                    return None
                return np.asarray(d.f["design/geometry/material"][int(w[0])])

            def draw(*_):
                if "t" not in state:
                    return
                i = int(rowv.get())
                i = max(0, min(i, len(state["rid"]) - 1))
                self.say("reading %s" % state["labels"][i])
                tr, pr = d.field(state["t"], state["s"], i)
                mat = material_for(int(state["gid"][i]))
                ttl = "%s, %s      %s" % (state["t"], state["s"],
                                          state["labels"][i])
                if view.get() == "compare":
                    draw_compare(p.fig, tr, pr, mat, ttl, cmap.get())
                elif view.get() == "predicted 3D":
                    draw_3d(p.fig, pr, mat, ttl + "   predicted",
                            cmap.get(), "fraction of the feed", cut=cut.get())
                elif view.get() == "simulated 3D":
                    draw_3d(p.fig, tr, mat, ttl + "   simulated",
                            cmap.get(), "fraction of the feed", cut=cut.get())
                else:
                    lim = float(np.nanmax(np.abs(pr - tr))) or 1e-9
                    draw_3d(p.fig, pr - tr, mat, ttl + "   difference",
                            "RdBu_r", "predicted minus simulated",
                            vmin=-lim, vmax=lim, cut=cut.get())
                lab.configure(text=state["labels"][i])
                p.redraw(); self.say("ready")

            c = p.controls
            cb = ttk.Combobox(c, textvariable=which, values=names, width=18,
                              state="readonly")
            cb.pack(side="left", padx=6)
            cb.bind("<<ComboboxSelected>>", reload_rows)
            ttk.Label(c, text="snapshot").pack(side="left", padx=(10, 2))
            sc = ttk.Scale(c, from_=0, to=1, orient="horizontal", length=240,
                           command=lambda v: (rowv.set(int(float(v))), draw()))
            sc.pack(side="left")
            lab = ttk.Label(c, text="", width=30)
            lab.pack(side="left", padx=6)
            cb2 = ttk.Combobox(c, textvariable=view, state="readonly", width=14,
                               values=["compare", "predicted 3D",
                                       "simulated 3D", "difference 3D"])
            cb2.pack(side="left", padx=6); cb2.bind("<<ComboboxSelected>>", draw)
            cb4 = ttk.Combobox(c, textvariable=cut, width=8, state="readonly",
                               values=["corner", "quarter", "half", "none"])
            cb4.pack(side="left", padx=4); cb4.bind("<<ComboboxSelected>>", draw)
            cb3 = ttk.Combobox(c, textvariable=cmap, values=CMAPS, width=11,
                               state="readonly")
            cb3.pack(side="left", padx=4); cb3.bind("<<ComboboxSelected>>", draw)
            reload_rows()

        # ------------------------------- simulations, from a results file
        def tab_simulations_results(self):
            """The same tab as for a dataset, fed from the stored volumes.

            A results file does not carry every run. It carries the simulated
            and predicted volumes of the held out snapshots, which is what the
            models were actually judged on, so that is what this shows.
            """
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Simulations")
            p = Panel(f); p.pack(fill="both", expand=True)
            pairs = d.field_pairs()
            if not pairs:
                p.fig.text(.5, .5,
                           "This file carries no volumes.\n\n"
                           "Simulated fields live under\n"
                           "results/<model>/<split>/fields.",
                           ha="center", va="center")
                p.redraw()
                return
            names = ["%s, %s" % x for x in pairs]
            which = tk.StringVar(value=names[0])
            view = tk.StringVar(value="slices")
            cmap = tk.StringVar(value="viridis")
            cut = tk.StringVar(value="corner")
            rowv = tk.IntVar(value=0)
            st = {}

            def reload_rows(*_):
                t, sp = pairs[names.index(which.get())]
                rid, gid, tn = d.field_rows(t, sp)
                st.update(t=t, sp=sp, rid=rid, gid=gid, tn=tn)
                sc.configure(to=len(rid) - 1)
                last = int(np.where(tn == tn.max())[0][0])
                rowv.set(last); sc.set(last)
                draw()

            def draw(*_):
                if "t" not in st:
                    return
                i = max(0, min(int(rowv.get()), len(st["rid"]) - 1))
                self.say("reading snapshot %d" % i)
                tr, _pr = d.field(st["t"], st["sp"], i)
                mat = d.material_for_gid(int(st["gid"][i]))
                ttl = ("simulated by CompLaB3D\n%s, %s   run %d  gid %d  "
                       "t/t_end %.2f" % (st["t"], st["sp"], st["rid"][i],
                                         st["gid"][i], st["tn"][i]))
                if view.get() == "slices":
                    draw_ortho(p.fig, tr, mat, ttl, cmap.get(),
                               "fraction of the feed")
                else:
                    draw_3d(p.fig, tr, mat, ttl, cmap.get(),
                            "fraction of the feed", cut=cut.get())
                lab.configure(text="run %d  gid %d  t/t_end %.2f"
                              % (st["rid"][i], st["gid"][i], st["tn"][i]))
                p.redraw(); self.say("ready")

            c = p.controls
            ttk.Label(c, text="model, split").pack(side="left")
            cb = ttk.Combobox(c, textvariable=which, values=names, width=18,
                              state="readonly")
            cb.pack(side="left", padx=6)
            cb.bind("<<ComboboxSelected>>", reload_rows)
            ttk.Label(c, text="snapshot").pack(side="left", padx=(10, 2))
            sc = ttk.Scale(c, from_=0, to=1, orient="horizontal", length=220,
                           command=lambda v: (rowv.set(int(float(v))), draw()))
            sc.pack(side="left")
            lab = ttk.Label(c, text="", width=30); lab.pack(side="left", padx=6)
            for v, t in (("slices", "slices"), ("3d", "3D")):
                ttk.Radiobutton(c, text=t, value=v, variable=view,
                                command=draw).pack(side="left", padx=4)
            cb2 = ttk.Combobox(c, textvariable=cut, width=8, state="readonly",
                               values=["corner", "quarter", "half", "none"])
            cb2.pack(side="left", padx=4); cb2.bind("<<ComboboxSelected>>", draw)
            cb3 = ttk.Combobox(c, textvariable=cmap, values=CMAPS, width=11,
                               state="readonly")
            cb3.pack(side="left", padx=4); cb3.bind("<<ComboboxSelected>>", draw)
            reload_rows()

        # ------------------------------- run inputs, from a results file
        def tab_inputs_results(self):
            """What each simulation was, and which side of the split it fell on.

            The dataset file keeps the literal input files of every run. A
            results file keeps the design instead: the two numbers that define
            each case, the pore space it used, how long it took, and which of
            train, validation or test its pore space belongs to.
            """
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Run inputs")
            top = ttk.Frame(f); top.pack(side="top", fill="x", padx=6, pady=4)
            tp = TextPanel(f); tp.pack(fill="both", expand=True)
            cols, tab = d.design_runs()
            models = d.models()
            which = tk.StringVar(value=(models[0] if models else ""))

            def side_of(gid_value, tag):
                out = []
                for sp in ("train", "val", "test"):
                    k = "models/%s/split/%s" % (tag, sp)
                    if k in d.f and gid_value in set(d.f[k][:].tolist()):
                        out.append(sp)
                return ", ".join(out) or "not in any split"

            def show(*_):
                if cols is None:
                    tp.set("This file has no run table at design/runs.")
                    return
                tag = which.get()
                n = len(tab[cols[0]])
                L = ["EVERY SIMULATION IN THE CAMPAIGN",
                     "=" * 33, "",
                     "One row per run. pe and da are the two numbers that "
                     "define a case;",
                     "gid names the pore space it was run on. The last column "
                     "says which",
                     "side of the split that pore space fell on for model "
                     "%s, and it is" % tag,
                     "a property of the pore space, never of the single run.",
                     ""]
                head = ["run", "gid", "pe", "da", "porosity",
                        "t_end s", "wall s", "split for %s" % tag]
                L.append("%-7s %-5s %-9s %-11s %-9s %-10s %-9s %s" % tuple(head))
                L.append("-" * 86)
                for i in range(n):
                    g = int(tab["gid"][i])
                    L.append("%-7d %-5d %-9.4g %-11.4g %-9.3f %-10.1f %-9.1f %s"
                             % (int(tab["run_id"][i]), g,
                                float(tab["pe"][i]), float(tab["da"][i]),
                                float(tab["porosity"][i]),
                                float(tab["t_end_seconds"][i]),
                                float(tab["wall_seconds"][i]),
                                side_of(g, tag)))
                tp.set("\n".join(L))

            ttk.Label(top, text="split shown for model").pack(side="left")
            cb = ttk.Combobox(top, textvariable=which, values=models, width=8,
                              state="readonly")
            cb.pack(side="left", padx=6); cb.bind("<<ComboboxSelected>>", show)
            ttk.Label(top, text="   every run, with the two numbers that "
                                "defined it").pack(side="left")
            show()

        # --------------------------------------------------------- settings
        def tab_settings(self):
            """Every input value, with a sentence saying what it does.

            This is the tab for a reader who wants to know what the run was
            told to do without opening any code. Each row is one setting, the
            value this file was made with, and what changing it would change.
            """
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Settings")
            top = ttk.Frame(f); top.pack(side="top", fill="x", padx=6, pady=4)
            tp = TextPanel(f); tp.pack(fill="both", expand=True)

            if d.kind == KIND_RESULTS:
                models = d.models()
                which = tk.StringVar(value=(models[0] if models else ""))

                def show(*_):
                    tag = which.get()
                    a = d.model_args(tag)
                    if not a:
                        tp.set("No settings recorded for this model.")
                        return
                    ck, why = d.model_checkpoint(tag)
                    lead = ("These are the values model %s was fitted with. "
                            "They are the complete record: running the same "
                            "command with the same campaign file and the same "
                            "seed reproduces this model exactly." % tag)
                    txt = settings_table(sorted(a.items()), PRT_SETTINGS,
                                         "SETTINGS FOR MODEL %s" % tag, lead)
                    if ck:
                        txt += "\n\n[ WHAT THE FITTING ACTUALLY DID ]\n\n"
                        order = ["epochs_run", "best_epoch", "best_val_loss",
                                 "compute_seconds"]
                        says = {
                            "epochs_run": "passes through the training set "
                                          "that were made before stopping",
                            "best_epoch": "the pass whose weights were kept. "
                                          "Every later pass was worse on the "
                                          "validation set",
                            "best_val_loss": "the validation error of that "
                                             "pass, which is the number "
                                             "stopping was decided on",
                            "compute_seconds": "time on the graphics card, in "
                                               "seconds"}
                        for k in order + [x for x in ck if x not in order]:
                            if k in ck:
                                txt += "%-28s  %-26s  %s\n" % (
                                    k, ("%g" % ck[k]), says.get(k, ""))
                        if why:
                            txt += "\nstopped because: %s\n" % why
                    tp.set(txt)

                ttk.Label(top, text="model").pack(side="left")
                cb = ttk.Combobox(top, textvariable=which, values=models,
                                  width=8, state="readonly")
                cb.pack(side="left", padx=6)
                cb.bind("<<ComboboxSelected>>", show)
                ttk.Label(top, text="   every value this model was fitted "
                                    "with, and what each one does"
                          ).pack(side="left")
                show()
            else:
                labels = [d.run_label(s) for s in range(d.n_runs)]
                which = tk.StringVar(value=labels[0] if labels else "")

                def show(*_):
                    try:
                        s_i = labels.index(which.get())
                    except ValueError:
                        s_i = 0
                    pairs = d.dataset_settings(s_i)
                    if not pairs:
                        tp.set("No settings recorded in this file.")
                        return
                    lead = ("These are the values this one simulation was run "
                            "with, read back from the input file the run "
                            "actually used. Every run in the campaign shares "
                            "all of them except the two that define the case.")
                    tp.set(settings_table(pairs, COMPLAB_SETTINGS,
                                          "SETTINGS FOR %s" % labels[s_i],
                                          lead))

                ttk.Label(top, text="run").pack(side="left")
                cb = ttk.Combobox(top, textvariable=which, values=labels,
                                  width=38, state="readonly")
                cb.pack(side="left", padx=6)
                cb.bind("<<ComboboxSelected>>", show)
                ttk.Label(top, text="   every value this run was given, and "
                                    "what each one does").pack(side="left")
                show()

        # -------------------------------------------------------- structure
        def tab_structure(self):
            """The layout the write-up draws, checked against this file."""
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Structure")
            tp = TextPanel(f); tp.pack(fill="both", expand=True)
            try:
                tp.set(structure_text(d))
            except Exception as e:
                tp.set("could not read the layout of this file\n\n%s" % e)

        # ------------------------------------------------------------ notes
        def tab_notes(self):
            """Whatever the reader wants to say, from a plain text file.

            This tab shows no notes of its own. It looks beside the .h5 for a
            text file and prints it, so what appears here is written by whoever
            is using the file rather than by whoever made it.
            """
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Notes")
            top = ttk.Frame(f); top.pack(side="top", fill="x", padx=6, pady=4)
            tp = TextPanel(f, mono=True); tp.pack(fill="both", expand=True)

            # Where to look, nearest first. The first two sit beside the .h5
            # and belong to that one file; the third is a folder note shared by
            # everything in the folder. The last two repeat that beside this
            # program, for the case where the .h5 is somewhere else entirely.
            stem = os.path.splitext(os.path.abspath(d.path))[0]
            here = os.path.dirname(os.path.abspath(d.path))
            mine = os.path.dirname(os.path.abspath(__file__))
            candidates = [stem + ".txt", stem + "_notes.txt",
                          os.path.join(here, "notes.txt"),
                          os.path.join(here, "NOTES.txt")]
            if os.path.normcase(mine) != os.path.normcase(here):
                candidates += [os.path.join(mine, "notes.txt"),
                               os.path.join(mine, "NOTES.txt")]
            seen, uniq = set(), []
            for c in candidates:
                k = os.path.normcase(c)
                if k not in seen:
                    seen.add(k); uniq.append(c)
            candidates = uniq
            current = {"path": None}

            def load(path=None):
                p = path
                if p is None:
                    p = next((c for c in candidates if os.path.isfile(c)), None)
                if p is None:
                    current["path"] = None
                    lbl.configure(text="no notes file yet")
                    tp.set(
                        "This tab prints a plain text file, so that anything "
                        "worth saying about\nthis .h5 can be written once and "
                        "read here by everyone who opens it.\n\n"
                        "Nothing is here yet. Make a text file at any of these "
                        "names and it will\nbe picked up, nearest first:\n\n"
                        + "\n".join("    " + c for c in candidates)
                        + "\n\nor press  Open a text file  and point at one "
                          "anywhere.\n\nPress  New  to create the first one and "
                          "start typing.\nPress  Reload  after editing it "
                          "outside this window.")
                    return
                try:
                    with open(p, "r", encoding="utf-8", errors="replace") as fh:
                        body = fh.read()
                except Exception as e:
                    tp.set("could not read %s\n\n%s" % (p, e)); return
                current["path"] = p
                lbl.configure(text=p)
                tp.set(body if body.strip() else
                       "%s is empty.\n\nWrite anything in it and press Reload."
                       % p)

            def choose():
                p = filedialog.askopenfilename(
                    title="Open a text file",
                    initialdir=here,
                    filetypes=[("text", "*.txt *.md"), ("all files", "*.*")])
                if p:
                    load(p)

            def create():
                p = current["path"] or candidates[0]
                if not os.path.isfile(p):
                    try:
                        with open(p, "w", encoding="utf-8") as fh:
                            fh.write("Notes on %s\n\n"
                                     % os.path.basename(d.path))
                    except Exception as e:
                        messagebox.showerror("Could not create the file",
                                             "%s\n\n%s" % (p, e))
                        return
                load(p)
                messagebox.showinfo(
                    "Notes file",
                    "Created\n\n%s\n\nOpen it in any text editor, write what "
                    "you want, then press Reload." % p)

            ttk.Button(top, text="Reload",
                       command=lambda: load(current["path"])).pack(side="left")
            ttk.Button(top, text="Open a text file",
                       command=choose).pack(side="left", padx=6)
            ttk.Button(top, text="New", command=create).pack(side="left")
            lbl = ttk.Label(top, text="", foreground="#555555")
            lbl.pack(side="left", padx=12)
            load()

        # ---------------------------------------------------------- figures
        def tab_figures(self):
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Figures")
            p = Panel(f); p.pack(fill="both", expand=True)
            names = d.figures()
            which = tk.StringVar(value=names[0])

            def draw(*_):
                import matplotlib.image as mpimg
                blob = d.f["figures"][which.get()][:].tobytes()
                img = mpimg.imread(io.BytesIO(blob), format="png")
                p.fig.clear()
                ax = p.fig.add_subplot(111)
                ax.imshow(img); ax.axis("off")
                ax.set_title(which.get(), fontsize=9)
                p.redraw()

            cb = ttk.Combobox(p.controls, textvariable=which, values=names,
                              width=42, state="readonly")
            cb.pack(side="left", padx=6); cb.bind("<<ComboboxSelected>>", draw)
            draw()

        # ----------------------------------------------------------- browse
        def tab_browse(self):
            d = self.doc
            f = ttk.Frame(self.nb); self.nb.add(f, text="Browse")
            pane = ttk.PanedWindow(f, orient="horizontal")
            pane.pack(fill="both", expand=True)
            left = ttk.Frame(pane); right = ttk.Frame(pane)
            pane.add(left, weight=1); pane.add(right, weight=2)

            tree = ttk.Treeview(left, columns=("shape",), height=28)
            tree.heading("#0", text="path")
            tree.heading("shape", text="shape and type")
            tree.column("shape", width=170)
            ys = ttk.Scrollbar(left, orient="vertical", command=tree.yview)
            tree.configure(yscrollcommand=ys.set)
            tree.pack(side="left", fill="both", expand=True)
            ys.pack(side="right", fill="y")
            tp = TextPanel(right); tp.pack(fill="both", expand=True)

            def add(parent, g, path):
                for k in g:
                    o = g[k]
                    full = "%s/%s" % (path, k) if path else k
                    if isinstance(o, h5py.Group):
                        n = tree.insert(parent, "end", iid=full, text=k,
                                        values=("group",))
                        add(n, o, full)
                    else:
                        tree.insert(parent, "end", iid=full, text=k,
                                    values=("%s %s" % (o.shape, o.dtype),))
            add("", d.f, "")

            def show(_evt=None):
                sel = tree.focus()
                if not sel or sel not in d.f:
                    return
                o = d.f[sel]
                L = ["PATH   /%s" % sel]
                for k, v in o.attrs.items():
                    s = _s(v)
                    if isinstance(s, str) and len(s) > 2000:
                        s = s[:2000] + " ..."
                    L.append("  @%-22s %s" % (k, s))
                if isinstance(o, h5py.Dataset):
                    L += ["", "SHAPE  %s   %s" % (o.shape, o.dtype)]
                    try:
                        if o.shape == ():
                            L.append("VALUE  %s" % _s(o[()]))
                        elif o.dtype.kind == "S":
                            L += ["", _s(o[0]) if o.shape[0] else ""]
                        elif o.size <= 4096:
                            L += ["", np.array2string(o[:], threshold=4096,
                                                      max_line_width=110)]
                        else:
                            a = np.asarray(o[tuple(slice(0, min(s, 4))
                                                   for s in o.shape)])
                            L += ["", "first corner:",
                                  np.array2string(a, max_line_width=110),
                                  "",
                                  "too large to print whole, %d values" % o.size]
                    except Exception as e:
                        L.append("could not read: %s" % e)
                tp.set("\n".join(L))

            tree.bind("<<TreeviewSelect>>", show)
            # Arrive on something rather than on a blank pane.
            tp.set(d.header() + "\n\nSelect anything on the left to see its "
                                "shape, its attributes and its values.")
            first = tree.get_children("")
            if first:
                tree.selection_set(first[0]); tree.focus(first[0])

    app = App()
    if initial:
        app.after(80, lambda: app.open(initial))
    if on_ready is not None:          # used by the screenshot harness
        app.after(400, lambda: on_ready(app))
    app.mainloop()
    return app


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    path = argv[0] if argv else None
    if path and not os.path.isfile(path):
        sys.exit("no such file: %s" % path)
    try:
        _run_gui(path)
    except Exception:
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
