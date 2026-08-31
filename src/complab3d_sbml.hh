/* This file is a part of the CompLaB program.
 *
 * The CompLaB3D software (3-D pore-scale extension) is developed since 2024
 * by the University of Georgia (Meile Lab, Department of Marine Sciences).
 *
 * Contact: Shahram Asgari, Christof Meile
 *          Department of Marine Sciences, University of Georgia
 *          shahram.asgari@uga.edu
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* ================================================================================================
 * complab3d_sbml.hh  --  READ A GENOME-SCALE MODEL DIRECTLY FROM SBML
 * ================================================================================================
 *
 *  WHAT THIS REPLACES
 *
 *  Until now, using a genome-scale model meant running extractMM.py first: a separate Python step,
 *  needing cobrapy, producing an intermediate matrix file that CompLB3D then read. That put a
 *  Python dependency in front of the GLPK path, which otherwise needs no Python at all, and it put
 *  a manual step between the user and a run.
 *
 *  This reads SBML Level 3 with the FBC package -- the format BiGG, ModelSEED and cobrapy all
 *  write -- straight into the S, b, c, lb, ub arrays the solver wants.  No Python, no intermediate
 *  file, no extra dependency: TinyXML is already compiled into every CompLB3D build, because
 *  CMakeLists.txt globs Palabos's externalLibraries/tinyxml.
 *
 *      <model_filename>iAF987</model_filename>
 *
 *  now finds iAF987.xml and works out for itself whether it is an SBML model or an extractMM
 *  matrix file, so existing input files keep working unchanged.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT IT READS, AND WHAT IT DELIBERATELY IGNORES
 *
 *  Read:
 *      listOfSpecies                     the metabolites, in document order
 *      listOfReactions                   the reactions, in document order
 *        listOfReactants / listOfProducts  stoichiometry, reactants negative
 *        fbc:lowerFluxBound / upperFluxBound   by parameter id, or as a literal number
 *      listOfParameters                  the bound values those ids point at
 *      fbc:listOfObjectives              the ACTIVE objective and its reaction
 *
 *  Ignored on purpose:
 *      gene product associations, annotations, notes, units, compartments, groups.
 *      A flux balance problem does not need them, and parsing them would triple the code for no
 *      numerical difference.
 *
 *  ORDER MATTERS AND IS PRESERVED.  <exchange_reaction_indices> in CompLaB.xml is positional, so
 *  the reaction order this reader produces must match what cobrapy would produce from the same
 *  file.  Both use document order, and that equivalence is checked against cobrapy in the test
 *  suite rather than assumed.
 *
 *  NAME PREFIXES.  SBML requires ids to start with a letter, so cobrapy writes metabolites as
 *  M_<id> and reactions as R_<id> and strips those prefixes when reading back.  This reader strips
 *  them too, so the names printed here match the ones in BiGG and in extractMM.py's index table.
 *
 *  ------------------------------------------------------------------------------------------------
 *  KNOWN LIMITS, stated plainly
 *
 *    * SBML Level 3 + FBC v2 only.  Level 2 models with COBRA-style <notes> bounds are not read;
 *      convert them once with cobrapy, or keep using extractMM.py.
 *    * .mat and .json models are not read here.  extractMM.py still handles those.
 *    * No unit checking.  FBC models are conventionally mmol/gDW/h and this assumes it, exactly as
 *      extractMM.py and cobrapy do.
 *    * Reversibility is taken from the bounds, not from the reversible="" attribute, because the
 *      two disagree in many published models and the bounds are what the LP actually uses.
 * ================================================================================================
 */
#ifndef COMPLAB3D_SBML_HH
#define COMPLAB3D_SBML_HH

#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "tinyxml.h"

namespace complab_sbml {

/* ------------------------------------------------------------------------------------------------
 *  The model, in exactly the shape load_metabolic_models3D() wants it.
 * ------------------------------------------------------------------------------------------------ */
struct SbmlModel {
    int nmet;
    int nrxn;
    int objLoc;                          // index of the objective reaction
    int objSense;                        // -1 maximize, +1 minimize

    std::vector<double> S;               // nmet*nrxn, row major, metabolite major
    std::vector<double> b;               // nmet, all zero for steady-state FBA
    std::vector<double> c;               // nrxn, objective coefficients
    std::vector<double> lb;              // nrxn
    std::vector<double> ub;              // nrxn

    std::vector<std::string> metNames;   // nmet
    std::vector<std::string> rxnNames;   // nrxn

    std::string modelId;
    std::string error;                   // empty on success

    SbmlModel() : nmet(0), nrxn(0), objLoc(0), objSense(-1) {}
    bool ok() const { return error.empty(); }
};

/* ------------------------------------------------------------------------------------------------ */
inline std::string stripPrefix(const std::string &id)
{
    if (id.size() > 2 && (id.compare(0, 2, "M_") == 0 || id.compare(0, 2, "R_") == 0))
        return id.substr(2);
    return id;
}

inline const char *attr(TiXmlElement *e, const char *name)
{
    const char *v = e->Attribute(name);
    return v ? v : 0;
}

/* Find a child element by local name, ignoring any namespace prefix.
 *
 * SBML files in the wild write <listOfReactants> and <sbml:listOfReactants> interchangeably, and
 * TinyXML does not resolve namespaces -- it hands back the raw tag text. Matching on the part after
 * the colon is what makes this reader work on files from more than one producer. */
inline bool localNameIs(const TiXmlElement *e, const char *want)
{
    std::string t = e->Value();
    const size_t colon = t.rfind(':');
    if (colon != std::string::npos) t = t.substr(colon + 1);
    return t == want;
}

inline TiXmlElement *firstChild(TiXmlElement *parent, const char *localName)
{
    if (!parent) return 0;
    for (TiXmlElement *e = parent->FirstChildElement(); e; e = e->NextSiblingElement())
        if (localNameIs(e, localName)) return e;
    return 0;
}

/* Attributes carry namespace prefixes too (fbc:lowerFluxBound). Try the prefixed spelling first,
 * because that is what conforming files use, then the bare one. */
inline const char *fbcAttr(TiXmlElement *e, const char *bare)
{
    std::string pref = std::string("fbc:") + bare;
    const char *v = e->Attribute(pref.c_str());
    if (v) return v;
    return e->Attribute(bare);
}

/* ------------------------------------------------------------------------------------------------
 *  Is this file SBML, or an extractMM matrix file?
 *
 *  Cheap and definite: look at the root element. <sbml> means SBML; <Metabolic_Model> means the
 *  matrix format. Anything else is neither, and saying so early is far better than failing later
 *  with a confusing message about array lengths.
 * ------------------------------------------------------------------------------------------------ */
enum ModelFileKind { MODEL_UNKNOWN = 0, MODEL_SBML = 1, MODEL_MATRIX = 2 };

inline ModelFileKind classify(const std::string &path)
{
    TiXmlDocument doc(path.c_str());
    if (!doc.LoadFile()) return MODEL_UNKNOWN;
    TiXmlElement *root = doc.RootElement();
    if (!root) return MODEL_UNKNOWN;
    if (localNameIs(root, "sbml")) return MODEL_SBML;
    if (localNameIs(root, "Metabolic_Model")) return MODEL_MATRIX;
    return MODEL_UNKNOWN;
}

/* ------------------------------------------------------------------------------------------------
 *  The reader.
 * ------------------------------------------------------------------------------------------------ */
inline SbmlModel readSbml(const std::string &path)
{
    SbmlModel M;

    TiXmlDocument doc(path.c_str());
    if (!doc.LoadFile()) {
        M.error = "cannot open or parse '" + path + "': " + doc.ErrorDesc();
        return M;
    }
    TiXmlElement *sbml = doc.RootElement();
    if (!sbml || !localNameIs(sbml, "sbml")) {
        M.error = "'" + path + "' has no <sbml> root element";
        return M;
    }
    TiXmlElement *model = firstChild(sbml, "model");
    if (!model) {
        M.error = "'" + path + "' has no <model> element";
        return M;
    }
    if (const char *id = attr(model, "id")) M.modelId = id;

    /* ---- parameters, so bound ids can be resolved to numbers ---- */
    std::map<std::string, double> param;
    if (TiXmlElement *lop = firstChild(model, "listOfParameters")) {
        for (TiXmlElement *p = lop->FirstChildElement(); p; p = p->NextSiblingElement()) {
            if (!localNameIs(p, "parameter")) continue;
            const char *pid = attr(p, "id");
            const char *val = attr(p, "value");
            if (pid && val) param[pid] = std::atof(val);
        }
    }

    /* ---- metabolites, in document order ---- */
    std::map<std::string, int> metIndex;
    if (TiXmlElement *los = firstChild(model, "listOfSpecies")) {
        for (TiXmlElement *s = los->FirstChildElement(); s; s = s->NextSiblingElement()) {
            if (!localNameIs(s, "species")) continue;
            const char *sid = attr(s, "id");
            if (!sid) continue;
            /* boundaryCondition species are sinks the modeller has already decided are infinite.
             * cobrapy drops them from the stoichiometric matrix, and so must we, or every exchange
             * reaction acquires a spurious row that makes the LP infeasible. */
            const char *bc = attr(s, "boundaryCondition");
            if (bc && (std::string(bc) == "true" || std::string(bc) == "1")) continue;
            metIndex[sid] = (int) M.metNames.size();
            M.metNames.push_back(stripPrefix(sid));
        }
    }
    M.nmet = (int) M.metNames.size();
    if (M.nmet == 0) {
        M.error = "'" + path + "' declares no non-boundary species";
        return M;
    }

    /* ---- reactions, in document order ---- */
    TiXmlElement *lor = firstChild(model, "listOfReactions");
    if (!lor) {
        M.error = "'" + path + "' has no <listOfReactions>";
        return M;
    }

    std::vector<TiXmlElement *> rxns;
    for (TiXmlElement *r = lor->FirstChildElement(); r; r = r->NextSiblingElement())
        if (localNameIs(r, "reaction")) rxns.push_back(r);
    M.nrxn = (int) rxns.size();
    if (M.nrxn == 0) {
        M.error = "'" + path + "' declares no reactions";
        return M;
    }

    M.S.assign((size_t) M.nmet * (size_t) M.nrxn, 0.0);
    M.b.assign(M.nmet, 0.0);
    M.c.assign(M.nrxn, 0.0);
    M.lb.assign(M.nrxn, 0.0);
    M.ub.assign(M.nrxn, 0.0);
    M.rxnNames.resize(M.nrxn);

    std::map<std::string, int> rxnIndex;

    for (int j = 0; j < M.nrxn; ++j) {
        TiXmlElement *r = rxns[j];
        const char *rid = attr(r, "id");
        M.rxnNames[j] = rid ? stripPrefix(rid) : "";
        if (rid) rxnIndex[rid] = j;

        /* bounds: an fbc attribute naming a parameter, or a literal number */
        const char *lbs = fbcAttr(r, "lowerFluxBound");
        const char *ubs = fbcAttr(r, "upperFluxBound");
        double lo = -1000.0, hi = 1000.0;
        if (lbs) {
            std::map<std::string, double>::const_iterator it = param.find(lbs);
            lo = (it != param.end()) ? it->second : std::atof(lbs);
        }
        if (ubs) {
            std::map<std::string, double>::const_iterator it = param.find(ubs);
            hi = (it != param.end()) ? it->second : std::atof(ubs);
        }
        M.lb[j] = lo;
        M.ub[j] = hi;

        /* stoichiometry.  Reactants negative, products positive, and ACCUMULATED rather than
         * assigned: a species may legitimately appear on both sides, and overwriting would silently
         * drop one of the two. */
        const char *sides[2] = { "listOfReactants", "listOfProducts" };
        const double sign[2] = { -1.0, +1.0 };
        for (int side = 0; side < 2; ++side) {
            TiXmlElement *lst = firstChild(r, sides[side]);
            if (!lst) continue;
            for (TiXmlElement *sr = lst->FirstChildElement(); sr; sr = sr->NextSiblingElement()) {
                if (!localNameIs(sr, "speciesReference")) continue;
                const char *sp = attr(sr, "species");
                if (!sp) continue;
                std::map<std::string, int>::const_iterator it = metIndex.find(sp);
                if (it == metIndex.end()) continue;      // a boundary species: correctly skipped
                const char *st = attr(sr, "stoichiometry");
                const double coeff = st ? std::atof(st) : 1.0;
                M.S[(size_t) it->second * (size_t) M.nrxn + (size_t) j] += sign[side] * coeff;
            }
        }
    }

    /* ---- objective ---- */
    TiXmlElement *loo = firstChild(model, "listOfObjectives");
    bool haveObj = false;
    if (loo) {
        const char *active = fbcAttr(loo, "activeObjective");
        for (TiXmlElement *o = loo->FirstChildElement(); o; o = o->NextSiblingElement()) {
            if (!localNameIs(o, "objective")) continue;
            const char *oid = fbcAttr(o, "id");
            if (active && oid && std::string(active) != std::string(oid)) continue;

            const char *type = fbcAttr(o, "type");
            M.objSense = (type && std::string(type) == "minimize") ? +1 : -1;

            TiXmlElement *lfo = firstChild(o, "listOfFluxObjectives");
            if (!lfo) continue;
            for (TiXmlElement *fo = lfo->FirstChildElement(); fo; fo = fo->NextSiblingElement()) {
                if (!localNameIs(fo, "fluxObjective")) continue;
                const char *rn = fbcAttr(fo, "reaction");
                const char *cf = fbcAttr(fo, "coefficient");
                if (!rn) continue;
                std::map<std::string, int>::const_iterator it = rxnIndex.find(rn);
                if (it == rxnIndex.end()) continue;
                M.c[it->second] = cf ? std::atof(cf) : 1.0;
                if (!haveObj) { M.objLoc = it->second; haveObj = true; }
            }
            break;
        }
    }
    if (!haveObj) {
        M.error = "'" + path + "' has no usable <fbc:listOfObjectives>. A flux balance problem "
                  "needs an objective; set one in the model, or use extractMM.py with an "
                  "explicit objective reaction.";
        return M;
    }

    return M;
}

/* ------------------------------------------------------------------------------------------------
 *  Which reactions are exchanges?
 *
 *  <exchange_reaction_indices> in CompLaB.xml is the one thing a user cannot look up anywhere else,
 *  so the loader prints this table at start-up. An exchange reaction touches exactly one metabolite:
 *  that is the structural definition and it does not depend on a naming convention, which matters
 *  because ModelSEED calls them EX_cpd..._e0 and BiGG calls them EX_..._e.
 * ------------------------------------------------------------------------------------------------ */
enum BoundaryKind { NOT_BOUNDARY = 0, BK_EXCHANGE = 1, BK_DEMAND = 2, BK_SINK = 3 };

inline const char *boundaryKindName(BoundaryKind k)
{
    switch (k) {
        case BK_EXCHANGE: return "exchange";
        case BK_DEMAND:   return "demand";
        case BK_SINK:     return "sink";
        default:          return "";
    }
}

/* Boundary reactions touch exactly one metabolite. That is the structural definition and it does
 * not depend on a naming convention, which matters because BiGG writes EX_glc__D_e and ModelSEED
 * writes EX_cpd00027_e0.
 *
 * cobrapy then SUBDIVIDES those into exchanges, demands and sinks by id prefix, and the split is
 * worth reproducing because it is what users see in BiGG. On iJO1366, for instance, 330 reactions
 * are structurally boundary: cobrapy calls 324 of them exchanges and 6 of them demands (the DM_
 * ones). Neither answer is wrong; they answer different questions.
 *
 * For <exchange_reaction_indices> what matters is "which reaction moves this metabolite across the
 * boundary", so ALL of them are candidates. The label is printed so you can tell which is which. */
inline BoundaryKind boundaryKind(const SbmlModel &M, int j)
{
    int touched = 0;
    for (int i = 0; i < M.nmet && touched < 2; ++i)
        if (M.S[(size_t) i * (size_t) M.nrxn + (size_t) j] != 0.0) ++touched;
    if (touched != 1) return NOT_BOUNDARY;

    const std::string &id = M.rxnNames[(size_t) j];
    if (id.compare(0, 3, "DM_") == 0) return BK_DEMAND;
    if (id.compare(0, 3, "SK_") == 0) return BK_SINK;
    return BK_EXCHANGE;
}

inline std::vector<int> exchangeReactions(const SbmlModel &M)
{
    std::vector<int> out;
    for (int j = 0; j < M.nrxn; ++j)
        if (boundaryKind(M, j) != NOT_BOUNDARY) out.push_back(j);
    return out;
}

/* Just the ones cobrapy would call exchanges, for when you want to compare like with like. */
inline std::vector<int> exchangeReactionsStrict(const SbmlModel &M)
{
    std::vector<int> out;
    for (int j = 0; j < M.nrxn; ++j)
        if (boundaryKind(M, j) == BK_EXCHANGE) out.push_back(j);
    return out;
}


/* ------------------------------------------------------------------------------------------------
 *  NAME EXCHANGES INSTEAD OF NUMBERING THEM
 *
 *  <exchange_reaction_indices>6 55 418 -1</exchange_reaction_indices> is positional. It depends on
 *  the model file having exactly the reaction order it had the day the number was copied. Publish a
 *  revision that inserts one reaction, or switch from BiGG to a ModelSEED build of the same
 *  organism, and every index after the insertion now points at a different reaction -- and the run
 *  proceeds, on the wrong exchanges, with no error at all. It is the worst kind of bug: silent, and
 *  it produces plausible numbers.
 *
 *      <exchange_reaction_names>EX_ac_e EX_fe3_e none</exchange_reaction_names>
 *
 *  cannot fail that way. A name either resolves or it does not, and if it does not the run stops
 *  with the name printed. Use "none", "-1" or "-" for a substrate the organism does not exchange.
 *
 *  Matching is exact first, then case-insensitive, then against the id with the compartment suffix
 *  removed -- so EX_ac_e finds EX_ac_e0 in a ModelSEED model. Anything looser would risk resolving
 *  to the wrong reaction quietly, which is the thing this exists to prevent.
 * ------------------------------------------------------------------------------------------------ */
inline std::string lowerCase(std::string s)
{
    for (size_t i = 0; i < s.size(); ++i)
        if (s[i] >= 'A' && s[i] <= 'Z') s[i] = (char) (s[i] - 'A' + 'a');
    return s;
}

/* "EX_ac_e0" and "EX_ac_e" both reduce to "ex_ac". Used only as a last resort. */
inline std::string reactionStem(const std::string &id)
{
    std::string s = lowerCase(id);
    size_t end = s.size();
    while (end > 0 && s[end - 1] >= '0' && s[end - 1] <= '9') --end;      // trailing digits
    if (end >= 2 && s[end - 1] == 'e' && s[end - 2] == '_') end -= 2;     // ..._e
    else if (end >= 2 && s[end - 2] == '_') end -= 2;                     // any one-letter suffix
    return s.substr(0, end);
}

inline int findReaction(const SbmlModel &M, const std::string &name)
{
    for (int j = 0; j < M.nrxn; ++j) if (M.rxnNames[(size_t) j] == name) return j;

    const std::string lc = lowerCase(name);
    for (int j = 0; j < M.nrxn; ++j) if (lowerCase(M.rxnNames[(size_t) j]) == lc) return j;

    const std::string stem = reactionStem(name);
    int hit = -1, hits = 0;
    for (int j = 0; j < M.nrxn; ++j)
        if (reactionStem(M.rxnNames[(size_t) j]) == stem) { hit = j; ++hits; }
    return (hits == 1) ? hit : -1;      // ambiguous is as good as not found
}

/* Resolve a whole <exchange_reaction_names> list. Returns false and fills `err` on the first name
 * that cannot be resolved, naming the closest candidates so the user has somewhere to go. */
inline bool resolveExchangeNames(const SbmlModel &M, const std::vector<std::string> &names,
                                 std::vector<int> &out, std::string &err)
{
    out.assign(names.size(), -1);
    for (size_t s = 0; s < names.size(); ++s) {
        const std::string &n = names[s];
        if (n == "none" || n == "-1" || n == "-" || n.empty()) continue;

        const int j = findReaction(M, n);
        if (j >= 0) { out[s] = j; continue; }

        err = "substrate " ;
        char b[32]; std::snprintf(b, sizeof(b), "%d", (int) s); err += b;
        err += ": no reaction called '" + n + "' in model '" + M.modelId + "'.\n";

        /* Suggest boundary reactions whose stem shares a prefix. A list of 2583 reaction names is
         * no help; three plausible ones is. */
        const std::string stem = reactionStem(n);
        std::string sug;
        int shown = 0;
        for (int j2 = 0; j2 < M.nrxn && shown < 6; ++j2) {
            if (boundaryKind(M, j2) == NOT_BOUNDARY) continue;
            const std::string s2 = reactionStem(M.rxnNames[(size_t) j2]);
            if (stem.size() >= 3 && s2.compare(0, stem.size() < s2.size() ? stem.size() : s2.size(),
                                               stem, 0, stem.size() < s2.size() ? stem.size() : s2.size()) == 0) {
                sug += "    " + M.rxnNames[(size_t) j2] + "\n";
                ++shown;
            }
        }
        if (!sug.empty()) err += "  Did you mean one of:\n" + sug;
        else err += "  Use 'none' if this organism does not exchange that substrate.\n";
        return false;
    }
    return true;
}

/* ------------------------------------------------------------------------------------------------
 *  A few structural checks worth making once, at start-up, rather than discovering per voxel.
 * ------------------------------------------------------------------------------------------------ */
inline std::string sanityReport(const SbmlModel &M)
{
    char buf[512];
    std::string out;

    size_t nz = 0;
    for (size_t k = 0; k < M.S.size(); ++k) if (M.S[k] != 0.0) ++nz;
    std::snprintf(buf, sizeof(buf),
                  "    %d metabolites, %d reactions, %lu nonzeros, objective '%s' (column %d, %s)\n",
                  M.nmet, M.nrxn, (unsigned long) nz,
                  M.objLoc < (int) M.rxnNames.size() ? M.rxnNames[M.objLoc].c_str() : "?",
                  M.objLoc, M.objSense < 0 ? "maximize" : "minimize");
    out += buf;

    int emptyRxn = 0, emptyMet = 0, badBound = 0;
    for (int j = 0; j < M.nrxn; ++j) {
        bool any = false;
        for (int i = 0; i < M.nmet && !any; ++i)
            if (M.S[(size_t) i * (size_t) M.nrxn + (size_t) j] != 0.0) any = true;
        if (!any) ++emptyRxn;
        if (M.lb[j] > M.ub[j]) ++badBound;
    }
    for (int i = 0; i < M.nmet; ++i) {
        bool any = false;
        for (int j = 0; j < M.nrxn && !any; ++j)
            if (M.S[(size_t) i * (size_t) M.nrxn + (size_t) j] != 0.0) any = true;
        if (!any) ++emptyMet;
    }
    if (emptyRxn) {
        std::snprintf(buf, sizeof(buf),
                      "    NOTE: %d reaction(s) touch no metabolite. Harmless, but they cannot "
                      "carry flux.\n", emptyRxn);
        out += buf;
    }
    if (emptyMet) {
        std::snprintf(buf, sizeof(buf),
                      "    NOTE: %d metabolite(s) appear in no reaction. Harmless; they make the "
                      "LP slightly larger than it needs to be.\n", emptyMet);
        out += buf;
    }
    if (badBound) {
        std::snprintf(buf, sizeof(buf),
                      "    WARNING: %d reaction(s) have lower bound above upper bound. The LP will "
                      "be infeasible.\n", badBound);
        out += buf;
    }
    return out;
}

}  // namespace complab_sbml

#endif  // COMPLAB3D_SBML_HH
