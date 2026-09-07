/* This file is part of the CompLB3D library.
 *
 * CompLB3D is developed since 2022 by the University of Georgia (United States)
 * and Chungnam National University (South Korea).
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more
 * details.
 * ============================================================================
 *
 *  complab3d_metabolic.hh
 *  --------------------------------------------------------------------------
 *  CONFIGURATION AND MODEL LOADING FOR THE OPTIONAL METABOLIC LAYER
 *
 *  CompLB3D shipped with hand-written Monod kinetics only.  This file, together
 *  with complab3d_processors_fba.hh and complab3d_processors_surrogate.hh, adds
 *  back the three metabolic solvers that CompLaB 2D v1.0.0 had:
 *
 *      * flux balance analysis through GLPK      (complab3d_glpkcpp.hh)
 *      * flux balance analysis through COBRApy   (complab3d_pythonAPI.hh)
 *      * a pre-trained surrogate model           (../surrogateModel.hh)
 *
 *  EVERYTHING HERE IS OPTIONAL.  Nothing in this file runs unless the user asks
 *  for it in CompLaB.xml.  A run with no <enable_fba_glpk>, no
 *  <enable_fba_cobrapy> and no <enable_surrogate> behaves exactly like the
 *  CompLB3D that came before this file existed.
 *
 *  There are TWO independent switches, and it helps to know why:
 *
 *    1. A BUILD switch, in CMakeLists.txt:  ENABLE_GLPK / ENABLE_COBRAPY.
 *       This decides whether the GLPK and Python code is put into the program
 *       at all when it is compiled.  If the machine you compile on does not
 *       have GLPK or Python installed, leave these OFF and CompLB3D still
 *       builds and runs -- it simply will not offer those two solvers.
 *
 *    2. A RUN switch, in CompLaB.xml:  <enable_fba_glpk> etc.
 *       This decides whether a program that WAS built with them actually uses
 *       them on this particular run.  This is the one you change day to day.
 *
 *       If you ask for a solver in the XML that was not built in, CompLB3D
 *       stops with a clear message telling you which CMake option to turn on,
 *       rather than silently ignoring you.
 *
 *  --------------------------------------------------------------------------
 *  UNITS -- read this before touching anything
 *
 *  CompLB3D keeps EVERY concentration and EVERY biomass pool in mol/L.  That is
 *  the convention defineKinetics.hh uses and the convention the 95-species
 *  equilibrium tableau uses.  Flux balance analysis, however, is universally
 *  posed in mmol per gram dry weight per hour, against a biomass in gDW/L.
 *
 *  We do NOT change CompLB3D's convention.  Instead every microbe that uses an
 *  FBA or surrogate solver carries one extra number in CompLaB.xml:
 *
 *      <biomass_molar_mass> 24.6 </biomass_molar_mass>     grams dry weight per
 *                                                          mole of cells
 *
 *  and the conversion happens at exactly one place, on the way into and out of
 *  the LP:
 *
 *      B[gDW/L]  =  B[mol/L] * biomass_molar_mass
 *
 *      dC[mol/L] =  v[mmol/gDW/h] * B[gDW/L] * dt[s] / 3.6e6
 *                       ( /1000 mmol->mol , /3600 h->s )
 *
 *      dB[mol/L] =  mu[1/h] * B[mol/L] * dt[s] / 3600
 *
 *  The helper functions fluxToDeltaC() and growthToDeltaB() below are the ONLY
 *  places these factors appear.  If a unit ever looks wrong, look there first
 *  and nowhere else.
 *
 *  The default 24.6 g/mol is the usual CH1.8O0.5N0.2 biomass formula weight.
 *  It is a default, not a truth; set it per organism if you know better.
 *
 *  --------------------------------------------------------------------------
 *  FREE versus TOTAL concentration
 *
 *  When the equilibrium solver is switched on, the substrate lattice for a
 *  component such as acetate holds only the FREE ion -- the rest of the
 *  acetate is sitting in CaAc+, FeAc+, MgAc and so on, each of which is its own
 *  substrate.  Feeding the free concentration straight to a Michaelis-Menten
 *  uptake bound therefore understates what the organism can actually reach,
 *  sometimes by a large factor.
 *
 *      <fba_concentration_basis> free </fba_concentration_basis>   (default)
 *      <fba_concentration_basis> total </fba_concentration_basis>
 *
 *  With "total", MetabolicTotals below rebuilds each component's total from the
 *  equilibrium tableau -- generically, by reading the stoichiometry matrix,
 *  rather than by hard-coding species indices the way defineKinetics.hh does --
 *  and consumption is drawn back down across the complexes in proportion to how
 *  much of the element each one holds, which keeps the element balance exact.
 *
 *  With equilibrium switched off the two settings are identical, because every
 *  substrate is then its own total.
 *
 *  "free" is the default only because it reproduces CompLaB 2D exactly.  If you
 *  are running with speciation on, "total" is almost certainly the physically
 *  correct choice.
 * ============================================================================
 */

#ifndef COMPLAB3D_METABOLIC_HH
#define COMPLAB3D_METABOLIC_HH

#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdlib>

/* [NEW] The native SBML reader.  Unconditional: it is header-only and depends only on tinyxml,
 * which Palabos already compiles into every build (externalLibraries/tinyxml).  Including it here
 * is what lets <model_filename> point straight at a genome-scale model instead of at the output of
 * extractMM.py. */
#include "complab3d_configpath.hh"   /* complab_input::configPath(), set from argv in main() */
#include "complab3d_sbml.hh"

#ifdef COMPLAB_ENABLE_GLPK
#include "complab3d_glpkcpp.hh"
#endif

/* The two refinements of the linear program itself.  Both are header-only and neither needs a
 * solver to be present -- run_glpk_lex() is the only part of either that touches GLPK, and it is
 * guarded inside its own file -- so they are included unconditionally and stay parseable in a
 * build with no metabolic back end at all. */
#include "complab3d_lexicographic.hh"
#include "complab3d_cybernetic.hh"

#ifdef COMPLAB_ENABLE_COBRAPY
#include "complab3d_pythonAPI.hh"
#endif


/* ============================================================================
 *  REACTION TYPE VOCABULARY
 *
 *  CompLaB 2D used bare integers scattered through the source (0 meant glpk, 1
 *  meant kinetics, 5 meant cobrapy ...) which made every comparison a puzzle.
 *  Here they are named, and there are predicates instead of magic numbers.
 *
 *  IMPORTANT -- the numbering is NOT the 2D numbering.  CompLB3D already
 *  shipped with 0 = none and 1 = kinetics, and existing input files rely on it,
 *  so those two keep their meaning and the new solvers are appended.  Never
 *  renumber 0 or 1.
 * ============================================================================
 */
namespace rxntype {
    enum Type {
        NONE       = 0,   // this microbe does not react at all (it can still be transported)
        KINETICS   = 1,   // defineKinetics.hh                       <-- CompLB3D original
        GLPK       = 2,   // flux balance analysis via GLPK
        SURROGATE  = 3,   // surrogateModel.hh
        COBRAPY    = 4,   // flux balance analysis via COBRApy
        GLPK_KNS   = 5,   // GLPK      and defineKinetics.hh together
        SRG_KNS    = 6,   // surrogate and defineKinetics.hh together
        CPY_KNS    = 7,   // COBRApy   and defineKinetics.hh together
        SYMBOLIC   = 8,   // an algebraic rate law read from a .sym file
        GRAPHNET   = 9,   // a message-passing network read from a .gnn file
        SYM_KNS    = 10,  // symbolic  and defineKinetics.hh together
        GNN_KNS    = 11   // graphnet  and defineKinetics.hh together
    };

    inline bool usesKinetics  (plint r) { return r==KINETICS || r==GLPK_KNS || r==SRG_KNS
                                              || r==CPY_KNS  || r==SYM_KNS  || r==GNN_KNS; }
    inline bool usesGlpk      (plint r) { return r==GLPK      || r==GLPK_KNS; }
    inline bool usesCobrapy   (plint r) { return r==COBRAPY   || r==CPY_KNS; }
    inline bool usesSurrogate (plint r) { return r==SURROGATE || r==SRG_KNS; }
    inline bool usesSymbolic  (plint r) { return r==SYMBOLIC  || r==SYM_KNS; }
    inline bool usesGraphnet  (plint r) { return r==GRAPHNET  || r==GNN_KNS; }
    inline bool usesFBA       (plint r) { return usesGlpk(r)  || usesCobrapy(r); }
    /* "learned" = fitted offline and read from a file at start-up, as opposed to solved in line.
     * The surrogate is deliberately NOT in here: it predates these two and is threaded through the
     * metabolic config, whereas these two carry their own files and their own runtime registry. */
    inline bool usesLearned   (plint r) { return usesSymbolic(r) || usesGraphnet(r); }
    inline bool usesMetabolic (plint r) { return usesFBA(r)   || usesSurrogate(r); }

    inline std::string name (plint r) {
        switch (r) {
            case NONE:      return "none";
            case KINETICS:  return "kinetics";
            case GLPK:      return "glpk";
            case SURROGATE: return "surrogate";
            case COBRAPY:   return "cobrapy";
            case GLPK_KNS:  return "glpk_and_kinetics";
            case SRG_KNS:   return "surrogate_and_kinetics";
            case CPY_KNS:   return "cobrapy_and_kinetics";
            case SYMBOLIC:  return "symbolic";
            case GRAPHNET:  return "graphnet";
            case SYM_KNS:   return "symbolic_and_kinetics";
            case GNN_KNS:   return "graphnet_and_kinetics";
            default:        return "UNKNOWN";
        }
    }

    /* Parse the <reaction_type> string.  Returns -1 if it is not recognised, so
     * the caller can print the offending text and stop.  Accepts the spellings
     * CompLaB 2D accepted, so 2D input files port unchanged. */
    inline plint parse (std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        /* [FIX-3D] "1" was missing, so an existing CompLB3D input file written as
         * <reaction_type>1</reaction_type> -- which the comment above promises
         * still works -- was rejected outright.  All eight numeric forms are
         * accepted now, matching the enum. */
        if (s=="none"      || s=="no"  || s=="0")                          return NONE;
        if (s=="kinetics"  || s=="kns" || s=="1")                          return KINETICS;
        if (s=="2") return GLPK;
        if (s=="3") return SURROGATE;
        if (s=="4") return COBRAPY;
        if (s=="5") return GLPK_KNS;
        if (s=="6") return SRG_KNS;
        if (s=="7") return CPY_KNS;
        if (s=="8") return SYMBOLIC;
        if (s=="9") return GRAPHNET;
        if (s=="10") return SYM_KNS;
        if (s=="11") return GNN_KNS;
        if (s=="glpk"      || s=="fba" || s=="fba_glpk")                   return GLPK;
        if (s=="surrogate" || s=="srg" || s=="ann")                        return SURROGATE;
        if (s=="cobrapy"   || s=="cpy" || s=="fba_cobrapy")                return COBRAPY;
        if (s=="symbolic"  || s=="sym" || s=="sym_law" || s=="rate_law")   return SYMBOLIC;
        if (s=="graphnet"  || s=="gnn" || s=="graph_network"
                           || s=="graph")                                  return GRAPHNET;
        if (s=="glpk_and_kinetics"      || s=="glpk+kinetics")             return GLPK_KNS;
        if (s=="surrogate_and_kinetics" || s=="surrogate+kinetics")        return SRG_KNS;
        if (s=="cobrapy_and_kinetics"   || s=="cobrapy+kinetics")          return CPY_KNS;
        if (s=="symbolic_and_kinetics"  || s=="symbolic+kinetics")         return SYM_KNS;
        if (s=="graphnet_and_kinetics"  || s=="graphnet+kinetics")         return GNN_KNS;
        return -1;
    }
}


/* ============================================================================
 *  MetabolicTotals
 *
 *  Turns the equilibrium tableau into a lookup that answers two questions:
 *
 *    "how much acetate is there in this voxel ALTOGETHER?"          total()
 *    "I need to remove X mol/L of acetate -- from which species,
 *     and how much from each?"                                      drawDown()
 *
 *  Built once at start-up from the same stoichiometry matrix the equilibrium
 *  solver uses, so it stays correct if the tableau is reordered.  That is the
 *  deliberate improvement over defineKinetics.hh, which hard-codes C[56],
 *  C[81]..C[86] and friends and silently corrupts itself if anyone edits the
 *  substrate order in CompLaB.xml.
 *
 *  When equilibrium is off, or the basis is "free", `active` stays false and
 *  every call degenerates to the plain free concentration.
 * ============================================================================
 */
struct MetabolicTotals {
    bool active;

    /* For each substrate index s that is an equilibrium COMPONENT, contrib[s]
     * lists every (species index, stoichiometric coefficient) pair that carries
     * that component -- including the component's own identity row.
     * For every other substrate, contrib[s] is empty. */
    std::vector< std::vector< std::pair<plint,T> > > contrib;

    MetabolicTotals() : active(false) {}

    /* eq_stoich[species][component], eq_names = component names,
     * subs_names = substrate names in CompLaB.xml order.
     * Returns the number of components successfully mapped onto a substrate. */
    plint build (const std::vector< std::vector<T> > &eq_stoich,
                 const std::vector<std::string>      &eq_names,
                 const std::vector<std::string>      &subs_names,
                 plint num_of_substrates)
    {
        active = false;
        contrib.assign(num_of_substrates, std::vector< std::pair<plint,T> >());
        if (eq_stoich.empty() || eq_names.empty()) return 0;

        const plint nSpecies = (plint) eq_stoich.size();
        plint mapped = 0;

        for (plint j = 0; j < (plint) eq_names.size(); ++j) {
            /* Which substrate carries this component?  Match by name. */
            plint sIdx = -1;
            for (plint s = 0; s < num_of_substrates && s < (plint) subs_names.size(); ++s) {
                if (subs_names[s] == eq_names[j]) { sIdx = s; break; }
            }
            if (sIdx < 0) continue;               // component has no transported substrate

            for (plint i = 0; i < nSpecies && i < num_of_substrates; ++i) {
                if (j >= (plint) eq_stoich[i].size()) continue;
                const T coeff = eq_stoich[i][j];
                if (coeff > (T) 1e-12 || coeff < (T) -1e-12) {
                    contrib[sIdx].push_back(std::make_pair(i, coeff));
                }
            }
            if (!contrib[sIdx].empty()) ++mapped;
        }

        active = (mapped > 0);
        return mapped;
    }

    /* Raw SIGNED component total: sum over the component itself and every
     * complex that carries it.  May legitimately be negative -- a proton
     * condition (TOTH) is negative in any alkaline water, because the OH-,
     * CO3-- and hydroxo-complex rows carry negative coefficients.
     *
     * [FIX-3D] The first version returned the FREE concentration whenever the
     * signed total came out <= 0, and drawDown then used that fallback as a
     * divisor.  For a real tableau with negative rows that silently inverted the
     * shares: asking to remove 1 nmol/L of the H component removed 10 umol/L of
     * ammonia, moving the element balance 10,000x the wrong way. */
    inline T rawTotal (const std::vector<T> &conc, plint s) const {
        T tot = T();
        for (size_t k = 0; k < contrib[s].size(); ++k) {
            const plint i = contrib[s][k].first;
            if (i < (plint) conc.size()) tot += contrib[s][k].second * conc[i];
        }
        return tot;
    }

    /* Sum of |coeff * C| over the carriers -- the scale against which a near
     * cancellation in rawTotal should be judged. */
    inline T absScale (const std::vector<T> &conc, plint s) const {
        T sc = T();
        for (size_t k = 0; k < contrib[s].size(); ++k) {
            const plint i = contrib[s][k].first;
            if (i < (plint) conc.size()) sc += std::fabs(contrib[s][k].second * conc[i]);
        }
        return sc;
    }

    /* True when substrate s has a usable component total: the table is active,
     * s is in range, it has carriers, and the signed total is not a near-exact
     * cancellation of much larger terms.
     *
     * [FIX-3D] The old guard was the absolute `tot > 1e-30`.  The amplification
     * factor is absScale/|tot|, which blows up long before 1e-30: a total of
     * 1e-9 built from two 1e-3 terms of opposite sign turned a 1e-10 removal
     * request into a 1e-4 change in two other substrates.  The test is relative
     * now. */
    inline bool usable (const std::vector<T> &conc, plint s) const {
        if (!active || s < 0 || s >= (plint) contrib.size()) return false;
        if (contrib[s].empty()) return false;
        const T tot = rawTotal(conc, s);
        return std::fabs(tot) > (T) 1e-6 * absScale(conc, s) && std::fabs(tot) > (T) 1e-30;
    }

    /* Total concentration of substrate s.  Falls back to the free value when the
     * tableau cannot give a trustworthy answer.
     *
     * [FIX-3D] The bounds test used to detect an out-of-range s and then index
     * conc[s] with it anyway -- a guard that read as protection and was not one.
     * Range is checked before any subscript now. */
    inline T total (const std::vector<T> &conc, plint s) const {
        if (s < 0 || s >= (plint) conc.size()) return T();
        if (!usable(conc, s)) return conc[s];
        return std::fabs(rawTotal(conc, s));
    }

    /* Remove `amount` (>0 consumed, <0 produced) of the element held by
     * substrate s, spread across its carriers in proportion to how much of the
     * element each one holds.  Accumulates into dC.
     *
     * `committed` carries the running dC of THIS voxel so the per-species
     * positivity clamp can see draws made for other components.
     *
     * [FIX-3D] Without that, two component draws that were each individually
     * legal could between them take a shared complex negative: CaAc+ carries
     * both acetate and calcium, and each draw was clamped only against its own
     * component total.  Now every species write is clamped against what is
     * actually left of that species. */
    inline void drawDown (const std::vector<T> &conc, plint s, T amount,
                          std::vector<T> &dC) const
    {
        if (s < 0 || s >= (plint) dC.size()) return;
        if (!usable(conc, s)) {
            /* [FIX-3D] The fallback branch used to subtract unconditionally, with no
             * positivity clamp, while the branch below had one -- so whether a
             * substrate could go negative depended on whether the tableau happened
             * to give it a usable total. */
            T take = amount;
            const T remaining = conc[s] + dC[s];
            if (take > remaining) take = (remaining > T()) ? remaining : T();
            dC[s] -= take;
            return;
        }

        const T tot = rawTotal(conc, s);

        for (size_t k = 0; k < contrib[s].size(); ++k) {
            const plint i     = contrib[s][k].first;
            const T     coeff = contrib[s][k].second;
            if (i >= (plint) conc.size() || i >= (plint) dC.size()) continue;
            if (!(conc[i] > T())) continue;

            /* share of the element held by species i; removing `amount` of the
             * element means removing amount*share/coeff of the species itself.
             * The coeff cancels exactly, so a stoichiometry of 2 correctly
             * removes half as many molecules per unit of element. */
            const T share = (coeff * conc[i]) / tot;
            T       take  = amount * share / coeff;

            /* per-species positivity clamp, aware of draws already committed */
            const T remaining = conc[i] + dC[i];
            if (take > remaining) take = (remaining > T()) ? remaining : T();

            dC[i] -= take;
        }
    }
};


/* ============================================================================
 *  MetabolicConfig
 *
 *  Everything the FBA / surrogate layer needs, in one object.
 *
 *  This is deliberately NOT bolted onto initialize_complab(), which already
 *  takes more than sixty arguments by reference and is the least maintainable
 *  thing in the code base.  initialize_metabolic() below re-opens CompLaB.xml
 *  and reads only the metabolic settings, so the existing parser is untouched
 *  and this whole feature can be lifted out again by deleting three lines from
 *  complab.cpp.
 * ============================================================================
 */
struct MetabolicConfig {

    /* ---- run switches, from <simulation_mode> ---------------------------- */
    bool enable_fba_glpk;
    bool enable_fba_cobrapy;
    bool enable_surrogate;
    /* The two learned paths.  Their files are read by integ::readConfig() and their networks live
     * in complab_sym::runtime() / complab_gnn::runtime(), so all this struct carries is whether
     * the switch is on and how many organisms asked for it -- which is what the dispatch in
     * complab.cpp gates on, exactly as it does for the three above. */
    bool enable_symbolic;
    bool enable_graphnet;
    /* The thermodynamic gate is not a rate path -- no microbe has reaction_type thermo -- so this
     * switch has no count beside it.  It is read here only so the summary below can say whether
     * the run is thermodynamically controlled; the file itself is loaded by integ::prepareThermo,
     * which runs later, once the substrate and microbe names are known. */
    bool enable_thermo;

    /* ---- derived counts -------------------------------------------------- */
    plint glpk_count;      // microbes whose reaction_type uses GLPK
    plint cpy_count;       // ... COBRApy
    plint srg_count;       // ... the surrogate
    plint sym_count;       // ... a .sym rate law
    plint gnn_count;       // ... a .gnn network
    plint mm_count;        // glpk_count + cpy_count  (microbes with a metabolic model)

    /* ---- per-microbe, indexed by GLOBAL microbe id ------------------------ */
    std::vector<std::string>            model_filename;     // "" when not an FBA microbe
    std::vector< std::vector<plint> >   subsLoc;            // [microbe][substrate] -> exchange column, -1 if unused
    /* [NEW] <exchange_reaction_names>: the same mapping written as reaction names rather than
     * column numbers.  Empty when the user gave indices instead.  Resolved into subsLoc by
     * load_metabolic_models3D(), because it needs the model's reaction list to do it -- which is
     * the whole point: a name survives a change of model revision, a column number does not. */
    std::vector< std::vector<std::string> > exchange_names;
    /* [FIX-3D] The Michaelis-Menten Vmax now comes from <maximum_uptake_flux>,
     * which is what the documentation, the 2-D code and every comment always
     * said it did.  The first version of this port used <substrate_lower_bounds>
     * as the Vmax, so a user who set <maximum_uptake_flux> saw it silently
     * ignored, and a user who set only <substrate_lower_bounds> got its
     * magnitude used as an uptake rate.  The two are different quantities:
     *     vmax       positive, mmol/gDW/h, the enzyme's maximum rate
     *     maxUptake  negative, mmol/gDW/h, a hard floor the LP may not pass
     * A vmax entry of 0 means "no enzyme limit" and selects the supply-limited
     * fallback, exactly as documented. */
    std::vector< std::vector<T> >       vmax;               // [microbe][substrate] Monod Vmax (>=0), 0 = unset
    std::vector< std::vector<T> >       maxUptake;          // [microbe][substrate] lower bound floor  (<=0)
    std::vector< std::vector<T> >       maxRelease;         // [microbe][substrate] upper bound        (>=0)
    std::vector<plint>                  objDir;             // -1 maximize, +1 minimize
    /* [v1.3] "the user wrote <objective_direction>" has to be distinguishable from "nobody
     * said", because both used to be stored as objDir = -1 and the SBML override below then
     * could not tell them apart: an explicit `maximize` was silently flipped to minimize by an
     * fbc:type="minimize" in the model file, which is the opposite of what its comment says. */
    std::vector<char>                   objDirExplicit;     // 1 where CompLaB.xml set it
    std::vector< std::vector<int> >     constraint_loc;     // extra bound overrides
    std::vector< std::vector<T> >       constraint_lb;
    std::vector< std::vector<T> >       constraint_ub;
    std::vector<bool>                   equate_bounds;      // force ub = lb + eps
    std::vector<T>                      biomass_molar_mass; // gDW per mol cells

    /* ---- per-substrate --------------------------------------------------- */
    std::vector<bool> fix_concentration;   // hold this substrate fixed (do not apply dC)
    std::vector<bool> fix_lower_bounds;    // clamp the uptake bound to -C exactly

    /* ---- concentration basis --------------------------------------------- */
    bool            useTotals;             // <fba_concentration_basis> total
    MetabolicTotals totals;

    /* ---- loaded metabolic models (FBA microbes only, indexed by GLOBAL id)  */
    std::vector< std::vector< std::vector<T> > > S3;     // dense stoichiometry
    std::vector< std::vector<T> >                S1;     // flattened, for COBRApy
    std::vector< std::vector<T> >                vec_b, vec_c, vec_lb, vec_ub;
    std::vector<plint>                           vec_objLoc;
    std::vector<int>                             vec_nmets, vec_nrxns;

    /* ---- refinements of the program, per microbe, both off by default ------ */
    /* <multi_step>: the lexicographic chain.  Empty or one stage == plain FBA. */
    std::vector<complab_lex::LexPlan>       lexPlan;
    /* <cybernetic>: alternative carbon sources and the weights between them.
     * Fewer than two sources == one growth option == the ordinary path. */
    std::vector<complab_cyb::CyberneticPlan> cybPlan;

    /* ---- solver handles --------------------------------------------------- */
#ifdef COMPLAB_ENABLE_GLPK
    std::vector<glp_prob*> vec_lp;
    std::vector<glp_smcp>  sParam;
    std::vector<glp_iocp>  iParam;
    std::vector<int>       method;     // 'S','I','T','E' as int
    std::vector<int>       isMIP;
    std::vector<char*>     ctype;
    std::vector<char*>     vtype;
#endif
#ifdef COMPLAB_ENABLE_COBRAPY
    std::vector<PyObject*> vec_model;
#endif
    int lpsolver;      // 1 simplex, 2 interior, 3 exact.  1 unless you know why.
    int save_pb;       // 0 off; 1 .lp, 2/3 .mps  -- debugging only
    bool released;     // set by release_metabolic_solvers so it is safe to call twice

    /* [FIX-3D] This struct owns raw glp_prob*, char* and PyObject* handles that
     * release_metabolic_solvers() frees.  With the implicit copy constructor a
     * single accidental `MetabolicConfig c2 = mmcfg;` would double-free them all
     * at shutdown.  The processors only ever hold a pointer, so nothing needs to
     * copy it. */
    MetabolicConfig(const MetabolicConfig&);            // not defined: copying is forbidden
    MetabolicConfig& operator=(const MetabolicConfig&); // not defined

    MetabolicConfig()
      : enable_fba_glpk(false), enable_fba_cobrapy(false), enable_surrogate(false),
        enable_symbolic(false), enable_graphnet(false), enable_thermo(false),
        glpk_count(0), cpy_count(0), srg_count(0), sym_count(0), gnn_count(0), mm_count(0),
        useTotals(false), lpsolver(1), save_pb(0), released(false)
    {}

    bool anyEnabled() const { return enable_fba_glpk || enable_fba_cobrapy || enable_surrogate
                                  || enable_symbolic || enable_graphnet; }
    /* The learned paths need vec_Kc / vec_mu sized and the totals tableau built, but they do NOT
     * need a metabolic model loaded, so they are separated out where that distinction matters. */
    bool anyLearned() const { return enable_symbolic || enable_graphnet; }

    /* ---- the two unit conversions, in one place only --------------------- */

    /* v [mmol/gDW/h] , B [mol/L] , dt [s]  ->  dC [mol/L] */
    inline T fluxToDeltaC (T v, T B_molL, plint microbeId, T dt) const {
        const T MW = biomass_molar_mass[microbeId];      // gDW per mol cells
        return v * (B_molL * MW) * dt / (T) 3.6e6;       // /1000 mmol->mol, /3600 h->s
    }
    /* mu [1/h] , B [mol/L] , dt [s]  ->  dB [mol/L]
     * For growth rates coming OUT of a metabolic model, which are per hour. */
    inline T growthToDeltaB (T mu, T B_molL, T dt) const {
        return mu * B_molL * dt / (T) 3600.0;
    }
    /* [FIX-3D] <decay_coefficient> is per SECOND everywhere in CompLB3D -- both
     * shipped XML files use 6.94e-8, which is 6e-3 per day expressed per second,
     * and defineKinetics.hh works in per-second units throughout.  The decay
     * fallback used to be pushed through growthToDeltaB, which divides by 3600,
     * so decay ran 3600x too slow: biomass that should have died away persisted
     * indefinitely, and the fallback fires at every voxel where the LP predicts
     * no growth -- i.e. constantly during the early transient.
     * kd [1/s] , B [mol/L] , dt [s]  ->  dB [mol/L] */
    inline T decayToDeltaB (T kd_per_s, T B_molL, T dt) const {
        return -kd_per_s * B_molL * dt;
    }
    /* The "cell-specific" fallback used when no <maximum_uptake_flux> is given:
     * the uptake flux that would exactly exhaust C in one time step.
     * Returns a NEGATIVE number (uptake), in mmol/gDW/h. */
    inline T exhaustionFlux (T C_molL, T B_molL, plint microbeId, T dt) const {
        if (microbeId < 0 || microbeId >= (plint) biomass_molar_mass.size()) return T();
        const T MW = biomass_molar_mass[microbeId];
        const T denom = B_molL * MW * dt;
        if (!(denom > (T) 1e-30)) return T();
        T v = -C_molL * (T) 3.6e6 / denom;
        /* [FIX-3D] A voxel holding barely more biomass than the participation
         * threshold (1e-12 mol/L) produced bounds around -1e17 mmol/gDW/h.  Handing
         * GLPK a column spanning twenty-odd orders of magnitude is a reliable way
         * to get GLP_ECOND or a meaningless optimum.  Cap at a flux no organism
         * approaches; 1e4 mmol/gDW/h is already three orders above anything
         * physiological. */
        const T CAP = (T) 1e4;
        if (v < -CAP) v = -CAP;
        return v;
    }
};


/* ============================================================================
 *  initialize_metabolic
 *
 *  Reads the metabolic settings out of CompLaB.xml.  Call it AFTER
 *  initialize_complab(), passing the reaction_type vector it produced.
 *
 *  Returns 0 on success, -1 on a fatal configuration error (message already
 *  printed).  Returns 0 and leaves everything switched off if the user asked
 *  for nothing -- that is not an error, it is the default CompLB3D.
 * ============================================================================
 */
inline int initialize_metabolic (MetabolicConfig &cfg,
                                 const std::vector<plint> &reaction_type,
                                 const std::vector<std::string> &subs_names,
                                 plint num_of_microbes, plint num_of_substrates,
                                 bool useEquilibrium,
                                 const std::vector< std::vector<T> > &eq_stoich,
                                 const std::vector<std::string>      &eq_component_names,
                                 const std::vector< std::vector<T> > &vec_Kc,
                                 const std::vector<T>                &vec_mu)
{
    std::string fin(complab_input::configPath());

    /* ---------------------------------------------------------------- switches */
    {
        XMLreader doc(fin);
        const char *keys[3] = { "enable_fba_glpk", "enable_fba_cobrapy", "enable_surrogate" };
        bool *dest[3]       = { &cfg.enable_fba_glpk, &cfg.enable_fba_cobrapy, &cfg.enable_surrogate };
        for (int k = 0; k < 3; ++k) {
            try {
                std::string tmp;
                doc["parameters"]["simulation_mode"][keys[k]].read(tmp);
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
                /* [FIX-3D] Anything unrecognised used to mean "off", so a typo such
                 * as <enable_fba_glpk>ture</enable_fba_glpk> silently disabled the
                 * feature the user had just asked for. */
                if      (tmp=="yes" || tmp=="true"  || tmp=="1" || tmp=="on")  *dest[k] = true;
                else if (tmp=="no"  || tmp=="false" || tmp=="0" || tmp=="off") *dest[k] = false;
                else {
                    pcout << "<" << keys[k] << "> is \"" << tmp
                          << "\", which is not a yes/no value. Use true or false. Terminating.\n";
                    return -1;
                }
            }
            catch (PlbIOException& exception) { *dest[k] = false; }
        }

        try {
            std::string tmp;
            doc["parameters"]["simulation_mode"]["fba_concentration_basis"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if      (tmp=="total") cfg.useTotals = true;
            else if (tmp=="free")  cfg.useTotals = false;
            else {
                pcout << "fba_concentration_basis \"" << tmp << "\" is not understood. Use free or total. "
                      << "Terminating the simulation.\n";
                return -1;
            }
        }
        catch (PlbIOException& exception) { cfg.useTotals = false; }

        try { doc["parameters"]["simulation_mode"]["glpk_lpsolver"].read(cfg.lpsolver); }
        catch (PlbIOException& exception) { cfg.lpsolver = 1; }
        try { doc["parameters"]["simulation_mode"]["glpk_save_problem"].read(cfg.save_pb); }
        catch (PlbIOException& exception) { cfg.save_pb = 0; }

        /* The two learned paths keep their switch inside their own block, beside the file they
         * read, because that is where a user looks for it.  Same yes/no discipline: a typo stops
         * the run rather than quietly turning the feature off. */
        {
            const char *blk[3]  = { "symbolic", "graphnet", "thermodynamics" };
            bool *ldest[3]      = { &cfg.enable_symbolic, &cfg.enable_graphnet,
                                    &cfg.enable_thermo };
            for (int k = 0; k < 3; ++k) {
                try {
                    std::string tmp;
                    doc["parameters"][blk[k]]["enabled"].read(tmp);
                    std::transform(tmp.begin(), tmp.end(), tmp.begin(),
                                   [](unsigned char c){ return std::tolower(c); });
                    if      (tmp=="yes" || tmp=="true"  || tmp=="1" || tmp=="on")  *ldest[k] = true;
                    else if (tmp=="no"  || tmp=="false" || tmp=="0" || tmp=="off") *ldest[k] = false;
                    else {
                        pcout << "<" << blk[k] << "><enabled> is \"" << tmp
                              << "\", which is not a yes/no value. Use true or false. Terminating.\n";
                        return -1;
                    }
                }
                catch (PlbIOException& exception) { *ldest[k] = false; }
            }
        }
    }

    /* ------------------------------------------------- count what is requested */
    cfg.glpk_count = cfg.cpy_count = cfg.srg_count = cfg.sym_count = cfg.gnn_count = 0;
    for (plint iM = 0; iM < num_of_microbes && iM < (plint) reaction_type.size(); ++iM) {
        if (rxntype::usesGlpk     (reaction_type[iM])) ++cfg.glpk_count;
        if (rxntype::usesCobrapy  (reaction_type[iM])) ++cfg.cpy_count;
        if (rxntype::usesSurrogate(reaction_type[iM])) ++cfg.srg_count;
        if (rxntype::usesSymbolic (reaction_type[iM])) ++cfg.sym_count;
        if (rxntype::usesGraphnet (reaction_type[iM])) ++cfg.gnn_count;
    }
    cfg.mm_count = cfg.glpk_count + cfg.cpy_count;

    /* ------------------------------------------------------ consistency checks */
    if (cfg.glpk_count > 0 && !cfg.enable_fba_glpk) {
        pcout << "A microbe asks for reaction_type glpk but <enable_fba_glpk> is false. "
              << "Set it to true, or change the microbe's reaction_type. Terminating.\n";
        return -1;
    }
    if (cfg.cpy_count > 0 && !cfg.enable_fba_cobrapy) {
        pcout << "A microbe asks for reaction_type cobrapy but <enable_fba_cobrapy> is false. "
              << "Set it to true, or change the microbe's reaction_type. Terminating.\n";
        return -1;
    }
    if (cfg.srg_count > 0 && !cfg.enable_surrogate) {
        pcout << "A microbe asks for reaction_type surrogate but <enable_surrogate> is false. "
              << "Set it to true, or change the microbe's reaction_type. Terminating.\n";
        return -1;
    }
    if (cfg.sym_count > 0 && !cfg.enable_symbolic) {
        pcout << "A microbe asks for reaction_type symbolic but <symbolic><enabled> is false. "
              << "Set it to true, or change the microbe's reaction_type. Terminating.\n";
        return -1;
    }
    if (cfg.gnn_count > 0 && !cfg.enable_graphnet) {
        pcout << "A microbe asks for reaction_type graphnet but <graphnet><enabled> is false. "
              << "Set it to true, or change the microbe's reaction_type. Terminating.\n";
        return -1;
    }
#ifndef COMPLAB_ENABLE_GLPK
    if (cfg.enable_fba_glpk) {
        pcout << "<enable_fba_glpk> is true but this executable was built without GLPK.\n"
              << "  Rebuild with:  cmake -DENABLE_GLPK=ON ..\n"
              << "  (and make sure libglpk and glpk.h are installed). Terminating.\n";
        return -1;
    }
#endif
#ifndef COMPLAB_ENABLE_COBRAPY
    if (cfg.enable_fba_cobrapy) {
        pcout << "<enable_fba_cobrapy> is true but this executable was built without Python.\n"
              << "  Rebuild with:  cmake -DENABLE_COBRAPY=ON ..\n"
              << "  (and make sure Python.h and cobrapy are installed). Terminating.\n";
        return -1;
    }
#endif
    /* CompLaB 2D refused to mix the two FBA back ends in one run.  The reason given
     * was that both write into the same increment lattices with different flux index
     * conventions.  That reason no longer holds: both back ends address substrates
     * through cfg.subsLoc[globalMicrobe][substrate] and convert through the one
     * shared fluxToDeltaC(), so a flux means the same thing on either path.
     *
     * What actually broke was narrower.  Both processors were handed ONE list of FBA
     * microbes, built over usesFBA(), so each tried to solve for the other back end's
     * organisms -- and cfg.vec_lp[gM] is null for a COBRApy microbe.  complab.cpp now
     * builds glpk_globalId and cpy_globalId separately and gives each processor its
     * own lattice vector, so neither can see the other's organisms.  The two then add
     * into the shared dC/dB increments in the same way the surrogate, symbolic and
     * graphnet paths already coexist, and update_rxnLattices applies the sum once.
     *
     * One honest caveat, the same one that applies to any two paths running side by
     * side: if a GLPK organism and a COBRApy organism draw on the SAME substrate in
     * the SAME voxel, each sizes its uptake without knowing about the other, and only
     * the positivity clamp in update_rxnLattices stops the pair from over-drawing.
     * The within-back-end exhaustion repair does not cross back ends.  Give competing
     * organisms the same back end if that matters for your problem. */
    if (cfg.glpk_count > 0 && cfg.cpy_count > 0) {
        pcout << "  note: this run mixes GLPK and COBRApy organisms ("
              << cfg.glpk_count << " GLPK, " << cfg.cpy_count << " COBRApy).\n"
              << "        Each back end solves only its own organisms and both add into the same\n"
              << "        increment lattices. Shared-substrate exhaustion is repaired within a back\n"
              << "        end, not across the two.\n";
    }

    if (!cfg.anyEnabled()) return 0;   // nothing asked for: leave everything alone

    /* [FIX-3D] vec_Kc and vec_mu are indexed [microbe][substrate] by all three
     * metabolic processors, and this layer is their FIRST consumer -- run_kinetics
     * stores them and never reads them.  The existing parser is lax: a wrong-length
     * <half_saturation_constants> prints a message but does NOT stop the run and
     * does not push a row, and a missing tag pushes a single-element {-99} row.
     * Either leaves the vectors ragged, and the processors would read past the end
     * of a one-element row for every substrate above the first.  Check here, before
     * anything can index them. */
    if ((plint) vec_Kc.size() != num_of_microbes) {
        pcout << "  <half_saturation_constants> is missing or the wrong length for at least one microbe: "
              << "found " << (plint) vec_Kc.size() << " rows for " << num_of_microbes << " microbes.\n"
              << "  Every microbe needs exactly " << num_of_substrates
              << " values (use 0 for substrates it does not use). Terminating.\n";
        return -1;
    }
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (!rxntype::usesMetabolic(reaction_type[iM])) continue;
        if ((plint) vec_Kc[iM].size() != num_of_substrates) {
            pcout << "  microbe" << iM << ": <half_saturation_constants> has "
                  << (plint) vec_Kc[iM].size() << " values but there are " << num_of_substrates
                  << " substrates. Terminating.\n";
            return -1;
        }
    }
    if ((plint) vec_mu.size() != num_of_microbes) {
        pcout << "  <decay_coefficient> is missing for at least one microbe ("
              << (plint) vec_mu.size() << " values for " << num_of_microbes << " microbes). Terminating.\n";
        return -1;
    }

    /* -------------------------------------------------- allocate per-microbe -- */
    cfg.model_filename.assign(num_of_microbes, std::string());
    cfg.subsLoc      .assign(num_of_microbes, std::vector<plint>(num_of_substrates, -1));
    cfg.exchange_names.assign(num_of_microbes, std::vector<std::string>());
    cfg.vmax         .assign(num_of_microbes, std::vector<T>(num_of_substrates, T()));
    cfg.maxUptake    .assign(num_of_microbes, std::vector<T>(num_of_substrates, (T) -1e30));
    cfg.maxRelease   .assign(num_of_microbes, std::vector<T>(num_of_substrates, (T)  1e30));
    cfg.objDir       .assign(num_of_microbes, -1);
    cfg.objDirExplicit.assign(num_of_microbes, (char) 0);
    cfg.constraint_loc.assign(num_of_microbes, std::vector<int>());
    cfg.constraint_lb .assign(num_of_microbes, std::vector<T>());
    cfg.constraint_ub .assign(num_of_microbes, std::vector<T>());
    cfg.equate_bounds .assign(num_of_microbes, false);
    cfg.biomass_molar_mass.assign(num_of_microbes, (T) 24.6);

    cfg.S3.assign(num_of_microbes, std::vector< std::vector<T> >());
    cfg.S1.assign(num_of_microbes, std::vector<T>());
    cfg.vec_b .assign(num_of_microbes, std::vector<T>());
    cfg.vec_c .assign(num_of_microbes, std::vector<T>());
    cfg.vec_lb.assign(num_of_microbes, std::vector<T>());
    cfg.vec_ub.assign(num_of_microbes, std::vector<T>());
    cfg.vec_objLoc.assign(num_of_microbes, 0);
    cfg.vec_nmets .assign(num_of_microbes, 0);
    cfg.vec_nrxns .assign(num_of_microbes, 0);
    cfg.lexPlan.assign(num_of_microbes, complab_lex::LexPlan());
    cfg.cybPlan.assign(num_of_microbes, complab_cyb::CyberneticPlan());

    /* ------------------------------------------------------- per-substrate --- */
    cfg.fix_concentration.assign(num_of_substrates, false);
    cfg.fix_lower_bounds .assign(num_of_substrates, false);
    {
        XMLreader doc(fin);
        /* [FIX-3D] CompLaB 2D writes these as words -- "no yes no yes" -- and the
         * first version of this parser read them as ints, so every 2D input file
         * threw and silently fell back to all-false.  Accept both spellings. */
        const char *flagKeys[2] = { "fix_concentration", "fix_lower_bounds" };
        std::vector<bool> *flagDest[2] = { &cfg.fix_concentration, &cfg.fix_lower_bounds };
        for (int fk = 0; fk < 2; ++fk) {
            try {
                std::vector<std::string> tmp;
                doc["parameters"]["chemistry"][flagKeys[fk]].read(tmp);
                if ((plint) tmp.size() != num_of_substrates) {
                    pcout << "<" << flagKeys[fk] << "> has " << (plint) tmp.size()
                          << " entries but there are " << num_of_substrates << " substrates. Terminating.\n";
                    return -1;
                }
                for (plint s = 0; s < num_of_substrates; ++s) {
                    std::string v = tmp[s];
                    std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return std::tolower(c); });
                    if      (v=="yes" || v=="true"  || v=="1") (*flagDest[fk])[s] = true;
                    else if (v=="no"  || v=="false" || v=="0") (*flagDest[fk])[s] = false;
                    else {
                        pcout << "Element " << s << " of <" << flagKeys[fk] << "> is \"" << tmp[s]
                              << "\", which is not a yes/no value. Terminating.\n";
                        return -1;
                    }
                }
            } catch (PlbIOException& exception) {}
        }
    }

    /* ------------------------------------------------- per-microbe settings --- */
    {
        XMLreader doc(fin);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            const plint rt = (iM < (plint) reaction_type.size()) ? reaction_type[iM] : (plint) rxntype::NONE;
            if (!rxntype::usesMetabolic(rt)) continue;

            const std::string bioname = "microbe" + std::to_string(iM);

            /* biomass molar mass -- needed by every metabolic solver */
            try { doc["parameters"]["microbiology"][bioname]["biomass_molar_mass"].read(cfg.biomass_molar_mass[iM]); }
            catch (PlbIOException& exception) {
                pcout << "  note: " << bioname << " has no <biomass_molar_mass>; using the default 24.6 g/mol "
                      << "(CH1.8O0.5N0.2).\n";
            }
            if (!(cfg.biomass_molar_mass[iM] > (T) 0)) {
                pcout << bioname << ": <biomass_molar_mass> must be positive. Terminating.\n";
                return -1;
            }

            /* which substrates this microbe exchanges, and where in its model.
             *
             * [NEW] Two spellings are accepted.  <exchange_reaction_indices> is the original:
             * column numbers, which are fast and exact but silently mean something different the
             * moment the model file changes revision.  <exchange_reaction_names> is the safer one:
             * reaction names, resolved against the model itself in load_metabolic_models3D().  If
             * both are given the names win, because they are the ones that can be checked. */
            bool haveNames = false;
            try {
                std::vector<std::string> nm;
                doc["parameters"]["microbiology"][bioname]["exchange_reaction_names"].read(nm);
                if ((plint) nm.size() != num_of_substrates) {
                    pcout << bioname << ": <exchange_reaction_names> has " << (plint) nm.size()
                          << " entries but there are " << num_of_substrates << " substrates. "
                          << "Use 'none' for a substrate this microbe does not exchange. Terminating.\n";
                    return -1;
                }
                cfg.exchange_names[iM] = nm;
                haveNames = true;
            }
            catch (PlbIOException& exception) {}

            try {
                std::vector<plint> loc;
                doc["parameters"]["microbiology"][bioname]["exchange_reaction_indices"].read(loc);
                if ((plint) loc.size() != num_of_substrates) {
                    pcout << bioname << ": <exchange_reaction_indices> has " << (plint) loc.size()
                          << " entries but there are " << num_of_substrates << " substrates. "
                          << "Use -1 for a substrate this microbe does not exchange. Terminating.\n";
                    return -1;
                }
                /* [FIX-3D] CompLaB 2D marks an unused substrate with -99; this code
                 * uses -1.  Accept both, so 2D input files port unchanged. */
                for (size_t q = 0; q < loc.size(); ++q) if (loc[q] < 0) loc[q] = -1;
                if (haveNames)
                    pcout << "  note: " << bioname << " gives both <exchange_reaction_names> and "
                          << "<exchange_reaction_indices>. The names are used; the indices are ignored.\n";
                else
                    cfg.subsLoc[iM] = loc;
            }
            catch (PlbIOException& exception) {
                if (rxntype::usesFBA(rt) && !haveNames) {
                    pcout << bioname << ": reaction_type " << rxntype::name(rt)
                          << " requires <exchange_reaction_indices> or <exchange_reaction_names>. "
                          << "Terminating.\n";
                    return -1;
                }
            }

            /* [FIX-3D] Michaelis-Menten Vmax.  This is the tag the documentation
             * has always named; it was not being read at all.  Length is checked
             * rather than silently truncated -- adding a substrate and forgetting
             * to extend this list used to leave the new substrate at "no enzyme
             * limit", consumed to exhaustion every step, with no message. */
            /* [FIX-3D] Tag collision.  <maximum_uptake_flux> is ALSO read by
             * complab_functions.hh into vec_Vmax for the kinetics path, where it is in
             * mol substrate / mol cells / s -- not mmol/gDW/h.  A microbe with a
             * combined reaction_type (glpk_and_kinetics and friends) would have to
             * satisfy both meanings with one number, which is impossible.  Prefer a
             * dedicated tag; fall back to the shared one with a warning so existing
             * FBA-only input files keep working. */
            bool vmaxFromSharedTag = false;
            try {
                std::vector<T> v;
                try {
                    doc["parameters"]["microbiology"][bioname]["fba_maximum_uptake_flux"].read(v);
                } catch (PlbIOException& exception) {
                    doc["parameters"]["microbiology"][bioname]["maximum_uptake_flux"].read(v);
                    vmaxFromSharedTag = true;
                }
                if ((plint) v.size() != num_of_substrates) {
                    pcout << bioname << ": <maximum_uptake_flux> has " << (plint) v.size()
                          << " entries but there are " << num_of_substrates
                          << " substrates. Use 0 where the organism has no enzyme limit. Terminating.\n";
                    return -1;
                }
                for (plint s = 0; s < num_of_substrates; ++s) cfg.vmax[iM][s] = std::fabs(v[s]);
                if (vmaxFromSharedTag && rxntype::usesKinetics(rt)) {
                    pcout << "  WARNING: " << bioname << " uses BOTH kinetics and a metabolic solver, and\n"
                          << "           takes its Vmax from <maximum_uptake_flux>, which defineKinetics.hh\n"
                          << "           reads in different units (mol/mol cells/s, not mmol/gDW/h).\n"
                          << "           Give the metabolic solver its own <fba_maximum_uptake_flux>.\n";
                }
            } catch (PlbIOException& exception) {
                pcout << "  note: " << bioname << " has no <maximum_uptake_flux>; every substrate falls back "
                      << "to the supply-limited bound (eat what is locally available in one step).\n";
            }

            /* uptake floor / release ceiling, per substrate */
            try {
                std::vector<T> v;
                doc["parameters"]["microbiology"][bioname]["substrate_lower_bounds"].read(v);
                if ((plint) v.size() != num_of_substrates) {
                    pcout << bioname << ": <substrate_lower_bounds> has " << (plint) v.size()
                          << " entries but there are " << num_of_substrates << " substrates. Terminating.\n";
                    return -1;
                }
                for (plint s = 0; s < num_of_substrates; ++s) cfg.maxUptake[iM][s] = v[s];
            } catch (PlbIOException& exception) {}
            try {
                std::vector<T> v;
                doc["parameters"]["microbiology"][bioname]["substrate_upper_bounds"].read(v);
                if ((plint) v.size() != num_of_substrates) {
                    pcout << bioname << ": <substrate_upper_bounds> has " << (plint) v.size()
                          << " entries but there are " << num_of_substrates << " substrates. Terminating.\n";
                    return -1;
                }
                for (plint s = 0; s < num_of_substrates; ++s) cfg.maxRelease[iM][s] = v[s];
            } catch (PlbIOException& exception) {}

            /* objective direction */
            try {
                std::string tmp;
                doc["parameters"]["microbiology"][bioname]["objective_direction"].read(tmp);
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
                if      (tmp=="maximize" || tmp=="max") { cfg.objDir[iM] = -1; cfg.objDirExplicit[iM] = 1; }
                else if (tmp=="minimize" || tmp=="min") { cfg.objDir[iM] =  1; cfg.objDirExplicit[iM] = 1; }
                else {
                    pcout << bioname << ": <objective_direction> must be maximize or minimize. Terminating.\n";
                    return -1;
                }
            } catch (PlbIOException& exception) { cfg.objDir[iM] = -1; }

            /* ================================================================================
             *  <multi_step> -- the lexicographic chain
             * --------------------------------------------------------------------------------
             *      <multi_step>
             *          <stage_reactions>  BIOMASS EX_pyr_e EX_ac_e </stage_reactions>
             *          <retain_fraction>  0.6721  0.6848   1.0     </retain_fraction>
             *          <stage_direction>  max     max      max     </stage_direction>
             *      </multi_step>
             *
             *  Reactions are named, not numbered, because a column index is meaningless the
             *  moment the model file is regenerated.  They are resolved against the model's own
             *  reaction names in load_metabolic_models3D(), which is the first point at which
             *  those names exist; here we only record the request.
             *
             *  <retain_fraction> is what makes this calibration rather than pure mathematics --
             *  see the header of complab3d_lexicographic.hh.  It defaults to 1.0 for every
             *  stage, which is strict lexicographic optimisation and costs no fitted parameter.
             * ============================================================================== */
            try {
                std::vector<std::string> names;
                doc["parameters"]["microbiology"][bioname]["multi_step"]["stage_reactions"].read(names);

                if (names.size() < 2) {
                    pcout << bioname << ": <multi_step> needs at least two <stage_reactions>; "
                          << "one stage is ordinary FBA and the block should be omitted. Terminating.\n";
                    return -1;
                }

                std::vector<T> alphas;
                try { doc["parameters"]["microbiology"][bioname]["multi_step"]["retain_fraction"].read(alphas); }
                catch (PlbIOException& exception) {}
                if (alphas.empty()) alphas.assign(names.size(), (T) 1.0);
                if (alphas.size() != names.size()) {
                    pcout << bioname << ": <retain_fraction> has " << (plint) alphas.size()
                          << " entries but <stage_reactions> has " << (plint) names.size()
                          << ". Terminating.\n";
                    return -1;
                }
                for (size_t s = 0; s < alphas.size(); ++s) {
                    if (alphas[s] <= (T) 0 || alphas[s] > (T) 1) {
                        pcout << bioname << ": <retain_fraction> entry " << (plint) s << " is "
                              << alphas[s] << "; it must lie in (0, 1]. Terminating.\n";
                        return -1;
                    }
                }

                std::vector<std::string> dirs;
                try { doc["parameters"]["microbiology"][bioname]["multi_step"]["stage_direction"].read(dirs); }
                catch (PlbIOException& exception) {}
                if (dirs.empty()) dirs.assign(names.size(), std::string("max"));
                if (dirs.size() != names.size()) {
                    pcout << bioname << ": <stage_direction> has " << (plint) dirs.size()
                          << " entries but <stage_reactions> has " << (plint) names.size()
                          << ". Terminating.\n";
                    return -1;
                }

                cfg.lexPlan[iM].stages.clear();
                for (size_t s = 0; s < names.size(); ++s) {
                    std::string d = dirs[s];
                    std::transform(d.begin(), d.end(), d.begin(),
                                   [](unsigned char c){ return std::tolower(c); });
                    if (d != "max" && d != "maximize" && d != "min" && d != "minimize") {
                        pcout << bioname << ": <stage_direction> entry " << (plint) s
                              << " is \"" << dirs[s] << "\"; it must be max or min. Terminating.\n";
                        return -1;
                    }
                    const bool mx = (d == "max" || d == "maximize");
                    /* column stays -1 until the model's reaction names are known */
                    cfg.lexPlan[iM].stages.push_back(
                        complab_lex::LexStage(-1, (double) alphas[s], mx, names[s]));
                }
            } catch (PlbIOException& exception) {}

            /* ================================================================================
             *  <cybernetic> -- competing growth options on alternative carbon sources
             * --------------------------------------------------------------------------------
             *      <cybernetic>
             *          <sources>       lactate pyruvate acetate </sources>
             *          <substrate_ids>    0        1       2    </substrate_ids>
             *          <carbon_number>    3        3       2    </carbon_number>
             *          <uptake_kmax>   22.1      8.19    4.39   </uptake_kmax>
             *          <uptake_half_saturation> 0.02 0.02 0.02  </uptake_half_saturation>
             *          <weight_floor>  1e-3 </weight_floor>
             *      </cybernetic>
             *
             *  <substrate_ids> are indices into this simulation's substrate list, which is what
             *  ties each growth option to a concentration on the lattice.  Every source is
             *  solved with the OTHER sources closed, so the exclusion list is derived here
             *  rather than being asked for: writing it by hand is an invitation to leave one
             *  open, and an open second carbon source silently destroys the switching.
             * ============================================================================== */
            try {
                std::vector<std::string> snames;
                doc["parameters"]["microbiology"][bioname]["cybernetic"]["sources"].read(snames);

                if (snames.size() < 2) {
                    pcout << bioname << ": <cybernetic> needs at least two <sources>; with one "
                          << "growth option there is nothing to switch between. Terminating.\n";
                    return -1;
                }

                std::vector<plint> sids;
                std::vector<T>     ncarbon, kmax, khalf;
                doc["parameters"]["microbiology"][bioname]["cybernetic"]["substrate_ids"]  .read(sids);
                doc["parameters"]["microbiology"][bioname]["cybernetic"]["carbon_number"]  .read(ncarbon);
                doc["parameters"]["microbiology"][bioname]["cybernetic"]["uptake_kmax"]    .read(kmax);
                doc["parameters"]["microbiology"][bioname]["cybernetic"]["uptake_half_saturation"].read(khalf);

                if (sids.size() != snames.size() || ncarbon.size() != snames.size() ||
                    kmax.size() != snames.size() || khalf.size() != snames.size()) {
                    pcout << bioname << ": every list inside <cybernetic> must have one entry per "
                          << "source (" << (plint) snames.size() << " expected). Terminating.\n";
                    return -1;
                }
                for (size_t s = 0; s < sids.size(); ++s) {
                    if (sids[s] < 0 || sids[s] >= num_of_substrates) {
                        pcout << bioname << ": <substrate_ids> entry " << (plint) s << " is "
                              << sids[s] << ", outside the " << num_of_substrates
                              << " substrates declared. Terminating.\n";
                        return -1;
                    }
                    if (ncarbon[s] <= (T) 0) {
                        pcout << bioname << ": <carbon_number> entry " << (plint) s
                              << " must be positive. Terminating.\n";
                        return -1;
                    }
                    if (khalf[s] <= (T) 0) {
                        pcout << bioname << ": <uptake_half_saturation> entry " << (plint) s
                              << " must be positive. Terminating.\n";
                        return -1;
                    }
                }

                T floorv = (T) 1e-3;
                try { doc["parameters"]["microbiology"][bioname]["cybernetic"]["weight_floor"].read(floorv); }
                catch (PlbIOException& exception) {}
                if (floorv < (T) 0 || floorv >= (T) 1) {
                    pcout << bioname << ": <weight_floor> must lie in [0, 1). Terminating.\n";
                    return -1;
                }

                cfg.cybPlan[iM].sources.clear();
                cfg.cybPlan[iM].weightFloor = (double) floorv;
                for (size_t s = 0; s < snames.size(); ++s) {
                    complab_cyb::CyberneticSource src;
                    src.name      = snames[s];
                    src.subsIndex = (int) sids[s];
                    src.carbon    = (double) ncarbon[s];
                    src.kmax      = (double) kmax[s];
                    src.Ks        = (double) khalf[s];
                    /* everything that is not this source is closed while this source is solved */
                    for (size_t o = 0; o < sids.size(); ++o)
                        if (o != s) src.exclusive.push_back((int) sids[o]);
                    cfg.cybPlan[iM].sources.push_back(src);
                }
            } catch (PlbIOException& exception) {}

            /* optional extra bound overrides applied to the model at load time */
            try { doc["parameters"]["microbiology"][bioname]["constraint_indices"]     .read(cfg.constraint_loc[iM]); } catch (PlbIOException& exception) {}
            try { doc["parameters"]["microbiology"][bioname]["constraint_lower_bounds"].read(cfg.constraint_lb [iM]); } catch (PlbIOException& exception) {}
            try { doc["parameters"]["microbiology"][bioname]["constraint_upper_bounds"].read(cfg.constraint_ub [iM]); } catch (PlbIOException& exception) {}
            if (cfg.constraint_loc[iM].size() != cfg.constraint_lb[iM].size() ||
                cfg.constraint_loc[iM].size() != cfg.constraint_ub[iM].size()) {
                pcout << bioname << ": <constraint_indices>, <constraint_lower_bounds> and "
                      << "<constraint_upper_bounds> must be the same length. Terminating.\n";
                return -1;
            }

            /* [FIX-3D] CompLaB 2D stored equate_bounds in a single scalar that the
             * LAST microbe parsed overwrote for everyone, even though the tag lives
             * under <microbe i>.  It is per-microbe here, as it always should have been. */
            try {
                std::string tmp;
                doc["parameters"]["microbiology"][bioname]["equate_bounds"].read(tmp);
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
                cfg.equate_bounds[iM] = (tmp=="yes" || tmp=="true" || tmp=="1" || tmp=="on");
            } catch (PlbIOException& exception) { cfg.equate_bounds[iM] = false; }

            /* the metabolic model file itself */
            if (rxntype::usesFBA(rt)) {
                try { doc["parameters"]["microbiology"][bioname]["model_filename"].read(cfg.model_filename[iM]); }
                catch (PlbIOException& exception) {
                    /* [NEW] No longer fatal here.  <model_source> can supply the file instead, and
                     * that block is read later, by the integration layer.  An FBA microbe with
                     * neither is still an error -- it is just reported by
                     * load_metabolic_models3D(), which is the first point that knows both. */
                    cfg.model_filename[iM].clear();
                }
            }
        }
    }

    /* ---------------------------------------------- free / total concentration */
    if (cfg.useTotals) {
        if (!useEquilibrium) {
            pcout << "  note: <fba_concentration_basis> is total but the equilibrium solver is off. "
                  << "Every substrate is already its own total, so this changes nothing.\n";
            cfg.useTotals = false;
        } else {
            const plint mapped = cfg.totals.build(eq_stoich, eq_component_names, subs_names, num_of_substrates);
            pcout << "  metabolic layer: using TOTAL concentrations; " << mapped
                  << " of " << (plint) eq_component_names.size()
                  << " equilibrium components mapped onto a transported substrate.\n";
            if (mapped == 0) {
                pcout << "  WARNING: no component name matched a substrate name, so totals fall back to "
                      << "free concentrations. Check that <components> uses the same names as "
                      << "<name_of_substrates>.\n";
            }
        }
    }

    /* ------------------------------------------------------------- summary --- */
    pcout << "\n  Metabolic layer\n"
          << "    GLPK FBA        : " << (cfg.enable_fba_glpk    ? "on " : "off") << "  (" << cfg.glpk_count << " microbe(s))\n"
          << "    COBRApy FBA     : " << (cfg.enable_fba_cobrapy ? "on " : "off") << "  (" << cfg.cpy_count  << " microbe(s))\n"
          << "    Surrogate model : " << (cfg.enable_surrogate   ? "on " : "off") << "  (" << cfg.srg_count  << " microbe(s))\n"
          << "    Symbolic law    : " << (cfg.enable_symbolic   ? "on " : "off") << "  (" << cfg.sym_count  << " microbe(s))\n"
          << "    Graph network   : " << (cfg.enable_graphnet   ? "on " : "off") << "  (" << cfg.gnn_count  << " microbe(s))\n"
          << "    Thermo gate     : " << (cfg.enable_thermo       ? "on " : "off")
          << "  (multiplies whichever path each microbe uses)\n"
          << "    Uptake bounds   : " << (cfg.useTotals ? "TOTAL (speciation-aware)" : "FREE ion") << "\n";

    return 0;
}


/* ============================================================================
 *  load_metabolic_models3D
 *
 *  Reads one metabolic model per FBA microbe and applies the <constraint_*>
 *  overrides.  TWO file formats are accepted, told apart by their root element
 *  rather than by their extension:
 *
 *    <sbml>            a genome-scale model as it is distributed -- BiGG,
 *                      ModelSEED, anything with the FBC v2 package.  Read
 *                      natively by complab3d_sbml.hh.  This is the path that
 *                      makes extractMM.py optional.
 *
 *    <Metabolic_Model> the flattened matrix file extractMM.py produces:
 *                      nmet nrxn objLoc S b c lb ub, with S row-major and
 *                      metabolite-major, S[i][j] = flat[i*nrxn + j].  Kept
 *                      because every existing input file uses it.
 *
 *  <model_filename> may name the file with or without its extension.  Without,
 *  ".xml" is appended, which is what every CompLaB 2D and 3D input file has
 *  always assumed.
 *
 *  Returns 0 on success, -1 on failure.
 * ============================================================================
 */

/* Does `s` end with `suffix`?  Written out because this header is compiled as
 * C++98 in some build configurations and std::string::ends_with is C++20. */
inline bool endsWithCI (const std::string &s, const std::string &suffix)
{
    if (s.size() < suffix.size()) return false;
    for (size_t i = 0; i < suffix.size(); ++i) {
        char a = s[s.size() - suffix.size() + i];
        char b = suffix[i];
        if (a >= 'A' && a <= 'Z') a = (char)(a - 'A' + 'a');
        if (b >= 'A' && b <= 'Z') b = (char)(b - 'A' + 'a');
        if (a != b) return false;
    }
    return true;
}

/* <model_filename> -> a path to open.
 *
 * An absolute path, or one that already looks like a file, is taken as given;
 * a bare stem gets input_path in front and ".xml" behind.  That keeps every
 * existing input file working while letting <model_source> hand back a path
 * such as "input/e_coli_core.xml" that must NOT have ".xml" appended twice. */
inline std::string resolveModelPath (const std::string &input_path, const std::string &name)
{
    if (name.empty()) return name;
    const bool absolute = (name[0] == '/');
    const bool hasDir   = (name.find('/') != std::string::npos);
    const bool hasExt   = endsWithCI(name, ".xml") || endsWithCI(name, ".sbml")
                       || endsWithCI(name, ".xml.gz");
    if (absolute) return name;
    if (hasExt && hasDir) return name;             // already a usable relative path
    if (hasExt)           return input_path + name;
    return input_path + name + ".xml";
}

/* The half of the load that is the same whichever format the file was in:
 * override bounds, check the exchange columns, reshape S, store, report.
 * Returns 0 or -1.  `nmet`/`nrxn` are taken from the caller's own parse. */
inline int storeMetabolicModel (MetabolicConfig &cfg, plint iM, const std::string &fname,
                                int nmet, int nrxn, plint objLoc,
                                std::vector<T> &s1, std::vector<T> &b, std::vector<T> &c,
                                std::vector<T> &lb, std::vector<T> &ub)
{
    if ((plint) s1.size() != (plint) nmet * (plint) nrxn) {
        pcout << "  " << fname << ": <S> has " << (plint) s1.size()
              << " entries but nmet*nrxn = " << (plint) nmet * (plint) nrxn << ". Terminating.\n";
        return -1;
    }
    if ((plint) b.size()  != nmet) { pcout << "  " << fname << ": <b> must have nmet entries.\n";  return -1; }
    if ((plint) c.size()  != nrxn) { pcout << "  " << fname << ": <c> must have nrxn entries.\n";  return -1; }
    if ((plint) lb.size() != nrxn) { pcout << "  " << fname << ": <lb> must have nrxn entries.\n"; return -1; }
    if ((plint) ub.size() != nrxn) { pcout << "  " << fname << ": <ub> must have nrxn entries.\n"; return -1; }

    /* apply the per-microbe constraint overrides */
    for (size_t k = 0; k < cfg.constraint_loc[iM].size(); ++k) {
        const int j = cfg.constraint_loc[iM][k];
        if (j < 0 || j >= nrxn) {
            pcout << "  microbe" << iM << ": <constraint_indices> entry " << j
                  << " is outside 0.." << nrxn-1 << ". Terminating.\n";
            return -1;
        }
        lb[j] = cfg.constraint_lb[iM][k];
        ub[j] = cfg.constraint_ub[iM][k];
    }
    /* and sanity-check the exchange indices while we have nrxn to hand */
    for (plint s = 0; s < (plint) cfg.subsLoc[iM].size(); ++s) {
        const plint j = cfg.subsLoc[iM][s];
        if (j >= (plint) nrxn) {
            pcout << "  microbe" << iM << ": exchange reaction " << j
                  << " for substrate " << s << " is outside 0.." << nrxn-1 << ". Terminating.\n";
            return -1;
        }
    }

    /* reshape */
    std::vector< std::vector<T> > s2(nmet, std::vector<T>(nrxn, T()));
    for (int i = 0; i < nmet; ++i)
        for (int j = 0; j < nrxn; ++j)
            s2[i][j] = s1[(size_t) i * (size_t) nrxn + (size_t) j];

    cfg.S3[iM] = s2;
    cfg.S1[iM] = s1;
    cfg.vec_b[iM]  = b;
    cfg.vec_c[iM]  = c;
    cfg.vec_lb[iM] = lb;
    cfg.vec_ub[iM] = ub;
    cfg.vec_objLoc[iM] = objLoc;
    cfg.vec_nmets[iM]  = nmet;
    cfg.vec_nrxns[iM]  = nrxn;

    plint nz = 0;
    for (size_t k = 0; k < s1.size(); ++k) if (s1[k] != T()) ++nz;
    pcout << "    microbe" << iM << "  " << fname
          << " : " << nmet << " metabolites, " << nrxn << " reactions, "
          << nz << " nonzeros, objective at column " << objLoc << "\n";
    return 0;
}

inline int load_metabolic_models3D (MetabolicConfig &cfg, const std::string &input_path,
                                    const std::vector<plint> &reaction_type,
                                    plint num_of_microbes)
{
    pcout << "  loading metabolic models ...\n";

    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (!rxntype::usesFBA(reaction_type[iM])) continue;

        if (cfg.model_filename[iM].empty()) {
            pcout << "  microbe" << iM << ": reaction_type " << rxntype::name(reaction_type[iM])
                  << " needs a metabolic model, but neither <model_filename> nor <model_source>\n"
                  << "  supplied one. Terminating.\n";
            return -1;
        }

        const std::string fname = resolveModelPath(input_path, cfg.model_filename[iM]);

        /* ------------------------------------------------------------------ *
         *  Which format is this?  Decided by the root element, not the name.  *
         * ------------------------------------------------------------------ */
        const complab_sbml::ModelFileKind kind = complab_sbml::classify(fname);

        if (kind == complab_sbml::MODEL_UNKNOWN) {
            pcout << "  could not read the metabolic model " << fname << "\n"
                  << "  Its root element is neither <sbml> nor <Metabolic_Model>, or the file is\n"
                  << "  missing or not well-formed XML.\n"
                  << "  A gzipped model (.xml.gz) must be decompressed first; <model_source> does\n"
                  << "  that for you. Terminating.\n";
            return -1;
        }

        /* ================================================================== *
         *  SBML, read natively.                                              *
         * ================================================================== */
        if (kind == complab_sbml::MODEL_SBML) {
            complab_sbml::SbmlModel M = complab_sbml::readSbml(fname);
            if (!M.ok()) {
                pcout << "  could not read the SBML model " << fname << "\n  " << M.error
                      << "\n  Terminating.\n";
                return -1;
            }

            /* Resolve <exchange_reaction_names> now that the reaction list exists.  This is the
             * reason names are better than indices: a wrong name is caught here, by name, with
             * suggestions -- a wrong index is not caught at all. */
            if (!cfg.exchange_names[iM].empty()) {
                std::vector<int> cols;
                std::string err;
                if (!complab_sbml::resolveExchangeNames(M, cfg.exchange_names[iM], cols, err)) {
                    pcout << "  microbe" << iM << ", model " << fname << ":\n  " << err
                          << "  Terminating.\n";
                    return -1;
                }
                cfg.subsLoc[iM].assign(cols.size(), -1);
                for (size_t s = 0; s < cols.size(); ++s) cfg.subsLoc[iM][s] = (plint) cols[s];
                pcout << "    microbe" << iM << "  resolved " << (plint) cols.size()
                      << " exchange reaction name(s) against " << M.modelId << "\n";
            }

            /* Resolve <multi_step><stage_reactions> the same way, and for the same reason.  A
             * mistyped stage is far more damaging than a mistyped exchange: the chain would
             * silently skip that stage at run time and the flux vector would go back to being
             * whatever the simplex chose, which is exactly the ambiguity the block was added to
             * remove -- and nothing in the output would say so.  So it is fatal here. */
            if (!cfg.lexPlan[iM].stages.empty()) {
                for (size_t s = 0; s < cfg.lexPlan[iM].stages.size(); ++s) {
                    complab_lex::LexStage &st = cfg.lexPlan[iM].stages[s];
                    const int j = complab_sbml::findReaction(M, st.name);
                    if (j < 0) {
                        pcout << "  microbe" << iM << ", model " << fname << ":\n"
                              << "  <multi_step> stage " << (plint) s << " names reaction '"
                              << st.name << "', which is not in model '" << M.modelId << "'.\n"
                              << "  Terminating.\n";
                        return -1;
                    }
                    st.column = j;
                }
                /* Two stages on the same column would pin a quantity against itself and make
                 * every later stage vacuous.  Cheap to check, impossible to debug from output. */
                for (size_t a = 0; a < cfg.lexPlan[iM].stages.size(); ++a)
                    for (size_t b = a + 1; b < cfg.lexPlan[iM].stages.size(); ++b)
                        if (cfg.lexPlan[iM].stages[a].column == cfg.lexPlan[iM].stages[b].column) {
                            pcout << "  microbe" << iM << ": <multi_step> stages " << (plint) a
                                  << " and " << (plint) b << " both resolve to reaction '"
                                  << cfg.lexPlan[iM].stages[a].name << "'. Terminating.\n";
                            return -1;
                        }
                pcout << "    microbe" << iM << "  multi-step: "
                      << cfg.lexPlan[iM].describe() << "\n";
            }

            /* The objective sense the file itself declares.  <objective_direction> in CompLaB.xml
             * still wins if the user set it, because a user who overrides it means to. */
            if (!cfg.objDirExplicit[iM] && M.objSense == 1) cfg.objDir[iM] = 1;

            pcout << complab_sbml::sanityReport(M);

            if (storeMetabolicModel(cfg, iM, fname, M.nmet, M.nrxn, (plint) M.objLoc,
                                    M.S, M.b, M.c, M.lb, M.ub) != 0) return -1;
            continue;
        }

        /* ================================================================== *
         *  The extractMM.py matrix format, exactly as before.                *
         * ================================================================== */
        if (!cfg.exchange_names[iM].empty()) {
            pcout << "  microbe" << iM << ": <exchange_reaction_names> needs an SBML model, because\n"
                  << "  the matrix format produced by extractMM.py does not carry reaction names.\n"
                  << "  Point <model_filename> at the SBML file, or use <exchange_reaction_indices>.\n"
                  << "  Terminating.\n";
            return -1;
        }
        if (!cfg.lexPlan[iM].stages.empty()) {
            pcout << "  microbe" << iM << ": <multi_step> names reactions, and the matrix format\n"
                  << "  produced by extractMM.py does not carry reaction names.  Point\n"
                  << "  <model_filename> at the SBML file instead.\n"
                  << "  Terminating.\n";
            return -1;
        }
        try {
            XMLreader doc(fname);

            int nmet = 0, nrxn = 0;
            doc["Metabolic_Model"]["nmet"].read(nmet);
            doc["Metabolic_Model"]["nrxn"].read(nrxn);
            if (nmet <= 0 || nrxn <= 0) {
                pcout << "  " << fname << ": nmet and nrxn must both be positive. Terminating.\n";
                return -1;
            }

            plint objLoc = 0;
            try { doc["Metabolic_Model"]["objLoc"].read(objLoc); }
            catch (PlbIOException& exception) { objLoc = 0; }

            std::vector<T> s1, b, c, lb, ub;
            doc["Metabolic_Model"]["S"] .read(s1);
            doc["Metabolic_Model"]["b"] .read(b);
            doc["Metabolic_Model"]["c"] .read(c);
            doc["Metabolic_Model"]["lb"].read(lb);
            doc["Metabolic_Model"]["ub"].read(ub);

            /* [FIX-3D] CompLaB 2D never checked these lengths, so a truncated or
             * mis-generated model file silently produced garbage fluxes.  Check now,
             * loudly, at start-up instead of quietly at every voxel. */
            if (storeMetabolicModel(cfg, iM, fname, nmet, nrxn, objLoc, s1, b, c, lb, ub) != 0)
                return -1;
        }
        catch (PlbIOException& exception) {
            pcout << "  could not read the metabolic model " << fname << "\n"
                  << "  " << exception.what() << "\n"
                  << "  (produce it with:  python3 extractMM.py <model.sbml>, or point\n"
                  << "   <model_filename> straight at the SBML file)  Terminating.\n";
            return -1;
        }
    }
    return 0;
}


/* ============================================================================
 *  setup_metabolic_solvers
 *
 *  Builds the persistent solver state: one glp_prob* per GLPK microbe, or the
 *  cobra model objects for the COBRApy path.  Call once, after
 *  load_metabolic_models3D().
 * ============================================================================
 */
inline int setup_metabolic_solvers (MetabolicConfig &cfg,
                                    const std::vector<plint> &reaction_type,
                                    plint num_of_microbes,
                                    char *pyFileName, char *src_path)
{
    /* Silence unused-parameter warnings in the build configurations where one or
     * both back ends are compiled out.  Everything below is guarded. */
    (void) cfg; (void) reaction_type; (void) num_of_microbes;
    (void) pyFileName; (void) src_path;

#ifdef COMPLAB_ENABLE_GLPK
    if (cfg.glpk_count > 0) {
        cfg.vec_lp .assign(num_of_microbes, (glp_prob*) 0);
        cfg.sParam .resize(num_of_microbes);
        cfg.iParam .resize(num_of_microbes);
        cfg.method .assign(num_of_microbes, (int) 'S');
        cfg.isMIP  .assign(num_of_microbes, 0);
        cfg.ctype  .assign(num_of_microbes, (char*) 0);
        cfg.vtype  .assign(num_of_microbes, (char*) 0);

        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesGlpk(reaction_type[iM])) continue;

            const int nmet = cfg.vec_nmets[iM];
            const int nrxn = cfg.vec_nrxns[iM];

            /* every row is an equality, Sv = b, with b = 0 : the FBA steady state */
            cfg.ctype[iM] = (char*) calloc((size_t) nmet + 1, sizeof(char));
            cfg.vtype[iM] = (char*) calloc((size_t) nrxn + 1, sizeof(char));
            if (cfg.ctype[iM] == 0 || cfg.vtype[iM] == 0) {
                pcout << "  out of memory allocating GLPK type arrays for microbe" << iM << ". Terminating.\n";
                return -1;
            }
            memset(cfg.ctype[iM], 'S', (size_t) nmet);
            memset(cfg.vtype[iM], 'C', (size_t) nrxn);

            cfg.vec_lp[iM] = glp_create_prob();
            if (cfg.vec_lp[iM] == 0) {
                pcout << "  glp_create_prob failed for microbe" << iM << ". Terminating.\n";
                return -1;
            }

            int method = (int) 'S';
            int isMIP  = 0;
            const int erck = initialize_glpk((int) iM, cfg.vec_lp[iM],
                                             cfg.S3[iM], cfg.vec_b[iM], cfg.vec_c[iM],
                                             cfg.vec_lb[iM], cfg.vec_ub[iM],
                                             cfg.ctype[iM], cfg.vtype[iM],
                                             (int) cfg.objDir[iM], cfg.lpsolver, cfg.save_pb,
                                             method, isMIP, cfg.sParam[iM], cfg.iParam[iM]);
            if (erck != 0) {
                pcout << "  initialize_glpk failed for microbe" << iM << " (code " << erck << "). Terminating.\n";
                return -1;
            }
            cfg.method[iM] = method;
            cfg.isMIP [iM] = isMIP;
        }
    }
#endif

#ifdef COMPLAB_ENABLE_COBRAPY
    if (cfg.cpy_count > 0) {
        /* Start the embedded interpreter here rather than unconditionally at the
         * top of main() the way CompLaB 2D did.  2D called Py_Initialize() on
         * every rank of every run even when COBRApy was never used, which cost
         * memory and start-up time for nothing. */
        if (!Py_IsInitialized()) Py_Initialize();

        /* [FIX-3D] CRITICAL.  vec_model was declared but never sized, and
         * prep_cobrapy writes vec_model[j] with operator[] -- a heap write past a
         * zero-capacity allocation on the very first model.  Every other
         * per-microbe vector is assigned in initialize_metabolic; this one was
         * missed.  Size it here, where the count is known. */
        cfg.vec_model.assign((size_t) cfg.cpy_count, (PyObject*) 0);

        /* prep_cobrapy wants compacted per-model vectors, in microbe order */
        std::vector< std::vector<T> > S1c, cc, lbc, ubc;
        std::vector<int>   nrxnc, nmetc, sensec;
        std::vector<plint> objc;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesCobrapy(reaction_type[iM])) continue;
            /* [FIX-3D] <objective_direction> reached GLPK (as `sense`) but never
             * reached COBRApy, so a microbe asking to MINIMISE was silently
             * maximised on that back end.  Carry it through. */
            sensec.push_back((int) cfg.objDir[iM]);
            S1c .push_back(cfg.S1[iM]);
            cc  .push_back(cfg.vec_c[iM]);
            lbc .push_back(cfg.vec_lb[iM]);
            ubc .push_back(cfg.vec_ub[iM]);
            nrxnc.push_back(cfg.vec_nrxns[iM]);
            nmetc.push_back(cfg.vec_nmets[iM]);
            objc .push_back(cfg.vec_objLoc[iM]);
        }
        const int erck = prep_cobrapy(pyFileName, src_path, cfg.vec_model,
                                      S1c, cc, lbc, ubc, nrxnc, nmetc, objc, sensec);
        if (erck != 0) {
            pcout << "  prep_cobrapy failed (code " << erck << "). Terminating.\n";
            return -1;
        }

        /* Ship the lexicographic chains, once, now that the module is imported.
         *
         * NOTE ON THE INDEX.  prep_cobrapy() was handed a COMPACTED list -- only the microbes
         * whose reaction_type uses COBRApy -- so the Python side numbers its models 0..n-1 in
         * that compacted order, not by global microbe id.  The counter below has to walk the
         * same filter in the same order or a chain lands on the wrong organism, which would be
         * silent and would look like a modelling error rather than a wiring one. */
        plint slot = 0;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesCobrapy(reaction_type[iM])) continue;

            if (cfg.lexPlan[iM].active()) {
                std::vector<double> flat;
                flat.reserve(3 * cfg.lexPlan[iM].stages.size());
                for (size_t s = 0; s < cfg.lexPlan[iM].stages.size(); ++s) {
                    const complab_lex::LexStage &st = cfg.lexPlan[iM].stages[s];
                    flat.push_back((double) st.column);
                    flat.push_back(st.alpha);
                    flat.push_back(st.maximise ? 1.0 : 0.0);
                }
                /* The growth column travels with the chain so both back ends read the growth
                 * rate out of the same reaction of the same final vector.  Without it the two
                 * could disagree by tens of percent at any retain fraction below 1, and the
                 * GLPK-versus-COBRApy cross-check -- the only independent check either has --
                 * would be comparing two different quantities. */
                const int erm = register_multistep_cobrapy((int) slot, flat,
                                                           (int) cfg.vec_objLoc[iM]);
                if (erm != 0) {
                    pcout << "  microbe" << iM << ": could not register the <multi_step> chain "
                          << "with COBRApy (code " << erm << ").\n"
                          << "  Continuing would silently solve plain FBA instead. Terminating.\n";
                    return -1;
                }
                pcout << "    microbe" << iM << " (cobrapy slot " << slot << ")  multi-step: "
                      << cfg.lexPlan[iM].describe() << "\n";
            }
            ++slot;
        }
    }
#endif

    return 0;
}


/* ============================================================================
 *  release_metabolic_solvers -- call once at the end of main()
 * ============================================================================
 */
inline void release_metabolic_solvers (MetabolicConfig &cfg)
{
    if (cfg.released) return;      // idempotent: the scope guard may fire after an explicit call
    cfg.released = true;

#ifdef COMPLAB_ENABLE_GLPK
    for (size_t iM = 0; iM < cfg.vec_lp.size(); ++iM) {
        if (cfg.vec_lp[iM] != 0) { glp_delete_prob(cfg.vec_lp[iM]); cfg.vec_lp[iM] = 0; }
        if (cfg.ctype [iM] != 0) { free(cfg.ctype[iM]); cfg.ctype[iM] = 0; }
        if (cfg.vtype [iM] != 0) { free(cfg.vtype[iM]); cfg.vtype[iM] = 0; }
    }
    if (cfg.glpk_count > 0) glp_free_env();
#endif
#ifdef COMPLAB_ENABLE_COBRAPY
    if (cfg.cpy_count > 0) {
        finalize_cobrapy(&cfg.vec_model);
        if (Py_IsInitialized()) Py_FinalizeEx();
    }
#endif
    (void) cfg;
}

/* ============================================================================
 *  MetabolicScopeGuard
 *
 *  [FIX-3D] main() has fourteen early returns after the solvers are built --
 *  stability-check failures, the CA deadlock dumps, the flow-only exit -- and
 *  only one of them released the GLPK problems and shut the Python interpreter
 *  down.  Rather than patch fourteen call sites and hope the fifteenth
 *  remembers, declare one of these next to the config and every exit path is
 *  covered, including an exception unwinding out of Palabos.
 *
 *      MetabolicConfig    mmcfg;
 *      MetabolicScopeGuard mmguard(mmcfg);
 * ============================================================================
 */
struct MetabolicScopeGuard {
    MetabolicConfig &cfg;
    explicit MetabolicScopeGuard (MetabolicConfig &c) : cfg(c) {}
    ~MetabolicScopeGuard () { release_metabolic_solvers(cfg); }
private:
    MetabolicScopeGuard(const MetabolicScopeGuard&);
    MetabolicScopeGuard& operator=(const MetabolicScopeGuard&);
};

#endif  // COMPLAB3D_METABOLIC_HH
