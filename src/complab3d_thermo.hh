/* ================================================================================================
 *  complab3d_thermo.hh  --  THERMODYNAMIC CONTROL OF REACTION RATES
 *
 *  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
 *  Meile Lab, University of Georgia.  shahram.asgari@uga.edu
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT THIS IS FOR
 *
 *  Every rate path in CompLaB3D -- defineKinetics.hh, GLPK, COBRApy, the surrogate, a .sym law, a
 *  .gnn network -- answers the question "how fast CAN this organism go, given what is here".  None
 *  of them asks the other question: "is there enough energy in this reaction, at THESE local
 *  concentrations, for it to run at all".
 *
 *  In a pore-scale domain that second question has a different answer in different voxels.  A
 *  reaction that is strongly exergonic where the substrate is delivered can be at, or past, its
 *  thermodynamic threshold four voxels further in, once its own products have built up.  A Monod
 *  law cannot see that: it takes the substrate concentration and returns a positive rate whatever
 *  the product concentration is doing.  The result is a simulation that keeps a reaction running in
 *  a place where no organism could conserve energy from it.
 *
 *  This header adds the missing factor.  It computes, per voxel, how much energy the reaction
 *  actually releases at the local composition, subtracts what the organism must conserve as ATP,
 *  and returns a number F_T between 0 and 1 that multiplies the rate the other path produced.
 *  Far from the threshold F_T is 1 and nothing changes.  At the threshold it is 0 and the reaction
 *  stops.  Between the two it falls off smoothly.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE EQUATIONS, AND WHERE THEY COME FROM
 *
 *  The thermodynamic factor is the form of Jin & Bethke (2003, 2005), used for anaerobic oxidation
 *  of methane by He et al. (2021) and by Craig (2024, Ch. 3, Eqs. 3.3-3.6):
 *
 *      (1)  Q      = PROD_i  a_i ^ nu_i                             reaction quotient
 *      (2)  dG     = dG0 + R T ln Q                                 Craig Eq. 3.6
 *      (3)  f      = dG + m dG_ATP                                  Craig Eq. 3.4
 *      (4)  F_T    = max(0, 1 - exp( f / (chi R T) ))               Craig Eq. 3.3
 *
 *  SIGN CONVENTION.  dG here is negative when the reaction releases energy, which is the usual
 *  chemical convention and the one Palabos-side concentrations make natural.  Craig writes Eq. 3.4
 *  as f_x = -dG_net - m dG_ATP, i.e. as a positive "energy available", and correspondingly carries
 *  the minus sign inside the exponential.  The two are the same expression written twice; this file
 *  keeps dG negative-when-favourable so that a reader can check a printed dG against a table.
 *
 *  m dG_ATP is the ENERGY THRESHOLD: the amount of free energy the reaction must release before the
 *  organism can make its ATP.  With m = 0.25 mol ATP per mol reaction and dG_ATP = 50 kJ/mol the
 *  threshold is 12.5 kJ/mol, which is inside the 10-20 kJ/mol band reported for methanogens and
 *  sulfate reducers (Schink 1997; Hoehler et al. 2001).  chi is the average stoichiometric number,
 *  the number of times the rate-limiting step turns per reaction turnover; chi = 1 unless there is
 *  a measurement saying otherwise (Jin & Bethke 2003).
 *
 *  The energy-based growth yield is Heijnen & van Dijken (1992) as used by Craig (Eqs. 3.9-3.10):
 *
 *      (5)  dG_dis = 200 + 18 (6 - NoC)^1.8
 *                        + exp[ ((3.8 - gamma)^2)^0.16 * (3.6 + 0.4 NoC) ]     kJ / C-mol biomass
 *      (6)  Y      = dG / -( dG_ana + dG_dis )                        C-mol biomass / mol reaction
 *
 *  NoC is the number of carbon atoms in the carbon source and gamma its degree of reduction.
 *  dG_dis is a positive dissipation, so the denominator is negative and Y comes out positive.
 *
 *  Finally, for the two FBA paths, Craig Eq. 3.20 turns a biomass flux into something comparable
 *  with a measured growth efficiency:
 *
 *      (7)  gamma_FBA = nu_C_per_BM * F_BM / F_substrate
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHY A FILE RATHER THAN XML
 *
 *  A .thm file is read at run time, the same way a .sym rate law or a .gnn network is.  The
 *  reasoning is the same one: the standard free energy of a reaction and its energy threshold are
 *  scientific claims that get revised, compared, and quoted in a paper.  They belong in something a
 *  reader can open, diff and cite, not in a rebuild.  CompLaB.xml carries one switch and one file
 *  name; everything else is in the file.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE FILE FORMAT
 *
 *      # anything after a hash is a comment
 *      provenance   fitted by ... / taken from ...
 *      version      1
 *      temperature  277.15                     kelvin; the dG0 below must be quoted at this T
 *      energy_units kJ                         kJ (default) or J
 *      gas_constant 8.314e-3                   optional; defaults to kJ/mol/K
 *
 *      reaction AOM                            one block per reaction; the name is for the log
 *        microbe   ANME                        a name from <name_of_microbes>, an index, or all
 *        dG0       -33.12                      standard free energy at the stated temperature
 *        stoich    CH4 -1  SO4 -1  HS 1  HCO3 1
 *        atp       0.25                        mol ATP conserved per mol reaction
 *        dGatp     50.0                        free energy of ATP synthesis in situ
 *        chi       1.0                         average stoichiometric number
 *        cmin      1e-9                        activity floor, mol/L; keeps ln Q finite
 *        yield     on                          optional; off by default
 *        dGana     -30.0                       free energy of anabolism, kJ per C-mol biomass
 *        carbon    1                           NoC of the carbon source
 *        reduction 8.0                         degree of reduction of the carbon source
 *        dGdis     0                           optional; > 0 overrides Heijnen Eq. (5)
 *        yieldmax  0.2                         cap, C-mol biomass per mol reaction
 *
 *  A block runs until the next `reaction` line or the end of the file.  `stoich` takes name/number
 *  pairs, negative for reactants and positive for products; the names must be substrates of this
 *  run, and a name that is not is a hard error rather than a silently dropped term -- a missing
 *  product is the difference between a gate that closes and one that never does.
 *
 *  CONCENTRATIONS ARE USED AS ACTIVITIES.  This is the standard pore-scale simplification and it is
 *  what He et al. (2021) do.  Below `cmin` the value is floored, because a concentration that has
 *  been driven to zero by the transport solver would otherwise send ln Q to minus infinity and make
 *  the gate report an energy the system does not have.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT THIS HEADER DOES NOT DO
 *
 *  It does not model electron transfer.  Craig's Eqs. 3.5, 3.7, 3.8, 3.17 and 3.18 subtract
 *  activation and ohmic losses along a conductive nanowire network, which needs a solver for the
 *  cytochrome redox state that CompLaB3D does not have.  Setting those losses to zero is the
 *  assumption of direct contact, which is what every one-organism case in this repository is.  The
 *  file format leaves room for the term -- `dGloss` is parsed and subtracted -- so a user who has
 *  computed a loss elsewhere can supply it as a constant.
 *
 *  It gates BIOTIC rate paths only: defineKinetics.hh, GLPK, COBRApy, the surrogate, a .sym law and
 *  a .gnn network.  The abiotic paths are left alone, because an abiotic reaction in CompLaB3D is
 *  either a kinetic law the user wrote, in which case it is theirs to gate, or an equilibrium, in
 *  which case the equilibrium solver has already answered the same question exactly.
 *
 *  ONE CAVEAT ON THE defineKinetics.hh PATH.  That path computes one combined substrate rate vector
 *  for every organism in a voxel at once, so there is no way to give two organisms two different
 *  gates there.  A .thm file with more than one reaction block is therefore refused when the
 *  kinetics path is in use.  The other five paths compute rates per organism and have no such
 *  restriction.
 *
 *  It needs no Palabos and no solver, so tests/test_thermo.cpp compiles it on its own.
 * ================================================================================================ */

#ifndef COMPLAB3D_THERMO_HH
#define COMPLAB3D_THERMO_HH

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace complab_thermo {

/* ================================================================================================
 *  ONE REACTION'S ENERGETICS
 * ================================================================================================ */
struct Spec {
    std::string name;                       /* label for the log */
    std::string microbeName;                /* as written in the file: a name, an index, or "all" */
    int         microbe;                    /* resolved global id; -1 means every organism */

    std::vector<std::string> species;       /* the stoich line, in file order */
    std::vector<double>      nu;            /* negative reactant, positive product */
    std::vector<int>         subsOfSpecies; /* resolved substrate index, filled by bind() */

    double dG0;        /* kJ/mol, at the file's temperature */
    double atp;        /* mol ATP per mol reaction        (m in Eq. 3) */
    double dGatp;      /* kJ/mol ATP                                    */
    double chi;        /* average stoichiometric number   (chi in Eq. 4) */
    double cmin;       /* activity floor, mol/L */
    double dGloss;     /* kJ/mol subtracted from dG0 + RT lnQ; see the note above. >= 0 */

    bool   yieldOn;
    double dGana;      /* kJ per C-mol biomass */
    double NoC;        /* carbon number of the carbon source */
    double reduction;  /* degree of reduction of the carbon source */
    double dGdisFixed; /* > 0 overrides Heijnen */
    double yieldMax;   /* cap on Y, C-mol biomass per mol reaction */

    Spec()
      : microbe(-1), dG0(0.0), atp(0.0), dGatp(50.0), chi(1.0), cmin(1e-9), dGloss(0.0),
        yieldOn(false), dGana(0.0), NoC(1.0), reduction(4.0), dGdisFixed(0.0), yieldMax(1.0) {}

    bool valid() const { return !species.empty() && species.size() == nu.size(); }
};

/* ================================================================================================
 *  THE WHOLE FILE
 * ================================================================================================ */
struct Model {
    std::string provenance;
    std::string unitName;      /* "kJ" or "J" as written */
    double      unitScale;     /* multiply every energy in the file by this to get kJ/mol */
    double      Tkelvin;
    double      Rgas;          /* kJ / mol / K */
    /* Activities are read straight off the lattices, and the lattices carry whatever unit the
     * case was written in.  ln Q is not scale-free, so a case posed in mmol/L and one posed in
     * mol/L give different energies from the same numbers.  This factor converts the solver's
     * value to mol/L and appears in the log, so the assumption is stated rather than assumed. */
    double      concScale;
    std::vector<Spec> specs;

    Model() : unitName("kJ"), unitScale(1.0), Tkelvin(298.15), Rgas(8.314462618e-3),
              concScale(1.0) {}

    bool valid() const {
        if (specs.empty()) return false;
        for (size_t i = 0; i < specs.size(); ++i) if (!specs[i].valid()) return false;
        return true;
    }
    double RT() const { return Rgas * Tkelvin; }
};

/* ================================================================================================
 *  THE ARITHMETIC
 *
 *  Each of the four steps is its own function so that a test can check it on its own and so that
 *  the log can print the intermediate a reader would want to see.  None of them touches the file
 *  parser or the runtime, so they are pure and hand-checkable.
 * ================================================================================================ */

/* Eq. (1), as a logarithm: ln Q = SUM_i nu_i ln a_i, with a_i floored at cmin. */
inline double lnQ(const Spec &s, const std::vector<double> &conc)
{
    double acc = 0.0;
    for (size_t i = 0; i < s.nu.size(); ++i) {
        const int iS = (i < s.subsOfSpecies.size()) ? s.subsOfSpecies[i] : -1;
        double a = (iS >= 0 && iS < (int) conc.size()) ? conc[(size_t) iS] : 0.0;
        if (!(a > s.cmin)) a = s.cmin;                 /* also catches NaN */
        acc += s.nu[i] * std::log(a);
    }
    return acc;
}

/* Eq. (2), with the optional constant electron-transfer loss of the note above. */
inline double deltaG(const Spec &s, const Model &M, const std::vector<double> &conc)
{
    return s.dG0 + M.RT() * lnQ(s, conc) + s.dGloss;
}

/* Eq. (3).  Negative means there is energy left over after the ATP has been paid for. */
inline double driving(const Spec &s, double dG) { return dG + s.atp * s.dGatp; }

/* Eq. (4).  Returns a number in [0, 1].  dGout, when given, receives the dG that produced it, so a
 * caller can report the energy alongside the factor without computing ln Q twice. */
inline double factor(const Spec &s, const Model &M, const std::vector<double> &conc,
                     double *dGout = 0)
{
    const double dG = deltaG(s, M, conc);
    if (dGout) *dGout = dG;

    const double f   = driving(s, dG);
    const double chi = (s.chi > 0.0) ? s.chi : 1.0;
    const double x   = f / (chi * M.RT());

    if (x >= 0.0) return 0.0;                 /* at or past the threshold: nothing runs */
    if (x < -700.0) return 1.0;               /* exp() would underflow; the answer is 1 */
    const double F = 1.0 - std::exp(x);
    if (F < 0.0) return 0.0;
    if (F > 1.0) return 1.0;
    return F;
}

/* Eq. (5), Heijnen & van Dijken (1992).  kJ per C-mol biomass, always positive.
 *
 * [FIX] The first term used to be pow(6 - NoC, 1.8) unguarded.  Above six carbons the base goes
 * negative, the exponent is fractional, and std::pow returns NaN -- so glucose (6 C) was the
 * last carbon source this worked for, and sucrose or a fatty acid printed `dG_dis nan` in the
 * log and silently returned a zero yield.  The reduction term two lines down squares its
 * argument for exactly this reason, so the omission was on the carbon term alone.
 *
 * The correlation is |6 - NoC| in the original: it is a distance from the six-carbon reference,
 * and Heijnen tabulates it for carbon sources on both sides of six. */
inline double dissipation(double NoC, double reduction)
{
    const double a = 200.0 + 18.0 * std::pow(std::fabs(6.0 - NoC), 1.8);
    const double d = (3.8 - reduction) * (3.8 - reduction);
    const double b = std::exp(std::pow(d, 0.16) * (3.6 + 0.4 * NoC));
    return a + b;
}

/* Eq. (6).  dG is the value from deltaG(); the result is C-mol biomass per mol reaction. */
inline double yieldOf(const Spec &s, double dG)
{
    if (!s.yieldOn) return 0.0;
    const double dis = (s.dGdisFixed > 0.0) ? s.dGdisFixed : dissipation(s.NoC, s.reduction);
    const double den = -(s.dGana + dis);
    if (den >= 0.0) return 0.0;               /* anabolism cheaper than free: not physical */
    double Y = dG / den;
    if (!(Y > 0.0)) return 0.0;
    if (Y > s.yieldMax) Y = s.yieldMax;
    return Y;
}

/* Eq. (7), Craig Eq. 3.20.  A reporting quantity for the FBA paths: the fraction of the substrate
 * carbon that ends up in biomass.  F_BM and F_sub are fluxes in the same units; nuC is the number
 * of carbon atoms per unit of biomass in the model's biomass reaction. */
inline double fbaGrowthEfficiency(double nuC, double F_BM, double F_sub)
{
    const double d = std::fabs(F_sub);
    if (!(d > 0.0)) return 0.0;
    return nuC * std::fabs(F_BM) / d;
}

/* ================================================================================================
 *  READING A .thm FILE
 *
 *  Strict on purpose.  An unknown keyword, a missing dG0, a stoich line with an odd number of
 *  tokens: all stop the run with the line number.  A gate that silently defaults to 1 everywhere is
 *  indistinguishable in the output from not having a gate at all, which is the worst way for this
 *  feature to fail.
 * ================================================================================================ */
namespace detail {

inline std::vector<std::string> tokens(const std::string &line)
{
    std::vector<std::string> t;
    std::string cur;
    for (size_t i = 0; i <= line.size(); ++i) {
        const char c = (i < line.size()) ? line[i] : ' ';
        if (c == ' ' || c == '\t' || c == '\r') { if (!cur.empty()) { t.push_back(cur); cur.clear(); } }
        else cur += c;
    }
    return t;
}

inline bool isNumber(const std::string &s, double &out)
{
    if (s.empty()) return false;
    char *end = 0;
    out = std::strtod(s.c_str(), &end);
    return end != 0 && *end == '\0';
}

inline bool truthy(std::string s)
{
    for (size_t i = 0; i < s.size(); ++i)
        if (s[i] >= 'A' && s[i] <= 'Z') s[i] = (char) (s[i] - 'A' + 'a');
    return s == "on" || s == "yes" || s == "true" || s == "1";
}

inline bool falsy(std::string s)
{
    for (size_t i = 0; i < s.size(); ++i)
        if (s[i] >= 'A' && s[i] <= 'Z') s[i] = (char) (s[i] - 'A' + 'a');
    return s == "off" || s == "no" || s == "false" || s == "0";
}

} /* namespace detail */

inline bool load(Model &M, const std::string &path, std::string *err = 0)
{
    M = Model();
    std::ifstream in(path.c_str());
    if (!in) { if (err) *err = "cannot open '" + path + "'"; return false; }

    /* The energies in the file are read as written and scaled at the END of the load, once the
     * energy_units line has certainly been seen.  Scaling as we go would make the file
     * order-dependent, and a units line at the bottom is not a user error worth punishing. */
    bool sawUnits = false;
    Spec *cur = 0;
    std::string line;
    int lineNo = 0;

    while (std::getline(in, line)) {
        ++lineNo;
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        std::vector<std::string> t = detail::tokens(line);
        if (t.empty()) continue;

        char where[64];
        std::sprintf(where, " (line %d)", lineNo);
        const std::string at(where);
        const std::string key = t[0];

        /* ---- file-level keys ------------------------------------------------------------- */
        if (key == "provenance") {
            std::string p;
            for (size_t i = 1; i < t.size(); ++i) { if (i > 1) p += " "; p += t[i]; }
            M.provenance = p;
            continue;
        }
        if (key == "version") continue;                       /* accepted, currently unused */

        if (key == "temperature") {
            double v;
            if (t.size() != 2 || !detail::isNumber(t[1], v) || !(v > 0.0)) {
                if (err) *err = "temperature needs one positive number in kelvin" + at;
                return false;
            }
            M.Tkelvin = v; continue;
        }
        if (key == "energy_units") {
            if (t.size() != 2) { if (err) *err = "energy_units needs one word" + at; return false; }
            if      (t[1] == "kJ" || t[1] == "kj") { M.unitScale = 1.0;   M.unitName = "kJ"; }
            else if (t[1] == "J"  || t[1] == "j")  { M.unitScale = 1e-3;  M.unitName = "J";  }
            else { if (err) *err = "energy_units '" + t[1] + "' is not kJ or J" + at; return false; }
            sawUnits = true; continue;
        }
        if (key == "concentration_scale") {
            double v;
            if (t.size() != 2 || !detail::isNumber(t[1], v) || !(v > 0.0)) {
                if (err) *err = "concentration_scale needs one positive number" + at;
                return false;
            }
            M.concScale = v; continue;
        }
        if (key == "gas_constant") {
            double v;
            if (t.size() != 2 || !detail::isNumber(t[1], v) || !(v > 0.0)) {
                if (err) *err = "gas_constant needs one positive number" + at;
                return false;
            }
            M.Rgas = v; continue;
        }

        /* ---- a new reaction block --------------------------------------------------------- */
        if (key == "reaction") {
            if (t.size() != 2) { if (err) *err = "reaction needs exactly one name" + at; return false; }
            for (size_t i = 0; i < M.specs.size(); ++i)
                if (M.specs[i].name == t[1]) {
                    if (err) *err = "two reactions are both named '" + t[1] + "'" + at;
                    return false;
                }
            Spec s; s.name = t[1];
            M.specs.push_back(s);
            cur = &M.specs.back();
            continue;
        }

        if (cur == 0) {
            if (err) *err = "'" + key + "' appears before any reaction block" + at;
            return false;
        }
        /* M.specs may have reallocated; re-point at the last block. */
        cur = &M.specs.back();

        /* ---- per-reaction keys ------------------------------------------------------------ */
        if (key == "stoich") {
            if (t.size() < 3 || ((t.size() - 1) % 2) != 0) {
                if (err) *err = "stoich needs name/number pairs" + at;
                return false;
            }
            for (size_t i = 1; i + 1 < t.size(); i += 2) {
                double v;
                if (!detail::isNumber(t[i+1], v)) {
                    if (err) *err = "'" + t[i+1] + "' is not a number in the stoich line" + at;
                    return false;
                }
                if (v == 0.0) {
                    if (err) *err = "species '" + t[i] + "' has a zero coefficient" + at;
                    return false;
                }
                cur->species.push_back(t[i]);
                cur->nu.push_back(v);
            }
            continue;
        }
        if (key == "microbe") {
            if (t.size() != 2) { if (err) *err = "microbe needs one name, index or 'all'" + at; return false; }
            cur->microbeName = t[1];
            continue;
        }
        if (key == "yield") {
            if (t.size() != 2) { if (err) *err = "yield needs on or off" + at; return false; }
            if      (detail::truthy(t[1])) cur->yieldOn = true;
            else if (detail::falsy (t[1])) cur->yieldOn = false;
            else { if (err) *err = "yield '" + t[1] + "' is not on or off" + at; return false; }
            continue;
        }

        /* every remaining key is one number */
        {
            const char *names[12] = { "dG0","atp","dGatp","chi","cmin","dGloss",
                                      "dGana","carbon","reduction","dGdis","yieldmax", 0 };
            double *slots[12]     = { &cur->dG0, &cur->atp, &cur->dGatp, &cur->chi, &cur->cmin,
                                      &cur->dGloss, &cur->dGana, &cur->NoC, &cur->reduction,
                                      &cur->dGdisFixed, &cur->yieldMax, 0 };
            int hit = -1;
            for (int k = 0; names[k]; ++k) if (key == names[k]) { hit = k; break; }
            if (hit < 0) {
                if (err) *err = "'" + key + "' is not a keyword this file understands" + at;
                return false;
            }
            double v;
            if (t.size() != 2 || !detail::isNumber(t[1], v)) {
                if (err) *err = key + " needs exactly one number" + at;
                return false;
            }
            *slots[hit] = v;
            continue;
        }
    }

    if (M.specs.empty()) { if (err) *err = "the file declares no reaction blocks"; return false; }

    for (size_t i = 0; i < M.specs.size(); ++i) {
        Spec &s = M.specs[i];
        if (s.species.empty()) {
            if (err) *err = "reaction '" + s.name + "' has no stoich line";
            return false;
        }
        if (!(s.chi > 0.0)) {
            if (err) *err = "reaction '" + s.name + "' has chi <= 0";
            return false;
        }
        if (!(s.cmin > 0.0)) {
            if (err) *err = "reaction '" + s.name + "' has cmin <= 0, so ln Q can be -inf";
            return false;
        }
        if (s.dGloss < 0.0) {
            if (err) *err = "reaction '" + s.name + "' has a negative dGloss; a transfer loss "
                            "makes a reaction less favourable, so it is a positive number";
            return false;
        }
        /* Scale the energies into kJ/mol.  cmin, chi, atp, carbon and reduction are not energies. */
        if (M.unitScale != 1.0) {
            s.dG0 *= M.unitScale;  s.dGatp *= M.unitScale;  s.dGloss *= M.unitScale;
            s.dGana *= M.unitScale; s.dGdisFixed *= M.unitScale;
        }
    }
    (void) sawUnits;
    return true;
}

/* ================================================================================================
 *  BINDING TO THIS RUN
 *
 *  A .thm file names species and organisms as strings; the solver addresses them by index.  A file
 *  written against a different <name_of_substrates> ordering must be refused rather than bound to
 *  the wrong lattice, exactly as in complab3d_symbolic.hh.
 *
 *  Returns an empty string on success, or the message to print before terminating.
 * ================================================================================================ */
inline std::string bindToRun(Model &M,
                             const std::vector<std::string> &subsNames,
                             const std::vector<std::string> &microbeNames)
{
    if (!M.valid()) return "the file declares no usable reaction blocks.";

    for (size_t i = 0; i < M.specs.size(); ++i) {
        Spec &s = M.specs[i];

        s.subsOfSpecies.assign(s.species.size(), -1);
        for (size_t k = 0; k < s.species.size(); ++k) {
            int hit = -1;
            for (size_t j = 0; j < subsNames.size(); ++j)
                if (subsNames[j] == s.species[k]) { hit = (int) j; break; }
            if (hit < 0) {
                std::string m = "reaction '" + s.name + "' uses species '" + s.species[k]
                              + "', which is not in <name_of_substrates>. This run has:";
                for (size_t j = 0; j < subsNames.size(); ++j) m += " " + subsNames[j];
                m += ".";
                return m;
            }
            s.subsOfSpecies[k] = hit;
        }

        /* microbe: a name, an index, "all", or absent (which means all) */
        if (s.microbeName.empty() || s.microbeName == "all" || s.microbeName == "*") {
            s.microbe = -1;
        } else {
            int hit = -1;
            for (size_t j = 0; j < microbeNames.size(); ++j)
                if (microbeNames[j] == s.microbeName) { hit = (int) j; break; }
            if (hit < 0) {
                double v;
                if (detail::isNumber(s.microbeName, v)
                    && v >= 0.0 && v < (double) microbeNames.size()
                    && v == std::floor(v)) {
                    hit = (int) v;
                } else {
                    std::string m = "reaction '" + s.name + "' names microbe '" + s.microbeName
                                  + "', which is not in <name_of_microbes>. This run has:";
                    for (size_t j = 0; j < microbeNames.size(); ++j) m += " " + microbeNames[j];
                    m += ".";
                    return m;
                }
            }
            s.microbe = hit;
        }
    }

    /* Two blocks gating the same organism would multiply, which is never what a user means: a
     * second reaction for the same organism is a modelling choice that needs its own rate path, not
     * a second gate on the first one. */
    for (size_t i = 0; i < M.specs.size(); ++i)
        for (size_t j = i + 1; j < M.specs.size(); ++j)
            if (M.specs[i].microbe == M.specs[j].microbe) {
                char b[64]; std::sprintf(b, "%d", M.specs[i].microbe);
                return "reactions '" + M.specs[i].name + "' and '" + M.specs[j].name
                     + "' both gate microbe " + std::string(M.specs[i].microbe < 0 ? "all" : b)
                     + ". One organism gets one gate.";
            }

    /* [FIX] The loop above compares RESOLVED ids, so `all` (-1) never collides with a specific
     * organism and a file carrying both passed.  registerModel() then wrote the blocks in file
     * order and let the `all` branch overwrite every slot, so a block written for one organism
     * was silently discarded whenever a catch-all appeared after it -- and that organism was
     * then gated by the wrong reaction's dG0, stoichiometry and threshold, in every voxel of
     * every step, with describe() still printing both blocks as though both applied.
     *
     * Mixing the two is ambiguous however it is resolved, so refuse it and say which blocks
     * collide.  A file that wants a default plus an exception should name every organism. */
    {
        int nAll = 0; size_t iAll = 0;
        for (size_t i = 0; i < M.specs.size(); ++i)
            if (M.specs[i].microbe < 0) { ++nAll; iAll = i; }
        if (nAll > 0 && M.specs.size() > (size_t) nAll) {
            for (size_t i = 0; i < M.specs.size(); ++i)
                if (M.specs[i].microbe >= 0)
                    return "reaction '" + M.specs[iAll].name + "' gates 'all' while reaction '"
                         + M.specs[i].name + "' gates one organism by name. Which one applies to "
                         "that organism would depend on the order the blocks appear in the file, "
                         "so this is refused: name every organism explicitly instead.";
        }
    }
    return std::string();
}

/* ================================================================================================
 *  THE RUNTIME
 *
 *  One model per run, plus a per-microbe lookup so the inner loop is an array index rather than a
 *  string comparison, plus the counters the end-of-run report needs.
 * ================================================================================================ */
struct Runtime {
    Model model;
    bool  active;
    std::vector<int> specOfMicrobe;   /* global microbe id -> index into model.specs, -1 = none */

    long   evaluations;
    long   blocked;                   /* F_T came out exactly 0 */
    double sumF;
    double minF;
    double maxF;

    Runtime() : active(false), evaluations(0), blocked(0), sumF(0.0), minF(1.0), maxF(0.0) {}
};

inline Runtime &runtime() { static Runtime R; return R; }

/* Fill specOfMicrobe from a bound model.  Call once, after bindToRun(). */
inline void registerModel(int nMicrobes)
{
    Runtime &R = runtime();
    R.specOfMicrobe.assign((size_t) (nMicrobes > 0 ? nMicrobes : 0), -1);
    for (size_t i = 0; i < R.model.specs.size(); ++i) {
        const int m = R.model.specs[i].microbe;
        if (m < 0) {
            for (size_t k = 0; k < R.specOfMicrobe.size(); ++k) R.specOfMicrobe[k] = (int) i;
        } else if (m < (int) R.specOfMicrobe.size()) {
            R.specOfMicrobe[(size_t) m] = (int) i;
        }
    }
    R.active = !R.model.specs.empty();
}

inline bool enabled() { return runtime().active; }

/* ------------------------------------------------------------------------------------------------
 *  THE ONE CALL THE RATE PROCESSORS MAKE
 *
 *  Returns the factor to multiply this organism's whole rate vector by.  An organism with no block
 *  in the file returns 1, so mixing gated and ungated organisms in one run is allowed and costs
 *  the ungated ones one array lookup.
 *
 *  `conc` is the same vector the rate path used -- free or total, whichever <fba_concentration_basis>
 *  selected -- so the gate sees exactly the chemistry the rate law saw.
 * ---------------------------------------------------------------------------------------------- */
/* Copy the solver's concentration vector into the double vector the gate wants, applying
 * concentration_scale.  Written as a template so it takes Palabos's T without this header having
 * to know what T is. */
template <class Vec>
inline void fillConc(const Vec &avail, int n, std::vector<double> &out)
{
    const double s = runtime().model.concScale;
    out.assign((size_t) (n > 0 ? n : 0), 0.0);
    for (int i = 0; i < n; ++i) out[(size_t) i] = (double) avail[(size_t) i] * s;
}

inline double gateFor(int globalMicrobe, const std::vector<double> &conc, double *dGout = 0)
{
    Runtime &R = runtime();
    if (!R.active) return 1.0;
    if (globalMicrobe < 0 || globalMicrobe >= (int) R.specOfMicrobe.size()) return 1.0;
    const int si = R.specOfMicrobe[(size_t) globalMicrobe];
    if (si < 0) return 1.0;

    const double F = factor(R.model.specs[(size_t) si], R.model, conc, dGout);

    ++R.evaluations;
    R.sumF += F;
    if (F < R.minF) R.minF = F;
    if (F > R.maxF) R.maxF = F;
    if (F <= 0.0) ++R.blocked;
    return F;
}

/* ------------------------------------------------------------------------------------------------
 *  [v1.3.2]  THE SAME ARITHMETIC, WITHOUT THE COUNTERS
 *
 *  The energy snapshot walks every fluid voxel once per output interval and asks for dG and F_T
 *  there.  It must not go through gateFor().  That call feeds the end-of-run line
 *
 *      [THM] N gate evaluations, mean F_T ..., range ...
 *
 *  which is meant to describe the reaction over the run.  A snapshot pass would add millions of
 *  evaluations that no rate law ever used, and the mean would then depend on how often the user
 *  asked for output -- a statistic about the write schedule wearing the name of a statistic about
 *  the chemistry.  So these two do the identical arithmetic and touch no counter.
 *
 *  specOf() also answers "does this organism have a block in the .thm file at all", which is what
 *  decides whether a pair of energy fields is allocated for it.
 * ---------------------------------------------------------------------------------------------- */
inline int specOf(int globalMicrobe)
{
    Runtime &R = runtime();
    if (!R.active) return -1;
    if (globalMicrobe < 0 || globalMicrobe >= (int) R.specOfMicrobe.size()) return -1;
    return R.specOfMicrobe[(size_t) globalMicrobe];
}

/* false when this organism is ungated, in which case dGout and ftOut are left alone. */
inline bool gateProbe(int globalMicrobe, const std::vector<double> &conc,
                      double &dGout, double &ftOut)
{
    const int si = specOf(globalMicrobe);
    if (si < 0) return false;
    Runtime &R = runtime();
    ftOut = factor(R.model.specs[(size_t) si], R.model, conc, &dGout);
    return true;
}

/* The energy yield for one organism, or a negative number when this organism has no yield block --
 * which the caller reads as "keep whatever yield you were going to use". */
inline double yieldFor(int globalMicrobe, const std::vector<double> &conc)
{
    Runtime &R = runtime();
    if (!R.active) return -1.0;
    if (globalMicrobe < 0 || globalMicrobe >= (int) R.specOfMicrobe.size()) return -1.0;
    const int si = R.specOfMicrobe[(size_t) globalMicrobe];
    if (si < 0) return -1.0;
    const Spec &s = R.model.specs[(size_t) si];
    if (!s.yieldOn) return -1.0;
    return yieldOf(s, deltaG(s, R.model, conc));
}

/* ------------------------------------------------------------------------------------------------
 *  WHAT THE LOG SAYS
 * ---------------------------------------------------------------------------------------------- */
inline std::string describe(const Model &M, const std::string &path)
{
    char b[256];
    std::string s = "  [THM] energetics from '" + path + "'\n";
    if (!M.provenance.empty()) s += "  [THM]   provenance: " + M.provenance + "\n";
    std::sprintf(b, "  [THM]   T = %.2f K, RT = %.4f kJ/mol, energies read in %s\n",
                 M.Tkelvin, M.RT(), M.unitName.c_str());
    s += b;
    std::sprintf(b, "  [THM]   activities = lattice concentration x %g, taken as mol/L\n",
                 M.concScale);
    s += b;
    for (size_t i = 0; i < M.specs.size(); ++i) {
        const Spec &sp = M.specs[i];
        s += "  [THM]   " + sp.name + ":";
        for (size_t k = 0; k < sp.species.size(); ++k) {
            std::sprintf(b, " %+g %s", sp.nu[k], sp.species[k].c_str());
            s += b;
        }
        s += "\n";
        std::sprintf(b, "  [THM]     dG0 %.4g kJ/mol, threshold %.4g kJ/mol (%g ATP x %g), chi %g\n",
                     sp.dG0, sp.atp * sp.dGatp, sp.atp, sp.dGatp, sp.chi);
        s += b;
        if (sp.dGloss > 0.0) {
            std::sprintf(b, "  [THM]     electron transfer loss %.4g kJ/mol added to dG\n", sp.dGloss);
            s += b;
        }
        if (sp.yieldOn) {
            const double dis = (sp.dGdisFixed > 0.0) ? sp.dGdisFixed
                                                     : dissipation(sp.NoC, sp.reduction);
            std::sprintf(b, "  [THM]     yield on: dG_ana %.4g, dG_dis %.4g kJ/C-mol, cap %.4g\n",
                         sp.dGana, dis, sp.yieldMax);
            s += b;
        }
        s += "  [THM]     applies to: ";
        if (sp.microbe < 0) s += "every organism\n";
        else { std::sprintf(b, "microbe %d (%s)\n", sp.microbe, sp.microbeName.c_str()); s += b; }
    }
    return s;
}

inline std::string runtimeReport()
{
    Runtime &R = runtime();
    if (R.evaluations == 0) return std::string();
    char b[512];
    const double mean = R.sumF / (double) R.evaluations;
    const double pct  = 100.0 * (double) R.blocked / (double) R.evaluations;
    std::string s;
    std::sprintf(b, "  [THM] %ld gate evaluations, mean F_T %.4f, range %.4f..%.4f\n",
                 R.evaluations, mean, R.minF, R.maxF);
    s = b;
    std::sprintf(b, "  [THM] %ld of them (%.2f%%) were at or past the energy threshold, "
                    "so the rate was set to zero\n", R.blocked, pct);
    s += b;
    if (R.maxF <= 0.0)
        s += "  [THM] WARNING: the gate was closed everywhere, every time. Either dG0 has the\n"
             "  [THM] wrong sign, or the stoich line has reactants and products the wrong way\n"
             "  [THM] round, or the energy threshold is larger than the reaction can ever give.\n";
    else if (R.minF >= 1.0)
        s += "  [THM] note: the gate never closed anywhere. The run is identical to one with\n"
             "  [THM] <thermodynamics> switched off, which is worth knowing before reporting it\n"
             "  [THM] as a thermodynamically controlled result.\n";
    return s;
}

}  /* namespace complab_thermo */

#endif  /* COMPLAB3D_THERMO_HH */
