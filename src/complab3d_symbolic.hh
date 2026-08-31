/* ================================================================================================
 *  complab3d_symbolic.hh  --  RATE LAWS READ FROM A FILE AT RUN TIME
 *
 *  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
 *  Meile Lab, University of Georgia.  shahram.asgari@uga.edu
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT THIS IS FOR
 *
 *  defineKinetics.hh holds rate laws as C++, so changing one means editing a header and rebuilding
 *  the solver.  That is the right home for a rate law you wrote and intend to keep.  It is the wrong
 *  home for a rate law that a fitting procedure produced, because those arrive in batches, get
 *  revised, and have to be compared against one another.
 *
 *  This header reads rate expressions from a plain text file instead.  Symbolic regression, whether
 *  by PySR, SINDy, or a person with a pencil, produces a formula; you paste the formula into a .sym
 *  file and the next run uses it.  No rebuild, and the file is a scientific result that can be read,
 *  diffed and quoted in a paper without a tool -- the same reasoning that put the surrogate weights
 *  in a .srg file rather than in a compiled header.
 *
 *  This header does NOT do the fitting.  Genetic programming is slow and its answer depends on the
 *  seed, neither of which belongs inside a simulation.  Fit offline, write the file, run.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE FILE FORMAT
 *
 *      # anything after a hash is a comment
 *      version 1
 *      units   per_second                        or per_hour; see below
 *      vars    glucose o2 acetate Ecoli          names the solver will bind, in any order
 *
 *      rate    growth  = 0.412 * glucose * o2 / (0.083 + o2)
 *      rate    acetate = -0.31 * growth
 *
 *      range   glucose 0.1 10.0                   the box the fit was made over
 *      range   o2      0.1 20.0
 *
 *  A "rate" line names an output and gives an expression.  The output name is either a substrate
 *  name, in which case the value is that substrate's rate of change, or the reserved word "growth",
 *  which is the specific growth rate handed back as bioR.
 *
 *  EXPRESSIONS MAY REFER TO EARLIER OUTPUTS.  "growth" is defined above and then used in the acetate
 *  line, so a stoichiometric ratio can be written as one.  Expressions are evaluated in file order
 *  and a forward reference is a hard error, not a zero.
 *
 *  UNITS.  defineKinetics.hh works in mol/L per SECOND, and a growth rate per second, so that is
 *  what this header assumes and what "units per_second" (the default) means.  A formula fitted to
 *  flux balance output almost always comes out per HOUR, because that is the convention metabolic
 *  models use.  Writing "units per_hour" divides every rate by 3600 on load.  Getting this wrong is
 *  a factor of 3600, which is large enough to look like a modelling result rather than a mistake,
 *  so the loaded value is printed in the log either way.
 *
 *  RANGE LINES ARE NOT ADVICE.  A formula fitted between 0.1 and 10 mM is no more valid at 50 mM
 *  than a neural network is.  Every evaluation is clamped to the box and the clamping is counted;
 *  the end of a run reports how often it happened, and says so plainly when it happened a lot.
 *  A variable with no range line is not clamped, which is the right default for something like a
 *  biomass density that has no meaningful upper bound.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE EXPRESSION LANGUAGE
 *
 *      numbers        1  1.5  1e-3  6.02E23
 *      variables      any name from the vars line, or an earlier rate output
 *      operators      + - * / ^        unary minus       parentheses
 *      functions      exp log log10 sqrt abs min max pow tanh
 *      constants      pi  e
 *
 *  Precedence is the usual one: ^ binds tighter than * and /, which bind tighter than + and -.
 *  ^ is right associative, so 2^3^2 is 2^(3^2).  Division by zero returns 0 rather than an
 *  infinity, because a rate law that divides by a concentration that has gone to zero is a
 *  situation the solver has to survive, not a reason to stop.
 *
 *  An expression is parsed ONCE at start-up into a small tree of nodes and evaluated per voxel by
 *  walking that tree.  Parsing per voxel would cost more than the linear program it replaces.
 * ================================================================================================ */

#ifndef COMPLAB3D_SYMBOLIC_HH
#define COMPLAB3D_SYMBOLIC_HH

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace complab_sym {

/* ================================================================================================
 *  The parsed expression
 *
 *  A tree of nodes held in one flat vector; children are indices into it rather than pointers, so a
 *  Program copies and moves without any ownership bookkeeping.
 * ================================================================================================ */
enum NodeKind {
    NK_CONST,      /* value                       */
    NK_VAR,        /* a: slot in the value vector */
    NK_NEG,        /* a                           */
    NK_ADD, NK_SUB, NK_MUL, NK_DIV, NK_POW,   /* a, b */
    NK_EXP, NK_LOG, NK_LOG10, NK_SQRT, NK_ABS, NK_TANH,  /* a */
    NK_MIN, NK_MAX  /* a, b */
};

struct Node {
    NodeKind kind;
    double   value;   /* NK_CONST */
    int      a, b;    /* child indices, -1 when unused */
    Node() : kind(NK_CONST), value(0.0), a(-1), b(-1) {}
};

/* ------------------------------------------------------------------------------------------------
 *  Recursive-descent parser.  Small enough to read in one sitting, which is the point: a rate law
 *  the solver evaluates should not depend on a parser nobody in the group can check.
 * ---------------------------------------------------------------------------------------------- */
class Parser {
public:
    /* names[i] is the variable bound to value slot i */
    Parser(const std::vector<std::string> &names, std::vector<Node> &out)
        : nm(names), n(out), p(0), ok(true) {}

    bool parse(const std::string &text, int &root, std::string &err) {
        s = text; p = 0; ok = true; msg.clear();
        skip();
        root = expr();
        if (ok) { skip(); if (p != s.size()) fail("unexpected text after the expression"); }
        if (!ok) { err = msg; return false; }
        return true;
    }

private:
    const std::vector<std::string> &nm;
    std::vector<Node> &n;
    std::string s, msg;
    size_t p;
    bool ok;

    void fail(const std::string &m) {
        if (ok) { char b[32]; std::sprintf(b, " at character %d", (int) p + 1); msg = m + b; ok = false; }
    }
    void skip() { while (p < s.size() && (s[p] == ' ' || s[p] == '\t')) ++p; }
    bool eat(char c) { skip(); if (p < s.size() && s[p] == c) { ++p; return true; } return false; }
    char peek() { skip(); return p < s.size() ? s[p] : '\0'; }

    int mk(NodeKind k, int a = -1, int b = -1, double v = 0.0) {
        Node nd; nd.kind = k; nd.a = a; nd.b = b; nd.value = v;
        n.push_back(nd); return (int) n.size() - 1;
    }

    /* expr := term (('+'|'-') term)* */
    int expr() {
        int lhs = term();
        for (;;) {
            if (!ok) return lhs;
            skip();
            if (p < s.size() && s[p] == '+')      { ++p; lhs = mk(NK_ADD, lhs, term()); }
            else if (p < s.size() && s[p] == '-') { ++p; lhs = mk(NK_SUB, lhs, term()); }
            else return lhs;
        }
    }
    /* term := unary (('*'|'/') unary)* */
    int term() {
        int lhs = unary();
        for (;;) {
            if (!ok) return lhs;
            skip();
            if (p < s.size() && s[p] == '*')      { ++p; lhs = mk(NK_MUL, lhs, unary()); }
            else if (p < s.size() && s[p] == '/') { ++p; lhs = mk(NK_DIV, lhs, unary()); }
            else return lhs;
        }
    }
    /* unary := '-' unary | '+' unary | power */
    int unary() {
        skip();
        if (p < s.size() && s[p] == '-') { ++p; return mk(NK_NEG, unary()); }
        if (p < s.size() && s[p] == '+') { ++p; return unary(); }
        return power();
    }
    /* power := atom ('^' unary)?      right associative, and the exponent may be negative */
    int power() {
        int base = atom();
        skip();
        if (p < s.size() && s[p] == '^') { ++p; return mk(NK_POW, base, unary()); }
        return base;
    }

    int atom() {
        skip();
        if (p >= s.size()) { fail("expression ended early"); return mk(NK_CONST); }

        if (s[p] == '(') {
            ++p; int e = expr();
            if (!eat(')')) fail("missing ')'");
            return e;
        }
        if (std::isdigit((unsigned char) s[p]) || s[p] == '.') {
            const char *b = s.c_str() + p; char *e = 0;
            double v = std::strtod(b, &e);
            if (e == b) { fail("could not read a number"); return mk(NK_CONST); }
            p += (size_t) (e - b);
            return mk(NK_CONST, -1, -1, v);
        }
        if (std::isalpha((unsigned char) s[p]) || s[p] == '_') {
            size_t st = p;
            while (p < s.size() && (std::isalnum((unsigned char) s[p]) || s[p] == '_')) ++p;
            std::string id = s.substr(st, p - st);

            if (peek() == '(') {                       /* a function call */
                ++p;
                int a = expr(), b = -1;
                bool two = (id == "min" || id == "max" || id == "pow");
                if (two) { if (!eat(',')) fail("'" + id + "' takes two arguments"); b = expr(); }
                if (!eat(')')) fail("missing ')' after '" + id + "'");
                if      (id == "exp")   return mk(NK_EXP,   a);
                else if (id == "log")   return mk(NK_LOG,   a);
                else if (id == "log10") return mk(NK_LOG10, a);
                else if (id == "sqrt")  return mk(NK_SQRT,  a);
                else if (id == "abs")   return mk(NK_ABS,   a);
                else if (id == "tanh")  return mk(NK_TANH,  a);
                else if (id == "min")   return mk(NK_MIN,   a, b);
                else if (id == "max")   return mk(NK_MAX,   a, b);
                else if (id == "pow")   return mk(NK_POW,   a, b);
                fail("unknown function '" + id + "'");
                return mk(NK_CONST);
            }
            if (id == "pi") return mk(NK_CONST, -1, -1, 3.14159265358979323846);
            if (id == "e")  return mk(NK_CONST, -1, -1, 2.71828182845904523536);
            for (size_t i = 0; i < nm.size(); ++i)
                if (nm[i] == id) return mk(NK_VAR, (int) i);
            fail("'" + id + "' is not one of the declared variables");
            return mk(NK_CONST);
        }
        fail("unexpected character");
        return mk(NK_CONST);
    }
};

/* ------------------------------------------------------------------------------------------------
 *  Evaluation.  No allocation, no recursion beyond the tree depth, no branching on strings.
 * ---------------------------------------------------------------------------------------------- */
inline double evalNode(const std::vector<Node> &n, int i, const std::vector<double> &v) {
    if (i < 0 || i >= (int) n.size()) return 0.0;
    const Node &d = n[(size_t) i];
    switch (d.kind) {
        case NK_CONST: return d.value;
        case NK_VAR:   return (d.a >= 0 && d.a < (int) v.size()) ? v[(size_t) d.a] : 0.0;
        case NK_NEG:   return -evalNode(n, d.a, v);
        case NK_ADD:   return evalNode(n, d.a, v) + evalNode(n, d.b, v);
        case NK_SUB:   return evalNode(n, d.a, v) - evalNode(n, d.b, v);
        case NK_MUL:   return evalNode(n, d.a, v) * evalNode(n, d.b, v);
        case NK_DIV: {
            const double q = evalNode(n, d.b, v);
            /* A rate law that divides by a concentration which has reached zero is a state the
             * solver must survive.  Returning zero keeps the field finite; an infinity would
             * propagate into the lattice and end the run several thousand steps later, somewhere
             * unrelated. */
            if (std::fabs(q) < 1e-300) return 0.0;
            return evalNode(n, d.a, v) / q;
        }
        case NK_POW: {
            const double b = evalNode(n, d.a, v), e = evalNode(n, d.b, v);
            if (b < 0.0 && e != std::floor(e)) return 0.0;   /* would be NaN */
            return std::pow(b, e);
        }
        case NK_EXP:   return std::exp(evalNode(n, d.a, v));
        case NK_LOG: { const double x = evalNode(n, d.a, v); return x > 0.0 ? std::log(x)   : 0.0; }
        case NK_LOG10:{ const double x = evalNode(n, d.a, v); return x > 0.0 ? std::log10(x): 0.0; }
        case NK_SQRT: { const double x = evalNode(n, d.a, v); return x > 0.0 ? std::sqrt(x) : 0.0; }
        case NK_ABS:   return std::fabs(evalNode(n, d.a, v));
        case NK_TANH:  return std::tanh(evalNode(n, d.a, v));
        case NK_MIN: { const double x = evalNode(n, d.a, v), y = evalNode(n, d.b, v); return x < y ? x : y; }
        case NK_MAX: { const double x = evalNode(n, d.a, v), y = evalNode(n, d.b, v); return x > y ? x : y; }
    }
    return 0.0;
}

/* ================================================================================================
 *  A whole .sym file
 * ================================================================================================ */
struct Rate {
    std::string name;    /* "growth", or a substrate name */
    std::string source;  /* the expression as written, kept for the log and for provenance */
    int         root;    /* index into Program::nodes */
    Rate() : root(-1) {}
};

struct Program {
    std::vector<std::string> vars;     /* declared variable names, in file order              */
    std::vector<Node>        nodes;    /* every expression's tree, sharing one arena          */
    std::vector<Rate>        rates;    /* evaluated in this order; later may use earlier      */

    std::vector<double> loRange, hiRange;   /* per variable; hi <= lo means "no range given"  */
    std::string provenance;                 /* the comment lines, kept verbatim               */
    double      unitScale;                  /* 1 for per_second, 1/3600 for per_hour          */
    std::string unitName;

    Program() : unitScale(1.0), unitName("per_second") {}

    bool valid() const { return !rates.empty(); }

    int indexOfVar(const std::string &s) const {
        for (size_t i = 0; i < vars.size(); ++i) if (vars[i] == s) return (int) i;
        return -1;
    }
    int indexOfRate(const std::string &s) const {
        for (size_t i = 0; i < rates.size(); ++i) if (rates[i].name == s) return (int) i;
        return -1;
    }
};

/* ------------------------------------------------------------------------------------------------
 *  Runtime state.  One Program per microbe, plus the clamp counters, held the way
 *  complab_srg::runtime() holds the surrogate bindings.
 * ---------------------------------------------------------------------------------------------- */
struct Binding {
    const Program *prog;
    std::vector<int> subsOfVar;   /* variable slot -> substrate index, or -1 for "growth"/biomass */
    std::vector<int> subsOfRate;  /* rate index    -> substrate index, or -1 for "growth"         */
    Binding() : prog(0) {}
};

struct Runtime {
    std::vector<Binding> byMicrobe;

    /* The abiotic law, if there is one.  It is held on its own and not in byMicrobe because it
     * belongs to no organism: it is swept over every fluid voxel, biomass or not, which is the
     * whole difference between abiotic and biotic chemistry.  prog == 0 means none was loaded. */
    Binding abiotic;

    long evaluations, clamped;
    Runtime() : evaluations(0), clamped(0) {}
};
inline Runtime &runtime() { static Runtime R; return R; }

/* ================================================================================================
 *  Reading the file
 * ================================================================================================ */
inline bool load(Program &P, const std::string &path, std::string *err = 0)
{
    std::FILE *f = std::fopen(path.c_str(), "r");
    if (!f) { if (err) *err = "cannot open '" + path + "'"; return false; }

    P = Program();
    char line[4096];
    int lineNo = 0;
    bool sawVars = false;

    while (std::fgets(line, sizeof line, f)) {
        ++lineNo;
        std::string L(line);
        /* [FIX] A '#' line is now a plain comment and is dropped, matching
         * complab3d_graphnet.hh.  It used to be swept into the provenance string, so a file
         * with a normal explanatory header printed that whole header back as one run-on line
         * in the log -- which buried the thing provenance is for.  Where a file came from now
         * goes on its own 'provenance' line, read below. */
        const size_t h = L.find('#');
        if (h != std::string::npos) L = L.substr(0, h);
        /* split off the leading keyword */
        size_t i = 0;
        while (i < L.size() && std::isspace((unsigned char) L[i])) ++i;
        if (i >= L.size()) continue;
        size_t j = i;
        while (j < L.size() && !std::isspace((unsigned char) L[j])) ++j;
        const std::string key = L.substr(i, j - i);
        std::string rest = L.substr(j);

        char no[32]; std::sprintf(no, "%d", lineNo);
        const std::string where = "'" + path + "' line " + no + ": ";

        if (key == "version") {
            continue;                                   /* accepted and ignored, as in .srg */
        }
        else if (key == "provenance") {
     /* Where this file came from: printed in the log so a result carries it. */
     std::string c = L.substr(j);
     size_t b = c.find_first_not_of(" \t");
     while (!c.empty() && (c[c.size()-1]=='\n' || c[c.size()-1]=='\r')) c.erase(c.size()-1);
     if (b != std::string::npos) {
         c = c.substr(b);
         if (!P.provenance.empty()) P.provenance += " | ";
         P.provenance += c;
     }
     continue;
 }
 if (key == "units") {
            char u[64] = {0};
            if (std::sscanf(rest.c_str(), "%63s", u) != 1) {
                std::fclose(f); if (err) *err = where + "'units' names nothing"; return false; }
            std::string un(u);
            for (size_t k = 0; k < un.size(); ++k) un[k] = (char) std::tolower((unsigned char) un[k]);
            if      (un == "per_second" || un == "per_sec" || un == "s"   ) { P.unitScale = 1.0;          P.unitName = "per_second"; }
            else if (un == "per_hour"   || un == "per_hr"  || un == "h"   ) { P.unitScale = 1.0 / 3600.0; P.unitName = "per_hour";   }
            else { std::fclose(f); if (err) *err = where + "'units' must be per_second or per_hour, not '" + un + "'"; return false; }
        }
        else if (key == "vars") {
            size_t k = 0;
            while (k < rest.size()) {
                while (k < rest.size() && std::isspace((unsigned char) rest[k])) ++k;
                size_t st = k;
                while (k < rest.size() && !std::isspace((unsigned char) rest[k])) ++k;
                if (k > st) P.vars.push_back(rest.substr(st, k - st));
            }
            if (P.vars.empty()) { std::fclose(f); if (err) *err = where + "'vars' names nothing"; return false; }
            P.loRange.assign(P.vars.size(), 0.0);
            P.hiRange.assign(P.vars.size(), 0.0);
            /* an unset range is marked by hi <= lo, which no real range satisfies */
            for (size_t v = 0; v < P.vars.size(); ++v) { P.loRange[v] = 0.0; P.hiRange[v] = -1.0; }
            sawVars = true;
        }
        else if (key == "range") {
            if (!sawVars) { std::fclose(f); if (err) *err = where + "'range' appears before 'vars'"; return false; }
            char nameBuf[256]; double lo = 0, hi = 0;
            if (std::sscanf(rest.c_str(), "%255s %lf %lf", nameBuf, &lo, &hi) != 3) {
                std::fclose(f); if (err) *err = where + "'range' needs a name, a low value and a high value"; return false; }
            const int vi = P.indexOfVar(nameBuf);
            if (vi < 0) { std::fclose(f); if (err) *err = where + "'" + nameBuf + "' is not in the vars line"; return false; }
            if (!(hi > lo)) { std::fclose(f); if (err) *err = where + "range high must exceed range low"; return false; }
            P.loRange[(size_t) vi] = lo; P.hiRange[(size_t) vi] = hi;
        }
        else if (key == "rate") {
            if (!sawVars) { std::fclose(f); if (err) *err = where + "'rate' appears before 'vars'"; return false; }
            const size_t eq = rest.find('=');
            if (eq == std::string::npos) { std::fclose(f); if (err) *err = where + "a rate line needs '='"; return false; }
            std::string nm = rest.substr(0, eq), ex = rest.substr(eq + 1);
            /* trim */
            while (!nm.empty() && std::isspace((unsigned char) nm[0])) nm.erase(0, 1);
            while (!nm.empty() && std::isspace((unsigned char) nm[nm.size() - 1])) nm.erase(nm.size() - 1);
            while (!ex.empty() && (ex[ex.size() - 1] == '\n' || ex[ex.size() - 1] == '\r')) ex.erase(ex.size() - 1);
            if (nm.empty()) { std::fclose(f); if (err) *err = where + "the rate has no name"; return false; }
            if (P.indexOfRate(nm) >= 0) { std::fclose(f); if (err) *err = where + "'" + nm + "' is defined twice"; return false; }

            /* An expression may use the declared variables AND any rate defined above it, so the
             * parser is given both, with the earlier rates appended after the variables. */
            std::vector<std::string> scope = P.vars;
            for (size_t r = 0; r < P.rates.size(); ++r) scope.push_back(P.rates[r].name);

            Rate R; R.name = nm; R.source = ex;
            Parser parser(scope, P.nodes);
            std::string perr;
            if (!parser.parse(ex, R.root, perr)) {
                std::fclose(f); if (err) *err = where + perr; return false; }
            P.rates.push_back(R);
        }
        else {
            std::fclose(f);
            if (err) *err = where + "unknown keyword '" + key + "'";
            return false;
        }
    }
    std::fclose(f);

    if (!sawVars)      { if (err) *err = "'" + path + "' has no 'vars' line"; return false; }
    if (P.rates.empty()){ if (err) *err = "'" + path + "' defines no rates"; return false; }
    return true;
}

/* ================================================================================================
 *  Evaluation
 *
 *  values[] is laid out as the parser expects: the declared variables first, then one slot per rate
 *  already computed, so an expression that refers to an earlier rate reads it from the same vector.
 * ================================================================================================ */
inline void evaluate(const Program &P, const std::vector<double> &varValues,
                     std::vector<double> &out, bool clamp = true, long *clampCount = 0)
{
    const size_t nv = P.vars.size();
    std::vector<double> v(nv + P.rates.size(), 0.0);

    for (size_t i = 0; i < nv; ++i) {
        double x = (i < varValues.size()) ? varValues[i] : 0.0;
        if (clamp && P.hiRange[i] > P.loRange[i]) {          /* a range was given for this one */
            if (x < P.loRange[i]) { x = P.loRange[i]; if (clampCount) ++*clampCount; }
            else if (x > P.hiRange[i]) { x = P.hiRange[i]; if (clampCount) ++*clampCount; }
        }
        v[i] = x;
    }
    out.assign(P.rates.size(), 0.0);
    for (size_t r = 0; r < P.rates.size(); ++r) {
        double y = evalNode(P.nodes, P.rates[r].root, v);
        if (y != y) y = 0.0;                                  /* NaN */
        if (y > 1e300) y = 1e300; else if (y < -1e300) y = -1e300;
        /* Chained expressions see the UNSCALED value, so a line like
         *     rate acetate = 0.31 * growth
         * stays a pure stoichiometric ratio whichever unit the file declares. The scale is
         * applied only on the way out. */
        v[nv + r] = y;
        out[r] = y * P.unitScale;
    }
}

/* ------------------------------------------------------------------------------------------------
 *  A one-line summary for the end of a run, in the same voice as complab_srg::runtimeReport().
 * ---------------------------------------------------------------------------------------------- */
inline std::string runtimeReport()
{
    Runtime &R = runtime();
    if (R.evaluations == 0) return std::string();
    char buf[512];
    const double pct = 100.0 * (double) R.clamped / (double) R.evaluations;
    std::sprintf(buf, "  [SYM] %ld expression evaluations, %ld clamped to the fitted range (%.2f%%)\n",
                 R.evaluations, R.clamped, pct);
    std::string s(buf);
    if (pct > 5.0)
        s += "  [SYM] WARNING: more than 5% of evaluations were outside the range these rate laws\n"
             "  [SYM] were fitted over. The results should not be used until the ranges are widened\n"
             "  [SYM] and the expressions refitted.\n";
    return s;
}

/* ------------------------------------------------------------------------------------------------
 *  What the log prints when a file is loaded.  A rate law is a scientific claim; the run should say
 *  out loud which one it is using.
 * ---------------------------------------------------------------------------------------------- */
inline std::string describe(const Program &P, const std::string &path)
{
    std::string s = "  [SYM] rate laws from '" + path + "'\n";
    if (!P.provenance.empty()) s += "  [SYM]   provenance: " + P.provenance + "\n";
    s += "  [SYM]   units: " + P.unitName;
    if (P.unitScale != 1.0) s += "  (every rate divided by 3600 on load)";
    s += "\n";
    for (size_t r = 0; r < P.rates.size(); ++r)
        s += "  [SYM]   " + P.rates[r].name + " = " + P.rates[r].source + "\n";
    bool any = false;
    for (size_t i = 0; i < P.vars.size(); ++i) if (P.hiRange[i] > P.loRange[i]) any = true;
    if (any) {
        s += "  [SYM]   valid over:";
        for (size_t i = 0; i < P.vars.size(); ++i) {
            if (P.hiRange[i] <= P.loRange[i]) continue;
            char b[128];
            std::sprintf(b, " %s %.6g..%.6g", P.vars[i].c_str(), P.loRange[i], P.hiRange[i]);
            s += b;
        }
        s += "\n";
    } else {
        s += "  [SYM]   no range lines: evaluations will not be clamped or counted.\n";
    }
    return s;
}

}  // namespace complab_sym

#endif  // COMPLAB3D_SYMBOLIC_HH
