#include <cstdio>
#include <cmath>
#include "complab3d_symbolic.hh"
using namespace complab_sym;

static int fails = 0;
static void ck(const char *what, double got, double want, double tol = 1e-9) {
    bool ok = std::fabs(got - want) <= tol * (1.0 + std::fabs(want));
    if (!ok) ++fails;
    std::printf("  %-46s got %-14.8g want %-14.8g %s\n", what, got, want, ok ? "ok" : "** FAIL **");
}

static double one(const char *expr, const std::vector<std::string> &vars,
                  const std::vector<double> &vals, bool *parsed = 0) {
    std::vector<Node> nodes; int root = -1; std::string err;
    Parser p(vars, nodes);
    bool good = p.parse(expr, root, err);
    if (parsed) *parsed = good;
    if (!good) { std::printf("     parse error: %s\n", err.c_str()); return 0.0; }
    return evalNode(nodes, root, vals);
}

int main() {
    std::vector<std::string> V; V.push_back("a"); V.push_back("b");
    std::vector<double> X; X.push_back(2.0); X.push_back(3.0);

    std::printf("--- arithmetic and precedence\n");
    ck("1+2*3",            one("1+2*3", V, X), 7);
    ck("(1+2)*3",          one("(1+2)*3", V, X), 9);
    ck("2^3^2 right assoc",one("2^3^2", V, X), 512);
    ck("-a^2 is -(a^2)",   one("-a^2", V, X), -4);
    ck("a*b/(a+b)",        one("a*b/(a+b)", V, X), 6.0/5.0);
    ck("1e-3*2",           one("1e-3*2", V, X), 2e-3);
    ck("2^-1",             one("2^-1", V, X), 0.5);

    std::printf("--- functions and constants\n");
    ck("exp(0)",           one("exp(0)", V, X), 1);
    ck("log(e)",           one("log(e)", V, X), 1);
    ck("sqrt(a*b*6)",      one("sqrt(a*b*6)", V, X), 6);
    ck("min(a,b)",         one("min(a,b)", V, X), 2);
    ck("max(a,b)",         one("max(a,b)", V, X), 3);
    ck("pow(a,b)",         one("pow(a,b)", V, X), 8);
    ck("abs(a-b)",         one("abs(a-b)", V, X), 1);
    ck("pi",               one("pi", V, X), 3.14159265358979323846);

    std::printf("--- the failure modes a rate law must survive\n");
    ck("a/0 returns 0 not inf",  one("a/0", V, X), 0);
    ck("log(0) returns 0",       one("log(0)", V, X), 0);
    ck("sqrt(-1) returns 0",     one("sqrt(-1)", V, X), 0);
    ck("(-2)^0.5 returns 0",     one("(-2)^0.5", V, X), 0);

    std::printf("--- malformed input is rejected, not guessed at\n");
    const char *bad[] = {"1+", "(1+2", "a*", "nosuchvar", "min(a)", "1 2", "exp(", ")"};
    for (int i = 0; i < 8; ++i) {
        bool parsed = true; one(bad[i], V, X, &parsed);
        std::printf("  %-46s %s\n", bad[i], parsed ? "** FAIL: accepted **" : "rejected ok");
        if (parsed) ++fails;
    }

    std::printf("--- a real Monod law against the closed form\n");
    std::vector<std::string> M; M.push_back("S"); M.push_back("O");
    std::vector<double> mv; mv.push_back(1.0); mv.push_back(1.0);
    double got  = one("10*S/(0.05+S) * O/(0.05+O)", M, mv);
    double want = 10.0*1.0/(0.05+1.0) * 1.0/(0.05+1.0);
    ck("dual Monod", got, want);

    std::printf("\n%s\n", fails ? "SOME CHECKS FAILED" : "all expression checks passed");
    return fails;
}
