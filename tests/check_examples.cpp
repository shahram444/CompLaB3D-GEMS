/* Every .sym and .gnn file shipped in an example or a pipeline must load.
 * A file that no longer parses is a broken example, and this catches it
 * without needing Palabos or a simulation run. */
#include <cstdio>
#include "complab3d_symbolic.hh"
#include "complab3d_graphnet.hh"

static const char* SYMS[] = {
    "../examples/17_symbolic_law/input/growth.sym",
    "../examples/17_symbolic_law/input/growth_stoich.sym",
    "../pipelines/B_offline_models/B3_symbolic_law/expected/ecoli.sym",
    "../pipelines/B_offline_models/B3_symbolic_law/expected/ecoli_discovered.sym",
    "../pipelines/B_offline_models/B3_symbolic_law/expected/abiotic.sym",
    0
};
static const char* GNNS[] = {
    "../pipelines/B_offline_models/B4_graph_network/expected/aom.gnn",
    0
};

int main()
{
    int bad = 0;
    for (int i = 0; SYMS[i]; ++i) {
        complab_sym::Program P; std::string err;
        if (complab_sym::load(P, SYMS[i], &err))
            std::printf("  %-70s loaded ok\n", SYMS[i]);
        else { std::printf("  %-70s REJECTED: %s\n", SYMS[i], err.c_str()); ++bad; }
    }
    for (int i = 0; GNNS[i]; ++i) {
        complab_gnn::Network N; std::string err;
        if (complab_gnn::load(N, GNNS[i], &err))
            std::printf("  %-70s loaded ok\n", GNNS[i]);
        else { std::printf("  %-70s REJECTED: %s\n", GNNS[i], err.c_str()); ++bad; }
    }
    if (bad) { std::printf("\n%d shipped file(s) do not load\n", bad); return 1; }
    std::printf("\nevery shipped rate-law file loads\n");
    return 0;
}
