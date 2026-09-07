/* This file is a part of the CompLaB program.  AGPL-3.0-or-later.
 * Meile Lab, University of Georgia.  shahram.asgari@uga.edu
*/

/* ================================================================================================
 * complab3d_integration.hh  --  ONE PLACE WHERE THE NEW XML BLOCKS ARE READ AND ACTED ON
 * ================================================================================================
 *
 *  The five modules added in this round -- SBML reading, model fetching, geometry generation,
 *  surrogate training and CSV diagnostics -- are each self-contained and each testable on their
 *  own. This header is the seam between them and complab.cpp.
 *
 *  IT EXISTS SO THAT complab.cpp CHANGES AS LITTLE AS POSSIBLE. That file is 150 kB and cannot be
 *  compiled outside Palabos, so every line added to it is a line that cannot be tested until the
 *  first cluster build. Everything that CAN live outside it, does: all the XML parsing, all the
 *  validation, all the messages. complab.cpp gains five short calls.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE FIVE CALL SITES, in the order complab.cpp reaches them
 *
 *    1. after initialize_complab()          integ::readConfig(cfg)
 *    2. before the geometry is read         integ::provideGeometry(cfg, ...)
 *    3. while loading metabolic models      integ::resolveModel(cfg, name, ...)
 *    4. after the models are loaded         integ::prepareSurrogate(cfg, ...)
 *    5. at each save interval               integ::recordDiagnostics(cfg, ...)
 *
 *  Each is a no-op unless the matching XML block is present, so an existing input file behaves
 *  exactly as it did before.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT IS NOT TESTED
 *
 *  Everything in the five modules is tested standalone. This header is NOT: it needs Palabos to
 *  compile, so the first time it is built is on the cluster. It is written to be as close to
 *  trivial as possible for that reason -- parsing, string building, and calls into code that has
 *  been exercised. There is no arithmetic here worth getting wrong.
 * ================================================================================================
 */
#ifndef COMPLAB3D_INTEGRATION_HH
#define COMPLAB3D_INTEGRATION_HH

#include <string>
#include <vector>

#include "complab3d_sbml.hh"
#include "complab3d_fetch.hh"
#include "complab3d_geometry.hh"
#include "complab3d_surrogate.hh"
#include "complab3d_symbolic.hh"
#include "complab3d_graphnet.hh"
#include "complab3d_thermo.hh"
#include "complab3d_diagnostics.hh"

namespace integ {

/* ------------------------------------------------------------------------------------------------
 *  Everything the new blocks configure, in one struct.
 * ------------------------------------------------------------------------------------------------ */
struct Config {
    /* <model_source> and friends */
    std::string modelSource;          // "bigg:iJO1366", or empty
    std::string modelCache;           // default "input"
    std::string modelBundle;          // default "models"
    bool allowDownload;

    /* <domain><generate> / <import_raw> */
    bool generateGeometry;
    complab_geom::GenOptions gen;
    std::string importRaw;
    std::string rawDtype;
    double rawThreshold;
    bool rawInvert;
    std::string writeGeometryTo;      // where to put the generated .dat, for the record

    /* <surrogate> */
    bool srgEnabled;
    std::string srgWeights;
    bool srgTrainIfMissing;
    /* Which microbe this network is for. -1, the default, means every microbe whose
     * <reaction_type> is a surrogate type -- the common case of one organism, where making the
     * user write an index they cannot get wrong is just an opportunity to get it wrong. */
    int srgMicrobe;
    complab_srg::TrainOptions srgTrain;

    /* <symbolic> and <graphnet>.  Each may carry a biotic file, an abiotic file, or both:
     * the biotic one belongs to an organism and fires only where its biomass is; the abiotic one
     * belongs to nobody and fires in every fluid voxel.  They are independent switches on one
     * <enabled>, because a run that reads a .sym file at all wants the [SYM] report either way. */
    bool symEnabled;
    std::string symFile;              // <symbolic><expressions_file>
    std::string symAbioticFile;       // <symbolic><abiotic_file>
    int symMicrobe;                   // -1 = every microbe whose reaction_type is symbolic

    bool gnnEnabled;
    std::string gnnFile;              // <graphnet><network_file>
    std::string gnnAbioticFile;       // <graphnet><abiotic_file>
    int gnnMicrobe;

    /* <thermodynamics>.  Not a rate path: a factor applied to whichever rate path each organism
     * already uses.  One switch and one file name, for the same reason the two above have one --
     * the standard free energy of a reaction is a scientific claim and belongs in a file a reader
     * can open and cite, not in a rebuild.  See complab3d_thermo.hh. */
    bool thmEnabled;
    std::string thmFile;              // <thermodynamics><energetics_file>

    /* <diagnostics> */
    bool diagEnabled;
    std::string diagCsv;
    long diagInterval;
    double diagTol;

    /* <upscaling>.  Turns one resolved aggregate into the effectiveness factor and Thiele modulus
     * a continuum model needs.  Off unless asked for; see complab3d_upscale.hh. */
    bool upsEnabled;
    std::vector<plint> upsAggregateMat;   // which materials are the aggregate
    int  upsSpecies;                      // whose rate defines the reaction
    double upsRadius;                     // micrometres; 0 = derive from the measured volume
    double upsDiffusivity;                // m2/s; 0 = take the substrate's in-biofilm value
    /* Hold biomass fixed while the concentration profile relaxes.
     *
     * The effectiveness factor is a property of the CONCENTRATION profile inside the aggregate,
     * and that profile relaxes on R^2/D. Biomass moves on its own, much slower, timescale, and it
     * does not move uniformly -- the rim grows faster than the core because it is better fed. So
     * eta keeps drifting after the profile has settled, for a reason that has nothing to do with
     * transport, and a sweep can never report a steady number. Measured on this case: every run
     * still reported NOT STEADY at five diffusion times.
     *
     * With this on, the biomass increments are discarded each step while the substrate increments
     * are applied as usual. The reaction still consumes and still produces; only the catalyst is
     * held. That is exactly the question an effectiveness factor asks -- how much of THIS
     * aggregate is working -- and it is a diagnostic mode, not a way to run a simulation. */
    bool upsFreezeBiomass;
    std::string upsCsv;
    std::vector<std::string> diagConserve;

    Config()
      : modelCache("input"), modelBundle("models"), allowDownload(false),
        generateGeometry(false), rawDtype("uint8"), rawThreshold(128.0), rawInvert(false),
        srgEnabled(false), srgTrainIfMissing(false), srgMicrobe(-1),
        symEnabled(false), symMicrobe(-1),
        gnnEnabled(false), gnnMicrobe(-1),
        thmEnabled(false),
        diagEnabled(false), diagInterval(0), diagTol(1e-6),
        upsEnabled(false), upsSpecies(0), upsRadius(0), upsDiffusivity(0),
        upsFreezeBiomass(false), upsCsv("upscaling.csv") {}
};

/* ------------------------------------------------------------------------------------------------
 *  Parsing helpers.
 *
 *  Written against Palabos's XMLreader through a tiny shim so that this header does not need
 *  Palabos to be READ by a human, and so the pattern "absent means keep the default" is written
 *  once rather than forty times. Every optional tag behaves the same way: absent is fine, present
 *  but malformed stops the run with the tag named.
 * ------------------------------------------------------------------------------------------------ */
template <class Reader, class T>
inline bool getOpt(Reader &doc, const char *a, const char *b, const char *c, T &out)
{
    try {
        if (c) doc[a][b][c].read(out);
        else if (b) doc[a][b].read(out);
        else doc[a].read(out);
        return true;
    } catch (...) { return false; }
}

inline bool truthy(std::string s)
{
    for (size_t i = 0; i < s.size(); ++i)
        if (s[i] >= 'A' && s[i] <= 'Z') s[i] = (char) (s[i] - 'A' + 'a');
    return s == "true" || s == "yes" || s == "1" || s == "on";
}

/* Split "acetate 1e-3 10 log" into an input specification. */
inline bool parseInputSpec(const std::string &line, complab_srg::InputSpec &s, std::string &err)
{
    std::vector<std::string> tok;
    std::string cur;
    for (size_t i = 0; i <= line.size(); ++i) {
        if (i == line.size() || line[i] == ' ' || line[i] == '\t') {
            if (!cur.empty()) { tok.push_back(cur); cur.clear(); }
        } else cur += line[i];
    }
    if (tok.size() < 3) {
        err = "'" + line + "' should read: <substrate name> <low> <high> [log]";
        return false;
    }
    s.name = tok[0];
    s.lo = std::atof(tok[1].c_str());
    s.hi = std::atof(tok[2].c_str());
    s.logScale = (tok.size() > 3 && (tok[3] == "log" || tok[3] == "logarithmic"));
    if (s.hi <= s.lo) { err = "'" + line + "': the high end must exceed the low end"; return false; }
    if (s.logScale && s.lo <= 0.0) {
        err = "'" + line + "': a logarithmic sweep needs a strictly positive low end";
        return false;
    }
    return true;
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 1: read the new blocks. Returns "" on success, or a message to print before stopping.
 * ------------------------------------------------------------------------------------------------ */
template <class Reader>
inline std::string readConfig(Reader &doc, Config &cfg)
{
    const char *P = "parameters";
    std::string tmp;

    /* ---- model source ---- */
    if (getOpt(doc, P, "model_source", (const char *) 0, tmp)) cfg.modelSource = tmp;
    if (getOpt(doc, P, "model_cache", (const char *) 0, tmp)) cfg.modelCache = tmp;
    if (getOpt(doc, P, "model_bundle", (const char *) 0, tmp)) cfg.modelBundle = tmp;
    if (getOpt(doc, P, "allow_download", (const char *) 0, tmp)) cfg.allowDownload = truthy(tmp);

    /* ---- geometry ---- */
    if (getOpt(doc, "parameters", "LB_numerics", "domain", tmp)) { /* presence only */ }
    try {
        std::string kind;
        doc[P]["LB_numerics"]["domain"]["generate"].read(kind);
        cfg.generateGeometry = true;
        cfg.gen.kind = kind;
        double d;
        int i;
        if (getOpt(doc, P, "LB_numerics", "domain", tmp)) { }
        try { doc[P]["LB_numerics"]["domain"]["porosity"].read(d); cfg.gen.porosity = d; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["grain_radius"].read(d); cfg.gen.radius = d; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["grain_count"].read(i); cfg.gen.count = i; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["aperture"].read(d); cfg.gen.aperture = d; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["roughness"].read(d); cfg.gen.roughness = d; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["layers"].read(i); cfg.gen.layers = i; } catch (...) {}
        try { doc[P]["LB_numerics"]["domain"]["seed"].read(i); cfg.gen.seed = (unsigned long long) i; } catch (...) {}
    } catch (...) { cfg.generateGeometry = false; }

    /* [v1.3] <walls> used to be read INSIDE the <generate> try-block, after <generate> itself.
     * On the import path <generate> is absent, so the exception aborted the block before this
     * line and cfg.gen.walls kept its "y" default -- which importRaw() acts on, quietly turning
     * two whole planes of a segmented CT volume into inert wall, changing the porosity and the
     * pore network the run reports. <walls>none</walls> could not switch it off because it was
     * never read. It belongs to both paths, so it is read on its own. */
    try { doc[P]["LB_numerics"]["domain"]["walls"].read(tmp); cfg.gen.walls = tmp; } catch (...) {}

    try { doc[P]["LB_numerics"]["domain"]["import_raw"].read(cfg.importRaw); } catch (...) {}
    try { doc[P]["LB_numerics"]["domain"]["raw_dtype"].read(cfg.rawDtype); } catch (...) {}
    try { doc[P]["LB_numerics"]["domain"]["threshold"].read(cfg.rawThreshold); } catch (...) {}
    try { std::string s; doc[P]["LB_numerics"]["domain"]["invert"].read(s); cfg.rawInvert = truthy(s); } catch (...) {}
    try { doc[P]["LB_numerics"]["domain"]["write_geometry"].read(cfg.writeGeometryTo); } catch (...) {}

    if (cfg.generateGeometry && !cfg.importRaw.empty())
        return "  [GEOM] both <generate> and <import_raw> are set. Pick one.\n";

    /* ---- symbolic ---- */
    try {
        std::string en;
        doc[P]["symbolic"]["enabled"].read(en);
        cfg.symEnabled = truthy(en);
    } catch (...) { cfg.symEnabled = false; }

    if (cfg.symEnabled) {
        try { doc[P]["symbolic"]["expressions_file"].read(cfg.symFile); } catch (...) {}
        try { doc[P]["symbolic"]["abiotic_file"].read(cfg.symAbioticFile); } catch (...) {}
        try { int m; doc[P]["symbolic"]["microbe"].read(m); cfg.symMicrobe = m; } catch (...) {}
        if (cfg.symFile.empty() && cfg.symAbioticFile.empty())
            return "  [SYM] <symbolic> is enabled but neither <expressions_file> nor "
                   "<abiotic_file> is set.\n";
    }

    /* ---- graphnet ---- */
    try {
        std::string en;
        doc[P]["graphnet"]["enabled"].read(en);
        cfg.gnnEnabled = truthy(en);
    } catch (...) { cfg.gnnEnabled = false; }

    if (cfg.gnnEnabled) {
        try { doc[P]["graphnet"]["network_file"].read(cfg.gnnFile); } catch (...) {}
        try { doc[P]["graphnet"]["abiotic_file"].read(cfg.gnnAbioticFile); } catch (...) {}
        try { int m; doc[P]["graphnet"]["microbe"].read(m); cfg.gnnMicrobe = m; } catch (...) {}
        if (cfg.gnnFile.empty() && cfg.gnnAbioticFile.empty())
            return "  [GNN] <graphnet> is enabled but neither <network_file> nor "
                   "<abiotic_file> is set.\n";
    }

    /* ---- thermodynamics ---- */
    try {
        std::string en;
        doc[P]["thermodynamics"]["enabled"].read(en);
        cfg.thmEnabled = truthy(en);
    } catch (...) { cfg.thmEnabled = false; }

    if (cfg.thmEnabled) {
        try { doc[P]["thermodynamics"]["energetics_file"].read(cfg.thmFile); } catch (...) {}
        if (cfg.thmFile.empty())
            return "  [THM] <thermodynamics> is enabled but <energetics_file> is not set.\n";
    }

    /* ---- surrogate ---- */
    try {
        std::string en;
        doc[P]["surrogate"]["enabled"].read(en);
        cfg.srgEnabled = truthy(en);
    } catch (...) { cfg.srgEnabled = false; }

    if (cfg.srgEnabled) {
        try { doc[P]["surrogate"]["weights_file"].read(cfg.srgWeights); } catch (...) {}
        try { std::string s; doc[P]["surrogate"]["train_if_missing"].read(s); cfg.srgTrainIfMissing = truthy(s); } catch (...) {}
        try { int m; doc[P]["surrogate"]["microbe"].read(m); cfg.srgMicrobe = m; } catch (...) {}
        if (cfg.srgWeights.empty())
            return "  [SRG] <surrogate> is enabled but <weights_file> is not set.\n";

        int i;
        try { doc[P]["surrogate"]["train"]["samples"].read(i); cfg.srgTrain.samples = i; } catch (...) {}
        try { doc[P]["surrogate"]["train"]["restarts"].read(i); cfg.srgTrain.restarts = i; } catch (...) {}
        try { doc[P]["surrogate"]["train"]["epochs"].read(i); cfg.srgTrain.epochs = i; } catch (...) {}
        try { doc[P]["surrogate"]["train"]["verify_points"].read(i); cfg.srgTrain.verifyPoints = i; } catch (...) {}
        try { doc[P]["surrogate"]["train"]["seed"].read(i); cfg.srgTrain.seed = (unsigned long long) i; } catch (...) {}
        try { std::string s; doc[P]["surrogate"]["train"]["log_output"].read(s); cfg.srgTrain.logOutput = truthy(s); } catch (...) {}
        try {
            std::vector<int> L;
            doc[P]["surrogate"]["train"]["layers"].read(L);
            if (!L.empty()) cfg.srgTrain.hidden = L;
        } catch (...) {}

        cfg.srgTrain.inputs.clear();
        for (int k = 0; k < 8; ++k) {
            char tag[32];
            std::snprintf(tag, sizeof(tag), "input%d", k);
            std::string line;
            try { doc[P]["surrogate"]["train"][tag].read(line); } catch (...) { break; }
            complab_srg::InputSpec s;
            std::string err;
            if (!parseInputSpec(line, s, err))
                return "  [SRG] <" + std::string(tag) + ">: " + err + "\n";
            cfg.srgTrain.inputs.push_back(s);
        }
        if (cfg.srgTrainIfMissing && cfg.srgTrain.inputs.empty())
            return "  [SRG] <train_if_missing> is true but no <inputN> ranges are given, so there\n"
                   "  [SRG] is nothing to sweep. Add one line per input: <input0>name lo hi log</input0>\n";
        if (cfg.srgTrain.inputs.size() > 3)
            return "  [SRG] more than 3 training inputs. The sample count needed grows faster than\n"
                   "  [SRG] the fit improves; use 1 to 3, or fit offline with surrogate_training/.\n";
    }

    /* ---- diagnostics ---- */
    try {
        std::string en;
        doc[P]["diagnostics"]["enabled"].read(en);
        cfg.diagEnabled = truthy(en);
    } catch (...) { cfg.diagEnabled = false; }

    if (cfg.diagEnabled) {
        cfg.diagCsv = "summary.csv";
        try { doc[P]["diagnostics"]["summary_csv"].read(cfg.diagCsv); } catch (...) {}
        try { int i; doc[P]["diagnostics"]["interval"].read(i); cfg.diagInterval = i; } catch (...) {}
        try { doc[P]["diagnostics"]["tolerance"].read(cfg.diagTol); } catch (...) {}

        /* <upscaling>: pore scale in, one continuum number out. */
        try {
            std::string on;
            doc[P]["upscaling"]["enabled"].read(on);
            std::transform(on.begin(), on.end(), on.begin(),
                           [](unsigned char c){ return (char) std::tolower(c); });
            cfg.upsEnabled = (on == "true" || on == "yes" || on == "1");
        } catch (...) { cfg.upsEnabled = false; }
        if (cfg.upsEnabled) {
            try { doc[P]["upscaling"]["aggregate_materials"].read(cfg.upsAggregateMat); } catch (...) {}
            try { doc[P]["upscaling"]["species"].read(cfg.upsSpecies); } catch (...) {}
            try { doc[P]["upscaling"]["radius"].read(cfg.upsRadius); } catch (...) {}
            try { doc[P]["upscaling"]["diffusivity"].read(cfg.upsDiffusivity); } catch (...) {}
            try { doc[P]["upscaling"]["record_csv"].read(cfg.upsCsv); } catch (...) {}
            try {
                std::string fz;
                doc[P]["upscaling"]["freeze_biomass"].read(fz);
                std::transform(fz.begin(), fz.end(), fz.begin(),
                               [](unsigned char c){ return (char) std::tolower(c); });
                cfg.upsFreezeBiomass = (fz == "true" || fz == "yes" || fz == "1");
            } catch (...) { cfg.upsFreezeBiomass = false; }
            if (cfg.upsAggregateMat.empty())
                return "  <upscaling> is on but <aggregate_materials> is empty. Name the material\n"
                       "  number(s) the aggregate is made of -- normally the biofilm materials.\n"
                       "  Terminating.\n";
        }
        for (int k = 0; k < 32; ++k) {
            std::string e;
            if (k == 0) { try { doc[P]["diagnostics"]["conserve"].read(e); } catch (...) { break; } }
            else {
                char tag[32];
                std::snprintf(tag, sizeof(tag), "conserve%d", k);
                try { doc[P]["diagnostics"][tag].read(e); } catch (...) { break; }
            }
            if (!e.empty()) cfg.diagConserve.push_back(e);
        }
    }
    return "";
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 2: geometry.
 *
 *  Returns "" if the run should read <filename> as before. Otherwise generates or imports, writes
 *  a .dat so the run is reproducible from its own output, and returns the path to read.
 * ------------------------------------------------------------------------------------------------ */
inline std::string provideGeometry(const Config &cfg, int nx, int ny, int nz,
                                   const std::string &inputDir, std::string &log, bool &fatal)
{
    fatal = false;
    if (!cfg.generateGeometry && cfg.importRaw.empty()) return "";

    complab_geom::Grid g;
    std::string err;

    if (cfg.generateGeometry) {
        g = complab_geom::generate(nx, ny, nz, cfg.gen, &err);
        if (g.size() == 0) { log += "  [GEOM] " + err + "\n"; fatal = true; return ""; }
        log += "  [GEOM] generated a '" + cfg.gen.kind + "' pore space\n";
    } else {
        const std::string src = inputDir + cfg.importRaw;
        if (!complab_geom::importRaw(g, src, nx, ny, nz, cfg.rawDtype, cfg.rawThreshold,
                                     cfg.rawInvert, cfg.gen.walls, &err)) {
            log += "  [GEOM] " + err + "\n";
            fatal = true;
            return "";
        }
        log += "  [GEOM] imported " + src + "\n";
    }

    const complab_geom::Inspection R = complab_geom::inspect(g);
    log += complab_geom::inspectionText(g, R);
    if (!R.percolates) fatal = true;      // a sealed domain is not a well-posed flow problem

    const std::string out = inputDir + (cfg.writeGeometryTo.empty()
                                        ? std::string("generated_geometry.dat")
                                        : cfg.writeGeometryTo);
    if (complab_geom::writeDat(g, out, &err))
        log += "  [GEOM] written to " + out + " so this run can be reproduced exactly\n";
    else
        log += "  [GEOM] note: " + err + "\n";
    return out;
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 3: find a model file, and say plainly whether it is the one expected.
 * ------------------------------------------------------------------------------------------------ */
inline std::string resolveModel(const Config &cfg, const std::string &fallbackName,
                                bool isMaster, std::string &log, bool &fatal)
{
    fatal = false;
    if (cfg.modelSource.empty())
        return cfg.modelCache + "/" + fallbackName + ".xml";

    complab_fetch::Result r = complab_fetch::ensureModel(cfg.modelSource, cfg.modelCache,
                                                         cfg.allowDownload, isMaster,
                                                         cfg.modelBundle);
    log += r.message;
    if (!r.ok) fatal = true;
    return r.path;
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 4: load the surrogate, training it first if it is missing and training is allowed.
 * ------------------------------------------------------------------------------------------------ */
inline bool prepareSurrogate(const Config &cfg, complab_srg::Network &net,
                             complab_srg::FbaFn fba, void *ctx,
                             bool isMaster, std::string &log)
{
    if (!cfg.srgEnabled) return true;

    std::string err;
    if (complab_srg::load(net, cfg.srgWeights, &err)) {
        log += "  [SRG] loaded " + cfg.srgWeights + "\n";

        /* [FIX] Does this file belong to THIS model?
         *
         * Nothing used to check. Change <model_source> from bigg:e_coli_core to
         * bigg:iJO1366, forget to delete the old .srg, and the run quietly keeps using
         * the network fitted to the first organism -- every growth rate wrong, every
         * number plausible, and no line in the log to suggest it.
         *
         * The file records what it was trained on. Compare, and say so if they differ.
         * A warning rather than a stop: a network trained offline against the same model
         * under a different name is legitimate, and refusing to run would be worse than
         * saying what was noticed. */
        if (!cfg.modelSource.empty() && !net.provenance.empty()
            && net.provenance.find(cfg.modelSource) == std::string::npos) {
            log += "  [SRG] WARNING: this weights file says it was '" + net.provenance + "',\n"
                   "  [SRG] but <model_source> is '" + cfg.modelSource + "'. If the metabolic\n"
                   "  [SRG] model changed, this network was fitted to the OTHER one and every\n"
                   "  [SRG] growth rate it returns is wrong. Delete " + cfg.srgWeights + "\n"
                   "  [SRG] to refit, or ignore this if you know the two are the same model.\n";
        }
        char b[256];
        for (int i = 0; i < net.nIn; ++i) {
            std::snprintf(b, sizeof(b), "  [SRG]   input %d valid over %.6g .. %.6g\n",
                          i, net.trainMin[(size_t) i], net.trainMax[(size_t) i]);
            log += b;
        }
        log += "  [SRG] outside that box the input is held at the edge and counted; the end of\n"
               "  [SRG] the run says how often that happened.\n";

        /* HOW MANY NUMBERS DOES IT RETURN, AND WHAT DOES THE SOLVER DO WITH THE REST?
         *
         * A network that predicts growth alone leaves the solver to work out substrate
         * consumption from a Monod term. That is exact only where the swept uptake bound was the
         * binding constraint, and everywhere else it draws too much; the organism also cannot
         * excrete anything, because there is no flux to excrete. Both are real modelling
         * consequences of a choice made in the trainer, hours earlier, in another program.
         *
         * Nothing used to say so. The run started, the growth rates were plausible, and the
         * substrate field was wrong in a way that looks like a boundary condition. One line at
         * start-up is the whole fix. */
        if (net.nOut() <= 1) {
            log += "  [SRG] this network returns GROWTH ONLY. Substrate consumption will be\n"
                   "  [SRG] computed from a Monod term instead of from the metabolic model, which\n"
                   "  [SRG] is exact only where the swept uptake bound was the binding constraint;\n"
                   "  [SRG] elsewhere it draws too much. The organism cannot excrete a product at\n"
                   "  [SRG] all on this network. Retrain with the flux columns to remove both\n"
                   "  [SRG] limits -- the sweep already records them.\n";
        } else {
            /* [v1.3] This used to say "consumption comes from the model rather than from a Monod
             * term", which is not true on THIS path and was the most misleading line in the
             * start-up log.
             *
             * The run-time registry evaluates a .srg through evalBound(), and evalBound() returns
             * one number: N.eval(x), the growth output. Network::evalAll() -- which does return
             * all nOut() values -- is called from nowhere in the solver, only from the tests. So a
             * multi-output network loaded through <weights_file> has its exchange fluxes read,
             * validated, range-checked and then discarded, and Fout keeps the Monod estimate the
             * solver pre-filled. The rates looked right; the substrate field was wrong by the same
             * margin as a growth-only network, and the log said the opposite.
             *
             * Saying so is the fix that belongs in a correctness pass. Consuming the extra outputs
             * needs an output-to-substrate binding that this file format does not yet carry: the
             * outputnames are exchange reaction ids (EX_ac_e), not substrate names, and nothing
             * maps one to the other. The compiled path in surrogateModel.hh can already use them,
             * because there a human writes that mapping out by hand. */
            char nb[512];
            std::snprintf(nb, sizeof(nb),
                          "  [SRG] this network returns growth and %d exchange flux(es), but the\n"
                          "  [SRG] run-time <weights_file> path USES ONLY THE GROWTH OUTPUT. The\n"
                          "  [SRG] flux outputs are read and checked and then dropped, so substrate\n"
                          "  [SRG] consumption still comes from a Monod term and this organism still\n"
                          "  [SRG] cannot excrete a product -- exactly as for a growth-only network.\n"
                          "  [SRG] To use the fluxes, paste the network into surrogateModel.hh with\n"
                          "  [SRG] tools/surrogate/exportSurrogateHeader (the compiled path calls\n"
                          "  [SRG] evalAll and writes Fout), and rebuild.\n",
                          net.nOut() - 1);
            log += nb;
        }
        return true;
    }

    if (!cfg.srgTrainIfMissing) {
        log += "  [SRG] " + err + "\n"
               "  [SRG] Set <train_if_missing>true</train_if_missing> to fit one now, or produce\n"
               "  [SRG] the file offline with surrogate_training/.\n";
        return false;
    }
    if (!fba) {
        log += "  [SRG] training needs a flux balance solver, and this executable was built\n"
               "  [SRG] without one. Rebuild with -DENABLE_GLPK=ON, or supply the weights file.\n";
        return false;
    }

    if (!isMaster) {
        /* Ranks other than 0 must not train: they would each fit a different network. The caller
         * barriers and then loads the file rank 0 wrote. */
        return true;
    }

    log += "  [SRG] no weights file; training one now. This happens once.\n";
    complab_srg::TrainReport rep;
    if (!complab_srg::train(net, cfg.srgTrain, fba, ctx, rep, &err)) {
        log += "  [SRG] training failed: " + err + "\n";
        return false;
    }
    net.provenance = "trained in-run from " + cfg.modelSource;
    log += complab_srg::reportText(net, rep);
    if (complab_srg::save(net, cfg.srgWeights, &err))
        log += "  [SRG] written to " + cfg.srgWeights + "; later runs will load it directly.\n";
    else
        log += "  [SRG] note: could not write " + cfg.srgWeights + " (" + err + ")\n";
    return true;
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 5: diagnostics. The caller has already reduced the totals across ranks.
 * ------------------------------------------------------------------------------------------------ */
inline void setupDiagnostics(const Config &cfg, complab_diag::Diagnostics &D, bool isMaster,
                             const std::string &outputDir = std::string())
{
    if (!cfg.diagEnabled) return;

    /* Put the CSV where the rest of the run's output goes. Without this it lands in whatever
     * directory the job was launched from, which on a cluster is neither where the user looks
     * nor, for an array job, unique. An absolute path in <summary_csv> is honoured as given. */
    std::string path = cfg.diagCsv;
    if (!outputDir.empty() && !path.empty() && path[0] != '/')
        path = outputDir + path;

    D.configure(true, path, cfg.diagTol, isMaster);
    for (size_t i = 0; i < cfg.diagConserve.size(); ++i) D.addConserve(cfg.diagConserve[i]);
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 4b: the two learned rate paths.
 *
 *  Loads whichever of the four files the configuration names, binds each to this run's substrate
 *  and microbe names, and registers the bindings in the runtime the processors read.  Returns false
 *  after appending to `log` if anything is wrong; the caller prints and stops.
 *
 *  The storage lives in the CALLER, because the registry keeps pointers into it and those pointers
 *  are followed at every voxel of every step.  Passing them in rather than owning them here is what
 *  makes that lifetime obvious at the call site.
 *
 *  Binding is per microbe rather than "register for all", so that each organism's own biomass
 *  variable resolves to its own name.  Two organisms may share one .sym file and each will still
 *  see its own biomass.
 * ---------------------------------------------------------------------------------------------- */
inline bool prepareLearned(const Config &cfg,
                           const std::vector<std::string> &subsNames,
                           const std::vector<std::string> &microbeNames,
                           int nMicrobes,
                           const std::vector<int> &symUsers,
                           const std::vector<int> &gnnUsers,
                           complab_sym::Program &symProg,
                           complab_sym::Program &symAbioticProg,
                           complab_gnn::Network &gnnNet,
                           complab_gnn::Network &gnnAbioticNet,
                           std::string &log)
{
    /* ---------------------------------------------------------------- symbolic, biotic */
    if (cfg.symEnabled && !cfg.symFile.empty()) {
        std::string err;
        if (!complab_sym::load(symProg, cfg.symFile, &err)) {
            log += "  [SYM] cannot read " + cfg.symFile + ": " + err + "\n";
            return false;
        }
        log += complab_sym::describe(symProg, cfg.symFile);

        if (symUsers.empty())
            log += "  [SYM] note: a file was loaded but no microbe has reaction_type symbolic, "
                   "so nothing will evaluate it.\n";

        for (size_t k = 0; k < symUsers.size(); ++k) {
            const int gM = symUsers[k];
            const std::string nm = (gM >= 0 && gM < (int) microbeNames.size())
                                 ? microbeNames[(size_t) gM] : std::string();
            complab_sym::Binding b;
            const std::string berr = complab_sym::bindToSubstrates(symProg, subsNames, nm, b, false);
            if (!berr.empty()) {
                log += "  [SYM] " + cfg.symFile + ", for microbe '" + nm + "': " + berr + "\n";
                return false;
            }
            complab_sym::registerProgram(gM, b, nMicrobes);
        }
    }

    /* ---------------------------------------------------------------- symbolic, abiotic */
    if (cfg.symEnabled && !cfg.symAbioticFile.empty()) {
        std::string err;
        if (!complab_sym::load(symAbioticProg, cfg.symAbioticFile, &err)) {
            log += "  [SYM] cannot read " + cfg.symAbioticFile + ": " + err + "\n";
            return false;
        }
        log += complab_sym::describe(symAbioticProg, cfg.symAbioticFile);
        complab_sym::Binding b;
        const std::string berr = complab_sym::bindToSubstrates(symAbioticProg, subsNames,
                                                               std::string(), b, true);
        if (!berr.empty()) {
            log += "  [SYM] " + cfg.symAbioticFile + " (abiotic): " + berr + "\n";
            return false;
        }
        complab_sym::registerAbiotic(b);
        log += "  [SYM] the abiotic law will be swept over every fluid voxel.\n";
    }

    /* ---------------------------------------------------------------- graphnet, biotic */
    if (cfg.gnnEnabled && !cfg.gnnFile.empty()) {
        std::string err;
        if (!complab_gnn::load(gnnNet, cfg.gnnFile, &err)) {
            log += "  [GNN] cannot read " + cfg.gnnFile + ": " + err + "\n";
            return false;
        }
        log += complab_gnn::describe(gnnNet, cfg.gnnFile);

        if (gnnUsers.empty())
            log += "  [GNN] note: a network was loaded but no microbe has reaction_type graphnet, "
                   "so nothing will evaluate it.\n";

        complab_gnn::Binding b;
        const std::string berr = complab_gnn::bindToSubstrates(gnnNet, subsNames, b, false);
        if (!berr.empty()) {
            log += "  [GNN] " + cfg.gnnFile + ": " + berr + "\n";
            return false;
        }
        for (size_t k = 0; k < gnnUsers.size(); ++k)
            complab_gnn::registerNetwork(gnnUsers[k], b, nMicrobes);
    }

    /* ---------------------------------------------------------------- graphnet, abiotic */
    if (cfg.gnnEnabled && !cfg.gnnAbioticFile.empty()) {
        std::string err;
        if (!complab_gnn::load(gnnAbioticNet, cfg.gnnAbioticFile, &err)) {
            log += "  [GNN] cannot read " + cfg.gnnAbioticFile + ": " + err + "\n";
            return false;
        }
        log += complab_gnn::describe(gnnAbioticNet, cfg.gnnAbioticFile);
        complab_gnn::Binding b;
        const std::string berr = complab_gnn::bindToSubstrates(gnnAbioticNet, subsNames, b, true);
        if (!berr.empty()) {
            log += "  [GNN] " + cfg.gnnAbioticFile + " (abiotic): " + berr + "\n";
            return false;
        }
        complab_gnn::registerAbiotic(b);
        log += "  [GNN] the abiotic network will be swept over every fluid voxel.\n";
    }

    return true;
}

/* ------------------------------------------------------------------------------------------------
 *  CALL SITE 4c: the thermodynamic gate.
 *
 *  Loads the .thm file, binds its species and organism names to this run, and registers it in the
 *  runtime the rate processors read.  Returns false after appending to `log` if anything is wrong.
 *
 *  Unlike the two above, the model is owned by the runtime rather than by the caller: it holds no
 *  pointers into anything, so there is no lifetime for the call site to keep track of.
 *
 *  `kineticsInUse` says whether any organism gets its rates from defineKinetics.hh.  That path
 *  produces one combined rate vector per voxel, so a file with two reaction blocks would have to
 *  pick one of them arbitrarily.  Refusing is better than picking.
 * ---------------------------------------------------------------------------------------------- */
inline bool prepareThermo(const Config &cfg,
                          const std::vector<std::string> &subsNames,
                          const std::vector<std::string> &microbeNames,
                          int nMicrobes,
                          bool kineticsInUse,
                          std::string &log)
{
    if (!cfg.thmEnabled || cfg.thmFile.empty()) return true;

    complab_thermo::Model M;
    std::string err;
    if (!complab_thermo::load(M, cfg.thmFile, &err)) {
        log += "  [THM] cannot read " + cfg.thmFile + ": " + err + "\n";
        return false;
    }

    const std::string berr = complab_thermo::bindToRun(M, subsNames, microbeNames);
    if (!berr.empty()) {
        log += "  [THM] " + cfg.thmFile + ": " + berr + "\n";
        return false;
    }

    /* [v1.3] The kinetics path asks for gateFor(0) by hard-coded index, because it computes one
     * combined rate vector for the whole voxel. A one-block .thm naming any OTHER organism passed
     * the check below, registered a specOfMicrobe with 0 unset, and gateFor(0) then returned 1.0
     * before incrementing the evaluation counter -- so the gate did nothing, runtimeReport()
     * returned the empty string, and the only trace was the start-up echo of a block that was not
     * being used. The run was bit-identical to <thermodynamics> being off. */
    if (kineticsInUse && M.specs.size() == 1 &&
        M.specs[0].microbe >= 0 && M.specs[0].microbe != 0) {
        log += "  [THM] " + cfg.thmFile + " declares its single reaction block for an organism\n"
               "  [THM] other than microbe0, but an organism in this run takes its rates from\n"
               "  [THM] defineKinetics.hh, which produces one combined rate vector for the whole\n"
               "  [THM] voxel and can only be gated as microbe0. As written the gate would never\n"
               "  [THM] apply and nothing would say so. Write the block as `microbe all`, or name\n"
               "  [THM] microbe0, or move the organisms onto a per-organism rate path.\n";
        return false;
    }

    if (kineticsInUse && M.specs.size() > 1) {
        log += "  [THM] " + cfg.thmFile + " declares more than one reaction block, but an organism\n"
               "  [THM] in this run takes its rates from defineKinetics.hh. That path computes one\n"
               "  [THM] combined rate vector for the whole voxel, so a second gate would have\n"
               "  [THM] nowhere to apply. Use one block, or move the organisms onto a per-organism\n"
               "  [THM] rate path (glpk, cobrapy, surrogate, symbolic, graphnet).\n";
        return false;
    }

    log += complab_thermo::describe(M, cfg.thmFile);
    complab_thermo::runtime().model = M;
    complab_thermo::registerModel(nMicrobes);
    log += "  [THM] the gate multiplies whichever rate path each organism already uses.\n";
    return true;
}

}  // namespace integ

#endif  // COMPLAB3D_INTEGRATION_HH
