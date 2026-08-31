/* ============================================================================
 * CompLaB3D - Three-Dimensional Biogeochemical Reactive Transport Solver
 * ============================================================================
 *
 * Author:      Shahram Asgari
 * Advisor:     Dr. Christof Meile
 * Laboratory:  Meile Lab
 * Institution: University of Georgia (UGA)
 *
 * ============================================================================
 * CALCULATION FLOW (10 PHASES):
 * ───────────────────────────────────────────────────────────────────────────
 * PHASE 1:  Load XML configuration and validate inputs
 * PHASE 2:  Geometry setup and preprocessing
 * PHASE 3:  Navier-Stokes flow field simulation
 *           └─ STEP 3.1: Initial pressure simulation → measure u₀
 *           └─ STEP 3.2: Calculate permeability: k = (u₀ × ν × L) / ΔP₀
 *           └─ STEP 3.3: Calculate target velocity: u_target = (Pe × D) / L
 *           └─ STEP 3.4: Corrected pressure: ΔP_new = (u_target × ν × L) / k
 *           └─ STEP 3.5: Second NS simulation → achieve target velocity
 *           └─ STEP 3.6: Stability checks (Ma, CFL, τ)
 * PHASE 4:  Reactive transport lattice setup (substrates + biomass)
 * PHASE 5:  NS-ADE velocity field coupling
 * PHASE 6:  Main simulation loop
 *           └─ STEP 6.1: Collision step (LBM)
 *           └─ STEP 6.2: Kinetics reactions (Monod, decay)
 *           └─ STEP 6.3: Equilibrium chemistry solver
 *           └─ STEP 6.4: Biomass expansion (CA/FD)
 *           └─ STEP 6.5: Flow field update (if biofilm changed)
 *           └─ STEP 6.6: Streaming step (LBM)
 * PHASE 7:  Output VTI/CHK files
 * PHASE 8:  Calculate moments and BTC analysis
 * PHASE 9:  Write summary files
 * PHASE 10: Finalize and cleanup
 * ───────────────────────────────────────────────────────────────────────────
 *
 * SIMULATION MODES:
 *   - biotic_mode: true/false (with/without microbes)
 *   - enable_kinetics: true/false (biotic kinetics reactions on/off)
 *   - enable_abiotic_kinetics: true/false (abiotic chemical reactions on/off)
 *   - enable_validation_diagnostics: true/false (detailed per-iteration output)
 *
 * OUTPUT FILES:
 *   - VTI: Concentration, biomass, velocity fields
 *   - CHK: Binary checkpoints for restart
 *   - CSV: BTC timeseries, domain properties, moments summary
 *
 * ============================================================================
 */

#include "complab_functions.hh"
#include "complab3d_processors.hh"
#include "../defineKinetics.hh"        // For KineticsStats namespace - in project root
#include "../defineAbioticKinetics.hh" // For abiotic kinetics (substrate-only reactions)
#include "complab3d_rates.hh"      // [RATE-OUT] optional reaction rate output, off unless CompLaB.xml asks
#include "precipitationVOP.hh"         // [PRECIP-VOP] surface precipitation + pore-clogging feedback
#include "dissolutionVOP.hh"           // [DISSOL-VOP] mineral dissolution + pore re-opening
#include "complab3d_processors_symbolic.hh"  // [SYM] rate laws read from a .sym file at run time
#include "complab3d_processors_graphnet.hh"  // [GNN] a graph network read from a .gnn file
#include "complab3d_processors_abiotic_learned.hh"  // [SYM][GNN] the same two, swept abiotically
#include <algorithm>
#include <cctype>

#include <chrono>
#include <string>
#include <iostream>
#include <cstring>
#include <vector>
#include <sys/stat.h>
#include <iomanip>
#include <cmath>

// ============================================================================
// STABILITY CHECK STRUCTURE  
// ============================================================================
struct StabilityReport {
    T Ma, CFL, tau_NS, tau_ADE, Pe_grid;
    bool Ma_ok, Ma_warning, CFL_ok, tau_NS_ok, tau_ADE_ok, Pe_grid_ok, all_ok, has_warnings;
};

StabilityReport performStabilityChecks(T u_max, T tau_NS, T tau_ADE, T D_lattice) {
    StabilityReport report;
    T cs = std::sqrt(1.0 / 3.0);
    report.Ma = u_max / cs;
    report.Ma_ok = (report.Ma < 1.0);
    report.Ma_warning = (report.Ma > 0.3);
    report.CFL = u_max;
    report.CFL_ok = (report.CFL < 1.0);
    report.tau_NS = tau_NS;
    report.tau_NS_ok = (tau_NS > 0.5 && tau_NS < 2.0);
    report.tau_ADE = tau_ADE;
    report.tau_ADE_ok = (tau_ADE > 0.5 && tau_ADE < 2.0);
    report.Pe_grid = (D_lattice > 1e-14) ? (u_max / D_lattice) : 0.0;
    report.Pe_grid_ok = (report.Pe_grid < 2.0);
    report.all_ok = report.Ma_ok && report.CFL_ok && report.tau_NS_ok && report.tau_ADE_ok;
    report.has_warnings = report.Ma_warning || !report.Pe_grid_ok;
    return report;
}

void printStabilityReport(const StabilityReport& report) {
    pcout << "\n╔════════════════════════════════════════════════════════════╗\n";
    pcout << "║              STABILITY CHECK REPORT                        ║\n";
    pcout << "╠════════════════════════════════════════════════════════════╣\n";
    pcout << "║ Ma = " << std::setprecision(4) << report.Ma << (report.Ma_ok ? " OK" : " FAIL");
    pcout << "   CFL = " << report.CFL << (report.CFL_ok ? " OK" : " FAIL") << "             ║\n";
    pcout << "║ tau_NS = " << report.tau_NS << (report.tau_NS_ok ? " OK" : " FAIL");
    pcout << "   tau_ADE = " << report.tau_ADE << (report.tau_ADE_ok ? " OK" : " FAIL") << "            ║\n";
    pcout << "║ Pe_grid = " << report.Pe_grid << (report.Pe_grid_ok ? " OK" : " WARN") << "                                       ║\n";
    pcout << "╚════════════════════════════════════════════════════════════╝\n\n";
}


int main(int argc, char **argv) {

    plbInit(&argc, &argv);
    /* [FLUSH-FIX 2026-07-24] On a SLURM cluster stdout->file is block-buffered, so pcout lines
     * that end in "\n" (Phase 3 setup, the per-iteration ITERATION block) sit unflushed for a long
     * time and the run LOOKS frozen even though it is progressing. unitbuf flushes after every write. */
    std::cout << std::unitbuf;
    global::timer("total").start();

    // ════════════════════════════════════════════════════════════════════════════
    // STARTUP BANNER
    // ════════════════════════════════════════════════════════════════════════════
    pcout << "\n";
    pcout << "╔══════════════════════════════════════════════════════════════════════════╗\n";
    pcout << "║                            CompLaB3D                                     ║\n";
    pcout << "║       Three-Dimensional Biogeochemical Reactive Transport Solver        ║\n";
    pcout << "║              Lattice Boltzmann Method (LBM) + Equilibrium                ║\n";
    pcout << "╠══════════════════════════════════════════════════════════════════════════╣\n";
    pcout << "║  Author:  Shahram Asgari                                                 ║\n";
    pcout << "║  Advisor: Dr. Christof Meile                                             ║\n";
    pcout << "║  Lab:     Meile Lab, University of Georgia                               ║\n";
    pcout << "╚══════════════════════════════════════════════════════════════════════════╝\n\n";

    ImageWriter<T> image("leeloo");

    // Diagnostic counters
    plint diag_ca_triggers = 0;
    plint diag_ca_redistributions = 0;
    T diag_initial_biomass = 0.0;

    // asserted variables
    plint kns_count=0, fd_count=0, lb_count=0, ca_count=0, bfilm_count=0, bfree_count=0;
    char *main_path = (char*)malloc(100 * sizeof(char));
    getcwd(main_path, 100 * sizeof(char));
    char *src_path = (char*)malloc(100 * sizeof(char));
    char *input_path = (char*)malloc(100 * sizeof(char));
    char *output_path = (char*)malloc(100 * sizeof(char));
    char *ns_filename = (char*)malloc(100 * sizeof(char));
    plint nx, ny, nz, num_of_microbes, num_of_substrates;
    T dx, dy, dz, deltaP, Pe, charcs_length;
    std::string geom_filename, mask_filename;
    std::vector<bool> vec_left_btype, vec_right_btype, bio_left_btype, bio_right_btype;
    std::vector<T> vec_c0, vec_b0_free, vec_left_bcondition, vec_right_bcondition, bio_left_bcondition, bio_right_bcondition, vec_permRatio;
    std::vector< std::vector<T> > vec_b0_all, vec_b0_film, vec_Kc_kns, vec_Vmax, vec_Vmax_kns;
    std::vector<T> vec_mu;
    std::vector< std::vector<T> > vec_Kc;

    // variables with default values
    std::string ade_filename, bio_filename;
    bool read_NS_file=0, read_ADE_file=0, soluteDindex=0, bmassDindex=0, track_performance=0., halfflag=0;
    plint no_dynamics=0, bounce_back=1, ns_rerun_iT0=0, ns_update_interval=1, ade_update_interval=1,
        ns_maxiTer_1, ns_maxiTer_2, ade_rerun_iT0=0, ade_maxiTer=10000000, ade_VTI_iTer=1000, ade_CHK_iTer=1000000;
    T tau=0.8, max_bMassRho=1., ns_converge_iT1=1e-8, ns_converge_iT2=1e-4, ade_converge_iT=1e-8, thrd_bFilmFrac;

    T DarcyOutletUx=0., permeability=0., u_target=0., deltaP_new=0., u_final=0., Pe_achieved=0.;
    T tau_ADE_fixed=0.8, D_lattice_fixed=0., tortuosity_factor=3.0, safety_factor=1.5;
    plint estimated_iterations=0;

    std::vector<bool> bmass_type;
    std::vector<plint> pore_dynamics, solver_type, reaction_type;
    std::vector<T> vec_solute_poreD, vec_solute_bFilmD, vec_bMass_poreD, vec_bMass_bFilmD, vec_mu_kns;
    std::vector<std::string> vec_subs_names, vec_microbes_names;
    std::vector< std::vector<plint> > bio_dynamics;

    // Equilibrium chemistry variables
    bool useEquilibrium = false;
    EquilibriumChemistry<T> eqSolver;
    T eqtime = 0.0;
    std::vector<std::string> eq_component_names;
    std::vector<T> eq_logK_values;
    std::vector<std::vector<T>> eq_stoich_matrix;

    // Biotic/Abiotic and Kinetics control
    bool biotic_mode = true;      // true = with microbes, false = abiotic transport only
    bool enable_kinetics = true;  // true = kinetics enabled, false = equilibrium only
    bool enable_abiotic_kinetics = false;  // true = abiotic reactions (no microbes)
    bool enable_validation_diagnostics = false;  // true = detailed per-iteration diagnostics

    std::string str_mainDir=main_path;
    if (std::to_string(str_mainDir.back()).compare("/")!=0) { str_mainDir+="/"; }
    std::srand(std::time(nullptr));

    // ════════════════════════════════════════════════════════════════════════════
    // PHASE 1: LOAD CONFIGURATION
    // ════════════════════════════════════════════════════════════════════════════
    pcout << "┌────────────────────────────────────────────────────────────────────────┐\n";
    pcout << "│ PHASE 1: LOADING CONFIGURATION                                        │\n";
    pcout << "└────────────────────────────────────────────────────────────────────────┘\n";
    
    int erck = 0;
    try {
        erck=initialize_complab( main_path, src_path, input_path, output_path, ns_filename, ade_filename, bio_filename, geom_filename, mask_filename,
        read_NS_file, ns_rerun_iT0, ns_converge_iT1, ns_converge_iT2, ns_maxiTer_1, ns_maxiTer_2, ns_update_interval, ade_update_interval,
        read_ADE_file, ade_rerun_iT0, ade_VTI_iTer, ade_CHK_iTer, ade_converge_iT, ade_maxiTer, nx, ny, nz, dx, dy, dz, deltaP, tau,
        Pe, charcs_length, vec_solute_poreD, vec_solute_bFilmD, vec_bMass_poreD, vec_bMass_bFilmD, soluteDindex, bmassDindex, thrd_bFilmFrac, vec_permRatio, max_bMassRho,
        pore_dynamics, bounce_back, no_dynamics, bio_dynamics, num_of_microbes, num_of_substrates, vec_subs_names, vec_microbes_names,
        solver_type, fd_count, lb_count, ca_count, bfilm_count, bfree_count, kns_count, reaction_type,
        vec_c0, vec_left_btype, vec_right_btype, vec_left_bcondition, vec_right_bcondition, vec_b0_all, bio_left_btype, bio_right_btype, bio_left_bcondition, bio_right_bcondition,
        vec_Kc, vec_Kc_kns, vec_mu, vec_mu_kns, bmass_type, vec_b0_free, vec_b0_film, vec_Vmax, vec_Vmax_kns, track_performance, halfflag,
        useEquilibrium, eq_component_names, eq_logK_values, eq_stoich_matrix,
        biotic_mode, enable_kinetics, enable_abiotic_kinetics,
        enable_validation_diagnostics);
    }
    catch (PlbIOException& exception) {
        pcout << "  [ERROR] " << exception.what() << "\n";
        return -1;
    }
    if (erck!=0) { return -1; }
    pcout << "  [OK] XML configuration loaded and validated\n";

    // ============================================================================
    // OPTIONAL METABOLIC LAYER -- flux balance analysis and surrogate models.
    //
    //   Reads only the metabolic settings, out of the same CompLaB.xml.  It is
    //   deliberately kept out of initialize_complab(), which already takes over
    //   sixty arguments by reference; keeping it separate means this whole
    //   feature can be removed again by deleting these few lines plus three
    //   #includes.  See complab3d_metabolic.hh for the units and XML reference.
    //
    //   If the user asked for nothing, initialize_metabolic() returns with every
    //   switch off and the rest of the program behaves exactly as before.
    // ============================================================================
    MetabolicConfig mmcfg;
    // Releases the GLPK problems and shuts the embedded interpreter down on EVERY
    // exit path from main(), including the early returns further down.
    MetabolicScopeGuard mmguard(mmcfg);
    {
        const int mmerr = initialize_metabolic(mmcfg, reaction_type, vec_subs_names,
                                               num_of_microbes, num_of_substrates,
                                               useEquilibrium, eq_stoich_matrix, eq_component_names,
                                               vec_Kc, vec_mu);
        if (mmerr != 0) return -1;
    }

    // ============================================================================
    // [IMMOBILE] Per-substrate <immobile>true</immobile> (default false).
    //   An immobile species never advects/diffuses/streams/collides -- it only
    //   receives its reaction source term (stable at any diffusivity, incl. D=0).
    // ============================================================================
    std::vector<bool> vec_immobile(num_of_substrates, false);
    try {
        XMLreader immdoc("CompLaB.xml");
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            std::string chemname = "substrate" + std::to_string(iS);
            try {
                std::string tmp;
                immdoc["parameters"]["chemistry"][chemname]["immobile"].read(tmp);
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp.compare("true")==0 || tmp.compare("1")==0 || tmp.compare("yes")==0) vec_immobile[iS]=true;
            } catch (PlbIOException&) { /* flag absent -> mobile (default) */ }
        }
    } catch (PlbIOException&) { /* no XML re-read -> all mobile */ }
    {
        plint nimm=0; for (plint iS=0; iS<num_of_substrates; ++iS) if (vec_immobile[iS]) ++nimm;
        pcout << "  [IMMOBILE] immobile (solid/precipitate) substrates: " << nimm << " of " << num_of_substrates << "\n";
        for (plint iS=0; iS<num_of_substrates; ++iS) if (vec_immobile[iS])
            pcout << "    - substrate" << iS << " (" << vec_subs_names[iS] << "): IMMOBILE (reaction source only)\n";
    }

    // ============================================================================
    // [PRECIP-VOP] Precipitation-induced pore clogging (Kang VOP + flow feedback).
    //   All fields optional; if <precipitation> is absent the solver is unchanged.
    // ============================================================================
    int   precip_enabled = 0;         plint precip_solidSub = -1;
    T     max_precipRho = 1e30;       int   precip_surfaceOnly = 1;
    T     precip_permRatio = 0.0;     plint precip_update_interval = 200;
    try {
        XMLreader pdoc("CompLaB.xml");
        try { std::string s; pdoc["parameters"]["precipitation"]["enabled"].read(s);
              std::transform(s.begin(),s.end(),s.begin(),[](unsigned char c){ return std::tolower(c); });
              if (s=="true"||s=="1"||s=="yes") precip_enabled = 1; } catch (PlbIOException&) {}
        if (precip_enabled) {
            try { pdoc["parameters"]["precipitation"]["solid_substrate"].read(precip_solidSub); } catch (PlbIOException&) {}
            try { pdoc["parameters"]["precipitation"]["max_precipRho"].read(max_precipRho); } catch (PlbIOException&) {}
            try { pdoc["parameters"]["precipitation"]["surface_only"].read(precip_surfaceOnly); } catch (PlbIOException&) {}
            try { pdoc["parameters"]["precipitation"]["perm_ratio"].read(precip_permRatio); } catch (PlbIOException&) {}
            try { pdoc["parameters"]["precipitation"]["update_interval"].read(precip_update_interval); } catch (PlbIOException&) {}
        }
    } catch (PlbIOException&) {}
    if (precip_enabled) {
        pcout << "  [PRECIP-VOP] enabled: solid substrate=" << precip_solidSub
              << " max_precipRho=" << max_precipRho
              << " surface_only=" << precip_surfaceOnly
              << " perm_ratio=" << precip_permRatio
              << " update_interval=" << precip_update_interval << "\n";
    }

    // ============================================================================
    // [DISSOL-VOP] Mineral dissolution and pore re-opening.  Fully optional and
    //   independent of <precipitation>; if the block is absent nothing changes.
    //
    //   A "phase" is a declared dissolvable solid: which material number it
    //   occupies, which immobile substrate holds its inventory, how much a full
    //   voxel holds, and how full it starts.  Precipitated mineral and original
    //   rock are the same mechanism, differing only in how the inventory is
    //   seeded.  Anything not declared here can never dissolve.
    // ============================================================================
    DissolutionConfig dissolCfg;
    plint precip_phase_id = 0;      // phase id stamped on a voxel that seals
    try {
        XMLreader ddoc("CompLaB.xml");
        try { std::string s2; ddoc["parameters"]["dissolution"]["enabled"].read(s2);
              std::transform(s2.begin(),s2.end(),s2.begin(),[](unsigned char c){ return std::tolower(c); });
              dissolCfg.enabled = (s2=="true"||s2=="1"||s2=="yes"||s2=="on"); } catch (PlbIOException&) {}
        if (dissolCfg.enabled) {
            try { ddoc["parameters"]["dissolution"]["reopen_fraction"].read(dissolCfg.reopen_fraction); } catch (PlbIOException&) {}
            try { ddoc["parameters"]["dissolution"]["surface_only"].read(dissolCfg.surface_only); } catch (PlbIOException&) {}
            try { ddoc["parameters"]["dissolution"]["update_interval"].read(dissolCfg.update_interval); } catch (PlbIOException&) {}
            for (plint k = 0; k < 64; ++k) {                 // phase0, phase1, ...
                const std::string pname = "phase" + std::to_string(k);
                SolidPhase ph;
                try { ddoc["parameters"]["dissolution"][pname]["material_number"].read(ph.material_number); }
                catch (PlbIOException&) { break; }           // no more phases
                try { ddoc["parameters"]["dissolution"][pname]["substrate"].read(ph.substrate); } catch (PlbIOException&) {}
                try { ddoc["parameters"]["dissolution"][pname]["full_density"].read(ph.full_density); } catch (PlbIOException&) {}
                try { ddoc["parameters"]["dissolution"][pname]["initial_fill"].read(ph.initial_fill); } catch (PlbIOException&) {}
                try { ddoc["parameters"]["dissolution"][pname]["name"].read(ph.name); } catch (PlbIOException&) { ph.name = pname; }
                try { std::string s3; ddoc["parameters"]["dissolution"][pname]["is_precipitate"].read(s3);
                      std::transform(s3.begin(),s3.end(),s3.begin(),[](unsigned char c){ return std::tolower(c); });
                      ph.is_precipitate = (s3=="true"||s3=="1"||s3=="yes"); } catch (PlbIOException&) {}
                ph.id = (plint) dissolCfg.phases.size() + 1;
                dissolCfg.phases.push_back(ph);
                if (ph.is_precipitate) precip_phase_id = ph.id;
            }
        }
    } catch (PlbIOException&) {}

    if (dissolCfg.enabled) {
        if (dissolCfg.phases.empty()) {
            pcout << "  [DISSOL-VOP] <dissolution> is enabled but no <phase0> block was given, so there is\n"
                  << "               nothing that can dissolve. Declare at least one phase. Terminating.\n";
            return -1;
        }
        if (!(dissolCfg.reopen_fraction > 0.0) || dissolCfg.reopen_fraction >= 1.0) {
            pcout << "  [DISSOL-VOP] <reopen_fraction> must be strictly between 0 and 1 (0.9 is a good\n"
                  << "               default). A value of 1 makes a voxel on the threshold flip every\n"
                  << "               update interval, and each flip re-solves the flow. Terminating.\n";
            return -1;
        }
        pcout << "  [DISSOL-VOP] enabled: " << (plint) dissolCfg.phases.size() << " dissolvable phase(s), "
              << "reopen below " << dissolCfg.reopen_fraction << " x full, "
              << "surface_only=" << dissolCfg.surface_only
              << ", update_interval=" << dissolCfg.update_interval << "\n";
        for (size_t k = 0; k < dissolCfg.phases.size(); ++k) {
            const SolidPhase &ph = dissolCfg.phases[k];
            pcout << "               phase " << ph.id << " \"" << ph.name << "\": material " << ph.material_number
                  << ", substrate " << ph.substrate << ", full " << ph.full_density
                  << " mol/L, starts at " << ph.initial_fill
                  << (ph.is_precipitate ? "  (this is the precipitate)" : "") << "\n";
            if (ph.substrate < 0 || ph.substrate >= num_of_substrates) {
                pcout << "               ERROR: substrate " << ph.substrate << " is outside 0.."
                      << num_of_substrates-1 << ". Terminating.\n";
                return -1;
            }
            if (!(ph.full_density > 0.0)) {
                pcout << "               ERROR: <full_density> must be positive (it is the mol/L of a\n"
                      << "                      completely full voxel, = 1000 / molar volume in cm3/mol).\n"
                      << "                      Terminating.\n";
                return -1;
            }
        }
        if (precip_enabled && precip_phase_id == 0) {
            pcout << "               NOTE: precipitation is on but no phase is marked\n"
                  << "                     <is_precipitate>true</is_precipitate>, so voxels that seal will\n"
                  << "                     not be dissolvable. That is a valid choice; say so deliberately.\n";
        }
    }

    plint rxn_count = kns_count;

    // ════════════════════════════════════════════════════════════════════════════
    // PRINT CONFIGURATION SUMMARY
    // ════════════════════════════════════════════════════════════════════════════
    pcout << "\n┌────────────────────────────────────────────────────────────────────────┐\n";
    pcout << "│ CONFIGURATION SUMMARY                                                 │\n";
    pcout << "├────────────────────────────────────────────────────────────────────────┤\n";
    pcout << "│ Domain: " << nx << " x " << ny << " x " << nz << " = " << nx*ny*nz << " voxels\n";
    pcout << "│ Resolution: dx = " << std::scientific << dx << " m\n";
    pcout << "│ Peclet: " << std::fixed << Pe << "\n";
    pcout << "├────────────────────────────────────────────────────────────────────────┤\n";
    pcout << "│ SUBSTRATES (" << num_of_substrates << "):\n";
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        pcout << "│   [" << iS << "] " << vec_subs_names[iS] << "  C0=" << std::scientific << vec_c0[iS] << " M\n";
    }
    pcout << "├────────────────────────────────────────────────────────────────────────┤\n";
    pcout << "│ MICROBES (" << num_of_microbes << "):\n";
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        pcout << "│   [" << iM << "] " << vec_microbes_names[iM];
        pcout << " type=" << (bmass_type[iM] ? "biofilm" : "planktonic");
        pcout << " solver=" << (solver_type[iM]==1 ? "FD" : (solver_type[iM]==2 ? "CA" : "LBM"));
        pcout << " rxn=" << rxntype::name(reaction_type[iM]) << "\n";
    }
    /* ============================ [SYM] symbolic rate laws ============================
     * Read once, by every rank: the file is small, and reading it everywhere avoids the
     * rank-0 barrier the surrogate needs (which exists only so training cannot happen twice).
     * A file that will not parse, or a name that does not resolve, stops the run.  An index
     * that silently points at the wrong species produces a well-posed simulation with a wrong
     * answer in every voxel, which is the failure mode <exchange_reaction_names> was added to
     * prevent, and the same argument applies here. */
    complab_sym::Program symProg;
    std::vector<plint>   sym_globalId;
    plint sym_count = 0;
    for (plint iM = 0; iM < num_of_microbes; ++iM)
        if (rxntype::usesSymbolic(reaction_type[iM])) ++sym_count;

    if (sym_count > 0) {
        std::string symFile;
        try {
            std::string on;
            XMLreader sdoc("CompLaB.xml");
            sdoc["parameters"]["symbolic"]["enabled"].read(on);
            std::transform(on.begin(), on.end(), on.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            if (on == "true" || on == "1" || on == "yes")
                sdoc["parameters"]["symbolic"]["expressions_file"].read(symFile);
        } catch (PlbIOException&) { }

        if (symFile.empty()) {
            pcout << "  [SYM] " << sym_count << " microbe(s) ask for <reaction_type>symbolic, but\n"
                  << "  [SYM] <symbolic><enabled>true</enabled><expressions_file>...</> is not set.\n";
            return -1;
        }
        std::string serr;
        if (!complab_sym::load(symProg, symFile, &serr)) {
            pcout << "  [SYM] " << serr << "\n"
                  << "  [SYM] A rate law that will not parse is a stop, not a warning.\n";
            return -1;
        }
        pcout << complab_sym::describe(symProg, symFile);

        complab_sym::Runtime &SR = complab_sym::runtime();
        SR.byMicrobe.assign((size_t) num_of_microbes, complab_sym::Binding());
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesSymbolic(reaction_type[iM])) continue;
            complab_sym::Binding &B = SR.byMicrobe[(size_t) iM];
            B.prog = &symProg;

            B.subsOfVar.assign(symProg.vars.size(), -1);
            for (size_t v = 0; v < symProg.vars.size(); ++v) {
                if (symProg.vars[v] == vec_microbes_names[iM]) continue;   /* -1 means biomass */
                for (plint iS = 0; iS < num_of_substrates; ++iS)
                    if (vec_subs_names[iS] == symProg.vars[v]) { B.subsOfVar[v] = (int) iS; break; }
                if (B.subsOfVar[v] < 0) {
                    pcout << "  [SYM] '" << symProg.vars[v] << "' in the vars line is neither a "
                          << "substrate nor the biomass of " << vec_microbes_names[iM] << ".\n";
                    return -1;
                }
            }
            B.subsOfRate.assign(symProg.rates.size(), -1);
            for (size_t r = 0; r < symProg.rates.size(); ++r) {
                if (symProg.rates[r].name == "growth") continue;            /* -1 means bioR */
                for (plint iS = 0; iS < num_of_substrates; ++iS)
                    if (vec_subs_names[iS] == symProg.rates[r].name) { B.subsOfRate[r] = (int) iS; break; }
                if (B.subsOfRate[r] < 0) {
                    pcout << "  [SYM] rate '" << symProg.rates[r].name
                          << "' names neither a substrate nor 'growth'.\n";
                    return -1;
                }
            }
        }
    }
    /* ================================================================================== */

    /* ============================ [GNN] graph network ================================
     * Same shape as the symbolic block above, and for the same reasons: read by every rank
     * because the file is small, and a name that does not resolve stops the run rather than
     * silently pointing at the wrong species.
     *
     * The one difference is what has to be bound.  A .sym file names its own variables, so a
     * variable may be a substrate OR this organism's biomass.  A .gnn file names the species
     * of a reaction network, and every one of them has to be a substrate the simulation
     * actually carries -- biomass enters through the growth output, not through the graph. */
    complab_gnn::Network gnnNet;
    std::vector<plint>   gnn_globalId;
    plint gnn_count = 0;
    for (plint iM = 0; iM < num_of_microbes; ++iM)
        if (rxntype::usesGraphnet(reaction_type[iM])) ++gnn_count;

    if (gnn_count > 0) {
        std::string gnnFile;
        try {
            std::string on;
            XMLreader gdoc("CompLaB.xml");
            gdoc["parameters"]["graphnet"]["enabled"].read(on);
            std::transform(on.begin(), on.end(), on.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            if (on == "true" || on == "1" || on == "yes")
                gdoc["parameters"]["graphnet"]["network_file"].read(gnnFile);
        } catch (PlbIOException&) { }

        if (gnnFile.empty()) {
            pcout << "  [GNN] " << gnn_count << " microbe(s) ask for <reaction_type>graphnet, but\n"
                  << "  [GNN] <graphnet><enabled>true</enabled><network_file>...</> is not set.\n";
            return -1;
        }
        std::string gerr;
        if (!complab_gnn::load(gnnNet, gnnFile, &gerr)) {
            pcout << "  [GNN] " << gerr << "\n"
                  << "  [GNN] A network that will not load is a stop, not a warning.\n";
            return -1;
        }
        pcout << complab_gnn::describe(gnnNet, gnnFile);

        complab_gnn::Runtime &GR = complab_gnn::runtime();
        GR.byMicrobe.assign((size_t) num_of_microbes, complab_gnn::Binding());
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesGraphnet(reaction_type[iM])) continue;
            complab_gnn::Binding &B = GR.byMicrobe[(size_t) iM];
            B.net = &gnnNet;

            B.subsOfSpecies.assign((size_t) gnnNet.nS, -1);
            for (int sp = 0; sp < gnnNet.nS; ++sp) {
                const std::string &nm = gnnNet.species[(size_t) sp];
                for (plint iS = 0; iS < num_of_substrates; ++iS)
                    if (vec_subs_names[iS] == nm) { B.subsOfSpecies[(size_t) sp] = (int) iS; break; }
                if (B.subsOfSpecies[(size_t) sp] < 0) {
                    pcout << "  [GNN] the network names species '" << nm << "', which is not one of\n"
                          << "  [GNN] this simulation's substrates. A network fitted to a different\n"
                          << "  [GNN] set of species cannot be applied to this one.\n";
                    return -1;
                }
            }
            /* The growth rate is the last output when the file declares one.  Without it the
             * network changes the chemistry but no biomass grows, which is a legitimate setup
             * (an abiotic network) and so is not an error -- but it is worth saying out loud. */
            B.growthSlot = gnnNet.hasGrowth ? gnnNet.nS : -1;
            if (!gnnNet.hasGrowth)
                pcout << "  [GNN] note: this network has no growth output, so " << vec_microbes_names[iM]
                      << " will\n  [GNN] change the chemistry but will not itself grow.\n";
        }
    }
    /* ================================================================================== */

    /* ==================== [SYM][GNN] the ABIOTIC half ================================
     * Everything above binds a learned rate law to an ORGANISM, and the sweep that uses it
     * skips any voxel with no biomass.  That is right for a microbe and wrong for abiotic
     * chemistry, which happens in open water whether or not anything is alive nearby.
     *
     * So an abiotic law is loaded separately, bound to no organism, and swept over every fluid
     * voxel by run_symbolic_abiotic3D / run_graphnet_abiotic3D.  It is the learned counterpart of
     * defineAbioticKinetics.hh and it switches on the same way: name a file, and it runs. */
    complab_sym::Program symAbioticProg;
    complab_gnn::Network gnnAbioticNet;
    bool haveSymAbiotic = false, haveGnnAbiotic = false;
    {
        std::string symAbFile, gnnAbFile;
        try {
            XMLreader adoc("CompLaB.xml");
            try { adoc["parameters"]["symbolic"]["abiotic_file"].read(symAbFile); }
            catch (PlbIOException&) { }
            try { adoc["parameters"]["graphnet"]["abiotic_file"].read(gnnAbFile); }
            catch (PlbIOException&) { }
        } catch (PlbIOException&) { }

        if (!symAbFile.empty()) {
            std::string aerr;
            if (!complab_sym::load(symAbioticProg, symAbFile, &aerr)) {
                pcout << "  [SYM-ABIOTIC] " << aerr << "\n";
                return -1;
            }
            pcout << "  [SYM-ABIOTIC] this law is swept over EVERY fluid voxel, biomass or not\n";
            pcout << complab_sym::describe(symAbioticProg, symAbFile);

            complab_sym::Binding &A = complab_sym::runtime().abiotic;
            A.prog = &symAbioticProg;
            A.subsOfVar.assign(symAbioticProg.vars.size(), -1);
            for (size_t v = 0; v < symAbioticProg.vars.size(); ++v) {
                for (plint iS = 0; iS < num_of_substrates; ++iS)
                    if (vec_subs_names[iS] == symAbioticProg.vars[v]) { A.subsOfVar[v] = (int) iS; break; }
                if (A.subsOfVar[v] < 0) {
                    pcout << "  [SYM-ABIOTIC] '" << symAbioticProg.vars[v] << "' is not a substrate.\n"
                          << "  [SYM-ABIOTIC] An abiotic law has no organism, so every variable in it\n"
                          << "  [SYM-ABIOTIC] must be a substrate -- a biomass name here means this\n"
                          << "  [SYM-ABIOTIC] file was meant for <expressions_file>, not this one.\n";
                    return -1;
                }
            }
            A.subsOfRate.assign(symAbioticProg.rates.size(), -1);
            for (size_t r = 0; r < symAbioticProg.rates.size(); ++r) {
                if (symAbioticProg.rates[r].name == "growth") {
                    pcout << "  [SYM-ABIOTIC] this file has a 'growth' rate, and nothing on the\n"
                          << "  [SYM-ABIOTIC] abiotic path grows. Dropping the line silently would\n"
                          << "  [SYM-ABIOTIC] hide the mistake, so this is a stop.\n";
                    return -1;
                }
                for (plint iS = 0; iS < num_of_substrates; ++iS)
                    if (vec_subs_names[iS] == symAbioticProg.rates[r].name) { A.subsOfRate[r] = (int) iS; break; }
                if (A.subsOfRate[r] < 0) {
                    pcout << "  [SYM-ABIOTIC] rate '" << symAbioticProg.rates[r].name
                          << "' does not name a substrate.\n";
                    return -1;
                }
            }
            haveSymAbiotic = true;
        }

        if (!gnnAbFile.empty()) {
            std::string aerr;
            if (!complab_gnn::load(gnnAbioticNet, gnnAbFile, &aerr)) {
                pcout << "  [GNN-ABIOTIC] " << aerr << "\n";
                return -1;
            }
            if (gnnAbioticNet.hasGrowth) {
                pcout << "  [GNN-ABIOTIC] this network has a growth output, and nothing on the\n"
                      << "  [GNN-ABIOTIC] abiotic path grows. Refit it without one, or bind it to an\n"
                      << "  [GNN-ABIOTIC] organism through <network_file> instead.\n";
                return -1;
            }
            pcout << "  [GNN-ABIOTIC] this network is swept over EVERY fluid voxel, biomass or not\n";
            pcout << complab_gnn::describe(gnnAbioticNet, gnnAbFile);

            complab_gnn::Binding &A = complab_gnn::runtime().abiotic;
            A.net = &gnnAbioticNet;
            A.growthSlot = -1;
            A.subsOfSpecies.assign((size_t) gnnAbioticNet.nS, -1);
            for (int sp = 0; sp < gnnAbioticNet.nS; ++sp) {
                const std::string &nm = gnnAbioticNet.species[(size_t) sp];
                for (plint iS = 0; iS < num_of_substrates; ++iS)
                    if (vec_subs_names[iS] == nm) { A.subsOfSpecies[(size_t) sp] = (int) iS; break; }
                if (A.subsOfSpecies[(size_t) sp] < 0) {
                    pcout << "  [GNN-ABIOTIC] the network names species '" << nm
                          << "', which is not a substrate here.\n";
                    return -1;
                }
            }
            haveGnnAbiotic = true;
        }
    }
    /* ================================================================================== */

    pcout << "├────────────────────────────────────────────────────────────────────────┤\n";
    pcout << "│ SOLVERS ENABLED:\n";
    pcout << "│   [" << (kns_count > 0 ? "X" : " ") << "] Kinetics      - " << kns_count << " model(s)\n";
    pcout << "│   [" << (useEquilibrium ? "X" : " ") << "] Equilibrium   - " << eq_component_names.size() << " component(s)\n";
    pcout << "│   [" << (ca_count > 0 ? "X" : " ") << "] CA            - " << ca_count << " microbe(s)\n";
    pcout << "│   [" << (fd_count > 0 ? "X" : " ") << "] FD            - " << fd_count << " microbe(s)\n";
    pcout << "│   [" << (lb_count > 0 ? "X" : " ") << "] LB Diffusion  - " << lb_count << " microbe(s)\n";
    pcout << "│   [" << (mmcfg.glpk_count > 0 ? "X" : " ") << "] FBA (GLPK)    - " << mmcfg.glpk_count << " microbe(s)\n";
    pcout << "│   [" << (mmcfg.cpy_count  > 0 ? "X" : " ") << "] FBA (COBRApy) - " << mmcfg.cpy_count  << " microbe(s)\n";
    pcout << "│   [" << (mmcfg.srg_count  > 0 ? "X" : " ") << "] Surrogate     - " << mmcfg.srg_count  << " microbe(s)\n";
    pcout << "│   [" << (sym_count > 0 ? "X" : " ") << "] Symbolic      - " << sym_count << " microbe(s)\n";
    pcout << "├────────────────────────────────────────────────────────────────────────┤\n";
    pcout << "│ BIOMASS: Bmax=" << max_bMassRho << " kg/m3, threshold=" << thrd_bFilmFrac << "\n";
    pcout << "│ SIMULATION: max_iter=" << ade_maxiTer << ", VTI=" << ade_VTI_iTer << ", CHK=" << ade_CHK_iTer << "\n";
    pcout << "└────────────────────────────────────────────────────────────────────────┘\n\n";

    // Equilibrium setup
    if (useEquilibrium) {
        pcout << "  [EQ] Setting up equilibrium chemistry solver...\n";
        eqSolver.setSpeciesNames(vec_subs_names);
        if (!eq_component_names.empty()) eqSolver.setComponentNames(eq_component_names);
        if (!eq_stoich_matrix.empty()) eqSolver.setStoichiometryMatrix(eq_stoich_matrix);
        if (!eq_logK_values.empty()) eqSolver.setLogK(eq_logK_values);
        eqSolver.setMaxIterations(200);
        eqSolver.setTolerance(1e-10);
        eqSolver.setAndersonDepth(4);
        pcout << "  [EQ] Solver configured: Anderson+PCF, tol=1e-10, maxiter=200\n\n";
        if (!vec_c0.empty()) {
            std::vector<T> _c0chk(vec_c0.begin(), vec_c0.end());
            for (size_t _i=0;_i<_c0chk.size();++_i) _c0chk[_i]=std::max(_c0chk[_i], EquilibriumChemistry<T>::MIN_CONC);
            std::vector<T> _eqc = eqSolver.calculate_species_concentrations(_c0chk);
            const T _MINC = EquilibriumChemistry<T>::MIN_CONC;
            // ---- solved pH (free H+ species is named "Hp") ----
            T _pH = -1.0;
            for (size_t i=0;i<vec_subs_names.size() && i<_eqc.size();++i) if (vec_subs_names[i]=="Hp") { _pH = -std::log10(std::max(_eqc[i], _MINC)); break; }
            if (!eqSolver.didConverge()) {
                // ---- Locate the EXACT culprit from the solver's own mass balance ----
                size_t _ncmp = eq_component_names.size();
                size_t _nsp  = vec_subs_names.size();
                std::vector<T> _Ttar(_ncmp, 0.0), _Tsol(_ncmp, 0.0);
                for (size_t j=0;j<_ncmp;++j){
                    for (size_t i=0;i<_nsp && i<eq_stoich_matrix.size();++i){
                        if (j>=eq_stoich_matrix[i].size()) continue;
                        T s = eq_stoich_matrix[i][j];
                        if (std::abs(s) <= 1e-10) continue;   // species i does not contain component j
                        T ci  = (i<_c0chk.size())? std::max(std::min(_c0chk[i], (T)10.0), _MINC) : _MINC;
                        T cei = (i<_eqc.size())?   std::max(std::min(_eqc[i],   (T)10.0), _MINC) : _MINC;
                        _Ttar[j] += s*ci;   // TOTAL you asked for (from CompLaB.xml)
                        _Tsol[j] += s*cei;  // TOTAL the solver could actually hold
                    }
                    if (_Ttar[j] < _MINC) _Ttar[j] = _MINC;
                    if (_Tsol[j] < _MINC) _Tsol[j] = _MINC;
                }
                // worst REAL component (absent/zero-total components are handled by the solver guard)
                int _wj = -1; T _wrel = -1.0;
                for (size_t j=0;j<_ncmp;++j){
                    if (_Ttar[j] <= _MINC*10.0) continue;              // absent from this water
                    T rel = std::abs(_Tsol[j]-_Ttar[j]) / _Ttar[j];   // relative mass-balance violation
                    if (rel > _wrel){ _wrel = rel; _wj = (int)j; }
                }
                std::string _cn = (_wj>=0)? eq_component_names[_wj] : std::string("(unknown)");
                int _csub = -1;
                for (size_t i=0;i<_nsp;++i) if (vec_subs_names[i]==_cn){ _csub=(int)i; break; }

                pcout << "\n╔═══════════════════════════════════════════════════════════════════════╗\n";
                if (_wj < 0 || _wrel < 1e-4) {
                    // Every REAL component balances -> the leftover residual is the harmless
                    // zero-total-component artifact. If you still see this, the RUNNING BINARY
                    // predates the zero-total fix in complab3d_processors_part4_eqsolver.hh.
                    pcout << "║  NOTE: your INITIAL water is actually FEASIBLE - every real component balances.\n";
                    pcout << "║  This warning is only the zero-total-component artifact, which means the BINARY\n";
                    pcout << "║  you are running was built BEFORE the equilibrium fix. Rebuild to clear it:\n";
                    pcout << "║      cd build && make clean && make\n";
                    if (_pH>=0.0) pcout << "║  (solved pH = " << std::fixed << _pH << ", which is correct - nothing to change in CompLaB.xml.)\n";
                } else {
                    pcout << "║  INPUT ERROR (non-fatal): the equilibrium chemistry in CompLaB.xml is INFEASIBLE.\n";
                    pcout << "║  The speciation solver could not satisfy mass balance for your STARTING water.\n";
                    pcout << "║  (global residual=" << std::scientific << eqSolver.getLastResidual() << ", used " << eqSolver.getLastIterations() << "/" << eqSolver.getMaxIterations() << " iters)\n";
                    pcout << "║\n";
                    pcout << "║  WHAT IS WRONG  (exact culprit, computed from the solver's own mass balance):\n";
                    pcout << "║    Component  '" << _cn << "'  (equilibrium component #" << _wj << ") does NOT conserve mass.\n";
                    pcout << "║    You asked for TOTAL " << _cn << " = " << std::scientific << _Ttar[_wj] << " M,\n";
                    pcout << "║    but the only self-consistent speciation the solver can reach holds " << std::scientific << _Tsol[_wj] << " M\n";
                    pcout << "║    -> mass-balance violation = " << std::fixed << (_wrel*100.0) << " %  (this is what blocks convergence).\n";
                    if (_pH>=0.0) pcout << "║    Solver best-effort pH for this water = " << std::fixed << _pH << " .\n";
                    pcout << "║\n";
                    pcout << "║  WHY  (physical cause):\n";
                    if (_Tsol[_wj] < _Ttar[_wj]) {
                        pcout << "║    You are DEMANDING more '" << _cn << "' than the logK/stoichiometry in CompLaB.xml can hold\n";
                        pcout << "║    in solution at the pH set by the other components. Its free ion hits the numerical floor,\n";
                        pcout << "║    so the requested total can never be stored -> the fixed point cannot close.\n";
                    } else {
                        pcout << "║    The species built from '" << _cn << "' already EXCEED the total you asked for at this pH;\n";
                        pcout << "║    the logK values make '" << _cn << "' too abundant for the total you set.\n";
                    }
                    pcout << "║\n";
                    pcout << "║  WHICH INPUT TO FIX  (file -> block -> parameter):\n";
                    pcout << "║    File:      CompLaB.xml\n";
                    if (_csub>=0) {
                        pcout << "║    Block:     <substrate" << _csub << ">   (name = " << _cn << ")\n";
                        pcout << "║    Parameter: <initial_concentration>  (currently " << std::scientific << vec_c0[_csub] << " M)\n";
                        if (_cn=="Hp") {
                            pcout << "║    NOTE: 'Hp' is the TOTAL proton (TOTH), NOT free H+. For pH-7 water it is ~1.9e-3, not 1e-7.\n";
                            pcout << "║    FIX:  set <initial_concentration> AND <left_boundary_condition> for <substrate" << _csub << "> to the\n";
                            pcout << "║          total proton of your recipe (protons carried by CO2/H2PO4/NH4/HS... ~1.9e-3 for this water).\n";
                        } else if (_Tsol[_wj] < _Ttar[_wj]) {
                            pcout << "║    FIX:  LOWER <initial_concentration> (and <left_boundary_condition>) for <substrate" << _csub << ">\n";
                            pcout << "║          to a value the chemistry can hold, OR add the missing product/mineral species for '" << _cn << "',\n";
                            pcout << "║          OR correct this species' logK. (This is pure XML equilibrium data - defineKinetics is NOT involved.)\n";
                        } else {
                            pcout << "║    FIX:  RAISE <initial_concentration> (and <left_boundary_condition>) for <substrate" << _csub << ">,\n";
                            pcout << "║          OR correct the logK of the '" << _cn << "' species so it is less abundant at this pH.\n";
                        }
                    } else {
                        pcout << "║    Parameter: the <initial_concentration> of the <substrate*> whose name is '" << _cn << "'.\n";
                    }
                }
                pcout << "║\n";
                pcout << "║  The run will CONTINUE; chemistry for the affected component is untrustworthy until fixed.\n";
                pcout << "╚═══════════════════════════════════════════════════════════════════════╝\n";
            } else {
                pcout << "  [EQ] Initial-composition feasibility: converged (iters=" << eqSolver.getLastIterations() << ", residual=" << std::scientific << eqSolver.getLastResidual();
                if (_pH>=0.0) pcout << ", pH=" << std::fixed << _pH;
                pcout << ").\n";
            }
        }
    }

    std::string str_inputDir=input_path, str_outputDir=output_path;
    if (std::to_string(str_inputDir.back()).compare("/")!=0) { str_inputDir+="/"; }
    if (std::to_string(str_outputDir.back()).compare("/")!=0) { str_outputDir+="/"; }

    // ---- OPTIONAL METABOLIC LAYER: load the models and build the solvers -------
    //   Now that the input directory is known we can read one metabolic model per
    //   FBA microbe (the XML that extractMM.py produces) and build the persistent
    //   solver state: one glp_prob* per GLPK microbe, or the cobra model objects.
    //   Both are created ONCE and reused at every voxel and every time step.
    char *pyFileName = (char*)"complab3d_cobrapy";
    if (mmcfg.mm_count > 0) {
        if (load_metabolic_models3D(mmcfg, str_inputDir, reaction_type, num_of_microbes) != 0) return -1;
    }
    if (mmcfg.anyEnabled()) {
        if (setup_metabolic_solvers(mmcfg, reaction_type, num_of_microbes, pyFileName, src_path) != 0) return -1;
    }

    // ════════════════════════════════════════════════════════════════════════════
    // PHASE 2: GEOMETRY AND FLOW SETUP
    // ════════════════════════════════════════════════════════════════════════════
    struct stat statStruct;
    stat(output_path, &statStruct);

    pcout << "┌────────────────────────────────────────────────────────────────────────┐\n";
    pcout << "│ PHASE 2: GEOMETRY AND FLOW SETUP                                      │\n";
    pcout << "└────────────────────────────────────────────────────────────────────────┘\n";
    pcout << "  Main:   " << str_mainDir << "\n";
    pcout << "  Input:  " << main_path << "/" << input_path << "\n";
    pcout << "  Output: " << main_path << "/" << output_path << "\n";
    
    if (S_ISDIR(statStruct.st_mode)) {} else { mkdir(output_path, 0777); }
    global::directories().setOutputDir(str_outputDir);

    T PoreMeanU=0, PoreMaxUx=0;
    plint iT = 0;
    T nsLatticeTau = tau;
    T nsLatticeOmega = 1 / nsLatticeTau;
    T nsLatticeNu = NSDES<T>::cs2*(nsLatticeTau-0.5);
    char *ns_read_filename = strcat(strdup(str_inputDir.c_str()),ns_filename);

    pcout << "  [GEOM] Reading " << geom_filename << "...\n";
    MultiScalarField3D<int> geometry(nx,ny,nz);
    readGeometry(str_inputDir+geom_filename, geometry);
    saveGeometry("inputGeom", geometry);
    pcout << "  [GEOM] Geometry loaded\n";

    MultiScalarField3D<int> distanceDomain(nx,ny,nz);
    distanceDomain = geometry;
    std::vector< std::vector< std::vector<plint> > > distVec(nx);
    for (plint iX=0; iX<nx; ++iX) {
        distVec[iX]=std::vector< std::vector<plint> > (ny);
        for (plint iY=0; iY<ny; ++iY) { distVec[iX][iY]=std::vector<plint> (nz); }
    }
    calculateDistanceFromSolid(distanceDomain, no_dynamics, bounce_back, distVec);
    applyProcessingFunctional(new createDistanceDomain3D<int> (distVec), distanceDomain.getBoundingBox(), distanceDomain);

    MultiScalarField3D<int> ageDomain(nx,ny,nz);
    ageDomain = geometry;
    applyProcessingFunctional(new createAgeDomain3D<int> (pore_dynamics, bounce_back, no_dynamics), ageDomain.getBoundingBox(), ageDomain);
    pcout << "  [GEOM] Distance and age fields ready\n";
    
    if (track_performance == 1) { pcout << "  [PERF] Performance tracking ON - VTI disabled\n"; }

    pcout << "  [NS] Initializing fluid lattice (deltaP=" << deltaP << ")...\n";
    MultiBlockLattice3D<T,NSDES> nsLattice(nx, ny, nz, new IncBGKdynamics<T,NSDES>(nsLatticeOmega));
    util::ValueTracer<T> ns_convg1(1.0,1000.0,ns_converge_iT1);
    NSdomainSetup(nsLattice, createLocalBoundaryCondition3D<T,NSDES>(), geometry, deltaP, nsLatticeOmega, pore_dynamics, bounce_back, no_dynamics, bio_dynamics, vec_permRatio);

    // NS main loop
    global::timer("NS").start();
    if (Pe == 0) { pcout << "  [NS] Pe=0, skipping flow solver\n"; }
    else {
        pcout << "  [NS] tau=" << nsLatticeTau << ", omega=" << nsLatticeOmega << ", nu=" << nsLatticeNu << "\n";
        if (read_NS_file == 1 && track_performance == 0) {
            pcout << "  [NS] Loading checkpoint...\n";
            try { loadBinaryBlock(nsLattice, strcat(ns_read_filename,".chk")); }
            catch (PlbIOException& exception) { pcout << "  [NS] ERROR: " << exception.what() << "\n"; return -1; }
            if (ns_rerun_iT0 > 0) {
                iT = ns_rerun_iT0;
                for (; iT < ns_maxiTer_1; ++iT) {
                    nsLattice.collideAndStream();
                    ns_convg1.takeValue(getStoredAverageEnergy(nsLattice),true);
                    if (ns_convg1.hasConverged()) break;
                }
            }
        }
        else {
            pcout << "  [NS] Running new simulation...\n";
            for (; iT < ns_maxiTer_1; ++iT) {
                nsLattice.collideAndStream();
                ns_convg1.takeValue(getStoredAverageEnergy(nsLattice),true);
                if (ns_convg1.hasConverged()) break;
            }
        }
        pcout << "  [NS] Converged at iter=" << iT << "\n";

        // Calculate velocities
        if (bfilm_count > 0) {
            plint totalCount = 0; T totalVel = 0;
            for (size_t iT = 0; iT < pore_dynamics.size(); ++iT) {
                plint poreCount = MaskedScalarCounts3D(Box3D(1,nx-2,0,ny-1,0,nz-1), geometry, pore_dynamics[iT]);
                totalCount += poreCount;
                totalVel += computeAverage(*computeVelocityNorm(nsLattice, Box3D(1,nx-2,0,ny-1,0,nz-1)), geometry, pore_dynamics[iT]) * poreCount;
            }
            for (plint iT0 = 0; iT0 < bfilm_count; ++iT0) {
                plint bFilmCount = 0;
                for (size_t iT1 = 0; iT1 < bio_dynamics[iT0].size(); ++iT1) {
                    bFilmCount += MaskedScalarCounts3D(Box3D(1,nx-2,0,ny-1,0,nz-1), geometry, bio_dynamics[iT0][iT1]);
                }
                totalCount += bFilmCount;
                totalVel += computeAverage(*computeVelocityNorm(nsLattice, Box3D(1,nx-2,0,ny-1,0,nz-1)), geometry, bio_dynamics[iT0][0]) * bFilmCount;
            }
            PoreMeanU = totalVel / totalCount;
        }
        else { PoreMeanU = computeAverage(*computeVelocityNorm(nsLattice, Box3D(1,nx-2,0,ny-1,0,nz-1))); }

        PoreMaxUx = computeMax(*computeVelocityComponent(nsLattice, Box3D(1,nx-2,0,ny-1,0,nz-1), 0));
        DarcyOutletUx = computeAverage(*computeVelocityComponent(nsLattice, Box3D(nx-2,nx-2, 0,ny-1, 0,nz-1), 0));
        
        D_lattice_fixed = RXNDES<T>::cs2 * (tau_ADE_fixed - 0.5);
        permeability = DarcyOutletUx * nsLatticeNu * charcs_length / deltaP;
        pcout << "  [NS] Permeability k=" << permeability << " (lattice)\n";
        
        u_target = Pe * D_lattice_fixed / charcs_length;
        deltaP_new = u_target * nsLatticeNu * charcs_length / permeability;
        
        if (std::abs(deltaP_new - deltaP) / deltaP > 0.01) {
            pcout << "  [NS] Re-running with corrected deltaP=" << deltaP_new << "\n";
            NSdomainSetup(nsLattice, createLocalBoundaryCondition3D<T,NSDES>(), geometry, deltaP_new, nsLatticeOmega, pore_dynamics, bounce_back, no_dynamics, bio_dynamics, vec_permRatio);
            ns_convg1.resetValues();
            for (plint iT2 = 0; iT2 < ns_maxiTer_1; ++iT2) {
                nsLattice.collideAndStream();
                ns_convg1.takeValue(getStoredAverageEnergy(nsLattice), true);
                if (ns_convg1.hasConverged()) break;
            }
            PoreMeanU = computeAverage(*computeVelocityNorm(nsLattice, Box3D(1,nx-2,0,ny-1,0,nz-1)));
            PoreMaxUx = computeMax(*computeVelocityComponent(nsLattice, Box3D(1,nx-2,0,ny-1,0,nz-1), 0));
            DarcyOutletUx = computeAverage(*computeVelocityComponent(nsLattice, Box3D(nx-2,nx-2,0,ny-1,0,nz-1), 0));
            deltaP = deltaP_new;
        }
        
        u_final = DarcyOutletUx;
        Pe_achieved = u_final * charcs_length / D_lattice_fixed;
        pcout << "  [NS] Pe achieved=" << Pe_achieved << " (target=" << Pe << ")\n";
        
        StabilityReport stability = performStabilityChecks(PoreMaxUx, nsLatticeTau, tau_ADE_fixed, D_lattice_fixed);
        printStabilityReport(stability);
        if (!stability.all_ok) {
            pcout << "\n╔══════════════════════════════════════════════════════════════════════════╗\n";
            pcout << "║  INPUT ERROR - stability check FAILED before the run started. Fix CompLaB.xml:\n";
            if (!stability.Ma_ok)      pcout << "║    Ma = " << stability.Ma << " (must be < 1; aim <= 0.02): lower deltaP or Pe.\n";
            if (!stability.CFL_ok)     pcout << "║    CFL = " << stability.CFL << " (must be < 1): lower deltaP or Pe.\n";
            if (!stability.tau_NS_ok)  pcout << "║    tau_NS = " << stability.tau_NS << " (must be 0.5-2): adjust tau / viscosity in CompLaB.xml.\n";
            if (!stability.tau_ADE_ok) pcout << "║    tau_ADE = " << stability.tau_ADE << " (must be 0.5-2): adjust solute diffusivity in CompLaB.xml.\n";
            pcout << "╚══════════════════════════════════════════════════════════════════════════╝\n";
            return -1;
        }
        else if (stability.has_warnings) {
            pcout << "  [STABILITY] INPUT WARNING: ";
            if (stability.Ma_warning)  pcout << "Ma=" << stability.Ma << ">0.3 (flow re-solves may diverge; lower deltaP/Pe in CompLaB.xml). ";
            if (!stability.Pe_grid_ok) pcout << "Pe_grid=" << stability.Pe_grid << ">2 (advection may overshoot into negatives; lower Pe or raise diffusivity in CompLaB.xml). ";
            pcout << "\n";
        }

        T Ma = PoreMaxUx/sqrt(RXNDES<T>::cs2);
        if (Ma > 1) { pcout << "  [NS] ERROR: Ma=" << Ma << " > 1\n"; return -1; }
    }
    global::timer("NS").stop();
    T nstime = global::timer("NS").getTime();

    // Flow-only run: stop here.  MetabolicScopeGuard releases the solvers.
    if (ade_maxiTer == 0) { pcout << "  [ADE] ade_maxiTer=0, done.\n"; return 0; }

    // ════════════════════════════════════════════════════════════════════════════
    // PHASE 3: REACTIVE TRANSPORT SETUP
    // ════════════════════════════════════════════════════════════════════════════
    pcout << "\n┌────────────────────────────────────────────────────────────────────────┐\n";
    pcout << "│ PHASE 3: REACTIVE TRANSPORT SETUP                                     │\n";
    pcout << "└────────────────────────────────────────────────────────────────────────┘\n";

    T refNu, refTau;
    if (Pe > thrd) {
        refNu = PoreMeanU * charcs_length / Pe;
        refTau = refNu * RXNDES<T>::invCs2 + 0.5;
        if (refTau > 2 || refTau <= 0.5) { pcout << "  [ADE] ERROR: tau=" << refTau << " invalid\n"; return -1; }
    }
    else { refTau = tau; refNu = RXNDES<T>::cs2 * (refTau - 0.5); }
    T refOmega = 1/refTau;
    T ade_dt = refNu * dx * dx / vec_solute_poreD[0];

    std::vector<T> substrNUinPore(num_of_substrates), substrTAUinPore(num_of_substrates), substrOMEGAinPore(num_of_substrates), substrOMEGAinbFilm(num_of_substrates);
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        if (iS == 0) { substrNUinPore[iS]=refNu; substrTAUinPore[iS]=refTau; substrOMEGAinPore[iS]=refOmega; }
        else {
            substrNUinPore[iS] = substrNUinPore[0]*vec_solute_poreD[iS]/vec_solute_poreD[0];
            substrTAUinPore[iS] = substrNUinPore[iS]*RXNDES<T>::invCs2+0.5;
            substrOMEGAinPore[iS] = 1/substrTAUinPore[iS];
        }
        substrOMEGAinbFilm[iS] = 1/(refNu*vec_solute_bFilmD[iS]/vec_solute_poreD[0]*RXNDES<T>::invCs2+0.5);
    }

    std::vector<T> bioNUinPore(num_of_microbes), bioTAUinPore(num_of_microbes), bioOMEGAinPore(num_of_microbes), bioOMEGAinbFilm(num_of_microbes), bioTAUinbFilm(num_of_microbes);
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (vec_bMass_poreD[iM] > 0) {
            bioNUinPore[iM] = refNu * vec_bMass_poreD[iM] / vec_solute_poreD[0];
            bioTAUinPore[iM] = bioNUinPore[iM] * RXNDES<T>::invCs2 + 0.5;
            bioOMEGAinPore[iM] = 1/bioTAUinPore[iM];
        }
        else { bioNUinPore[iM] = 0.; bioTAUinPore[iM] = 0.; bioOMEGAinPore[iM] = 0.; }
        if (vec_bMass_bFilmD[iM] > 0) {
            bioOMEGAinbFilm[iM] = 1/(refNu*vec_bMass_bFilmD[iM]/vec_bMass_poreD[iM]*RXNDES<T>::invCs2+0.5);
            bioTAUinbFilm[iM] = 1/bioOMEGAinbFilm[iM];
        }
        else { bioOMEGAinbFilm[iM] = 0.; bioTAUinbFilm[iM] = 0.; }
    }

    pcout << "  [ADE] dt=" << ade_dt << " s/iter, total=" << ade_maxiTer*ade_dt << " s\n";

    // Create substrate lattices
    pcout << "  [ADE] Creating " << num_of_substrates << " substrate lattices...\n";
    MultiBlockLattice3D<T,RXNDES> substrLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(refOmega));
    std::vector< MultiBlockLattice3D<T,RXNDES> > vec_substr_lattices(num_of_substrates, substrLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > dC(num_of_substrates, substrLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > dC0(num_of_substrates, substrLattice);
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        soluteDomainSetup(vec_substr_lattices[iS], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry,
                          substrOMEGAinbFilm[iS], substrOMEGAinPore[iS], pore_dynamics, bounce_back, no_dynamics, bio_dynamics,
                          vec_c0[iS], vec_left_btype[iS], vec_right_btype[iS], vec_left_bcondition[iS], vec_right_bcondition[iS],
                          vec_immobile[iS]);   /* [FIX] an immobile mineral must be readable on solid voxels */
        soluteDomainSetup(dC[iS], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry,
                          substrOMEGAinbFilm[iS], substrOMEGAinPore[iS], pore_dynamics, bounce_back, no_dynamics, bio_dynamics,
                          0., vec_left_btype[iS], vec_right_btype[iS], vec_left_bcondition[iS], vec_right_bcondition[iS]);
    }
    dC0=dC;

    // Create biomass lattices
    pcout << "  [ADE] Creating " << bfilm_count << " biofilm + " << bfree_count << " planktonic lattices...\n";
    MultiBlockLattice3D<T,RXNDES> initbFilmLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice3D<T,RXNDES> copybFilmLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice3D<T,RXNDES> initbFreeLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice3D<T,RXNDES> copybFreeLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    std::vector< MultiBlockLattice3D<T,RXNDES> > vec_bFilm_lattices(bfilm_count, initbFilmLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > vec_bFcopy_lattices(bfilm_count, copybFilmLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > vec_bFree_lattices(bfree_count, initbFreeLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > vec_bPcopy_lattices(bfree_count, copybFreeLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > dBf(bfilm_count, initbFilmLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > dBp(bfree_count, initbFreeLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > dBf0(bfilm_count, initbFilmLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > dBp0(bfree_count, initbFreeLattice);

    plint tmpIT0=0, tmpIT1=0;
    std::vector<plint> loctrack;
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (bmass_type[iM]==1) {
            bmassDomainSetup(vec_bFilm_lattices[tmpIT0], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry, bioOMEGAinPore[iM], bioOMEGAinbFilm[iM],
                             pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM]);
            bmassDomainSetup(vec_bFcopy_lattices[tmpIT0], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry, 0., 0.,
                             pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM]);
            bmassDomainSetup(dBf[tmpIT0], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry, 0., 0.,
                             pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM]);
            loctrack.push_back(tmpIT0); ++tmpIT0;
        }
        else {
            if (solver_type[iM]==3) {
                soluteDomainSetup(vec_bFree_lattices[tmpIT1], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry, bioOMEGAinbFilm[iM], bioOMEGAinPore[iM],
                                  pore_dynamics, bounce_back, no_dynamics, bio_dynamics, vec_b0_free[tmpIT1], bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM]);
                bmassDomainSetup(vec_bPcopy_lattices[tmpIT1], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry, 0., 0.,
                                 pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM]);
                bmassDomainSetup(dBp[tmpIT1], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry, 0., 0.,
                                 pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[iM], bio_right_btype[iM], bio_left_bcondition[iM], bio_right_bcondition[iM]);
            }
            else if (solver_type[iM]==1) { pcout << "  [ADE] ERROR: FD not implemented\n"; return -1; }
            loctrack.push_back(tmpIT1); ++tmpIT1;
        }
    }
    dBp0=dBp; dBf0=dBf;
    
    MultiBlockLattice3D<T,RXNDES> totalbFilmLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    // Only setup biomass lattices if we have microbes (avoid out-of-bounds access in abiotic mode)
    if (num_of_microbes > 0) {
        bmassDomainSetup(totalbFilmLattice, createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry, bioOMEGAinPore[0], bioOMEGAinbFilm[0],
                         pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[0], bio_right_btype[0], bio_left_bcondition[0], bio_right_bcondition[0]);
        bmassDomainSetup(copybFilmLattice, createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry, 0., 0.,
                         pore_dynamics, bounce_back, no_dynamics, bio_dynamics, bio_left_btype[0], bio_right_btype[0], bio_left_bcondition[0], bio_right_bcondition[0]);
    } else {
        // Abiotic mode: initialize lattices with zero density (no biomass)
        Array<T,3> zeroVelocity(0., 0., 0.);
        initializeAtEquilibrium(totalbFilmLattice, totalbFilmLattice.getBoundingBox(), (T)0.0, zeroVelocity);
        initializeAtEquilibrium(copybFilmLattice, copybFilmLattice.getBoundingBox(), (T)0.0, zeroVelocity);
    }

    // Initialize biomass
    for (plint iM = 0; iM < bfilm_count; ++iM) {
        applyProcessingFunctional(new initializeScalarLattice3D<T,RXNDES,int>(vec_b0_film[iM], bio_dynamics[iM]), vec_bFilm_lattices[iM].getBoundingBox(), vec_bFilm_lattices[iM], geometry);
        std::vector<T> vec_b1(vec_b0_film[iM].size(),0.);
        applyProcessingFunctional(new initializeScalarLattice3D<T,RXNDES,int>(vec_b1, bio_dynamics[iM]), vec_bFcopy_lattices[iM].getBoundingBox(), vec_bFcopy_lattices[iM], geometry);
        initTotalbFilmLatticeDensity(vec_bFilm_lattices[iM], totalbFilmLattice);
    }

    if (bfilm_count > 0) {
        diag_initial_biomass = computeMax(*computeDensity(totalbFilmLattice));
        pcout << "  [ADE] Initial max biomass: " << diag_initial_biomass << " kg/m3\n";
    }

    // Mask and distance lattices
    MultiBlockLattice3D<T,RXNDES> maskLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice3D<T,RXNDES> ageLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    MultiBlockLattice3D<T,RXNDES> distLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    defineMaskLatticeDynamics(totalbFilmLattice, maskLattice, thrd_bFilmFrac);
    applyProcessingFunctional(new CopyGeometryScalar2maskLattice3D<T,RXNDES,int>(bio_dynamics), maskLattice.getBoundingBox(), maskLattice, geometry);
    applyProcessingFunctional(new CopyGeometryScalar2ageLattice3D<T,RXNDES,int>(), ageLattice.getBoundingBox(), ageLattice, ageDomain);
    applyProcessingFunctional(new CopyGeometryScalar2distLattice3D<T,RXNDES,int>(), distLattice.getBoundingBox(), distLattice, distanceDomain);

    // [DISSOL-VOP] The phase lattice records WHICH dissolvable solid occupies each
    //   voxel (0 = none).  It sits alongside maskLattice / ageLattice / distLattice
    //   and uses the same carry-an-integer-on-a-D3Q7-lattice pattern.  It is
    //   allocated unconditionally -- one lattice is cheap next to the substrates --
    //   but stays all zeros and is never read when dissolution is off.
    MultiBlockLattice3D<T,RXNDES> phaseLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));
    {
        std::vector<plint> phMat, phId;
        for (size_t k = 0; k < dissolCfg.phases.size(); ++k) {
            // A phase that IS the precipitate starts nowhere: its voxels are stamped
            // by precipNodeConversion3D as they seal.  Grain phases are seeded here.
            if (dissolCfg.phases[k].is_precipitate) continue;
            phMat.push_back(dissolCfg.phases[k].material_number);
            phId .push_back(dissolCfg.phases[k].id);
        }
        applyProcessingFunctional(new initPhaseLattice3D<T,RXNDES,int>(phMat, phId),
                                  phaseLattice.getBoundingBox(), phaseLattice, geometry);

        // Seed each grain phase's mineral inventory into its immobile substrate.
        // <initial_concentration> is applied uniformly over the whole domain, so it
        // cannot put calcite only inside the grains -- this can.
        for (size_t k = 0; k < dissolCfg.phases.size(); ++k) {
            const SolidPhase &ph = dissolCfg.phases[k];
            if (ph.is_precipitate || !(ph.initial_fill > 0.0)) continue;
            std::vector<T>     fill(1, ph.initial_fill);
            std::vector<plint> mat (1, ph.material_number);
            applyProcessingFunctional(new initializeScalarLattice3D<T,RXNDES,int>(fill, mat),
                                      vec_substr_lattices[ph.substrate].getBoundingBox(),
                                      vec_substr_lattices[ph.substrate], geometry);
            pcout << "  [DISSOL-VOP] seeded " << ph.initial_fill << " mol/L of \"" << ph.name
                  << "\" into substrate " << ph.substrate << " on material " << ph.material_number << "\n";
            
            /* [FIX] Read the seed back.  Two different mistakes used to pass silently here and then
             * destroy the run on the first iteration.
             *   1. If the phase substrate is not declared <immobile>, its solid voxels carry
             *      NoDynamics, whose density reads back as zero whatever was written.
             *   2. If initial_fill is below the reopen threshold, the phase is under-filled by
             *      construction and every voxel of it is reopened as pore at iT = 0.
             * Either way the declared grain pack vanishes on the first step, which looks like a
             * physical result and is not one. */
            {
                const T seen = computeMax(*computeDensity(vec_substr_lattices[ph.substrate]));
                if (!(seen > 0.5 * ph.initial_fill)) {
                    pcout << "  [DISSOL-VOP] ERROR: \"" << ph.name << "\" was seeded with "
                          << ph.initial_fill << " mol/L but reads back at most " << seen
                          << ".\n  [DISSOL-VOP] Substrate " << ph.substrate << " is almost certainly"
                          << " not declared <immobile>true</immobile>,\n"
                          << "  [DISSOL-VOP] so its solid voxels cannot hold a value. The phase could"
                          << " never dissolve, and every\n  [DISSOL-VOP] voxel of it would be"
                          << " reopened as pore on the first step.\n";
                    return -1;
                }
            }
            if (!(ph.initial_fill >= dissolCfg.reopen_fraction * ph.full_density)) {
                pcout << "  [DISSOL-VOP] ERROR: \"" << ph.name << "\" starts at " << ph.initial_fill
                      << " mol/L, below its own reopen threshold of "
                      << (dissolCfg.reopen_fraction * ph.full_density) << " mol/L.\n"
                      << "  [DISSOL-VOP] Every voxel of this phase would be turned to pore at"
                      << " iteration 0. Raise <initial_fill>,\n  [DISSOL-VOP] lower"
                      << " <full_density>, or lower <reopen_fraction>.\n";
                return -1;
            }
        }
    }

    // Lattice vectors for the dissolution processors.
    //   ptr_dissol       = [C.., dC.., mask, phase]      length 2*S + 2
    //   ptr_dissol_conv  = [C.., mask, phase]            length S + 2
    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_dissol, ptr_dissol_conv;
    for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_dissol.push_back(&vec_substr_lattices[iS]);
    for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_dissol.push_back(&dC[iS]);
    ptr_dissol.push_back(&maskLattice);
    ptr_dissol.push_back(&phaseLattice);
    for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_dissol_conv.push_back(&vec_substr_lattices[iS]);
    ptr_dissol_conv.push_back(&maskLattice);
    ptr_dissol_conv.push_back(&phaseLattice);

    pcout << "  [ADE] All lattices created\n";

    // Pointer vectors
    std::vector< MultiBlockLattice3D<T, RXNDES>* > substrate_lattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) { substrate_lattices.push_back(&vec_substr_lattices[iS]); }
    substrate_lattices.push_back(&maskLattice);

    std::vector< MultiBlockLattice3D<T, RXNDES>* > planktonic_lattices;
    for (size_t iP = 0; iP < vec_bFree_lattices.size(); ++iP) { planktonic_lattices.push_back(&vec_bFree_lattices[iP]); }
    planktonic_lattices.push_back(&maskLattice);

    // Kinetics lattices.  The filter is rxntype::usesKinetics() rather than
    // "== 1" so that the combined reaction types (glpk_and_kinetics and friends)
    // also get their kinetics term.  Layout: [C.., B.., dC.., dB.., mask].
    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_kns_lattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) { ptr_kns_lattices.push_back(&vec_substr_lattices[iS]); }
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (!rxntype::usesKinetics(reaction_type[iM])) continue;
        if (bmass_type[iM]==1) ptr_kns_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
        else                   ptr_kns_lattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
    }
    for (plint iS = 0; iS < num_of_substrates; ++iS) { ptr_kns_lattices.push_back(&dC[iS]); }
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (!rxntype::usesKinetics(reaction_type[iM])) continue;
        if (bmass_type[iM]==1) ptr_kns_lattices.push_back(&dBf[loctrack[iM]]);
        else                   ptr_kns_lattices.push_back(&dBp[loctrack[iM]]);
    }
    ptr_kns_lattices.push_back(&maskLattice);

    // ========================================================================
    // [RATE-OUT] OPTIONAL REACTION RATE OUTPUT
    //   Off unless CompLaB.xml carries <IO><save_reaction_rates>true</...>.
    //   When it is off, nothing here allocates and nothing here runs, so an
    //   input file written before this feature existed behaves exactly as it
    //   did before.
    // ========================================================================
    bool        rate_enabled  = false;
    plint       rate_interval = 0;
    std::string rate_filename = "rateLattice", rate_unit;
    readRateSettings(rate_enabled, rate_interval, rate_filename, rate_unit);

    // The biomass vector handed to defineKinetics.hh is COMPACTED over the
    // microbes that actually use kinetics, so the microbe channels must be
    // named from the same compacted list or every label would be off by one.
    std::vector<std::string> kns_microbe_names;
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (!rxntype::usesKinetics(reaction_type[iM])) continue;
        kns_microbe_names.push_back(vec_microbes_names[iM]);
    }

    RateOutput rateOut = buildRateOutput(
        rate_enabled, rate_interval, rate_filename, rate_unit,
        (enable_kinetics && kns_count > 0),
        (enable_abiotic_kinetics && num_of_substrates > 0),
        (precip_enabled && precip_surfaceOnly) ? 1 : 0,
        vec_subs_names, kns_microbe_names);
    if (rateOut.interval <= 0) rateOut.interval = ade_VTI_iTer;

    std::vector< MultiBlockLattice3D<T,RXNDES> >  rateLat;
    std::vector< MultiBlockLattice3D<T,RXNDES>* > ptr_rate_lattices;
    if (rateOut.enabled && rateOut.nch() > 0) {
        double _rateMB = (double) rateOut.nch() * (double) nx * (double) ny * (double) nz
                       * 7.0 * (double) sizeof(T) / 1048576.0;
        pcout << "  [RATE] Reaction rate output ON: " << rateOut.nch() << " channels, every "
              << rateOut.interval << " iterations, about " << (long) _rateMB << " MB held.\n";
        if (ade_VTI_iTer > 0 && rateOut.interval % ade_VTI_iTer != 0) {
            pcout << "  [RATE] WARNING <rate_interval> is not a multiple of <save_VTK_interval>. "
                  << "Rates are sampled on VTK output steps only, so this will write less often "
                  << "than asked, or not at all.\n";
        }
        // Plain copies of the substrate lattice: uniform advection-diffusion
        // dynamics everywhere, so computeDensity() is the sum of the populations
        // at every voxel including the solids, which is what we store into.
        rateLat.assign((size_t) rateOut.nch(), substrLattice);
        // Layout the sampler expects: [C.., B.., rate.., mask]
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_rate_lattices.push_back(&vec_substr_lattices[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesKinetics(reaction_type[iM])) continue;
            if (bmass_type[iM]==1) ptr_rate_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
            else                   ptr_rate_lattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
        }
        for (plint k = 0; k < rateOut.nch(); ++k) ptr_rate_lattices.push_back(&rateLat[(size_t) k]);
        ptr_rate_lattices.push_back(&maskLattice);
        writeRateChannelList(rateOut, str_outputDir);
    }
    else if (rate_enabled) {
        pcout << "  [RATE] <save_reaction_rates> is on, but no reaction block is active. "
              << "Nothing to write.\n";
    }

    // [FIX] defineKinetics.hh receives the biomass vector COMPACTED over the microbes
    //   that use kinetics, not the full list.  That was already true before the
    //   metabolic layer existed, but it only became REACHABLE once a microbe could be
    //   something other than kinetics.  The shipped defineKinetics.hh indexes
    //   B[0..11] by global microbe number and returns early if B.size() < 12, so
    //   mixing one FBA microbe into an otherwise all-kinetics run would silently
    //   switch every rate law off.  Say so rather than let it pass.
    if (kns_count > 0 && kns_count != num_of_microbes) {
        pcout << "\n  WARNING: " << kns_count << " of " << num_of_microbes
              << " microbes use kinetics, so defineRxnKinetics() receives a biomass\n"
              << "           vector of length " << kns_count << ", compacted over those microbes only --\n"
              << "           B[0] is the FIRST KINETICS microbe, not microbe0.  If defineKinetics.hh\n"
              << "           indexes B[] by the microbe numbers in CompLaB.xml (the shipped one does),\n"
              << "           either adjust it or give every microbe a kinetics reaction_type.\n\n";
    }

    // ---- OPTIONAL METABOLIC LAYER: its two lattice vectors --------------------
    //   Same layout as the kinetics vector, but each holds only the microbes that
    //   use that particular solver.  mm_globalId / srg_globalId translate a
    //   position in the compacted biomass block back to the microbe's number in
    //   CompLaB.xml, so per-microbe parameters are never indexed by position.
    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_mm_lattices, ptr_srg_lattices,
                                                  ptr_sym_lattices, ptr_gnn_lattices;
    std::vector<plint> mm_globalId, srg_globalId, mm_modelSlot;
    {
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_mm_lattices.push_back(&vec_substr_lattices[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesFBA(reaction_type[iM])) continue;
            mm_globalId.push_back(iM);
            if (bmass_type[iM]==1) ptr_mm_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
            else                   ptr_mm_lattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
        }
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_mm_lattices.push_back(&dC[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesFBA(reaction_type[iM])) continue;
            if (bmass_type[iM]==1) ptr_mm_lattices.push_back(&dBf[loctrack[iM]]);
            else                   ptr_mm_lattices.push_back(&dBp[loctrack[iM]]);
        }
        ptr_mm_lattices.push_back(&maskLattice);

        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_srg_lattices.push_back(&vec_substr_lattices[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesSurrogate(reaction_type[iM])) continue;
            srg_globalId.push_back(iM);
            if (bmass_type[iM]==1) ptr_srg_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
            else                   ptr_srg_lattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
        }
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_srg_lattices.push_back(&dC[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesSurrogate(reaction_type[iM])) continue;
            if (bmass_type[iM]==1) ptr_srg_lattices.push_back(&dBf[loctrack[iM]]);
            else                   ptr_srg_lattices.push_back(&dBp[loctrack[iM]]);
        }
        ptr_srg_lattices.push_back(&maskLattice);

        /* [SYM] same layout as the surrogate list above: substrates, biomass of the organisms
         * on this path, then dC, then dB, then the mask. */
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_sym_lattices.push_back(&vec_substr_lattices[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesSymbolic(reaction_type[iM])) continue;
            sym_globalId.push_back(iM);
            if (bmass_type[iM]==1) ptr_sym_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
            else                   ptr_sym_lattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
        }
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_sym_lattices.push_back(&dC[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesSymbolic(reaction_type[iM])) continue;
            if (bmass_type[iM]==1) ptr_sym_lattices.push_back(&dBf[loctrack[iM]]);
            else                   ptr_sym_lattices.push_back(&dBp[loctrack[iM]]);
        }
        ptr_sym_lattices.push_back(&maskLattice);

        /* [GNN] identical layout again: substrates, biomass of the organisms on this path,
         * then dC, then dB, then the mask. */
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_gnn_lattices.push_back(&vec_substr_lattices[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesGraphnet(reaction_type[iM])) continue;
            gnn_globalId.push_back(iM);
            if (bmass_type[iM]==1) ptr_gnn_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
            else                   ptr_gnn_lattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
        }
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_gnn_lattices.push_back(&dC[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesGraphnet(reaction_type[iM])) continue;
            if (bmass_type[iM]==1) ptr_gnn_lattices.push_back(&dBf[loctrack[iM]]);
            else                   ptr_gnn_lattices.push_back(&dBp[loctrack[iM]]);
        }
        ptr_gnn_lattices.push_back(&maskLattice);

        // cfg.vec_model is packed over COBRApy microbes only, in ascending
        // global order -- that is the order prep_cobrapy() was handed.  This map
        // turns a global microbe id into its position in that packed list, so
        // the processor never has to recompute it per voxel.
        mm_modelSlot.assign(num_of_microbes, 0);
        {
            plint slot = 0;
            for (plint iM = 0; iM < num_of_microbes; ++iM) {
                mm_modelSlot[iM] = slot;
                if (rxntype::usesCobrapy(reaction_type[iM])) ++slot;
            }
        }

        // Length assertions.  The layout IS the interface here, so check it once
        // at start-up rather than discovering a mis-sized vector as garbage
        // physics a thousand iterations in.
        if ((plint) ptr_mm_lattices.size()  != 2*(num_of_substrates + (plint) mm_globalId.size())  + 1 ||
            (plint) ptr_srg_lattices.size() != 2*(num_of_substrates + (plint) srg_globalId.size()) + 1 ||
            (plint) ptr_sym_lattices.size() != 2*(num_of_substrates + (plint) sym_globalId.size()) + 1 ||
            (plint) ptr_gnn_lattices.size() != 2*(num_of_substrates + (plint) gnn_globalId.size()) + 1) {
            pcout << "  [ERROR] metabolic lattice vector has the wrong length. This is a bug; please report it.\n";
            return -1;
        }
    }

    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_update_rxnLattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) { ptr_update_rxnLattices.push_back(&vec_substr_lattices[iS]); }
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (bmass_type[iM]==1) ptr_update_rxnLattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
        else ptr_update_rxnLattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
    }
    for (plint iS = 0; iS < num_of_substrates; ++iS) { ptr_update_rxnLattices.push_back(&dC[iS]); }
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (bmass_type[iM]==1) ptr_update_rxnLattices.push_back(&dBf[loctrack[iM]]);
        else ptr_update_rxnLattices.push_back(&dBp[loctrack[iM]]);
    }
    ptr_update_rxnLattices.push_back(&maskLattice);

    // Abiotic kinetics lattices (substrates only, no biomass)
    // Order: [C0, C1, ..., dC0, dC1, ..., mask]
    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_abiotic_kns_lattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) { ptr_abiotic_kns_lattices.push_back(&vec_substr_lattices[iS]); }
    for (plint iS = 0; iS < num_of_substrates; ++iS) { ptr_abiotic_kns_lattices.push_back(&dC[iS]); }
    ptr_abiotic_kns_lattices.push_back(&maskLattice);

    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_ca_lattices;
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (solver_type[iM]==2) {
            if (bmass_type[iM]==1) ptr_ca_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
            else { pcout << "  [CA] ERROR: CA only for biofilm\n"; return -1; }
        }
    }
    for (plint iM = 0; iM < num_of_microbes; ++iM) { if (solver_type[iM]==2) ptr_ca_lattices.push_back(&vec_bFcopy_lattices[loctrack[iM]]); }
    ptr_ca_lattices.push_back(&totalbFilmLattice);
    ptr_ca_lattices.push_back(&maskLattice);
    ptr_ca_lattices.push_back(&ageLattice);
    plint caLlen = ptr_ca_lattices.size();

    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_fd_lattices;
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (solver_type[iM]==1) {
            if (bmass_type[iM]==1) ptr_fd_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
            else ptr_fd_lattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
        }
    }
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (solver_type[iM]==1) {
            if (bmass_type[iM]==1) ptr_fd_lattices.push_back(&vec_bFcopy_lattices[loctrack[iM]]);
            else ptr_fd_lattices.push_back(&vec_bPcopy_lattices[loctrack[iM]]);
        }
    }
    ptr_fd_lattices.push_back(&maskLattice);
    plint fdLlen = ptr_fd_lattices.size();

    // [FIX] Same FD microbes, but in the layout updateLocalMaskNtotalLattices3D
    //   actually expects: [B.., copies.., totalBiomass, mask, age].
    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_fd_mask;
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (solver_type[iM]!=1) continue;
        if (bmass_type[iM]==1) ptr_fd_mask.push_back(&vec_bFilm_lattices[loctrack[iM]]);
        else                   ptr_fd_mask.push_back(&vec_bFree_lattices[loctrack[iM]]);
    }
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (solver_type[iM]!=1) continue;
        if (bmass_type[iM]==1) ptr_fd_mask.push_back(&vec_bFcopy_lattices[loctrack[iM]]);
        else                   ptr_fd_mask.push_back(&vec_bPcopy_lattices[loctrack[iM]]);
    }
    ptr_fd_mask.push_back(&totalbFilmLattice);
    ptr_fd_mask.push_back(&maskLattice);
    ptr_fd_mask.push_back(&ageLattice);
    plint fdMaskLen = ptr_fd_mask.size();

    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_eq_lattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) { ptr_eq_lattices.push_back(&vec_substr_lattices[iS]); }
    ptr_eq_lattices.push_back(&maskLattice);

    std::vector< MultiBlockLattice3D<T, RXNDES>* > ageNdistance_lattices;
    ageNdistance_lattices.push_back(&ageLattice);
    ageNdistance_lattices.push_back(&distLattice);
    ageNdistance_lattices.push_back(&totalbFilmLattice);

    // Initial mask update
    if (track_performance == 1) { global::timer("NS").restart(); }
    plint old_totMask = util::roundToInt(computeAverage(*computeDensity(maskLattice))*nx*ny*nz);
    if (bfilm_count > 0) {
        applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, caLlen, bounce_back, no_dynamics, bio_dynamics, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
    }
    plint new_totMask = util::roundToInt(computeAverage(*computeDensity(maskLattice))*nx*ny*nz);
    if (std::abs(old_totMask-new_totMask)>0) {
        old_totMask = new_totMask;
        if (soluteDindex == 1) applyProcessingFunctional(new updateSoluteDynamics3D<T,RXNDES>(num_of_substrates, bounce_back, no_dynamics, pore_dynamics, substrOMEGAinbFilm, substrOMEGAinPore), vec_substr_lattices[0].getBoundingBox(), substrate_lattices);
        if (bmassDindex == 1) applyProcessingFunctional(new updateBiomassDynamics3D<T,RXNDES>((plint)vec_bFree_lattices.size(), bounce_back, no_dynamics, pore_dynamics, bioOMEGAinbFilm, bioOMEGAinPore), vec_bFree_lattices[0].getBoundingBox(), planktonic_lattices);
        applyProcessingFunctional(new updateNsLatticesDynamics3D<T,NSDES,T,RXNDES>(nsLatticeOmega, vec_permRatio[0], pore_dynamics, no_dynamics, bounce_back), nsLattice.getBoundingBox(), nsLattice, maskLattice);
        pcout << "  [ADE] Biofilm mask changed geometry -> re-solving flow (up to " << ns_maxiTer_1 << " steps)...\n";
        for (plint iT2 = 0; iT2 < ns_maxiTer_1; ++iT2) {
            nsLattice.collideAndStream();
            ns_convg1.takeValue(getStoredAverageEnergy(nsLattice),false);
            if (ns_convg1.hasConverged()) { pcout << "  [ADE] flow re-solve converged at " << iT2 << "\n"; break; }
            if (iT2 % 5000 == 0 && iT2 > 0) pcout << "  [ADE] flow re-solve ... " << iT2 << "/" << ns_maxiTer_1 << "\n";
        }
    }
    if (read_NS_file==0 || (read_NS_file==1 && ns_rerun_iT0>0)) {
        if (track_performance == 0) {
            writeNsVTI(nsLattice,ns_maxiTer_1,"nsLatticeFinal1_");
            saveBinaryBlock(nsLattice, str_outputDir+ns_filename+".chk");
        }
    }
    if (track_performance == 1) { nstime += global::timer("NS").getTime(); global::timer("NS").stop(); }

    // Couple NS and ADE
    if (Pe > thrd) {
        pcout << "  [ADE] Coupling NS-ADE lattices...\n";
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            if (vec_immobile[iS]) continue;   // [immobile] no advection coupling
            latticeToPassiveAdvDiff(nsLattice, vec_substr_lattices[iS], vec_substr_lattices[iS].getBoundingBox());
        }
        tmpIT0=0;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (solver_type[iM] == 3) {
                latticeToPassiveAdvDiff(nsLattice, vec_bFree_lattices[tmpIT0], vec_bFree_lattices[tmpIT0].getBoundingBox());
                ++tmpIT0;
            }
        }
        pcout << "  [ADE] Stabilizing (10000 iter)...\n";
        for (plint iT=0; iT<10000; ++iT) {
            if (iT % 1000 == 0) pcout << "  [ADE] stabilizing ... " << iT << "/10000  (" << global::timer("total").getTime() << " s elapsed)\n";
            for (plint iS = 0; iS < num_of_substrates; ++iS) if (!vec_immobile[iS]) vec_substr_lattices[iS].collideAndStream();
            for (size_t iM = 0; iM < vec_bFree_lattices.size(); ++iM) vec_bFree_lattices[iM].collideAndStream();
        }
        pcout << "  [ADE] Stabilization done.\n";
        for (plint iS = 0; iS < num_of_substrates; ++iS) if (!vec_immobile[iS]) applyProcessingFunctional(new stabilizeADElattice3D<T,RXNDES,int>(vec_c0[iS], pore_dynamics, bio_dynamics), vec_substr_lattices[iS].getBoundingBox(), vec_substr_lattices[iS], geometry);
        for (size_t iM = 0; iM < vec_bFree_lattices.size(); ++iM) applyProcessingFunctional(new stabilizeADElattice3D<T,RXNDES,int>(vec_b0_free[iM], pore_dynamics, bio_dynamics), vec_bFree_lattices[iM].getBoundingBox(), vec_bFree_lattices[iM], geometry);
    }

    // Load checkpoints if needed
    iT = 0;
    if (read_ADE_file==1 && ade_rerun_iT0>0) {
        pcout << "  [ADE] Loading checkpoints...\n";
        for (plint iS = 0; iS < num_of_substrates; ++iS) loadBinaryBlock(vec_substr_lattices[iS], str_outputDir+ade_filename+"_"+std::to_string(iS));
        tmpIT0=0; tmpIT1=0;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (bmass_type[iM]==1) { loadBinaryBlock(vec_bFilm_lattices[tmpIT0], str_outputDir+bio_filename+"_"+std::to_string(iM)); ++tmpIT0; }
            else { loadBinaryBlock(vec_bFree_lattices[tmpIT1], str_outputDir+bio_filename+"_"+std::to_string(iM)); ++tmpIT1; }
        }
        iT = ade_rerun_iT0;
    }
    T catime = 0, adetime = 0, knstime = 0, cnstime = 0;
    T fbatime = 0, srgtime = 0, symtime = 0, gnntime = 0;   // optional metabolic layer

    // ════════════════════════════════════════════════════════════════════════════
    // PHASE 4: MAIN SIMULATION LOOP
    // ════════════════════════════════════════════════════════════════════════════
    pcout << "\n┌────────────────────────────────────────────────────────────────────────┐\n";
    pcout << "│ PHASE 4: MAIN SIMULATION LOOP                                         │\n";
    pcout << "├────────────────────────────────────────────────────────────────────────┤\n";
    pcout << "│ Max iterations: " << ade_maxiTer << "   VTI interval: " << ade_VTI_iTer << "\n";
    pcout << "│ Output files will use SPECIES NAMES from XML:\n";
    for (plint iS = 0; iS < num_of_substrates; ++iS) pcout << "│   " << vec_subs_names[iS] << "_*.vti\n";
    for (plint iM = 0; iM < num_of_microbes; ++iM) pcout << "│   " << vec_microbes_names[iM] << "_*.vti\n";
    pcout << "└────────────────────────────────────────────────────────────────────────┘\n\n";

    global::timer("ade").restart();
    util::ValueTracer<T> ns_convg2(1.0,1000.0,ns_converge_iT2);
    bool ns_saturate=0, percolationFlag=0;

    for (; iT < ade_maxiTer; ++iT) {
        /* [HEARTBEAT 2026-07-24] lightweight liveness line between the (every-VTI) ITERATION blocks,
         * explicitly flushed so you can see the reactive loop is advancing even with buffered stdout. */
        if (iT > 0 && iT % 50 == 0 && (ade_VTI_iTer <= 0 || iT % ade_VTI_iTer != 0))
            pcout << "  [ade] iter " << iT << "/" << ade_maxiTer << "  (" << global::timer("total").getTime() << " s elapsed)\n";
        // ════════════════════════════════════════════════════════════════════════
        // VTI OUTPUT AND DIAGNOSTICS
        // ════════════════════════════════════════════════════════════════════════
        if (ade_VTI_iTer > 0 && iT % ade_VTI_iTer == 0) {
            pcout << "\n╔════════════════════════════════════════════════════════════════════════╗\n";
            pcout << "║ ITERATION " << iT << "  |  Time: " << std::scientific << iT*ade_dt << " s" << std::fixed << "\n";
            pcout << "╠════════════════════════════════════════════════════════════════════════╣\n";
            
            // Substrate status
            pcout << "║ SUBSTRATES:\n";
            for (plint iS = 0; iS < num_of_substrates; ++iS) {
                T sMin = computeMin(*computeDensity(vec_substr_lattices[iS]));
                T sMax = computeMax(*computeDensity(vec_substr_lattices[iS]));
                T sAvg = computeAverage(*computeDensity(vec_substr_lattices[iS]));
                pcout << "║   " << vec_subs_names[iS] << ": min=" << std::scientific << sMin 
                      << " avg=" << sAvg << " max=" << sMax;
                if (sMin < 0) pcout << " [NEG!]";
                pcout << std::fixed << "\n";
            }
            
            // Biomass status
            if (bfilm_count > 0) {
                pcout << "║ BIOMASS:\n";
                for (plint iM = 0; iM < bfilm_count; ++iM) {
                    T bMin = computeMin(*computeDensity(vec_bFilm_lattices[iM]));
                    T bMax = computeMax(*computeDensity(vec_bFilm_lattices[iM]));
                    T bAvg = computeAverage(*computeDensity(vec_bFilm_lattices[iM]));
                    T growth = (diag_initial_biomass > 0) ? ((bMax - diag_initial_biomass) / diag_initial_biomass * 100.0) : 0.0;
                    pcout << "║   " << vec_microbes_names[iM] << ": max=" << std::scientific << bMax 
                          << "/" << max_bMassRho << std::fixed << " (" << growth << "% growth)";
                    if (bMax > max_bMassRho) pcout << " [>Bmax!]";
                    pcout << "\n";
                }
                pcout << "║   CA: triggers=" << diag_ca_triggers << " redistributions=" << diag_ca_redistributions << "\n";
            }
            
            // ════════════════════════════════════════════════════════════════════
            // KINETICS DEBUG STATS (efficient - only summary)
            // ════════════════════════════════════════════════════════════════════
            if (kns_count > 0) {
                long cells_bio, cells_grow;
                double sum_dB, max_B, max_dB, min_DOC;
                KineticsStats::getStats(cells_bio, cells_grow, sum_dB, max_B, max_dB, min_DOC);
                pcout << "║ KINETICS (last " << ade_VTI_iTer << " iters):\n";
                pcout << "║   Active cells: " << cells_bio << " (growing: " << cells_grow << ")\n";
                if (cells_bio > 0) {
                    pcout << "║   Sum dB/dt: " << std::scientific << sum_dB << " kg/m³/s\n";
                    pcout << "║   Max dB/dt: " << max_dB << " kg/m³/s\n";
                    pcout << "║   Min DOC in biofilm: " << min_DOC << " mol/L\n";
                    // Estimate iterations to Bmax
                    if (max_dB > 0 && max_B < max_bMassRho) {
                        double time_to_bmax = (max_bMassRho - max_B) / max_dB;
                        double iters_to_bmax = time_to_bmax / ade_dt;
                        pcout << "║   Est. iters to Bmax: " << std::fixed << (long)iters_to_bmax << "\n";
                    }
                    pcout << std::fixed;
                }
                // Reset stats for next interval
                KineticsStats::resetIteration();
            }
            
            pcout << "╚════════════════════════════════════════════════════════════════════════╝\n";
            
            // Write VTI with SPECIES NAMES
            if (track_performance == 0) {
                for (plint iS = 0; iS < num_of_substrates; ++iS) {
                    writeAdvVTI(vec_substr_lattices[iS], iT, vec_subs_names[iS]+"_");
                }
                tmpIT0=0; tmpIT1=0;
                for (plint iM = 0; iM < num_of_microbes; ++iM) {
                    if (bmass_type[iM]==1) { writeAdvVTI(vec_bFilm_lattices[tmpIT0], iT, vec_microbes_names[iM]+"_"); ++tmpIT0; }
                    else { writeAdvVTI(vec_bFree_lattices[tmpIT1], iT, vec_microbes_names[iM]+"_"); ++tmpIT1; }
                }
                if (Pe > thrd) writeNsVTI(nsLattice, iT, "nsLattice_");

                // [RATE-OUT] Sample the reaction rates, then write them.
                //   The sampler calls the SAME two kinetics functions the solver
                //   is about to call further down this same iteration, on the
                //   SAME concentrations (collision conserves density, so nothing
                //   in between changes them). These are the rates the solver
                //   used, not an estimate reconstructed from the output.
                //   Calling them twice also counts them twice in KineticsStats,
                //   so the counters are cleared again here. They were reset a few
                //   lines above, so nothing real is thrown away.
                if (rateOut.enabled && rateOut.nch() > 0 && iT % rateOut.interval == 0) {
                    applyProcessingFunctional(
                        new sample_reaction_rates<T,RXNDES>(
                            nx, ny, nz, num_of_substrates, kns_count, rateOut.nch(),
                            no_dynamics, bounce_back,
                            rateOut.doBio, rateOut.doAbio, rateOut.surfaceOnly,
                            rateOut.nrb, rateOut.nra,
                            rateOut.offBioC(), rateOut.offBioM(), rateOut.offBioR(),
                            rateOut.offAbioC(), rateOut.offAbioR()),
                        vec_substr_lattices[0].getBoundingBox(), ptr_rate_lattices);
                    KineticsStats::resetIteration();
                    writeRateVTI(rateLat, rateOut, iT, nx, ny, nz);
                }
            }
            adetime += global::timer("ade").getTime();
            pcout << "  Wall clock: " << global::timer("ade").getTime() << " s\n";
            global::timer("ade").restart();
        }
        
        // Checkpoint
        if (ade_CHK_iTer > 0 && iT % ade_CHK_iTer == 0 && iT > 0 && track_performance == 0) {
            pcout << "  [CHK] Saving checkpoint at iter=" << iT << "\n";
            for (plint iS = 0; iS < num_of_substrates; ++iS)
                saveBinaryBlock(vec_substr_lattices[iS], str_outputDir + ade_filename + std::to_string(iS) + "_" + std::to_string(iT) + ".chk");
            tmpIT0=0; tmpIT1=0;
            for (plint iM = 0; iM < num_of_microbes; ++iM) {
                if (bmass_type[iM]==1) { saveBinaryBlock(vec_bFilm_lattices[tmpIT0], str_outputDir + bio_filename + std::to_string(iM) + "_" + std::to_string(iT) + ".chk"); ++tmpIT0; }
                else { saveBinaryBlock(vec_bFree_lattices[tmpIT1], str_outputDir + bio_filename + std::to_string(iM) + "_" + std::to_string(iT) + ".chk"); ++tmpIT1; }
            }
        }

        if (track_performance == 1) global::timer("cns").restart();
        
        // Collision
        for (plint iS = 0; iS < num_of_substrates; ++iS) if (!vec_immobile[iS]) vec_substr_lattices[iS].collide();
        if (lb_count > 0) {
            for (plint iM = 0; iM < num_of_microbes; ++iM) {
                if (solver_type[iM]==3) {
                    if (bmass_type[iM]==1) vec_bFilm_lattices[loctrack[iM]].collide();
                    else vec_bFree_lattices[loctrack[iM]].collide();
                }
            }
        }
        if (track_performance == 1) { cnstime += global::timer("cns").getTime(); global::timer("cns").stop(); }

        // Kinetics (biotic - only if enable_kinetics is true and biotic_mode)
        dC=dC0; dBp=dBp0; dBf=dBf0;
        if (enable_kinetics && kns_count > 0) {
            if (track_performance == 1) global::timer("kns").restart();
            applyProcessingFunctional(new run_kinetics<T,RXNDES>(nx, num_of_substrates, kns_count, ade_dt, vec_Kc_kns, vec_mu_kns, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_kns_lattices);
            if (track_performance == 1) { knstime+=global::timer("kns").getTime(); global::timer("kns").stop(); }
        }
        // ---- OPTIONAL METABOLIC LAYER -------------------------------------------
        //   These run BEFORE update_rxnLattices, into the same dC/dB increment
        //   lattices the kinetics term used, so a microbe with a combined
        //   reaction_type (glpk_and_kinetics, surrogate_and_kinetics, ...) gets
        //   both contributions added and applied together, exactly once.
#ifdef COMPLAB_ENABLE_GLPK
        if (mmcfg.enable_fba_glpk && mmcfg.glpk_count > 0) {
            if (track_performance == 1) global::timer("fba").restart();
            applyProcessingFunctional(new runFBA_glpk3D<T,RXNDES>(nx, num_of_substrates, (plint) mm_globalId.size(),
                                          ade_dt, mm_globalId, &mmcfg, vec_Kc, vec_mu, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_mm_lattices);
            if (track_performance == 1) { fbatime += global::timer("fba").getTime(); global::timer("fba").stop(); }
        }
#endif
#ifdef COMPLAB_ENABLE_COBRAPY
        if (mmcfg.enable_fba_cobrapy && mmcfg.cpy_count > 0) {
            if (track_performance == 1) global::timer("fba").restart();
            applyProcessingFunctional(new runFBA_cobrapy3D<T,RXNDES>(nx, num_of_substrates, (plint) mm_globalId.size(),
                                          ade_dt, mm_globalId, &mmcfg, pyFileName, mm_modelSlot, vec_Kc, vec_mu, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_mm_lattices);
            if (track_performance == 1) { fbatime += global::timer("fba").getTime(); global::timer("fba").stop(); }
        }
#endif
        if (mmcfg.enable_surrogate && mmcfg.srg_count > 0) {
            if (track_performance == 1) global::timer("srg").restart();
            applyProcessingFunctional(new run_surrogate3D<T,RXNDES>(nx, num_of_substrates, (plint) srg_globalId.size(),
                                          num_of_microbes, ade_dt, srg_globalId, &mmcfg, vec_Kc, vec_mu,
                                          no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_srg_lattices);
            if (track_performance == 1) { srgtime += global::timer("srg").getTime(); global::timer("srg").stop(); }
        }

        if (sym_count > 0) {
            if (track_performance == 1) global::timer("sym").restart();
            applyProcessingFunctional(new run_symbolic3D<T,RXNDES>(nx, num_of_substrates,
                                          (plint) sym_globalId.size(), num_of_microbes, ade_dt,
                                          sym_globalId, &mmcfg, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_sym_lattices);
            if (track_performance == 1) { symtime += global::timer("sym").getTime(); global::timer("sym").stop(); }
        }

        if (gnn_count > 0) {
            if (track_performance == 1) global::timer("gnn").restart();
            applyProcessingFunctional(new run_graphnet3D<T,RXNDES>(nx, num_of_substrates,
                                          (plint) gnn_globalId.size(), num_of_microbes, ade_dt,
                                          gnn_globalId, &mmcfg, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_gnn_lattices);
            if (track_performance == 1) { gnntime += global::timer("gnn").getTime(); global::timer("gnn").stop(); }
        }

        // [FIX] Gate on the actual COUNTS, not on the XML switches.  anyEnabled()
        //   is true as soon as a switch is set, so <enable_surrogate>true</> with no
        //   surrogate microbe swept the whole domain over 2*(S+M)+1 lattices every
        //   step applying identically zero increments.
        if ((enable_kinetics && kns_count > 0) ||
            (mmcfg.enable_fba_glpk    && mmcfg.glpk_count > 0) ||
            (mmcfg.enable_fba_cobrapy && mmcfg.cpy_count  > 0) ||
            (mmcfg.enable_surrogate   && mmcfg.srg_count  > 0) ||
            (sym_count > 0) || (gnn_count > 0)) {
            if (track_performance == 1) global::timer("rxn").restart();
            applyProcessingFunctional(new update_rxnLattices<T,RXNDES>(nx, num_of_substrates, num_of_microbes, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_update_rxnLattices);
            if (track_performance == 1) { T rxntime=global::timer("rxn").getTime(); global::timer("rxn").stop(); if (kns_count>0) knstime+=rxntime; }
        }

        /* [SYM-ABIOTIC][GNN-ABIOTIC] The learned abiotic laws, on the same lattice vector the
         * hand-written abiotic kinetics uses.  They accumulate into the same dC[], so
         * update_abiotic_rxnLattices below applies all of it in one pass. */
        if (haveSymAbiotic) {
            if (track_performance == 1) global::timer("sym").restart();
            applyProcessingFunctional(new run_symbolic_abiotic3D<T,RXNDES>(nx, num_of_substrates,
                                          ade_dt, &mmcfg, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_abiotic_kns_lattices);
            if (track_performance == 1) { symtime += global::timer("sym").getTime(); global::timer("sym").stop(); }
        }
        if (haveGnnAbiotic) {
            if (track_performance == 1) global::timer("gnn").restart();
            applyProcessingFunctional(new run_graphnet_abiotic3D<T,RXNDES>(nx, num_of_substrates,
                                          ade_dt, &mmcfg, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_abiotic_kns_lattices);
            if (track_performance == 1) { gnntime += global::timer("gnn").getTime(); global::timer("gnn").stop(); }
        }
        if ((haveSymAbiotic || haveGnnAbiotic) && !(enable_abiotic_kinetics && num_of_substrates > 0)
            && !dissolCfg.enabled) {
            /* Nothing else is going to apply dC[] on the abiotic path this step, so do it here.
             * When the hand-written abiotic block or dissolution IS on, they already call this
             * once and calling it twice would apply the increments twice. */
            applyProcessingFunctional(new update_abiotic_rxnLattices<T,RXNDES>(nx, num_of_substrates, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_abiotic_kns_lattices);
        }

        // Abiotic kinetics (substrate-only reactions without microbes)
        if (enable_abiotic_kinetics && num_of_substrates > 0) {
            if (track_performance == 1) global::timer("abiotic_kns").restart();
            // [FIX] The abiotic block accumulates into the SAME dC[] lattices the
            //   biotic block just used, and those were reset only once, above the
            //   biotic block.  Without this reset the biotic increments are still
            //   sitting in dC[] when run_abiotic_kinetics adds to them, and
            //   update_abiotic_rxnLattices then applies the biotic contribution a
            //   SECOND time -- silently doubling every biotic reaction in any run
            //   with both enable_kinetics and enable_abiotic_kinetics on, which is
            //   the shipped configuration.  Reset first.
            dC = dC0;
            // Calculate abiotic reaction rates
            // [PRECIP-VOP] surface-gate the reaction so A->P fires only on interface (wall-adjacent) voxels
            if (precip_enabled && precip_surfaceOnly) {
                applyProcessingFunctional(new surfaceAbioticKinetics3D<T,RXNDES>(nx, ny, nz, num_of_substrates, ade_dt, no_dynamics, bounce_back, 1),
                                          vec_substr_lattices[0].getBoundingBox(), ptr_abiotic_kns_lattices);
            } else {
                applyProcessingFunctional(new run_abiotic_kinetics<T,RXNDES>(nx, num_of_substrates, ade_dt, no_dynamics, bounce_back),
                                          vec_substr_lattices[0].getBoundingBox(), ptr_abiotic_kns_lattices);
            }
            // [DISSOL-VOP] Mineral dissolution, into the SAME dC lattices, before they
            //   are applied.  It runs here rather than in its own block so that a
            //   substrate both produced by dissolution and consumed by an abiotic
            //   reaction gets one combined increment, applied once.
            if (dissolCfg.enabled) {
                applyProcessingFunctional(new surfaceDissolutionKinetics3D<T,RXNDES>(
                                              nx, ny, nz, num_of_substrates, ade_dt,
                                              no_dynamics, bounce_back, pore_dynamics, &dissolCfg),
                                          vec_substr_lattices[0].getBoundingBox(), ptr_dissol);
            }

            // Apply concentration changes
            applyProcessingFunctional(new update_abiotic_rxnLattices<T,RXNDES>(nx, num_of_substrates, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_abiotic_kns_lattices);
            if (track_performance == 1) { knstime += global::timer("abiotic_kns").getTime(); global::timer("abiotic_kns").stop(); }
        }

        // [DISSOL-VOP] If the abiotic block above is switched off, dissolution still
        //   has to run somewhere.  Give it its own reset/apply pair so it works
        //   independently of <enable_abiotic_kinetics>.
        if (dissolCfg.enabled && !(enable_abiotic_kinetics && num_of_substrates > 0)) {
            dC = dC0;
            applyProcessingFunctional(new surfaceDissolutionKinetics3D<T,RXNDES>(
                                          nx, ny, nz, num_of_substrates, ade_dt,
                                          no_dynamics, bounce_back, pore_dynamics, &dissolCfg),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_dissol);
            applyProcessingFunctional(new update_abiotic_rxnLattices<T,RXNDES>(nx, num_of_substrates, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_abiotic_kns_lattices);
        }

        // Equilibrium chemistry (runs regardless of enable_kinetics - controlled separately)
        if (useEquilibrium) {
            if (track_performance == 1) global::timer("eq").restart();
            applyProcessingFunctional(new run_equilibrium_biotic<T, RXNDES>(nx, num_of_substrates, eqSolver, no_dynamics, bounce_back),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_eq_lattices);
            if (track_performance == 1) { eqtime += global::timer("eq").getTime(); global::timer("eq").stop(); }
        }

        // ════════════════════════════════════════════════════════════════════════════
        // VALIDATION DIAGNOSTICS (per-iteration detailed output)
        // ════════════════════════════════════════════════════════════════════════════
        if (enable_validation_diagnostics && (iT % 100 == 0 || iT < 10)) {
            pcout << "\n┌─────────────────────────────────────────────────────────────────────────┐\n";
            pcout << "│ VALIDATION DIAGNOSTICS - Iteration " << iT << "                              │\n";
            pcout << "├─────────────────────────────────────────────────────────────────────────┤\n";
            pcout << "│ Time: " << std::scientific << std::setprecision(4) << iT*ade_dt << " s" << std::fixed << "\n";

            // Step-by-step data flow verification
            pcout << "├─────────────────────────────────────────────────────────────────────────┤\n";
            pcout << "│ STEP 6.1 [COLLISION]: LBM collision completed                           │\n";

            // Sample concentration at center of domain
            plint midX = nx/2, midY = ny/2, midZ = nz/2;
            pcout << "│ STEP 6.2 [KINETICS]: ";
            if (enable_kinetics && kns_count > 0) {
                pcout << "ACTIVE - " << kns_count << " reaction(s)\n";
                // Show sample values
                for (plint iS = 0; iS < std::min((plint)2, num_of_substrates); ++iS) {
                    T cMid = vec_substr_lattices[iS].get(midX, midY, midZ).computeDensity();
                    T dC_mid = dC[iS].get(midX, midY, midZ).computeDensity();
                    pcout << "│   " << vec_subs_names[iS] << " @center: C=" << std::scientific
                          << cMid << ", dC=" << dC_mid << std::fixed << "\n";
                }
                if (bfilm_count > 0) {
                    T bMid = vec_bFilm_lattices[0].get(midX, midY, midZ).computeDensity();
                    T dB_mid = dBf[0].get(midX, midY, midZ).computeDensity();
                    pcout << "│   Biomass @center: B=" << std::scientific << bMid
                          << ", dB=" << dB_mid << std::fixed << "\n";
                }
            } else {
                pcout << "DISABLED (enable_kinetics=" << enable_kinetics << ", kns_count=" << kns_count << ")\n";
            }

            pcout << "│ STEP 6.2b [ABIOTIC KINETICS]: ";
            if (enable_abiotic_kinetics) {
                pcout << "ACTIVE (substrate-only reactions)\n";
                for (plint iS = 0; iS < std::min((plint)2, num_of_substrates); ++iS) {
                    T cMid = vec_substr_lattices[iS].get(midX, midY, midZ).computeDensity();
                    pcout << "│   " << vec_subs_names[iS] << " @center: C=" << std::scientific
                          << cMid << std::fixed << "\n";
                }
            } else {
                pcout << "DISABLED\n";
            }

            pcout << "│ STEP 6.3 [EQUILIBRIUM]: ";
            if (useEquilibrium) {
                pcout << "ACTIVE\n";
                // Show sample equilibrium-adjusted values
                for (plint iS = 0; iS < std::min((plint)2, num_of_substrates); ++iS) {
                    T cMin = computeMin(*computeDensity(vec_substr_lattices[iS]));
                    T cMax = computeMax(*computeDensity(vec_substr_lattices[iS]));
                    pcout << "│   " << vec_subs_names[iS] << ": min=" << std::scientific
                          << cMin << ", max=" << cMax << std::fixed << "\n";
                }
            } else {
                pcout << "DISABLED\n";
            }

            // Mass balance check
            pcout << "├─────────────────────────────────────────────────────────────────────────┤\n";
            pcout << "│ MASS BALANCE CHECK:                                                     │\n";
            for (plint iS = 0; iS < std::min((plint)2, num_of_substrates); ++iS) {
                T totalMass = computeSum(*computeDensity(vec_substr_lattices[iS]));
                pcout << "│   " << vec_subs_names[iS] << " total: " << std::scientific << totalMass << std::fixed << "\n";
            }
            if (bfilm_count > 0) {
                T totalBiomass = computeSum(*computeDensity(totalbFilmLattice));
                pcout << "│   Total biomass: " << std::scientific << totalBiomass << std::fixed << "\n";
            }

            pcout << "└─────────────────────────────────────────────────────────────────────────┘\n";
        }

        // CA biomass expansion
        if (ca_count > 0) {
            applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, caLlen, bounce_back, no_dynamics, bio_dynamics, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
            T globalBmax = computeMax(*computeDensity(totalbFilmLattice));
            if (std::isnan(globalBmax) || std::isinf(globalBmax)) { pcout << "\n  [CA] ERROR: non-finite biomass (NaN/Inf) at iter=" << iT << " -- stopping cleanly\n"; percolationFlag = 1; }
            // [CA-FAIL] 2D-style hard stop: on an unresolvable CA state, dump every field and terminate with an explicit reason report.
            auto dumpAllFields = [&](const std::string& reason, plint pushSweeps, plint ageSweeps) {
                T _bmax = computeMax(*computeDensity(totalbFilmLattice));
                T _bavg = computeAverage(*computeDensity(totalbFilmLattice));
                T _bsum = computeSum(*computeDensity(totalbFilmLattice));
                // ---- measure the true state at failure so the cause is DETERMINED, not guessed ----
                T _cs = std::sqrt(1.0/3.0);
                T _umax = 0.0, _MaNow = 0.0; bool _flowFinite = true;
                if (Pe > thrd) {
                    _umax = computeMax(*computeVelocityNorm(nsLattice, Box3D(1,nx-2,0,ny-1,0,nz-1)));
                    _MaNow = _umax / _cs;
                    if (!std::isfinite(_umax)) _flowFinite = false;
                }
                T _worstMin = 0.0; std::string _worstSp = "(none)"; bool _chemFinite = true; plint _worstIdx = -1;
                for (plint iS = 0; iS < num_of_substrates; ++iS) {
                    T _mn = computeMin(*computeDensity(vec_substr_lattices[iS]));
                    if (!std::isfinite(_mn)) _chemFinite = false;
                    if (_mn < _worstMin) { _worstMin = _mn; _worstSp = vec_subs_names[iS]; _worstIdx = iS; }
                }
                bool _bmassFinite = std::isfinite(_bmax);
                T _capRatio = (max_bMassRho > 0.0) ? (_bmax / max_bMassRho) : 0.0;
                // ---- decide the ONE root cause. Ordering matters: a diverged flow corrupts advection,
                // which then poisons chemistry and biomass, so flow-divergence is checked first as the upstream root. ----
                std::string _CAT, _WHY1, _WHY2, _WHY3, _SOL1, _SOL2, _PLAIN, _PFIX;
                if (!_flowFinite || _MaNow > 0.3) {
                    _CAT  = "FLOW / PRESSURE-DRIVE  (LBM Mach instability from deltaP + geometry)";
                    _WHY1 = "The Navier-Stokes flow field went unstable: peak velocity pushes the Mach number past the LBM limit.";
                    _WHY2 = "LBM is only stable for Ma < ~0.1 and blows up past ~0.3. The pressure drop deltaP drives fluid through";
                    _WHY3 = "the tight pore throats too fast; the biofilm mask-flip flow re-solve then amplified it into overflow.";
                    _SOL1 = "Lower the drive so Ma <= 0.02: scale deltaP_new = deltaP * (0.02 / Ma_now), or reduce the target Pe.";
                    _SOL2 = "Also helps: raise tau_NS toward 1.0 (more numerical viscosity), or widen throats (higher permeability).";
                    _PLAIN = "The water was pushed through the narrow pore channels faster than this simulation method can handle, so the flow calculation blew up: velocities jumped to impossible values and wrecked everything downstream.";
                    _PFIX  = "Push the water more gently - lower the pressure difference (or the target flow speed) until the flow is slow enough to stay stable, or use a geometry with wider channels.";
                } else if (!_bmassFinite || _capRatio > 10.0) {
                    _CAT  = "BIOMASS  (kinetics growth blow-up)";
                    _WHY1 = "Biomass integrated far beyond its physical carrying cap (max/cap ratio is shown in the evidence below).";
                    _WHY2 = "The flow is still finite, so this is a reaction-integration blow-up, not a hydraulic one: the growth";
                    _WHY3 = "increment per step is too large (per-day rate constants applied per-second ~ 86,400x), overshooting the cap.";
                    _SOL1 = "Shrink the reactive step: smaller dx, convert the per-day rate constants to per-second, or cap per-step growth.";
                    _SOL2 = "Confirm the f_cap carrying-capacity limiter is active for every biofilm pool in defineKinetics.hh.";
                    _PLAIN = "The microbes were told to grow by a huge amount in a single step, far more than is physically possible, so the biomass number exploded. This normally means the growth rates are applied too fast (daily rates used as if per-second).";
                    _PFIX  = "Slow the growth calculation down - use smaller time steps or convert the daily rate constants to per-second, and make sure the biomass-cap limiter is switched on.";
                } else if (!_chemFinite || _worstMin < -1.0e-4) {
                    _CAT  = "CHEMISTRY  (equilibrium speciation produced negative concentrations)";
                    _WHY1 = "Speciation produced unphysical negative mass (the most-negative species is named in the evidence below).";
                    _WHY2 = "The flow is finite and biomass is near cap, so the fault is chemical: the reaction/advection step removed";
                    _WHY3 = "more of a species than was locally present, or the bulk composition handed to the EQ solver is infeasible.";
                    _SOL1 = "Shrink the reactive step (smaller dx / lower rate constants) or tighten the equilibrium tolerance.";
                    _SOL2 = "Clamp solutes to >= 0 after the reaction step, and verify the initial composition is charge-balanced.";
                    _PLAIN = "A dissolved species (" + _worstSp + ") ended up with a NEGATIVE amount here, which is physically impossible. The proof section below pins down which step actually made it negative - the chemistry solver, a reaction, or the transport that carries species through the pores.";
                    _PFIX  = "Stop values dropping below zero right after each transport step, take smaller time steps so nothing overshoots, and keep concentration gradients gentle. The proof section below names the exact step to fix.";
                } else {
                    _CAT  = "GEOMETRY  (pore-throat clogging deadlock)";
                    _WHY1 = "Biofilm filled pore voxels to the cap, but the push/pull CA cannot move the excess anywhere: those";
                    _WHY2 = "voxels are boxed in by solid grain / bounce-back wall / other capped biofilm, so no open pore neighbour";
                    _WHY3 = "can receive it. Flow is finite and chemistry is physical, so the pore throats themselves are too narrow.";
                    _SOL1 = "Use a looser/coarser geometry with wider throats (fewer dead-end pores), or raise resolution (smaller dx)";
                    _SOL2 = "so the CA has more pore voxels to spread into; or lower growth so biofilm does not saturate the throats.";
                    _PLAIN = "The microbes filled a pore completely and had nowhere left to expand, because that pore is boxed in by solid grains on every side. The program kept trying to push the extra biomass somewhere and never could, so it stopped.";
                    _PFIX  = "Give the biomass more room - use a looser grain packing with wider gaps, use a finer grid so there are more pore cells, or slow the growth so pores do not fill so fast.";
                }
                pcout << "\n";
                // ---- pinpoint the exact voxel of the offending quantity (MPI-safe: each Box3D reduction is global) ----
                auto _locate = [&](MultiScalarField3D<T>& fld, bool findMax, plint& gx, plint& gy, plint& gz)->T {
                    T tgt = findMax ? computeMax(fld) : computeMin(fld);
                    gx=0; gy=0; gz=0;
                    for (plint x=0;x<nx;++x){ T v=findMax?computeMax(fld,Box3D(x,x,0,ny-1,0,nz-1)):computeMin(fld,Box3D(x,x,0,ny-1,0,nz-1)); if(v==tgt){gx=x;break;} }
                    for (plint y=0;y<ny;++y){ T v=findMax?computeMax(fld,Box3D(gx,gx,y,y,0,nz-1)):computeMin(fld,Box3D(gx,gx,y,y,0,nz-1)); if(v==tgt){gy=y;break;} }
                    for (plint z=0;z<nz;++z){ T v=findMax?computeMax(fld,Box3D(gx,gx,gy,gy,z,z)):computeMin(fld,Box3D(gx,gx,gy,gy,z,z)); if(v==tgt){gz=z;break;} }
                    return tgt;
                };
                auto _valAt = [&](const std::string& nm, plint x, plint y, plint z)->T {
                    for (plint iS=0;iS<num_of_substrates;++iS) if (vec_subs_names[iS]==nm) return computeMax(*computeDensity(vec_substr_lattices[iS]), Box3D(x,x,y,y,z,z));
                    return 0.0;
                };
                plint _gx=0,_gy=0,_gz=0; T _locVal=0.0; std::string _locWhat="(n/a)";
                if (_CAT[0]=='C') { auto _f=computeDensity(vec_substr_lattices[_worstIdx>=0?_worstIdx:0]); _locVal=_locate(*_f,false,_gx,_gy,_gz); _locWhat=_worstSp+" (most negative)"; }
                else if (_CAT[0]=='F') { auto _f=computeVelocityNorm(nsLattice); _locVal=_locate(*_f,true,_gx,_gy,_gz); _locWhat="peak velocity |u|"; }
                else { auto _f=computeDensity(totalbFilmLattice); _locVal=_locate(*_f,true,_gx,_gy,_gz); _locWhat="peak biofilm biomass"; }
                T _mcode = computeMax(*computeDensity(maskLattice), Box3D(_gx,_gx,_gy,_gy,_gz,_gz));
                T _bHere = computeMax(*computeDensity(totalbFilmLattice), Box3D(_gx,_gx,_gy,_gy,_gz,_gz));
                T _uHere = (Pe>thrd) ? computeMax(*computeVelocityNorm(nsLattice, Box3D(_gx,_gx,_gy,_gy,_gz,_gz))) : 0.0;
                std::string _face = (_gx<=2) ? "INLET face (x low)" : ((_gx>=nx-3) ? "OUTLET face (x high)" : "interior (mid-domain)");
                // ---- evidence for WHY a concentration is negative: reaction increment, pre-reaction value, and an equilibrium re-solve AT the voxel ----
                bool _isPrimary=false; T _dChere=0.0, _preC=0.0;
                bool _eqConv=true, _eqDone=false; plint _eqIters=0; T _eqResid=0.0, _eqOut=0.0;
                if (_CAT[0]=='C' && _worstIdx>=0) {
                    for (size_t _c=0;_c<eq_component_names.size();++_c) if (eq_component_names[_c]==_worstSp) _isPrimary=true;
                    _dChere = computeMax(*computeDensity(dC[_worstIdx]), Box3D(_gx,_gx,_gy,_gy,_gz,_gz));
                    _preC   = _locVal - _dChere;
                    if (useEquilibrium) {
                        std::vector<T> _cc(num_of_substrates);
                        for (plint iS=0; iS<num_of_substrates; ++iS) { T _c0=computeMax(*computeDensity(vec_substr_lattices[iS]), Box3D(_gx,_gx,_gy,_gy,_gz,_gz)); _cc[iS]=std::max(_c0, EquilibriumChemistry<T>::MIN_CONC); }
                        std::vector<T> _eqc = eqSolver.calculate_species_concentrations(_cc);
                        _eqConv = eqSolver.didConverge(); _eqIters = eqSolver.getLastIterations(); _eqResid = eqSolver.getLastResidual();
                        _eqOut = (_worstIdx < (plint)_eqc.size()) ? _eqc[_worstIdx] : 0.0; _eqDone=true;
                    }
                }
                // ---- input-level evidence: cell-Peclet, gradient steepness, suggested deltaP ----
                T _peCell = (D_lattice_fixed>1e-30) ? (_uHere / D_lattice_fixed) : 0.0;
                T _dPfix  = (_MaNow>1e-30) ? (deltaP * 0.02 / _MaNow) : deltaP;
                T _spMax=0.0, _spAvg=0.0, _gradRatio=0.0;
                if (_CAT[0]=='C' && _worstIdx>=0) {
                    _spMax = computeMax(*computeDensity(vec_substr_lattices[_worstIdx]));
                    _spAvg = computeAverage(*computeDensity(vec_substr_lattices[_worstIdx]));
                    _gradRatio = (std::abs(_spAvg)>1e-30) ? (_spMax/std::abs(_spAvg)) : 0.0;
                }
                pcout << "╔══════════════════════════════════════════════════════════════════════════╗\n";
                pcout << "║  SIMULATION FAILED — the solver DETERMINED the cause from the state below   \n";
                pcout << "╠══════════════════════════════════════════════════════════════════════════╣\n";
                pcout << "║  ROOT CAUSE : " << _CAT << "\n";
                pcout << "║  Trigger    : " << reason << "\n";
                pcout << "║  Iteration  : " << iT << "   |   time = " << std::scientific << std::setprecision(4) << iT*ade_dt << std::fixed << " s\n";
                pcout << "║  Sweeps     : push-pull " << pushSweeps << "  |  age " << ageSweeps << "\n";
                pcout << "╟──────────────────────────────────────────────────────────────────────────╢\n";
                pcout << "║  IN PLAIN ENGLISH - what triggered the failure:\n";
                pcout << "║    " << _PLAIN << "\n";
                pcout << "║  IN PLAIN ENGLISH - the solution:\n";
                pcout << "║    " << _PFIX << "\n";
                pcout << "╟──────────────────────────────────────────────────────────────────────────╢\n";
                pcout << "║  WHY IT FAILED (read directly from the fields at failure):\n";
                pcout << "║    " << _WHY1 << "\n";
                pcout << "║    " << _WHY2 << "\n";
                pcout << "║    " << _WHY3 << "\n";
                pcout << "║  HOW TO FIX:\n";
                pcout << "║    " << _SOL1 << "\n";
                pcout << "║    " << _SOL2 << "\n";
                pcout << "╟──────────────────────────────────────────────────────────────────────────╢\n";
                pcout << "║  WHERE (the exact voxel the solver pinpointed):\n";
                pcout << "║    Worst voxel  : (x=" << _gx << ", y=" << _gy << ", z=" << _gz << ")  in domain " << nx << "x" << ny << "x" << nz << "\n";
                pcout << "║    Position     : " << _face << "\n";
                pcout << "║    Offender     : " << _locWhat << " = " << std::scientific << std::setprecision(4) << _locVal << "\n";
                pcout << "║    Local state  : mask-code=" << _mcode << " ; biomass=" << _bHere << " kg/m3" << (_bHere>0.0?" [biofilm voxel]":" [open pore]") << " ; |u|=" << _uHere << "\n";
                pcout << "║    S-system here: HS=" << _valAt("HS",_gx,_gy,_gz) << " SO4=" << _valAt("SO4",_gx,_gy,_gz) << " H2S=" << _valAt("H2S",_gx,_gy,_gz) << " Hp=" << _valAt("Hp",_gx,_gy,_gz) << " Fe2=" << _valAt("Fe2",_gx,_gy,_gz) << std::fixed << "\n";
                if (_CAT[0]=='C' && _worstIdx>=0) {
                    pcout << "╟──────────────────────────────────────────────────────────────────────────╢\n";
                    pcout << "║  WHY THE VALUE IS NEGATIVE (deduced from the run, with proof):\n";
                    pcout << "║    " << _worstSp << (_isPrimary?" is a PRIMARY component":" is a SECONDARY (equilibrium-computed) species") << ".\n";
                    pcout << "║    Stored value here = " << std::scientific << std::setprecision(4) << _locVal << " ; reaction increment dC = " << _dChere << " ; value just before the reaction ~ " << _preC << ".\n";
                    if (_eqDone) {
                        pcout << "║    Equilibrium re-solve AT this voxel: converged=" << (_eqConv?"YES":"NO") << " ; iters=" << _eqIters << "/" << eqSolver.getMaxIterations() << " ; residual=" << _eqResid << " (tol=" << eqSolver.getTolerance() << ").\n";
                        pcout << "║    Solver clamped output for " << _worstSp << " here = " << _eqOut << "  (the solver clamps its output to >= 1e-30, so it CANNOT emit a negative).\n";
                    }
                    if (_eqDone && !_eqConv) {
                        pcout << "║    PROOF: the equilibrium solver did NOT converge here (residual=" << std::scientific << std::setprecision(4) << _eqResid << " >> tol " << eqSolver.getTolerance() << "),\n";
                        pcout << "║    and it also failed to converge on your INITIAL CompLaB.xml water at startup. The speciation is\n";
                        pcout << "║    therefore unreliable and IS the source of the negative. Cell-Peclet here is " << std::fixed << std::setprecision(2) << _peCell << " (advection is\n";
                        pcout << "║    negligible below 2), so this is a CHEMISTRY-INPUT problem, not transport and not kinetics.\n";
                        pcout << "║    EASY FIX: make the CompLaB.xml equilibrium composition feasible (initial concentrations, logK,\n";
                        pcout << "║    stoichiometry, charge balance) so the solver converges - start from the startup INPUT WARNING above.\n";
                    } else if (_eqDone && _dChere < 0.0 && _preC >= 0.0) {
                        pcout << "║    PROOF: the equilibrium solver converged to a non-negative value here, so it is not the source.\n";
                        pcout << "║      -> a KINETICS reaction removed more than was present (before reaction ~" << std::scientific << std::setprecision(4) << _preC << ", dC=" << _dChere << "),\n";
                        pcout << "║         driving it below zero. EASY FIX: smaller reactive step, or cap consumption to the amount available.\n";
                    } else if (_eqDone && _eqOut >= 0.0 && _locVal < 0.0) {
                        pcout << "║    PROOF: the equilibrium solver converged to a non-negative value and no reaction consumed it (dC~0),\n";
                        pcout << "║    so the negative came from the TRANSPORT (LBM advection-diffusion) step.\n";
                        if (_peCell > 2.0) {
                            pcout << "║      -> cell-Peclet=" << std::fixed << std::setprecision(2) << _peCell << " (>2): advection overshoot across a steep gradient.\n";
                            pcout << "║      EASY FIX: clamp densities >= 0 after the stream step, or lower Pe so advection stops overshooting.\n";
                        } else {
                            pcout << "║      -> cell-Peclet=" << std::fixed << std::setprecision(2) << _peCell << " (<2, advection weak): a diffusion / boundary transport artifact.\n";
                            pcout << "║      EASY FIX: clamp densities >= 0 after the stream step; check boundary / initial values in CompLaB.xml.\n";
                        }
                    }
                    pcout << std::fixed;
                }
                pcout << "║  EVIDENCE (why this category and not the others):\n";
                pcout << "║    FLOW      : Ma_now = " << std::scientific << std::setprecision(4) << _MaNow << "  (LBM limit ~0.1) ; u_max = " << _umax << " ; deltaP = " << deltaP << " ; k = " << permeability << "\n";
                pcout << "║    BIOMASS   : max = " << _bmax << " kg/m3 = " << std::fixed << std::setprecision(2) << _capRatio << "x cap ; mean = " << std::scientific << std::setprecision(4) << _bavg << " ; total = " << _bsum << "\n";
                pcout << "║    CHEMISTRY : most-negative solute = " << _worstSp << " at " << _worstMin << " M  (physical floor = 0)\n";
                pcout << "║    FINITE?   : flow=" << (_flowFinite?"yes":"NO") << "  biomass=" << (_bmassFinite?"yes":"NO") << "  chem=" << (_chemFinite?"yes":"NO") << std::fixed << "\n";
                pcout << "╟──────────────────────────────────────────────────────────────────────────╢\n";
                pcout << "║  WHICH INPUT TO FIX (this is a configuration issue, not a code bug):\n";
                if (_CAT[0]=='F') {
                    pcout << "║    File      : CompLaB.xml\n";
                    pcout << "║    Parameter : deltaP (currently " << std::scientific << std::setprecision(4) << deltaP << ")\n";
                    pcout << "║    Change to : " << _dPfix << "   [= deltaP * 0.02 / Ma_now ; targets Ma ~ 0.02]\n";
                    pcout << "║    Or        : use a looser geometry.dat (wider throats raise permeability, now k=" << permeability << ").\n";
                }
                else if (_CAT[0]=='C') {
                    pcout << "║    Cell-Peclet at this voxel = " << std::fixed << std::setprecision(2) << _peCell << " (advection-dominated if > 2) ; gradient(" << _worstSp << ") max/avg = " << _gradRatio << "x\n";
                    if (_eqDone && !_eqConv) {
                        pcout << "║    Cause     : the equilibrium solver did NOT converge here - infeasible chemistry (it also failed on the initial water at startup).\n";
                        pcout << "║    File      : CompLaB.xml (equilibrium block)\n";
                        pcout << "║    Fix       : make initial concentrations / logK / stoichiometry feasible and charge-balanced so the solver converges.\n";
                    } else if (_dChere < 0.0 && _preC >= 0.0) {
                        pcout << "║    Cause     : a KINETICS reaction removed more " << _worstSp << " than was present (dC=" << std::scientific << std::setprecision(4) << _dChere << ").\n";
                        pcout << "║    File      : defineKinetics.hh (biotic)  or  defineAbioticKinetics.hh (abiotic)\n";
                        pcout << "║    Fix #1    : lower the rate constant (Vmax/k) of the reaction consuming " << _worstSp << ".\n";
                        pcout << "║    Fix #2    : reduce dx in CompLaB.xml (smaller reactive step ; now dx=" << dx << ").\n";
                    } else {
                        pcout << "║    Cause     : TRANSPORT of a steep gradient (solver +ve, value <0 before reaction ; cell-Peclet shown above).\n";
                        pcout << "║    File      : CompLaB.xml\n";
                        pcout << "║    Fix #1    : lower Pe from " << std::scientific << std::setprecision(4) << Pe << " to <= " << (Pe*0.3) << " (or clamp densities >=0 after the stream step)\n";
                        pcout << "║    Fix #2    : reduce the initial-concentration contrast of " << _worstSp << " (smooth the " << std::fixed << std::setprecision(1) << _gradRatio << "x cliff)\n";
                    }
                }
                else if (_CAT[0]=='B') {
                    pcout << "║    File      : defineKinetics.hh\n";
                    pcout << "║    Parameter : growth rate (Vmax) of the fastest-growing microbe (biomass hit " << std::fixed << std::setprecision(2) << _capRatio << "x cap)\n";
                    pcout << "║    Fix #1    : lower that Vmax ; confirm the f_cap limiter is active in defineKinetics.hh.\n";
                    pcout << "║    Fix #2    : reduce dx in CompLaB.xml (smaller step ; now dx=" << std::scientific << std::setprecision(4) << dx << ").\n";
                }
                else {
                    pcout << "║    File      : geometry.dat (pore throats too tight)  or  CompLaB.xml\n";
                    pcout << "║    Fix #1    : use a looser packing / wider gaps in geometry.dat.\n";
                    pcout << "║    Fix #2    : raise Bmax (max_bMassRho=" << std::scientific << std::setprecision(4) << max_bMassRho << ") in CompLaB.xml, or lower growth in defineKinetics.hh.\n";
                }
                pcout << std::fixed;
                pcout << "╚══════════════════════════════════════════════════════════════════════════╝\n";
                pcout << "  [CA-FAIL] Writing full VTI snapshot (all species, microbes, mask, age, flow) at iter=" << iT << "...\n";
                plint _a=0, _b=0;
                for (plint iS = 0; iS < num_of_substrates; ++iS) { writeAdvVTI(vec_substr_lattices[iS], iT, vec_subs_names[iS]+"_"); }
                for (plint iM = 0; iM < num_of_microbes; ++iM) {
                    if (bmass_type[iM]==1) { writeAdvVTI(vec_bFilm_lattices[_a], iT, vec_microbes_names[iM]+"_"); ++_a; }
                    else { writeAdvVTI(vec_bFree_lattices[_b], iT, vec_microbes_names[iM]+"_"); ++_b; }
                }
                if (Pe > thrd) { writeNsVTI(nsLattice, iT, "nsLattice_"); }
                writeAdvVTI(maskLattice, iT, mask_filename+"_");
                writeAdvVTI(ageLattice, iT, "ageLattice_");
                pcout << "  [CA-FAIL] All snapshot files written at iter=" << iT << ". Terminating simulation.\n";
            };
            plint whilecount=0;
            T prevBmax = globalBmax; plint stall = 0;   // [CA-ROBUST] progress tracker: stop early if biomass genuinely cannot drain (buried / precipitated-over / biofilm-boxed cells)
            if (!percolationFlag && globalBmax - max_bMassRho > thrd) {
                diag_ca_triggers++;
                if (track_performance == 1) global::timer("ca").restart();
                while (globalBmax - max_bMassRho > thrd) {
                    for (plint iM=0; iM<bfilm_count; ++iM) vec_bFcopy_lattices[iM]=copybFilmLattice;
                    if (halfflag == 0) applyProcessingFunctional(new pushExcessBiomass3D<T,RXNDES>(max_bMassRho, nx, ny, nz, 1, caLlen, no_dynamics, bounce_back, pore_dynamics), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
                    else applyProcessingFunctional(new halfPushExcessBiomass3D<T,RXNDES>(max_bMassRho, nx, ny, nz, 1, caLlen, no_dynamics, bounce_back, pore_dynamics), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
                    applyProcessingFunctional(new pullExcessBiomass3D<T,RXNDES>(nx, ny, nz, 1, caLlen), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
                    applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, caLlen, bounce_back, no_dynamics, bio_dynamics, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
                    globalBmax = computeMax(*computeDensity(totalbFilmLattice));
                    diag_ca_redistributions++;
                    // [CA-ROBUST] if a sweep no longer lowers the biofilm max, the remaining excess is boxed in
                    // (wall corner, FeS-sealed throat, or surrounded by full biofilm). Leave those cells full and
                    // stop early instead of grinding the full 2000 sweeps. Threshold spans the age-update interval (50).
                    if (globalBmax > prevBmax - thrd) { if (++stall >= 100) { dumpAllFields("PUSH-PULL stalled: biofilm cannot spread (buried/clogged cells at cap; no drainage in 100 sweeps)", whilecount, 0); return -1; } }
                    else stall = 0;
                    prevBmax = globalBmax;
                    if (whilecount%50 == 0) {
                        plint diff = 1, whilecount1 = 0;
                        while (diff != 0) {
                            plint old_totAge = util::roundToInt(computeAverage(*computeDensity(ageLattice))*nx*ny*nz);
                            applyProcessingFunctional(new updateAgeDistance3D<T,RXNDES>(max_bMassRho, nx, ny, nz), ageLattice.getBoundingBox(), ageNdistance_lattices);
                            plint new_totAge = util::roundToInt(computeAverage(*computeDensity(ageLattice))*nx*ny*nz);
                            diff = new_totAge-old_totAge;
                            ++whilecount1;
                            if (whilecount1 > 1000) { dumpAllFields("AGE-DISTANCE relaxation did not settle (>1000 sub-sweeps)", whilecount, whilecount1); return -1; }
                        }
                    }
                    if (whilecount > 1000) { dumpAllFields("PUSH-PULL redistribution did not settle (>1000 sweeps)", whilecount, 0); return -1; }
                    ++whilecount;
                }
                if (track_performance == 1) { catime+=global::timer("ca").getTime(); global::timer("ca").stop(); }
            }
        }
        if (fd_count > 0) {
            // [FIX] updateLocalMaskNtotalLattices3D assumes the CA layout
            //   [B.., copies.., totalBiomass, mask, age]  (mask at length-2).
            //   ptr_fd_lattices is [B.., copies.., mask] with mask at length-1, so
            //   calling it on that vector made it read the last COPY lattice as the
            //   mask and the one before as the total biomass -- and WRITE to both.
            //   Every finite-difference run silently corrupted its own scratch
            //   lattices and mis-classified the biofilm.  Use a correctly shaped
            //   vector instead; ptr_fd_lattices is still right for fdDiffusion3D,
            //   which expects mask at length-1.
            applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, fdMaskLen, bounce_back, no_dynamics, bio_dynamics, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_fd_mask);
            for (plint iM=0; iM<bfilm_count; ++iM) vec_bFcopy_lattices[iM]=vec_bFilm_lattices[iM];
            for (plint iP=0; iP<bfree_count; ++iP) vec_bPcopy_lattices[iP]=vec_bFree_lattices[iP];
            applyProcessingFunctional(new fdDiffusion3D<T,RXNDES>(nx, ny, nz, fdLlen, 1, bioNUinPore[0]), vec_bFilm_lattices[0].getBoundingBox(), ptr_fd_lattices);
            applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, fdMaskLen, bounce_back, no_dynamics, bio_dynamics, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_fd_mask);
        }

        // Update flow and dynamics
        if (ca_count > 0 || fd_count > 0) {
            if (track_performance == 1) global::timer("ca").restart();
            new_totMask = util::roundToInt(computeAverage(*computeDensity(maskLattice))*nx*ny*nz);
            if (std::abs(old_totMask-new_totMask)>0) {
                old_totMask = new_totMask;
                applyProcessingFunctional(new updateAgeDistance3D<T,RXNDES>(max_bMassRho, nx, ny, nz), ageLattice.getBoundingBox(), ageNdistance_lattices);
                if (iT % ade_update_interval == 0) {
                    if (soluteDindex == 1) applyProcessingFunctional(new updateSoluteDynamics3D<T,RXNDES>(num_of_substrates, bounce_back, no_dynamics, pore_dynamics, substrOMEGAinbFilm, substrOMEGAinPore), vec_substr_lattices[0].getBoundingBox(), substrate_lattices);
                    if (bmassDindex == 1) applyProcessingFunctional(new updateBiomassDynamics3D<T,RXNDES>((plint)vec_bFree_lattices.size(), bounce_back, no_dynamics, pore_dynamics, bioOMEGAinbFilm, bioOMEGAinPore), vec_bFree_lattices[0].getBoundingBox(), planktonic_lattices);
                }
                if (track_performance == 1) { catime+=global::timer("ca").getTime(); global::timer("ca").stop(); }
                if (iT % ns_update_interval == 0 && Pe > thrd && ns_saturate == 0) {
                    if (track_performance == 1) global::timer("NS").restart();
                    applyProcessingFunctional(new updateNsLatticesDynamics3D<T,NSDES,T,RXNDES>(nsLatticeOmega, vec_permRatio[0], pore_dynamics, no_dynamics, bounce_back), nsLattice.getBoundingBox(), nsLattice, maskLattice);
                    for (plint iT2 = 0; iT2 < ns_maxiTer_2; ++iT2) {
                        nsLattice.collideAndStream();
                        ns_convg2.takeValue(getStoredAverageEnergy(nsLattice),false);
                        if (ns_convg2.hasConverged()) break;
                        if (iT2 == (ns_maxiTer_2-1)) ns_saturate = 1;
                    }
                    if (ns_saturate == 1) {
                        T outletvel = computeAverage(*computeVelocityComponent(nsLattice, Box3D(nx-2,nx-2, 0,ny-1, 0,nz-1), 0));
                        if (outletvel > thrd) ns_saturate = 0;
                        else { pcout << "\n  [NS] Percolation limit reached at iter=" << iT << "\n"; percolationFlag = 1; }
                    }
                    // [NS-ROBUST] catch flow divergence (NaN/Inf energy) from near-complete clogging: stop cleanly
                    // instead of coupling a non-finite velocity into the solutes/biomass (which then segfaults the CA).
                    { T nsE = getStoredAverageEnergy(nsLattice); if (std::isnan(nsE) || std::isinf(nsE)) { pcout << "\n  [NS] Flow solve diverged at iter=" << iT << " (pore clogged / Ma runaway); stopping cleanly with last valid data.\n"; percolationFlag = 1; } }
                    if (!percolationFlag) {
                    for (plint iS = 0; iS < num_of_substrates; ++iS) if (!vec_immobile[iS]) latticeToPassiveAdvDiff(nsLattice, vec_substr_lattices[iS], vec_substr_lattices[iS].getBoundingBox());
                    if (lb_count > 0) {
                        for (plint iM = 0; iM < num_of_microbes; ++iM) {
                            if (solver_type[iM]==3) {
                                if (bmass_type[iM]==1) latticeToPassiveAdvDiff(nsLattice, vec_bFilm_lattices[loctrack[iM]], vec_bFilm_lattices[loctrack[iM]].getBoundingBox());
                                else latticeToPassiveAdvDiff(nsLattice, vec_bFree_lattices[loctrack[iM]], vec_bFree_lattices[loctrack[iM]].getBoundingBox());
                            }
                        }
                    }
                    }
                    if (track_performance == 1) { nstime += global::timer("NS").getTime(); global::timer("NS").stop(); }
                }
            }
            else { if (track_performance == 1) { catime+=global::timer("ca").getTime(); global::timer("ca").stop(); } }
        }

        // ============================================================================
        // [PRECIP-VOP] Node conversion (VOP) + flow feedback: full immobile-solid voxels
        //   (P >= max_precipRho) flip to solid; if geometry changed, re-solve the flow.
        // ============================================================================
        if (precip_enabled && precip_solidSub >= 0 && (iT % precip_update_interval == 0)) {
            std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_precip_conv;
            ptr_precip_conv.push_back(&vec_substr_lattices[precip_solidSub]);
            ptr_precip_conv.push_back(&maskLattice);
            ptr_precip_conv.push_back(&phaseLattice);   // [DISSOL-VOP] so a sealing voxel is stamped
            plint pre_totMask  = util::roundToInt(computeAverage(*computeDensity(maskLattice))*nx*ny*nz);
            applyProcessingFunctional(new precipNodeConversion3D<T,RXNDES>(nx, ny, nz, max_precipRho, no_dynamics, bounce_back, no_dynamics, precip_phase_id),
                                      vec_substr_lattices[precip_solidSub].getBoundingBox(), ptr_precip_conv);
            plint post_totMask = util::roundToInt(computeAverage(*computeDensity(maskLattice))*nx*ny*nz);
            if (std::abs(pre_totMask - post_totMask) > 0 && Pe > thrd && ns_saturate == 0) {
                applyProcessingFunctional(new updateNsLatticesDynamics3D<T,NSDES,T,RXNDES>(nsLatticeOmega, precip_permRatio, pore_dynamics, no_dynamics, bounce_back),
                                          nsLattice.getBoundingBox(), nsLattice, maskLattice);
                for (plint iT2 = 0; iT2 < ns_maxiTer_2; ++iT2) {
                    nsLattice.collideAndStream();
                    ns_convg2.takeValue(getStoredAverageEnergy(nsLattice),false);
                    if (ns_convg2.hasConverged()) break;
                    if (iT2 == (ns_maxiTer_2-1)) ns_saturate = 1;
                }
                if (ns_saturate == 1) {
                    T outletvel = computeAverage(*computeVelocityComponent(nsLattice, Box3D(nx-2,nx-2, 0,ny-1, 0,nz-1), 0));
                    if (outletvel > thrd) ns_saturate = 0;
                    else { pcout << "\n  [PRECIP-VOP] Percolation limit reached (clogged) at iter=" << iT << "\n"; percolationFlag = 1; }
                }
                for (plint iS = 0; iS < num_of_substrates; ++iS) if (!vec_immobile[iS]) latticeToPassiveAdvDiff(nsLattice, vec_substr_lattices[iS], vec_substr_lattices[iS].getBoundingBox());
            }
        }

        // ============================================================================
        // [DISSOL-VOP] Reverse node conversion: a phase-carrying voxel whose mineral
        //   has fallen below reopen_fraction x full_density becomes pore again, and
        //   its phase id is cleared.  If any voxel reopened, the pore space changed
        //   and the flow has to be re-solved -- the mirror of the block above.
        //
        //   The reopen threshold is deliberately BELOW the seal threshold: with a
        //   single threshold a voxel sitting on it would flip every update interval,
        //   and every flip costs a full Navier-Stokes re-solve.
        // ============================================================================
        if (dissolCfg.enabled && (iT % dissolCfg.update_interval == 0)) {
            const plint pre_totMask = util::roundToInt(computeAverage(*computeDensity(maskLattice))*nx*ny*nz);
            applyProcessingFunctional(new dissolNodeConversion3D<T,RXNDES>(
                                          num_of_substrates,
                                          pore_dynamics.empty() ? (plint) 2 : pore_dynamics[0],
                                          no_dynamics, bounce_back, &dissolCfg),
                                      vec_substr_lattices[0].getBoundingBox(), ptr_dissol_conv);
            const plint post_totMask = util::roundToInt(computeAverage(*computeDensity(maskLattice))*nx*ny*nz);

            if (std::abs(pre_totMask - post_totMask) > 0 && Pe > thrd) {
                pcout << "  [DISSOL-VOP] pore space reopened at iter=" << iT << "; re-solving the flow\n";
                // A reopened voxel is ordinary pore again, so it takes the normal pore
                // relaxation, not the reduced precipitate one.
                applyProcessingFunctional(new updateNsLatticesDynamics3D<T,NSDES,T,RXNDES>(nsLatticeOmega, vec_permRatio.empty() ? (T)1 : vec_permRatio[0], pore_dynamics, no_dynamics, bounce_back),
                                          nsLattice.getBoundingBox(), nsLattice, maskLattice);
                ns_saturate = 0;                       // flow may percolate again
                for (plint iT2 = 0; iT2 < ns_maxiTer_2; ++iT2) {
                    nsLattice.collideAndStream();
                    ns_convg2.takeValue(getStoredAverageEnergy(nsLattice),false);
                    if (ns_convg2.hasConverged()) break;
                }
                for (plint iS = 0; iS < num_of_substrates; ++iS)
                    if (!vec_immobile[iS]) latticeToPassiveAdvDiff(nsLattice, vec_substr_lattices[iS], vec_substr_lattices[iS].getBoundingBox());
                percolationFlag = 0;                   // reopening can undo a clog
            }
        }

        // Streaming
        if (track_performance == 1) global::timer("cns").restart();
        for (plint iS = 0; iS < num_of_substrates; ++iS) if (!vec_immobile[iS]) vec_substr_lattices[iS].stream();
        if (lb_count > 0) {
            for (plint iM = 0; iM < num_of_microbes; ++iM) {
                if (solver_type[iM]==3) {
                    if (bmass_type[iM]==1) vec_bFilm_lattices[loctrack[iM]].stream();
                    else vec_bFree_lattices[loctrack[iM]].stream();
                }
            }
        }
        if (track_performance == 1) { nstime += global::timer("cns").getTime(); global::timer("cns").stop(); }
        if (percolationFlag == 1) break;
    }
    // ════════════════════════════════════════════════════════════════════════════
    // PHASE 7: FINAL OUTPUT FILES
    // ════════════════════════════════════════════════════════════════════════════
    pcout << "\n┌────────────────────────────────────────────────────────────────────────┐\n";
    pcout << "│ PHASE 7: WRITING FINAL OUTPUT FILES                                   │\n";
    pcout << "└────────────────────────────────────────────────────────────────────────┘\n";

    // Final output
    if (track_performance == 0) {
        pcout << "  Saving VTI and CHK files...\n";
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            writeAdvVTI(vec_substr_lattices[iS], iT, vec_subs_names[iS]+"_");
            saveBinaryBlock(vec_substr_lattices[iS], str_outputDir+ade_filename+std::to_string(iS)+"_"+std::to_string(iT)+".chk");
            pcout << "    [OK] " << vec_subs_names[iS] << " saved\n";
        }
        tmpIT0=0; tmpIT1=0;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (bmass_type[iM]==1) {
                writeAdvVTI(vec_bFilm_lattices[tmpIT0], iT, vec_microbes_names[iM]+"_");
                saveBinaryBlock(vec_bFilm_lattices[tmpIT0], str_outputDir+bio_filename+std::to_string(iM)+"_"+std::to_string(iT)+".chk");
                pcout << "    [OK] " << vec_microbes_names[iM] << " saved\n";
                ++tmpIT0;
            }
            else {
                writeAdvVTI(vec_bFree_lattices[tmpIT1], iT, vec_microbes_names[iM]+"_");
                saveBinaryBlock(vec_bFree_lattices[tmpIT1], str_outputDir+bio_filename+std::to_string(iM)+"_"+std::to_string(iT)+".chk");
                pcout << "    [OK] " << vec_microbes_names[iM] << " saved\n";
                ++tmpIT1;
            }
        }
        writeAdvVTI(maskLattice, iT, mask_filename+"_");
        saveBinaryBlock(maskLattice, str_outputDir+mask_filename+"_"+std::to_string(iT)+".chk");
        pcout << "    [OK] Mask lattice saved\n";
        if (Pe > thrd) {
            writeNsVTI(nsLattice, iT, "nsLattice_");
            saveBinaryBlock(nsLattice, str_outputDir+ns_filename+".chk");
            pcout << "    [OK] Flow field saved\n";
        }

        // [RATE-OUT] The same sample and write for the final state, so the last
        //   snapshot has rates beside its concentrations.
        if (rateOut.enabled && rateOut.nch() > 0) {
            applyProcessingFunctional(
                new sample_reaction_rates<T,RXNDES>(
                    nx, ny, nz, num_of_substrates, kns_count, rateOut.nch(),
                    no_dynamics, bounce_back,
                    rateOut.doBio, rateOut.doAbio, rateOut.surfaceOnly,
                    rateOut.nrb, rateOut.nra,
                    rateOut.offBioC(), rateOut.offBioM(), rateOut.offBioR(),
                    rateOut.offAbioC(), rateOut.offAbioR()),
                vec_substr_lattices[0].getBoundingBox(), ptr_rate_lattices);
            KineticsStats::resetIteration();
            writeRateVTI(rateLat, rateOut, iT, nx, ny, nz);
            pcout << "    [OK] Reaction rates saved (" << rateOut.nch() << " channels)\n";
        }
    }

    // ════════════════════════════════════════════════════════════════════════════
    // PHASE 8-9: SUMMARY AND STATISTICS
    // ════════════════════════════════════════════════════════════════════════════
    T TET = global::timer("total").getTime(); global::timer("total").stop();

    pcout << "\n╔══════════════════════════════════════════════════════════════════════════╗\n";
    pcout << "║                         SIMULATION COMPLETE                              ║\n";
    pcout << "╠══════════════════════════════════════════════════════════════════════════╣\n";
    pcout << "║ TIMING:                                                                  ║\n";
    pcout << "║   Total iterations: " << iT << "\n";
    pcout << "║   Simulated time:   " << std::scientific << iT*ade_dt << " s\n" << std::fixed;
    pcout << "║   Wall clock:       " << TET << " s (" << TET/60 << " min)\n";
    pcout << "╠══════════════════════════════════════════════════════════════════════════╣\n";
    pcout << "║ SIMULATION MODE:                                                         ║\n";
    pcout << "║   Biotic mode:      " << (biotic_mode ? "YES (with microbes)" : "NO (abiotic)") << "\n";
    pcout << "║   Kinetics (biotic):" << (enable_kinetics ? " ENABLED" : " DISABLED") << "\n";
    pcout << "║   Kinetics (abiotic):" << (enable_abiotic_kinetics ? "ENABLED" : "DISABLED") << "\n";
    pcout << "║   Equilibrium:      " << (useEquilibrium ? "ENABLED" : "DISABLED") << "\n";
    pcout << "║   Validation diag:  " << (enable_validation_diagnostics ? "ENABLED" : "DISABLED") << "\n";
    if (bfilm_count > 0) {
        T finalBmax = computeMax(*computeDensity(totalbFilmLattice));
        T totalGrowth = (diag_initial_biomass > 0) ? ((finalBmax - diag_initial_biomass) / diag_initial_biomass * 100.0) : 0.0;
        pcout << "╠══════════════════════════════════════════════════════════════════════════╣\n";
        pcout << "║ BIOMASS RESULTS:                                                         ║\n";
        pcout << "║   Initial max:      " << std::scientific << diag_initial_biomass << " kg/m³\n";
        pcout << "║   Final max:        " << finalBmax << " kg/m³\n" << std::fixed;
        pcout << "║   Growth:           " << totalGrowth << "%\n";
        pcout << "║   CA triggers:      " << diag_ca_triggers << "\n";
        pcout << "║   Redistributions:  " << diag_ca_redistributions << "\n";
    }
    pcout << "╠══════════════════════════════════════════════════════════════════════════╣\n";
    pcout << "║ FINAL CONCENTRATIONS:                                                    ║\n";
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        T sMin = computeMin(*computeDensity(vec_substr_lattices[iS]));
        T sMax = computeMax(*computeDensity(vec_substr_lattices[iS]));
        T sAvg = computeAverage(*computeDensity(vec_substr_lattices[iS]));
        pcout << "║   " << vec_subs_names[iS] << ": min=" << std::scientific << sMin
              << " avg=" << sAvg << " max=" << sMax << std::fixed << "\n";
    }
    pcout << "╚══════════════════════════════════════════════════════════════════════════╝\n";

    if (track_performance == 1) {
        pcout << "\n┌────────────────────────────────────────────────────────────────────────┐\n";
        pcout << "│ PERFORMANCE TIMING BREAKDOWN                                           │\n";
        pcout << "├────────────────────────────────────────────────────────────────────────┤\n";
        pcout << "│   NS (flow):         " << nstime << " s\n";
        pcout << "│   ADE (transport):   " << adetime << " s\n";
        pcout << "│   Collide+Stream:    " << cnstime << " s\n";
        if (ca_count > 0) pcout << "│   CA (biomass):      " << catime << " s\n";
        if (kns_count > 0) pcout << "│   Kinetics:          " << knstime << " s\n";
        if (useEquilibrium) pcout << "│   Equilibrium:       " << eqtime << " s\n";
        if (mmcfg.mm_count  > 0) pcout << "│   FBA:               " << fbatime << " s\n";
        if (mmcfg.srg_count > 0) pcout << "│   Surrogate:         " << srgtime << " s\n";
        if (sym_count > 0)       pcout << "│   Symbolic:          " << symtime << " s\n";
        if (gnn_count > 0)       pcout << "│   Graph network:     " << gnntime << " s\n";
        pcout << "└────────────────────────────────────────────────────────────────────────┘\n";
    }

    if (useEquilibrium) eqSolver.printStatistics();

    /* [SYM][GNN] How often a fitted rate law or network was asked for a number outside the box it
     * was fitted over.  Neither one fails out there -- it returns a confident answer -- so this
     * count is the only thing standing between an extrapolated run and a result.  It is printed
     * unconditionally, next to the equilibrium statistics, and warns above 5%. */
    pcout << complab_sym::runtimeReport();
    pcout << complab_gnn::runtimeReport();

    // Free allocated memory
    // Optional metabolic layer: delete the GLPK problems / release the cobra
    // models and shut the interpreter down.  Safe to call when nothing was on.
    release_metabolic_solvers(mmcfg);

    free(main_path);
    free(src_path);
    free(input_path);
    free(output_path);
    free(ns_filename);

    pcout << "\n╔══════════════════════════════════════════════════════════════════════════╗\n";
    pcout << "║                       Simulation Finished!                               ║\n";
    pcout << "║                                                                          ║\n";
    pcout << "║  Author:  Shahram Asgari                                                 ║\n";
    pcout << "║  Advisor: Dr. Christof Meile                                             ║\n";
    pcout << "║  Lab:     Meile Lab, University of Georgia                               ║\n";
    pcout << "╚══════════════════════════════════════════════════════════════════════════╝\n\n";

    return 0;
}
