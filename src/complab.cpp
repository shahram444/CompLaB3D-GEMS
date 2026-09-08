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
#include "precipitationVOP.hh"         // [PRECIP-VOP] surface precipitation + pore-clogging feedback
#include "dissolutionVOP.hh"           // [DISSOL-VOP] mineral dissolution + pore re-opening
// [NEW] The pipeline blocks: <model_source>, geometry generation, <surrogate>, <diagnostics>.
//   Header-only, and every one of them is inert unless its XML block is present, so a run with an
//   unchanged CompLaB.xml takes exactly the path it took before these lines existed.
#include "complab3d_integration.hh"
#include "complab3d_outerfaces.hh"     // the four faces nobody gave a boundary condition
#include "complab3d_upscale.hh"        // one aggregate -> an effectiveness factor
#include "complab3d_energyfield.hh"   // dG and F_T as fields, when <thermodynamics> is on
#include "complab3d_srgtrain_glpk.hh"  // [NEW] the LP callback that in-run surrogate training uses
#include <algorithm>
#include <cctype>
#include <sstream>

#include <chrono>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
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


/* ==================================================================================================
 *  complab_upscale_sample  --  one diagnostic interval's worth of the upscaling record
 *
 *  Lives here rather than in complab3d_upscale.hh because it needs two things that belong to the
 *  program and not to a header: defineRxnKinetics(), which is the user's own compiled rate law, and
 *  the thermodynamic gate, which is bound to this run's microbes. Everything reusable -- the
 *  reductions, the classical curve, the record and the report -- is in the header.
 *
 *  The numerator is measured from the increment lattices, which hold the rate every path just
 *  computed. The denominator is the SAME rate law evaluated once at the bulk composition, which is
 *  exactly what a continuum model would do with one concentration per grid block. Their ratio is
 *  the error that upscaling would make, measured rather than argued.
 * ================================================================================================== */
/* The biomass a run is actually holding, summed over the per-microbe lattices.
 *
 * NOT over totalbFilmLattice, which is the combined field the cellular automaton works on. Its
 * wall voxels are never parked at zero the way the per-microbe ones are, so reading its
 * populations counts one unit of nothing for every wall voxel in the domain -- 2000 of them on
 * example 07, against a real biomass of 126.
 *
 * The per-microbe lattices ARE parked, by the stabilisation block in main(), so the population
 * sum over them is the honest total: what is in the open voxels plus what is in transit at a
 * wall, which computeDensity() alone cannot see. */
static T complab_total_biomass(
        std::vector< MultiBlockLattice3D<T,RXNDES> > &bFilm,
        std::vector< MultiBlockLattice3D<T,RXNDES> > &bFree,
        T *seenOut = 0)
{
    T total = T(), seen = T();
    for (size_t i = 0; i < bFilm.size(); ++i) {
        seen  += computeSum(*computeDensity(bFilm[i]));
        total += (T) PopulationSum3D(bFilm[i].getBoundingBox(), bFilm[i]);
    }
    for (size_t i = 0; i < bFree.size(); ++i) {
        seen  += computeSum(*computeDensity(bFree[i]));
        total += (T) PopulationSum3D(bFree[i].getBoundingBox(), bFree[i]);
    }
    if (seenOut) *seenOut = seen;
    return total;
}

static void complab_upscale_sample(
        plint iT, const integ::Config &icfg, T dt, T dx, plint nx, plint ny, plint nz,
        plint num_of_substrates, plint num_of_microbes,
        std::vector< MultiBlockLattice3D<T,RXNDES> > &subs,
        std::vector< MultiBlockLattice3D<T,RXNDES> > &dC,
        MultiBlockLattice3D<T,RXNDES> &maskLattice,
        std::vector< MultiBlockLattice3D<T,RXNDES> > &bFilm,
        std::vector< MultiBlockLattice3D<T,RXNDES> > &bFree,
        const std::vector<bool> &bmass_type, const std::vector<plint> &loctrack,
        const std::vector<plint> &pore_dynamics,
        const std::vector<T> &vec_solute_bFilmD)
{
    (void) ny; (void) nz;
    const plint iS = (plint) icfg.upsSpecies;
    if (iS < 0 || iS >= num_of_substrates) return;

    const Box3D box(1, nx-2, 0, ny-1, 0, nz-1);
    std::unique_ptr< MultiScalarField3D<T> > live = computeDensity(maskLattice);

    /* The bulk is open pore OUTSIDE the aggregate. Named as a material list rather than
     * "everything that is not aggregate", so a wall or a grain can never be counted as water. */
    std::vector<plint> bulkMat;
    for (size_t k = 0; k < pore_dynamics.size(); ++k) {
        bool isAgg = false;
        for (size_t j = 0; j < icfg.upsAggregateMat.size(); ++j)
            if (pore_dynamics[k] == icfg.upsAggregateMat[j]) { isAgg = true; break; }
        if (!isAgg) bulkMat.push_back(pore_dynamics[k]);
    }

    complab_upscale::Point p;
    p.iteration = (long) iT;

    /* < r > inside the aggregate, straight off the increment lattice. */
    {
        complab_upscale::RegionRateFunctional3D<T,RXNDES,T> f(icfg.upsAggregateMat, dt);
        applyProcessingFunctional(f, box, dC[iS], *live);
        p.aggregateVoxels = f.getCount();
        p.rateMean = (p.aggregateVoxels > 0) ? f.getSum() / p.aggregateVoxels : 0.0;
    }
    /* Every substrate at its bulk mean, so the rate law sees a composition a continuum model
     * could actually have had. */
    std::vector<double> Cbulk((size_t) num_of_substrates, 0.0);
    for (plint k = 0; k < num_of_substrates; ++k) {
        complab_upscale::RegionMeanFunctional3D<T,RXNDES,T> f(bulkMat);
        applyProcessingFunctional(f, box, subs[k], *live);
        if (k == iS) { p.bulkVoxels = f.getCount(); }
        Cbulk[(size_t) k] = (f.getCount() > 0) ? f.getSum() / f.getCount() : 0.0;
    }
    p.bulkConc = Cbulk[(size_t) iS];

    /* The mean biomass inside the aggregate. r(C_bulk) has to be evaluated at the same amount of
     * catalyst the aggregate actually holds, or the ratio measures the biomass rather than the
     * transport limitation it is supposed to measure. */
    std::vector<double> Bbulk((size_t) num_of_microbes, 0.0);
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        MultiBlockLattice3D<T,RXNDES> &bl = (bmass_type[iM] == 1) ? bFilm[loctrack[iM]]
                                                                 : bFree[loctrack[iM]];
        complab_upscale::RegionMeanFunctional3D<T,RXNDES,T> f(icfg.upsAggregateMat);
        applyProcessingFunctional(f, box, bl, *live);
        Bbulk[(size_t) iM] = (f.getCount() > 0) ? f.getSum() / f.getCount() : 0.0;
    }

    /* r(C_bulk): the user's own rate law, once, at the bulk composition. mask 2 so the law's own
     * "no biology in a wall" guard lets it through. */
    std::vector<double> subsR((size_t) num_of_substrates, 0.0), bioR((size_t) num_of_microbes, 0.0);
    defineRxnKinetics(Bbulk, Cbulk, subsR, bioR, (plint) 2);
    p.rateAtBulk = -subsR[(size_t) iS];

    /* The gate at the bulk composition. This is the number somebody upscaling by hand would reach
     * for, and reporting it beside eta is what shows it is not a substitute: r(C_bulk) already has
     * it in, so whatever eta departs from 1 is exactly the error that shortcut would make. */
    if (complab_thermo::enabled()) {
        double g = 1.0;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            const double gm = complab_thermo::gateFor((int) iM, Cbulk);
            if (gm < g) g = gm;
        }
        p.gateAtBulk = g;
        p.rateAtBulk *= g;
    }

    /* R from the measured aggregate volume when it was not declared: the radius of the sphere of
     * the same volume, which is the length the Thiele modulus is built on. dx is in metres. */
    double R = icfg.upsRadius * 1e-6;
    if (!(R > 0) && p.aggregateVoxels > 0)
        R = std::pow(3.0 * p.aggregateVoxels / (4.0 * 3.14159265358979323846), 1.0/3.0) * (double) dx;

    double D = icfg.upsDiffusivity;
    if (!(D > 0) && iS < (plint) vec_solute_bFilmD.size()) D = (double) vec_solute_bFilmD[(size_t) iS];

    p.eta = (std::fabs(p.rateAtBulk) > 0) ? p.rateMean / p.rateAtBulk : 0.0;
    const double k = (std::fabs(p.bulkConc) > 0) ? p.rateAtBulk / p.bulkConc : 0.0;   /* 1/s */
    p.thiele = (k > 0 && D > 0 && R > 0) ? R * std::sqrt(k / D) : 0.0;
    p.etaClassical = complab_upscale::classicalEta(p.thiele);

    complab_upscale::record().push_back(p);
}

int main(int argc, char **argv) {

    plbInit(&argc, &argv);

    /* The configuration file, from the command line when one is given.  Every
     * reader in the program goes through complab_input::configPath(); before
     * this line existed, argv[1] was accepted and silently ignored, so
     * `./complab variant.xml` ran CompLaB.xml instead. */
    if (argc > 1 && argv[1] && argv[1][0] != '\0') complab_input::configPath() = argv[1];
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
    /* The conserved starting total, so the closing report can say what actually grew
     * rather than what happened to the peak. */
    T diag_initial_total_biomass = 0.0;

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
    /* 0 = Dirichlet (held at a value), 1 = Neumann (zero gradient, an outflow),
     * 2 = closed (no flux). See the parser for why the third one had to exist. */
    std::vector<plint> vec_left_btype, vec_right_btype, bio_left_btype, bio_right_btype;
    std::vector<T> vec_c0, vec_b0_free, vec_left_bcondition, vec_right_bcondition, bio_left_bcondition, bio_right_bcondition, vec_permRatio;
    std::vector< std::vector<T> > vec_b0_all, vec_b0_film, vec_Kc_kns, vec_Vmax, vec_Vmax_kns;
    std::vector<T> vec_mu;
    std::vector< std::vector<T> > vec_Kc;

    // variables with default values
    std::string ade_filename, bio_filename;
    bool read_NS_file=0, read_ADE_file=0, soluteDindex=0, bmassDindex=0, track_performance=0., halfflag=0;
    plint no_dynamics=0, bounce_back=1, ns_rerun_iT0=0, ns_update_interval=1, ade_update_interval=1,
        ns_maxiTer_1, ns_maxiTer_2, ade_rerun_iT0=0, ade_maxiTer=10000000, ade_VTI_iTer=1000, ade_CHK_iTer=1000000;
    /* [v1.3] thrd_bFilmFrac had no initialiser. The parser only assigns it when
     * <thrd_biofilm_fraction> is present, and only DEMANDS it when a microbe is on the CA, so an
     * abiotic run or a finite-difference/LBM biofilm run that omits the tag read an indeterminate
     * double -- printed in the configuration summary, and used as the pore/biofilm reclassification
     * threshold in updateLocalMaskNtotalLattices3D, which decides voxel identity, solute omega and
     * the flow geometry. Zero is the value the [FIX-3D] note in processors_part2.hh already
     * assumed was the default. */
    T tau=0.8, max_bMassRho=1., ns_converge_iT1=1e-8, ns_converge_iT2=1e-4, ade_converge_iT=1e-8, thrd_bFilmFrac=0.;

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
    // [NEW] THE PIPELINE BLOCKS
    //
    //   <model_source>, <surrogate>, <diagnostics>, and geometry generation inside
    //   <domain>.  These are the settings that replace the Python and MATLAB steps a
    //   user used to run by hand before and after a simulation.
    //
    //   Read here, from a fresh XMLreader rather than threaded through
    //   initialize_complab(), which already takes over sixty arguments by reference.
    //   Every block is optional: with none of them present icfg keeps its defaults and
    //   nothing below this point behaves differently.
    // ============================================================================
    integ::Config icfg;
    {
        XMLreader idoc(complab_input::configPath());
        const std::string imsg = integ::readConfig(idoc, icfg);
        if (!imsg.empty()) { pcout << imsg; return -1; }
    }

    //   <diagnostics>: the run's own scalar record.  Until now CompLB3D wrote VTI volumes and
    //   nothing else, so answering "did mass balance?" or "how did porosity change?" meant
    //   post-processing the volumes in Python.  With this block on, the run writes a summary
    //   CSV as it goes and checks the conserved sums the user names.
    //   Only rank 0 writes the file; the numbers are Palabos reductions and are already global.
    //   Configured further down, once the output directory is known, so that summary.csv
    //   lands beside the VTI files rather than in whatever directory the job was launched
    //   from.  Declared here so it is in scope for the whole run.
    complab_diag::Diagnostics diag;

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
        XMLreader immdoc(complab_input::configPath());
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
        XMLreader pdoc(complab_input::configPath());
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
        XMLreader ddoc(complab_input::configPath());
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
    //   [NEW] Now that the output directory is known, point the summary CSV at it.
    integ::setupDiagnostics(icfg, diag, global::mpi().isMainProcessor(), str_outputDir);

    char *pyFileName = (char*)"complab3d_cobrapy";

    // ---- [NEW] <model_source>: find the genome-scale model before loading it -----
    //   Bundled with the code, already in the cache, or downloaded -- in that order, and
    //   only downloading if <allow_download> says so.  The file is checked against
    //   models/manifest.txt and the run is told plainly if it is not the revision the
    //   manifest describes.
    //
    //   Only rank 0 may write into the cache; the others wait at the barrier and then
    //   find the file already there.
    if (!icfg.modelSource.empty()) {
        std::string mlog;
        bool mfatal = false;
        const std::string mpath = integ::resolveModel(icfg, "", global::mpi().isMainProcessor(),
                                                      mlog, mfatal);
        pcout << mlog;
#ifdef PLB_MPI_PARALLEL
        global::mpi().barrier();
#endif
        if (mfatal) return -1;

        //   Give the resolved path to every FBA microbe that did not name its own file.
        //   A microbe with an explicit <model_filename> keeps it, so a two-organism run
        //   can mix a bundled model with a local one.
        //   usesMetabolic, not usesFBA: a `surrogate` microbe is not an FBA microbe at run
        //   time, but it still needs the model file if the surrogate is to be TRAINED from it.
        plint adopted = 0;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesMetabolic(reaction_type[iM])) continue;
            if (!mmcfg.model_filename[iM].empty()) continue;
            mmcfg.model_filename[iM] = mpath;
            ++adopted;
        }
        if (adopted > 0)
            pcout << "  [MODEL] " << adopted << " microbe(s) will use " << mpath << "\n";
    }

    if (mmcfg.mm_count > 0) {
        if (load_metabolic_models3D(mmcfg, str_inputDir, reaction_type, num_of_microbes) != 0) return -1;
    }
    if (mmcfg.anyEnabled()) {
        if (setup_metabolic_solvers(mmcfg, reaction_type, num_of_microbes, pyFileName, src_path) != 0) return -1;
    }

    // ---- [NEW] <surrogate>: load the network, training it here if it is missing ----
    //   srgNet must outlive the simulation, because complab_srg::registerNetwork() below
    //   stores a pointer to it that defineSurrogateModel() follows at every voxel.  It is
    //   declared in main()'s own scope for exactly that reason -- do not move it into the
    //   block.
    complab_srg::Network srgNet;
#ifdef COMPLAB_ENABLE_GLPK
    complab_srgtrain::GlpkContext srgCtx;
    complab_srgtrain::StandaloneLp srgLp;      // released automatically after training
#endif
    if (icfg.srgEnabled) {
        std::string slog;
        complab_srg::FbaFn fba = 0;
        void *fctx = 0;

#ifdef COMPLAB_ENABLE_GLPK
        //   Training needs a solver.  Two ways to get one:
        //
        //   a) the run already has a GLPK microbe.  Use ITS problem -- the same persistent
        //      glp_prob the simulation will solve, so the network is fitted to exactly the
        //      linear program that would otherwise run, <constraint_indices> and all.
        //
        //   b) it does not, which is the usual case for a surrogate-only run.  Build a
        //      throwaway problem from the surrogate microbe's own model and release it as
        //      soon as training is done.
        //   Only go to the trouble of building an LP if there is nothing to load.
        //   Without this test a run whose weights file already exists still parsed a
        //   genome-scale model and built a linear program it then threw away.
        bool haveWeights = false;
        {
            std::ifstream wf(icfg.srgWeights.c_str());
            haveWeights = wf.good();
        }

        if (icfg.srgTrainIfMissing && !haveWeights) {
            plint donor = -1;
            for (plint iM = 0; iM < num_of_microbes; ++iM)
                if (rxntype::usesGlpk(reaction_type[iM]) && mmcfg.vec_lp[iM] != 0) { donor = iM; break; }

            if (donor >= 0) {
                srgCtx.cfg = &mmcfg;
                srgCtx.microbe = donor;
                pcout << "  [SRG] training will sweep microbe" << donor
                      << "'s metabolic model, the one this run already solves\n";
            } else {
                //   Which microbe are we training FOR?  <microbe> if the user said, otherwise
                //   the first one whose reaction type is a surrogate type.
                plint owner = icfg.srgMicrobe;
                if (owner < 0)
                    for (plint iM = 0; iM < num_of_microbes; ++iM)
                        if (rxntype::usesSurrogate(reaction_type[iM])) { owner = iM; break; }

                if (owner < 0 || owner >= num_of_microbes) {
                    pcout << "  [SRG] <train_if_missing> is on but no microbe has "
                          << "<reaction_type>surrogate</reaction_type>,\n"
                          << "  [SRG] so there is nothing to train for. Set one, or name the microbe\n"
                          << "  [SRG] with <surrogate><microbe>N</microbe>.\n";
                    return -1;
                }

                pcout << "  [SRG] no GLPK microbe in this run; building a temporary linear program\n"
                      << "  [SRG] from microbe" << owner << "'s own model, for training only.\n";
                const std::string lerr = complab_srgtrain::buildTrainingLp(srgLp, mmcfg, owner,
                                                                           str_inputDir,
                                                                           num_of_substrates);
                if (!lerr.empty()) { pcout << "  [SRG] " << lerr << "\n"; return -1; }
                srgCtx.cfg = &srgLp.cfg;
                srgCtx.microbe = 0;
            }

            const std::string berr = complab_srgtrain::bindInputs(srgCtx, icfg.srgTrain.inputs,
                                                                  vec_subs_names);
            if (!berr.empty()) { pcout << "  [SRG] " << berr << "\n"; return -1; }
            fba  = &complab_srgtrain::solveWithUptake;
            fctx = &srgCtx;
        }
#endif

        const bool sok = integ::prepareSurrogate(icfg, srgNet, fba, fctx,
                                                 global::mpi().isMainProcessor(), slog);
        pcout << slog;
        if (!sok) return -1;

#ifdef PLB_MPI_PARALLEL
        //   Only rank 0 trains, so every other rank now reads the file it wrote.  Training
        //   on every rank would fit a different network per rank from a different random
        //   start, and the domain would grow at a different rate in each subdomain.
        global::mpi().barrier();

        //   [FIX] THE FAILURE HAS TO BE COLLECTIVE.
        //
        //   The first version of this returned -1 from whichever rank could not read the
        //   file. That is a deadlock: the failing ranks leave, the rest walk on into the
        //   next MPI collective in the geometry setup and block there until the wall clock
        //   kills the job. And because the message went through pcout, which prints on rank
        //   0 only, a non-master failure said nothing at all. A job that hangs silently
        //   after "written to output/..." is the worst possible way to report a missing file.
        //
        //   So: every rank reports, the answer is reduced, and either all of them continue
        //   or all of them stop.
        //
        //   The retry is for a shared filesystem. A barrier synchronises PROCESSES, not
        //   filesystem metadata: on Lustre or NFS a file another node closed a microsecond
        //   ago may not be visible yet. Three tries over three seconds costs nothing on the
        //   run that does not need it.
        {
            int rankOk = 1;
            if (!global::mpi().isMainProcessor()) {
                std::string serr;
                rankOk = 0;
                for (int attempt = 0; attempt < 3 && !rankOk; ++attempt) {
                    if (attempt) sleep(1);
                    if (complab_srg::load(srgNet, icfg.srgWeights, &serr)) rankOk = 1;
                }
                if (!rankOk)
                    std::cerr << "  [SRG] rank " << global::mpi().getRank()
                              << " could not read " << icfg.srgWeights << ": " << serr << std::endl;
            }
            int allOk = rankOk;
            global::mpi().reduceAndBcast(allOk, MPI_MIN);
            if (!allOk) {
                pcout << "  [SRG] at least one rank could not read " << icfg.srgWeights << ".\n"
                      << "  [SRG] See the .err file for which. Stopping on every rank.\n";
                return -1;
            }
        }
#endif
#ifdef COMPLAB_ENABLE_GLPK
        if (srgCtx.calls > 0)
            pcout << "  [SRG] " << srgCtx.calls << " linear program(s) solved while training, "
                  << srgCtx.infeasible << " infeasible\n";
#endif

        //   Work out which entry of the per-substrate flux vector feeds each network input.
        //   The network records the substrate names it was trained on, so this is a name
        //   match rather than an assumption about order -- see bindToSubstrates().
        std::vector<int> srgSubs;
        std::string sbwarn;
        const std::string sberr = complab_srg::bindToSubstrates(srgNet, vec_subs_names,
                                                                srgSubs, &sbwarn);
        if (!sberr.empty()) {
            pcout << "  [SRG] " << sberr << "\n";
            return -1;
        }
        if (!sbwarn.empty()) pcout << sbwarn;

        //   This is the line that makes <weights_file> do something: surrogateModel.hh asks
        //   the registry for a network before it uses its own compiled-in weights.
        complab_srg::registerNetwork((int) icfg.srgMicrobe, &srgNet, (int) num_of_microbes, srgSubs);
    }

    // ---- [NEW] <symbolic> and <graphnet>: load the learned rate laws ----------------------
    //   Same lifetime rule as srgNet above, and for the same reason: complab_sym::runtime()
    //   and complab_gnn::runtime() hold POINTERS into these four objects and follow them at
    //   every voxel of every step, so they are declared in main()'s scope and never moved.
    //
    //   Everything else -- reading the files, matching their names against
    //   <name_of_substrates>, refusing a mismatch -- is in integ::prepareLearned(), which
    //   compiles without Palabos and is tested on its own.
    complab_sym::Program symProg, symAbioticProg;
    complab_gnn::Network gnnNet, gnnAbioticNet;
    if (icfg.symEnabled || icfg.gnnEnabled) {
        //   Which organisms asked for each path.  Bound per microbe rather than to all, so
        //   that each organism's own biomass variable resolves to its own name.
        std::vector<int> symUsers, gnnUsers;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (rxntype::usesSymbolic(reaction_type[iM])) symUsers.push_back((int) iM);
            if (rxntype::usesGraphnet(reaction_type[iM])) gnnUsers.push_back((int) iM);
        }

        std::string llog;
        const bool lok = integ::prepareLearned(icfg, vec_subs_names, vec_microbes_names,
                                               (int) num_of_microbes, symUsers, gnnUsers,
                                               symProg, symAbioticProg, gnnNet, gnnAbioticNet,
                                               llog);
        pcout << llog;
        if (!lok) return -1;
    }

    //   The thermodynamic gate.  Not a rate path of its own: it multiplies whichever rate path
    //   each organism already uses, so it is prepared here, after every path is known, and read
    //   from inside all six of them.  Everything except this call -- reading the file, matching
    //   its species and organism names against <name_of_substrates> and <name_of_microbes>,
    //   refusing a mismatch -- is in integ::prepareThermo(), which is tested standalone by
    //   tests/test_thermo.cpp.
    if (icfg.thmEnabled) {
        //   The defineKinetics.hh path returns one combined rate vector per voxel, so it can
        //   carry only one gate.  prepareThermo() needs to know whether anyone is on it.
        bool kineticsInUse = false;
        for (plint iM = 0; iM < num_of_microbes; ++iM)
            if (rxntype::usesKinetics(reaction_type[iM])) kineticsInUse = true;

        std::string tlog;
        const bool tok = integ::prepareThermo(icfg, vec_subs_names, vec_microbes_names,
                                              (int) num_of_microbes, kineticsInUse, tlog);
        pcout << tlog;
        if (!tok) return -1;
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
    /* [v1.3] This used to be strcat(strdup(str_inputDir.c_str()), ns_filename), which appends
     * ns_filename past the end of an allocation sized for str_inputDir alone, and then line
     * ~1090 appends ".chk" past the end of that. Two heap overflows on every run, silent with
     * short paths and a malloc abort somewhere unrelated with long ones. A std::string owns its
     * own growth, so build the name here and hand the checkpoint reader a c_str(). */
    const std::string ns_read_base = str_inputDir + std::string(ns_filename ? ns_filename : "");
    const std::string ns_read_chk  = ns_read_base + ".chk";

    // ---- [NEW] generate or import the pore space, if the XML asks for it ---------
    //   <generate> builds a packing, a fracture or a layered medium in place; <import_raw>
    //   thresholds a binary volume from imaging.  Either way a .dat is written next to the
    //   other inputs, so the run is reproducible from its own output and the generated
    //   geometry can be inspected with the same tools as any other.
    //
    //   With neither tag present provideGeometry() returns an empty string and the file
    //   named by <filename> is read exactly as before.
    {
        std::string glog;
        bool gfatal = false;

        //   THE GEOMETRY FILE IS nx-2 SLICES WIDE, NOT nx.
        //   complab_functions.hh does `nx += 2` when it reads <nx>, and readGeometry() below
        //   reads x = 1 .. nx-2 from the file and DUPLICATES the first and last slice into the
        //   two ghost columns.  So the file the generator writes must have the width the user
        //   asked for in the XML, which is nx-2 here.
        //
        //   Getting this wrong does not fail: readGeometry() would simply stop after nx-2
        //   slices and quietly ignore the rest of a too-wide file, leaving the run with a
        //   truncated pore space that still looks plausible.  The first real build caught it
        //   as a four-voxel disagreement between the porosity the generator reported and the
        //   porosity the run counted.
        const std::string gpath = integ::provideGeometry(icfg, (int) (nx - 2), (int) ny, (int) nz,
                                                         str_inputDir, glog, gfatal);
        pcout << glog;
        if (gfatal) {
            //   A sealed pore space is not a well-posed problem: a pressure drop across it has
            //   no solution, and the flow solver would spend its whole iteration budget failing
            //   to converge to one.  Better to say so now.
            pcout << "  [GEOM] the geometry is unusable; stopping before the flow solver.\n";
            return -1;
        }
        if (!gpath.empty() && gpath.size() > str_inputDir.size())
            geom_filename = gpath.substr(str_inputDir.size());
    }

    pcout << "  [GEOM] Reading " << geom_filename << "...\n";
    MultiScalarField3D<int> geometry(nx,ny,nz);
    readGeometry(str_inputDir+geom_filename, geometry);

    // ════════════════════════════════════════════════════════════════════════════
    // THE FOUR FACES NOBODY GAVE A BOUNDARY CONDITION
    // ════════════════════════════════════════════════════════════════════════════
    //   west and east get an inlet and an outlet.  y=0, y=ny-1, z=0 and z=nz-1 get
    //   nothing at all -- so wherever the geometry leaves them open, the ADE lattices
    //   stream off the block and read back an envelope nothing updates.  Every example
    //   in this repository leaves at least one of them open, and nothing said so.
    //   complab3d_outerfaces.hh has the measurement, the closure it offers, and the
    //   reason the closure is not the default.  Done here, on the geometry field,
    //   before any lattice is built from it, so that every lattice, the porosity
    //   count and the mask see one consistent domain.
    {
        std::string faceMode = "open";
        try {
            XMLreader fdoc(complab_input::configPath());
            fdoc["parameters"]["LB_numerics"]["domain"]["outer_faces"].read(faceMode);
        } catch (PlbIOException &) { faceMode = "open"; }
        if (!complab_faces::apply(geometry, nx, ny, nz, bounce_back, no_dynamics, faceMode))
            return -1;
    }

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
            try { loadBinaryBlock(nsLattice, ns_read_chk); }
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
        /* [FIX] This line used to print "Converged" whether the loop had converged or run out
         * of iterations, so the one thing a reader checks first in the log could not be
         * believed.  A flow field that never settled is coupled into every solute below. */
        if (iT >= ns_maxiTer_1)
            pcout << "  [NS] WARNING: stopped at the iteration cap ns_max_iT1=" << ns_maxiTer_1
                  << " WITHOUT converging. The velocity field below is not steady; raise the cap\n"
                  << "       or loosen ns_converge_iT1 before believing anything downstream.\n";
        else
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
    /* [FIX] Every solute's relaxation time is scaled against substrate 0's pore diffusivity, so
     * a zero there divides the whole transport setup by zero.  It is an easy mistake to make:
     * an immobile species is declared with <in_pore>0.</in_pore>, and putting one in
     * <substrate0> rather than further down the list used to give ade_dt = inf, every
     * rate * dt increment inf or NaN, and a silent all-NaN field -- the tau screen below tests
     * `< TAU_REJECT` and `> 2.0`, both false for NaN, so nothing caught it. */
    if (!(vec_solute_poreD[0] > 0)) {
        pcout << "  [ADE] ERROR: substrate 0 has a pore diffusivity of " << vec_solute_poreD[0]
              << ". Every other solute's relaxation time is scaled against it, so it must be\n"
              << "        positive. If substrate 0 is meant to be immobile (a mineral, say),\n"
              << "        move it further down <name_of_substrates> and put a mobile species\n"
              << "        first. Terminating.\n";
        return -1;
    }
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
            /* [v1.3] The divisor was vec_bMass_poreD[iM], not vec_solute_poreD[0].
             *
             * The unit system is fixed once, four lines above the substrate loop:
             * ade_dt = refNu * dx^2 / vec_solute_poreD[0], so dt/dx^2 = refNu/vec_solute_poreD[0]
             * and the lattice viscosity of ANY physical D is refNu * D / vec_solute_poreD[0].
             * That is what the substrate-in-biofilm line does, and what bioNUinPore does two
             * lines up. Dividing by the microbe's own pore diffusivity instead applies the
             * biomass-to-solute conversion a second time, as a ratio to itself, and is right
             * only by accident when vec_bMass_poreD[iM] == vec_solute_poreD[0].
             *
             * Shipped example 07 (biomass 1e-10 both places, substrate0 5e-10) got
             * tau_pore = 0.56 but tau_biofilm = 0.8: biomass diffusing five times too fast
             * inside the biofilm, which is exactly where biofilm biomass lives. Example 08's
             * LBM microbe had the same factor. The tau screen below could not catch it because
             * it recomputes bioTAUinbFilm from this same expression. */
            bioOMEGAinbFilm[iM] = 1/(refNu*vec_bMass_bFilmD[iM]/vec_solute_poreD[0]*RXNDES<T>::invCs2+0.5);
            bioTAUinbFilm[iM] = 1/bioOMEGAinbFilm[iM];
        }
        else { bioOMEGAinbFilm[iM] = 0.; bioTAUinbFilm[iM] = 0.; }
    }

    pcout << "  [ADE] dt=" << ade_dt << " s/iter, total=" << ade_maxiTer*ade_dt << " s\n";

    // ════════════════════════════════════════════════════════════════════════════
    // EVERY RELAXATION TIME THAT WILL ACTUALLY RELAX, CHECKED BEFORE IT RUNS
    // ════════════════════════════════════════════════════════════════════════════
    //   One diffusion coefficient in CompLaB.xml is the reference; every other one
    //   is carried onto the lattice as a ratio to it, and each ratio becomes its own
    //   relaxation time tau = D/D_ref * (tau_ref - 0.5) + 0.5. The check above tests
    //   only tau_ref. A species a few thousand times less mobile than the reference
    //   lands within a ten-thousandth of 0.5, and BGK there does not diffuse slowly:
    //   it rings.
    //
    //   Measured on example 07 -- one biomass patch, every reaction switched off, a
    //   fully walled box, so the only correct answer is that the patch spreads and
    //   the total never changes:
    //
    //        tau        worst negative biomass, as a fraction of the peak
    //        0.50018        25 %          <- what this example shipped with
    //        0.510         1.6 %
    //        0.520        0.33 %
    //        0.550        0.02 %, and none at all by iteration 1000
    //        0.800           0
    //
    //   Hence the two thresholds below. They are the measurement, not a convention.
    //
    //   Only lattices that actually collide and stream are checked. An immobile
    //   species has D = 0 and tau exactly 0.5, and never streams, so its tau means
    //   nothing; the same is true of a biomass field whose <solver_type> is CA or
    //   FD, which is transported by its own solver and not by this lattice. Checking
    //   those would reject examples 13, 14 and 15 for a number that is never used.
    const T TAU_REJECT = 0.51;   // below this the ringing is a large part of the signal
    const T TAU_WARN   = 0.55;   // below this it is visible but small
    {
        std::vector<std::string> label;
        std::vector<T>           tauv;
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            if (vec_immobile[iS]) continue;               // never streams: tau is not used
            label.push_back("substrate" + std::to_string(iS) + " in pore");
            tauv .push_back(substrTAUinPore[iS]);
            label.push_back("substrate" + std::to_string(iS) + " in biofilm");
            tauv .push_back(1 / substrOMEGAinbFilm[iS]);
        }
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            //   biofilm biomass relaxes only under <solver_type>LBM; free (planktonic)
            //   biomass is collided and streamed whatever its solver says.
            const bool streams = (solver_type[iM] == 3) || (bmass_type[iM] != 1);
            if (!streams || bioTAUinPore[iM] <= 0) continue;
            label.push_back("microbe" + std::to_string(iM) + " biomass in pore");
            tauv .push_back(bioTAUinPore[iM]);
            if (bioTAUinbFilm[iM] > 0) {
                label.push_back("microbe" + std::to_string(iM) + " biomass in biofilm");
                tauv .push_back(bioTAUinbFilm[iM]);
            }
        }

        std::string tooSmall, marginal, tooLarge;
        for (size_t k = 0; k < tauv.size(); ++k) {
            std::ostringstream os;
            os << "\n           " << label[k] << ": tau = " << tauv[k];
            if      (tauv[k] <  TAU_REJECT) tooSmall += os.str();
            else if (tauv[k] <  TAU_WARN)   marginal += os.str();
            else if (tauv[k] >  2.0)        tooLarge += os.str();
        }
        if (!tooLarge.empty())
            pcout << "  [ADE] WARNING: a relaxation time above 2 is over-diffusive and inaccurate:"
                  << tooLarge << "\n         Lower <tau>, or raise the reference diffusion "
                  << "coefficient this one is scaled against.\n";
        if (!marginal.empty())
            pcout << "  [ADE] NOTE: relaxation times close to 0.5 ring. These are inside the "
                  << "usable range but not\n         comfortably so, and small negative values "
                  << "may appear in the field:" << marginal << "\n";
        if (!tooSmall.empty()) {
            pcout << "  [ADE] ERROR: a lattice-Boltzmann relaxation time below 0.51"
                  << " does not diffuse slowly, it\n         oscillates, and the field goes "
                  << "negative by a large fraction of its own peak:" << tooSmall << "\n\n"
                  << "         tau = (D / D_reference) * (tau_reference - 0.5) + 0.5, so a species "
                  << "thousands of times\n"
                  << "         less mobile than the reference cannot be carried on this lattice at "
                  << "this time step.\n"
                  << "         Three ways out, in the order they are usually right:\n"
                  << "           1. transport that field with <solver_type>CA</solver_type> or "
                  << "<solver_type>FD</solver_type>,\n"
                  << "              which have no relaxation time -- this is the answer for slow "
                  << "biomass;\n"
                  << "           2. raise <tau> so every ratio lands further from 0.5 (this also "
                  << "raises the time step);\n"
                  << "           3. give the species a diffusion coefficient the lattice can "
                  << "represent, and say in the\n"
                  << "              case notes that the number was chosen for the lattice rather "
                  << "than measured.\n"
                  << "         Terminating.\n";
            return -1;
        }
    }

    // Create substrate lattices
    pcout << "  [ADE] Creating " << num_of_substrates << " substrate lattices...\n";
    MultiBlockLattice3D<T,RXNDES> substrLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T,RXNDES>(refOmega));
    std::vector< MultiBlockLattice3D<T,RXNDES> > vec_substr_lattices(num_of_substrates, substrLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > dC(num_of_substrates, substrLattice);
    std::vector< MultiBlockLattice3D<T,RXNDES> > dC0(num_of_substrates, substrLattice);
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        soluteDomainSetup(vec_substr_lattices[iS], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry,
                          substrOMEGAinbFilm[iS], substrOMEGAinPore[iS], pore_dynamics, bounce_back, no_dynamics, bio_dynamics,
                          vec_c0[iS], vec_left_btype[iS], vec_right_btype[iS], vec_left_bcondition[iS], vec_right_bcondition[iS]);
        soluteDomainSetup(dC[iS], createLocalAdvectionDiffusionBoundaryCondition3D<T,RXNDES>(), geometry,
                          substrOMEGAinbFilm[iS], substrOMEGAinPore[iS], pore_dynamics, bounce_back, no_dynamics, bio_dynamics,
                          0., vec_left_btype[iS], vec_right_btype[iS], vec_left_bcondition[iS], vec_right_bcondition[iS]);
    }
    dC0=dC;

    /* [v1.3.1] Somewhere to keep the reaction rate so it can be written as a field.
     *
     *  dC[] is the increment lattice every rate path writes into, but it is RESET
     *  twice per step -- once before the biotic block and again before the abiotic
     *  one -- so at no single moment does it hold the whole step's reaction. These
     *  accumulate across both, are zeroed once per step, and are written on the VTI
     *  interval. One scalar field per substrate: 8 bytes a voxel, which is a fifth
     *  of what the D3Q7 lattice beside it already costs. */
    std::vector< MultiScalarField3D<T> > rateField(num_of_substrates,
                                                   MultiScalarField3D<T>(nx, ny, nz, (T)0.));

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

    // ════════════════════════════════════════════════════════════════════════════
    // BOUNCE-BACK NODES START AT THE BACKGROUND CONCENTRATION, NOT AT ONE
    // ════════════════════════════════════════════════════════════════════════════
    //   Palabos stores an advection-diffusion population as a DEVIATION from
    //   equilibrium at density one, so an all-zero cell reads back as C = 1, not
    //   C = 0.  initializeAtEquilibrium() cannot fix that on a bounce-back node,
    //   because BounceBack::computeEquilibrium returns zero whatever density it is
    //   handed: every wall voxel therefore starts at C = 1 mol/L.
    //
    //   Bounce-back conserves mass -- it swaps opposite populations -- so nothing
    //   ever removes that. Each wall voxel streams C = 1 fluid into its neighbours
    //   from the first step, and in a case whose real concentrations are
    //   millimolar that is a source three orders of magnitude larger than the
    //   chemistry. It shows up as every field rising several-fold over the first
    //   hundred steps and then relaxing, which looks like a transient and is not:
    //   it is the walls filling the domain.
    //
    //   The existing 10000-step "stabilization" loop further down was written
    //   against this symptom. It resets the FLUID to c0 afterwards but leaves the
    //   walls alone, so the source is still there when the run proper starts.
    //
    //   Three lines, applied once, before anything streams: set every bounce-back
    //   voxel to the same background its neighbours hold. From then on the swap
    //   keeps it consistent and the wall is the zero-flux boundary it was meant
    //   to be.
    //
    //   Diagnosed on examples/19_thermodynamic_gate: with bounce-back walls the
    //   closed-domain methane total rose 10.75 -> 65.13 with no reaction running;
    //   with the walls removed it held to five digits. This also explains the
    //   field divergence noted in examples 17 and 18.
    if (bounce_back >= 0 || no_dynamics >= 0) {
        /* Both kinds of non-fluid voxel need this. The walls always did; the solid grains need it
         * from the moment they became bounce-back a few lines below, for exactly the same reason. */
        std::vector<plint> wallOnly;                 /* inert walls: every lattice */
        if (bounce_back >= 0) wallOnly.push_back(bounce_back);
        std::vector<plint> wallAndSolid = wallOnly;  /* walls AND grains: mobile species only */
        if (no_dynamics  >= 0) wallAndSolid.push_back(no_dynamics);
        std::vector< std::vector<plint> > noBio;
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            /* THE CONCENTRATION lattice, in every non-fluid voxel including the grains.
             *
             * Safe for an immobile species too, even though a grain is exactly where such a
             * species keeps its inventory: the dissolution setup seeds that inventory further
             * down this function, so it writes last and wins. Leaving the grains out instead --
             * to protect an inventory that has not been written yet -- leaves them holding the
             * all-zero populations initializeAtEquilibrium could not reach, which an
             * advection-diffusion lattice reads back as 1 mol/L. Example 13's FeS starts at zero
             * and reported a total of 72 before the first step, which is 1.0 in each of its 72
             * grain voxels. */
            if (!wallAndSolid.empty())
                applyProcessingFunctional(new stabilizeADElattice3D<T,RXNDES,int>(vec_c0[iS], wallAndSolid, noBio),
                                          vec_substr_lattices[iS].getBoundingBox(),
                                          vec_substr_lattices[iS], geometry);

            /* THE INCREMENT lattices, both of them, in EVERY non-fluid voxel including the grains
             * -- immobile or not. An increment is not an inventory: it is zero at the start of
             * every step by definition, and there is nothing to protect.
             *
             * This is not cosmetic. initializeAtEquilibrium() cannot reach a cell whose dynamics
             * returns a zero equilibrium, and both NoDynamics and BounceBack do, so those cells
             * were left holding all-zero populations -- which an advection-diffusion lattice reads
             * back as density ONE, not zero. The increment applier then added +1 mol/L to every
             * grain voxel on every step. Measured on example 14 before this line existed: the
             * calcite inventory grew from 2023 to 131575 over 1800 steps, 72 per step, which is
             * exactly the +1 times the 72 grain voxels in that geometry. A dissolving mineral was
             * accumulating. */
            if (!wallAndSolid.empty()) {
                applyProcessingFunctional(new stabilizeADElattice3D<T,RXNDES,int>((T) 0., wallAndSolid, noBio),
                                          dC[iS].getBoundingBox(), dC[iS], geometry);
                applyProcessingFunctional(new stabilizeADElattice3D<T,RXNDES,int>((T) 0., wallAndSolid, noBio),
                                          dC0[iS].getBoundingBox(), dC0[iS], geometry);
            }
        }
        //   Biomass is zero outside its patches, so its walls belong at zero too.
        for (size_t iM = 0; iM < vec_bFilm_lattices.size(); ++iM) {
            applyProcessingFunctional(new stabilizeADElattice3D<T,RXNDES,int>((T) 0., wallOnly, noBio),
                                      vec_bFilm_lattices[iM].getBoundingBox(), vec_bFilm_lattices[iM], geometry);
            applyProcessingFunctional(new stabilizeADElattice3D<T,RXNDES,int>((T) 0., wallOnly, noBio),
                                      dBf[iM].getBoundingBox(), dBf[iM], geometry);
        }
        for (size_t iM = 0; iM < vec_bFree_lattices.size(); ++iM) {
            applyProcessingFunctional(new stabilizeADElattice3D<T,RXNDES,int>((T) 0., wallOnly, noBio),
                                      vec_bFree_lattices[iM].getBoundingBox(), vec_bFree_lattices[iM], geometry);
            applyProcessingFunctional(new stabilizeADElattice3D<T,RXNDES,int>((T) 0., wallOnly, noBio),
                                      dBp[iM].getBoundingBox(), dBp[iM], geometry);
        }
        pcout << "  [ADE] Bounce-back voxels set to the background concentration.\n";
    }

    // ════════════════════════════════════════════════════════════════════════════
    // AN IMMOBILE SPECIES INSIDE A SOLID VOXEL MUST REPORT WHAT IS IN IT
    // ════════════════════════════════════════════════════════════════════════════
    //   Solid voxels are given Palabos's NoDynamics, and NoDynamics::computeDensity
    //   returns its own stored rho -- 1.0 by default -- whatever the populations
    //   hold. That is the same mistake as the bounce-back one above, in a second
    //   place, and here it is worse: an IMMOBILE species exists precisely to sit in
    //   a solid voxel and hold an inventory there.
    //
    //   So every mineral inventory read back as 1.0 mol/L no matter what was seeded
    //   into it. Measured on example 14: 27.1 mol/L of calcite is seeded, the
    //   reopening test reads 1.0, finds it below 0.9 x 27.1, and converts every
    //   grain to pore on iteration ZERO -- before a single molecule has dissolved.
    //   The case then had no mineral left, dissolved nothing for the rest of the
    //   run, and reported concentrations from voxels that had been solid a moment
    //   earlier. The whole dissolution example was measuring nothing.
    //
    //   The fix is to give those voxels a dynamics whose computeDensity sums the
    //   populations. Omega is irrelevant: an immobile species never collides and
    //   never streams -- complab.cpp skips collideAndStream for it -- so the
    //   dynamics object is consulted for nothing else.
    //
    //   Only immobile species are changed. A mobile species has no business holding
    //   anything inside a solid voxel, and NoDynamics is the right choice there.
    if (no_dynamics >= 0) {
        plint nimm2 = 0, nmob2 = 0;
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            if (vec_immobile[iS]) {
                ++nimm2;
                defineDynamics(vec_substr_lattices[iS], geometry,
                               new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.), no_dynamics);
                defineDynamics(dC[iS],  geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.), no_dynamics);
                defineDynamics(dC0[iS], geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.), no_dynamics);
                continue;
            }

            // ────────────────────────────────────────────────────────────────
            // A SOLID GRAIN IS A NO-FLUX BOUNDARY, NOT A HOLE IN THE LATTICE
            // ────────────────────────────────────────────────────────────────
            //   Solid voxels were given Palabos's NoDynamics on every substrate
            //   lattice. NoDynamics does not collide -- but collideAndStream
            //   STREAMS the whole lattice regardless, so those cells hand their
            //   populations to their fluid neighbours every step and never get
            //   any back. They are a one-way drain.
            //
            //   Worse, of the wrong sign. Palabos stores an advection-diffusion
            //   population as a deviation from equilibrium at density one, so a
            //   cell holding concentration c stores populations summing to c-1
            //   -- close to MINUS ONE for any millimolar chemistry. Every solid
            //   voxel therefore poured a large negative deviation into the water
            //   around it, step after step.
            //
            //   Measured on example 14, with every reaction switched off and the
            //   only thing running being transport: Ca starts at 0 in a closed
            //   domain, so it cannot physically change, and it reached
            //   -1.00 mol/L. H, fed at 0.01, reached -1.38. Replacing the grains
            //   with pore in the same geometry gives exactly 0.000000 for both,
            //   which is what identified this.
            //
            //   A dissolved species meets a mineral grain at a no-flux boundary,
            //   and the lattice-Boltzmann spelling of no-flux is bounce-back:
            //   what arrives is reflected, nothing is created or destroyed. That
            //   is already what the inert-wall material gets. Solid voxels now
            //   get it too, and are seeded with the background concentration by
            //   the block above for the same reason the walls are.
            ++nmob2;
            defineDynamics(vec_substr_lattices[iS], geometry, new BounceBack<T,RXNDES>(), no_dynamics);

            /* The INCREMENT lattices keep a dynamics that reads its populations back.
             *
             * Bounce-back is the right answer for the concentration lattice, which streams: it
             * makes a grain a no-flux boundary. The increment lattices never stream -- they are
             * accumulators, zeroed at the start of every step and applied by
             * update_*_rxnLattices -- so no-flux means nothing to them, and bounce-back would
             * actively hurt: BounceBack::computeDensity returns its own stored density and
             * ignores the populations, so anything written into a grain's increment slot would
             * read back as zero.
             *
             * Something does write there. surfaceDissolutionKinetics3D parks each product share
             * in the mineral voxel's own increment slot for dissolutionGather3D to collect, and
             * with bounce-back on this lattice every one of 816000 parked shares read back as
             * zero and the water collected nothing. */
            defineDynamics(dC[iS],  geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.), no_dynamics);
            defineDynamics(dC0[iS], geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.), no_dynamics);
        }
        if (nimm2 > 0)
            pcout << "  [ADE] " << nimm2 << " immobile species: solid voxels now report their\n"
                  << "  [ADE] own contents rather than a fixed 1.0 (see the note in complab.cpp).\n";
        if (nmob2 > 0)
            pcout << "  [ADE] " << nmob2 << " mobile species: solid voxels are now a no-flux\n"
                  << "  [ADE] boundary rather than a one-way drain (see the note in complab.cpp).\n";
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

    // ---- OPTIONAL METABOLIC LAYER: one lattice vector per solver ---------------
    //   Same layout as the kinetics vector, but each holds only the microbes that
    //   use that particular solver.  glpk_globalId / cpy_globalId / srg_globalId
    //   translate a position in the compacted biomass block back to the microbe's
    //   number in CompLaB.xml, so per-microbe parameters are never indexed by
    //   position.
    //
    //   The two FBA back ends get SEPARATE lists, not one shared FBA list.  This is
    //   what lets a GLPK organism and a COBRApy organism live in the same run.  The
    //   older code built one list over usesFBA() and handed it to both processors,
    //   so each would have tried to solve for the other's organisms -- for a COBRApy
    //   microbe, cfg.vec_lp[gM] is null -- and that, not any flux-index convention,
    //   is what the ban in complab3d_metabolic.hh was really protecting against.
    //   With one list per back end each processor only ever sees its own organisms,
    //   and the two dispatch into the shared dC/dB increments exactly the way the
    //   surrogate, symbolic and graphnet paths already do alongside each other.
    std::vector< MultiBlockLattice3D<T, RXNDES>* > ptr_glpk_lattices, ptr_cpy_lattices,
                                                   ptr_srg_lattices,
                                                   ptr_sym_lattices, ptr_gnn_lattices;
    std::vector<plint> glpk_globalId, cpy_globalId, srg_globalId, sym_globalId,
                       gnn_globalId, mm_modelSlot;
    {
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_glpk_lattices.push_back(&vec_substr_lattices[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesGlpk(reaction_type[iM])) continue;
            glpk_globalId.push_back(iM);
            if (bmass_type[iM]==1) ptr_glpk_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
            else                   ptr_glpk_lattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
        }
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_glpk_lattices.push_back(&dC[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesGlpk(reaction_type[iM])) continue;
            if (bmass_type[iM]==1) ptr_glpk_lattices.push_back(&dBf[loctrack[iM]]);
            else                   ptr_glpk_lattices.push_back(&dBp[loctrack[iM]]);
        }
        ptr_glpk_lattices.push_back(&maskLattice);

        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_cpy_lattices.push_back(&vec_substr_lattices[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesCobrapy(reaction_type[iM])) continue;
            cpy_globalId.push_back(iM);
            if (bmass_type[iM]==1) ptr_cpy_lattices.push_back(&vec_bFilm_lattices[loctrack[iM]]);
            else                   ptr_cpy_lattices.push_back(&vec_bFree_lattices[loctrack[iM]]);
        }
        for (plint iS = 0; iS < num_of_substrates; ++iS) ptr_cpy_lattices.push_back(&dC[iS]);
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (!rxntype::usesCobrapy(reaction_type[iM])) continue;
            if (bmass_type[iM]==1) ptr_cpy_lattices.push_back(&dBf[loctrack[iM]]);
            else                   ptr_cpy_lattices.push_back(&dBp[loctrack[iM]]);
        }
        ptr_cpy_lattices.push_back(&maskLattice);

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

        // The two learned paths use the identical layout -- substrates, their biomass
        // lattices, the dC increments, the dB increments, the mask -- because
        // run_symbolic3D and run_graphnet3D are run_surrogate3D with the evaluator
        // swapped.  Written out longhand rather than folded into a loop, so that these
        // read the same as the two blocks above them: the layout IS the interface, and
        // the place it is built should be the easiest thing in this file to check.
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
        if ((plint) ptr_glpk_lattices.size() != 2*(num_of_substrates + (plint) glpk_globalId.size()) + 1 ||
            (plint) ptr_cpy_lattices.size()  != 2*(num_of_substrates + (plint) cpy_globalId.size())  + 1 ||
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

    /* [FIX] updateLocalMaskNtotalLattices3D loops over the `bio` vector it is handed and reads
     * lattices[iM] for each entry.  It used to be handed bio_dynamics, which has one row per
     * BIOFILM microbe -- but ptr_ca_lattices and ptr_fd_mask carry only the microbes on THAT
     * solver.  In a run where the two counts differ, and shipped example 08 is exactly such a
     * run (microbe0 on the CA, microbe1 on the LBM, both biofilm), the loop walked past the
     * biomass block into the copy lattices: in-transit biomass was counted twice, the second
     * microbe's biomass was never summed into totalbFilmLattice at all, and with enough
     * microbes it would have indexed past the end of the vector.
     *
     * These two vectors hold the material numbers of the microbes each list actually contains,
     * in the same order the lattices were pushed, so the loop and the lattice list agree. */
    std::vector< std::vector<plint> > bio_ca, bio_fd;
    for (plint iM = 0; iM < num_of_microbes; ++iM) {
        if (bmass_type[iM] != 1) continue;                    /* biofilm rows only, as before */
        if (solver_type[iM] == 2) bio_ca.push_back(bio_dynamics[loctrack[iM]]);
        else if (solver_type[iM] == 1) bio_fd.push_back(bio_dynamics[loctrack[iM]]);
    }

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
        applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, caLlen, bounce_back, no_dynamics, bio_ca, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
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
        // [FIX-3D] SEGFAULT.  This loop indexed vec_bFree_lattices with a counter that
        //   advanced once per LBM microbe, but that vector is sized bfree_count and holds
        //   only the PLANKTONIC microbes.  A run with a biofilm organism whose
        //   <solver_type> is LBM -- example 08 is exactly that -- indexed past the end of an
        //   often EMPTY vector and crashed here, after several minutes of flow solving.
        //
        //   Two things were wrong and both are fixed by walking the free microbes instead:
        //   the index now matches how the vector was filled (the same tmpIT0/tmpIT1 pairing
        //   used by the checkpoint loader below), and biofilm biomass is skipped, which is
        //   correct on its own terms -- attached biomass does not advect with the flow.
        tmpIT1=0;
        for (plint iM = 0; iM < num_of_microbes; ++iM) {
            if (bmass_type[iM]==1) continue;              // biofilm: fixed in place, never advected
            if (solver_type[iM] == 3)
                latticeToPassiveAdvDiff(nsLattice, vec_bFree_lattices[tmpIT1], vec_bFree_lattices[tmpIT1].getBoundingBox());
            ++tmpIT1;                                     // advances for every free microbe, coupled or not
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

    /* [v1.3.1] Whether this run has a rate worth mapping.
     *
     *  "If applicable" is the whole point: a diffusion-only case has substrates and
     *  no reaction, and writing a file of zeros for it every interval would double
     *  the output volume to say nothing. So the field is written when at least one
     *  rate path is actually running -- biotic or abiotic, hand-written or learned,
     *  a linear program or a dissolving mineral. All of them accumulate into the
     *  same dC[] lattices, which is exactly why one accumulator catches every one
     *  of them and no path had to be modified to be included.
     *
     *  <track_performance> already suppresses VTI output altogether; this follows it. */
    bool rateFieldsOn = (num_of_substrates > 0) && (track_performance == 0) &&
        (   (enable_kinetics && kns_count > 0)
         || enable_abiotic_kinetics
         || dissolCfg.enabled
         || complab_sym::haveAbiotic() || complab_gnn::haveAbiotic()
         || (mmcfg.enable_symbolic && mmcfg.sym_count > 0)
         || (mmcfg.enable_graphnet && mmcfg.gnn_count > 0)
         || (mmcfg.enable_surrogate && mmcfg.srg_count > 0)
#ifdef COMPLAB_ENABLE_GLPK
         || (mmcfg.enable_fba_glpk && mmcfg.glpk_count > 0)
#endif
#ifdef COMPLAB_ENABLE_COBRAPY
         || (mmcfg.enable_fba_cobrapy && mmcfg.cpy_count > 0)
#endif
        );

    /* Add whatever is sitting in dC[] right now into the step's accumulator. Called
     * after each of the three apply processors, because dC[] is reset between them
     * and the file is meant to carry the whole step's reaction, not its last third. */
    auto accumulateRates = [&]() {
        if (!rateFieldsOn) return;
        for (plint iS = 0; iS < num_of_substrates; ++iS)
            addInPlace(rateField[iS], *computeDensity(dC[iS], rateField[iS].getBoundingBox()));
    };

    if (rateFieldsOn) {
        pcout << "│ Rate fields ON: rate_<species>_*.vti, mol/L/s, one per VTI interval\n";
        pcout << "│   positive = produced, negative = consumed, same box as every other field\n";
    }

    /* [v1.3.2] THE ENERGY FIELDS.
     *
     *  Allocated only for organisms that have a block in the .thm file, and only when VTI
     *  output is on at all.  Two fields per gated organism, recomputed from the standing
     *  concentrations at each write, so the steps in between cost nothing.  See
     *  complab3d_energyfield.hh for why this is a snapshot and the rate field is not. */
    const std::vector<int> gatedMic = complab_efield::gatedMicrobes((int) num_of_microbes);
    const bool energyFieldsOn = !gatedMic.empty() && (track_performance == 0)
                             && (num_of_substrates > 0) && (ade_VTI_iTer > 0);

    std::vector< MultiScalarField3D<T> > dGField(
        energyFieldsOn ? gatedMic.size() : 0, MultiScalarField3D<T>(nx, ny, nz, (T)0.));
    std::vector< MultiScalarField3D<T> > ftField(
        energyFieldsOn ? gatedMic.size() : 0, MultiScalarField3D<T>(nx, ny, nz, (T)0.));

    /* One block list per gated organism: [C.., mask, dG, F_T], the order
     * computeEnergyField3D::processGenericBlocks() expects. */
    std::vector< std::vector<MultiBlock3D*> > ptr_energy;
    if (energyFieldsOn) {
        for (size_t k = 0; k < gatedMic.size(); ++k) {
            std::vector<MultiBlock3D*> v;
            for (plint iS = 0; iS < num_of_substrates; ++iS) v.push_back(&vec_substr_lattices[iS]);
            v.push_back(&maskLattice);
            v.push_back(&dGField[k]);
            v.push_back(&ftField[k]);
            ptr_energy.push_back(v);
        }
        pcout << "│ Energy fields ON: dG_<microbe>_*.vti (kJ/mol) and FT_<microbe>_*.vti (0..1)\n";
        pcout << "│   dG is the energy the reaction yields at the LOCAL composition, negative\n";
        pcout << "│   when it yields anything; F_T is the factor the rate was multiplied by.\n";
        pcout << "│   Zero in either file means solid, bounce-back or an outer column, not a\n";
        pcout << "│   reaction sitting exactly at equilibrium.\n";
        pcout << "│   Basis: lattice concentrations scaled by <concentration_scale>. An FBA path\n";
        pcout << "│   set to a TOTAL basis gated on speciated totals, so its field and its gate\n";
        pcout << "│   agree only where complexation is weak.\n";
        for (size_t k = 0; k < gatedMic.size(); ++k)
            pcout << "│   [" << gatedMic[k] << "] " << vec_microbes_names[gatedMic[k]] << "\n";
    }

    /* Recompute both fields from the concentrations standing right now, then write them. */
    auto writeEnergyFields = [&](plint iter) {
        if (!energyFieldsOn) return;
        for (size_t k = 0; k < gatedMic.size(); ++k) {
            applyProcessingFunctional(
                new complab_efield::computeEnergyField3D<T,RXNDES>(
                        nx, num_of_substrates, gatedMic[k], no_dynamics, bounce_back),
                dGField[k].getBoundingBox(), ptr_energy[k]);
            const std::string nm = vec_microbes_names[gatedMic[k]];
            writeFieldVTI(dGField[k], iter, "dG_" + nm + "_", (T)1., "DeltaG");
            writeFieldVTI(ftField[k], iter, "FT_" + nm + "_", (T)1., "F_T");
        }
    };

    global::timer("ade").restart();
    util::ValueTracer<T> ns_convg2(1.0,1000.0,ns_converge_iT2);
    bool ns_saturate=0, percolationFlag=0;
    bool ns_warned_unconverged=false;   /* [v1.3] the unconverged-re-solve warning, once per run */

    /* The conserved biomass the run STARTS from, taken here rather than earlier because the
     * per-microbe lattices are parked at their background between the two points, and a sum
     * taken before that parking counts one unit of nothing for every wall voxel -- 2000 of
     * them on example 07, against a real biomass of 126. Taken once, immediately before the
     * first step, so the closing report compares like with like. */
    if (bfilm_count > 0 || bfree_count > 0)
        diag_initial_total_biomass = complab_total_biomass(vec_bFilm_lattices, vec_bFree_lattices);

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
                /* [v1.3.1] THE PEAK IS NOT GROWTH, and this line said it was.
                 *
                 * The closing report was corrected in v1.3 -- see the note beside
                 * finalBmax in PHASE 7 -- and this per-iteration line, which computes the
                 * same quantity from the same variable, was left behind. So example 06
                 * printed "(-3.2699% growth)" every interval for a population whose total
                 * biomass was rising the whole time: the patch spreads from 108 voxels to
                 * 876, which lowers the peak and raises the total. A reader watching the
                 * run would have concluded the organism was dying.
                 *
                 * Both numbers are printed now, each labelled for what it is, and the
                 * average comes with them -- it was already being computed and thrown
                 * away, which is also what the two unused-variable warnings on this block
                 * were about. */
                for (plint iM = 0; iM < bfilm_count; ++iM) {
                    T bMax = computeMax(*computeDensity(vec_bFilm_lattices[iM]));
                    T bAvg = computeAverage(*computeDensity(vec_bFilm_lattices[iM]));
                    T peakChange = (diag_initial_biomass > 0)
                                 ? ((bMax - diag_initial_biomass) / diag_initial_biomass * 100.0) : 0.0;
                    pcout << "║   " << vec_microbes_names[iM] << ": peak=" << std::scientific << bMax
                          << "/" << max_bMassRho << " avg=" << bAvg
                          << std::fixed << " (peak " << peakChange << "%)";
                    if (bMax > max_bMassRho) pcout << " [>Bmax!]";
                    pcout << "\n";
                }
                {
                    const T runningTotal = complab_total_biomass(vec_bFilm_lattices, vec_bFree_lattices);
                    const T totalChange = (diag_initial_total_biomass > 0)
                                        ? ((runningTotal - diag_initial_total_biomass) / diag_initial_total_biomass * 100.0) : 0.0;
                    pcout << "║   TOTAL biomass=" << std::scientific << runningTotal
                          << std::fixed << " (" << totalChange << "%, THIS is the growth)\n";
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
                /* [v1.3.2] Written HERE, beside the concentrations, and not down in the rate
                 * block. dG and F_T are functions of the concentrations at one instant, so
                 * they belong with the file that holds those concentrations: a voxel's dG in
                 * this file is exactly the number its own species files at the same iteration
                 * produce. The rate field is an interval average and is written after the
                 * reaction for the opposite reason. */
                writeEnergyFields(iT);
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

        // ════════════════════════════════════════════════════════════════════
        // WHERE CHEMISTRY IS ALLOWED TO HAPPEN
        // ════════════════════════════════════════════════════════════════════
        //   Not the whole bounding box.  readGeometry() duplicates the first and
        //   last slice of the geometry file into x = 0 and x = nx-1, and those two
        //   planes are the boundary condition: Dirichlet overwrites them every
        //   step, Neumann copies them from the layer inside, `closed` bounces off
        //   them.  They are not part of the domain, and every reduction in this
        //   program already excludes them.
        //
        //   Running reactions there is not harmless.  A grain in the duplicated
        //   plane dissolves like any other: the mineral it loses is outside every
        //   total, and the products it parks are collected by a voxel at x = 1,
        //   which is inside.  Mass appears in the water with no matching loss
        //   anywhere in the books.  On example 14 that manufactured about 0.7 of
        //   5.4 mol/L of calcium and hid about 0.3 of the calcite loss -- the last
        //   piece of a balance that was failing by 26%.
        //
        //   Same box as the diagnostics, the porosity count and saveGeometry().
        const Box3D reactionBox(1, nx-2, 0, ny-1, 0, nz-1);

        // Kinetics (biotic - only if enable_kinetics is true and biotic_mode)
        dC=dC0; dBp=dBp0; dBf=dBf0;
        /* [v1.3.1] The rate accumulator is zeroed HERE, once, and not again until the
         * next step -- unlike dC[] below it, which is reset a second time before the
         * abiotic block. That is the whole difference between the two, and the reason
         * the accumulator exists. */
        if (rateFieldsOn)
            for (plint iS = 0; iS < num_of_substrates; ++iS)
                setToConstant(rateField[iS], rateField[iS].getBoundingBox(), (T)0.);
        if (enable_kinetics && kns_count > 0) {
            if (track_performance == 1) global::timer("kns").restart();
            applyProcessingFunctional(new run_kinetics<T,RXNDES>(nx, num_of_substrates, kns_count, ade_dt, vec_Kc_kns, vec_mu_kns, no_dynamics, bounce_back),
                                      reactionBox, ptr_kns_lattices);
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
            applyProcessingFunctional(new runFBA_glpk3D<T,RXNDES>(nx, num_of_substrates, (plint) glpk_globalId.size(),
                                          ade_dt, glpk_globalId, &mmcfg, vec_Kc, vec_mu, no_dynamics, bounce_back),
                                      reactionBox, ptr_glpk_lattices);
            if (track_performance == 1) { fbatime += global::timer("fba").getTime(); global::timer("fba").stop(); }
        }
#endif
#ifdef COMPLAB_ENABLE_COBRAPY
        if (mmcfg.enable_fba_cobrapy && mmcfg.cpy_count > 0) {
            if (track_performance == 1) global::timer("fba").restart();
            applyProcessingFunctional(new runFBA_cobrapy3D<T,RXNDES>(nx, num_of_substrates, (plint) cpy_globalId.size(),
                                          ade_dt, cpy_globalId, &mmcfg, pyFileName, mm_modelSlot, vec_Kc, vec_mu, no_dynamics, bounce_back),
                                      reactionBox, ptr_cpy_lattices);
            if (track_performance == 1) { fbatime += global::timer("fba").getTime(); global::timer("fba").stop(); }
        }
#endif
        if (mmcfg.enable_surrogate && mmcfg.srg_count > 0) {
            if (track_performance == 1) global::timer("srg").restart();
            applyProcessingFunctional(new run_surrogate3D<T,RXNDES>(nx, num_of_substrates, (plint) srg_globalId.size(),
                                          num_of_microbes, ade_dt, srg_globalId, &mmcfg, vec_Kc, vec_mu,
                                          no_dynamics, bounce_back),
                                      reactionBox, ptr_srg_lattices);
            if (track_performance == 1) { srgtime += global::timer("srg").getTime(); global::timer("srg").stop(); }
        }
        //   The two learned paths.  Same slot, same dC/dB increment lattices, same shared
        //   substrate budget, so a symbolic organism and a surrogate organism in one voxel
        //   obey one rule rather than two.  Gated on the counts for the reason given below:
        //   a switch that is on with no organism behind it would sweep the whole domain
        //   applying identically zero increments.
        if (mmcfg.enable_symbolic && mmcfg.sym_count > 0) {
            if (track_performance == 1) global::timer("sym").restart();
            applyProcessingFunctional(new run_symbolic3D<T,RXNDES>(nx, num_of_substrates,
                                          (plint) sym_globalId.size(), num_of_microbes, ade_dt,
                                          sym_globalId, &mmcfg, no_dynamics, bounce_back),
                                      reactionBox, ptr_sym_lattices);
            if (track_performance == 1) { symtime += global::timer("sym").getTime(); global::timer("sym").stop(); }
        }
        if (mmcfg.enable_graphnet && mmcfg.gnn_count > 0) {
            if (track_performance == 1) global::timer("gnn").restart();
            applyProcessingFunctional(new run_graphnet3D<T,RXNDES>(nx, num_of_substrates,
                                          (plint) gnn_globalId.size(), num_of_microbes, ade_dt,
                                          gnn_globalId, &mmcfg, no_dynamics, bounce_back),
                                      reactionBox, ptr_gnn_lattices);
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
            (mmcfg.enable_symbolic    && mmcfg.sym_count  > 0) ||
            (mmcfg.enable_graphnet    && mmcfg.gnn_count  > 0)) {
            // ════════════════════════════════════════════════════════════════
            // PORE SCALE OUT, CONTINUUM SCALE IN: SAMPLED HERE AND NOWHERE ELSE
            // ════════════════════════════════════════════════════════════════
            //   The effectiveness factor is the ratio of the rate the aggregate
            //   actually achieves to the rate it would achieve if every point in
            //   it saw the bulk.  The numerator is the per-voxel reaction rate,
            //   and this is the only moment in the step when it exists: every
            //   rate path has written its increment into dC, and
            //   update_rxnLattices below is about to consume them.  Recomputing
            //   the rate anywhere else would be averaging a different model from
            //   the one that ran.
            if (icfg.upsEnabled && diag.active()) {
                const plint upsEvery = (icfg.diagInterval > 0) ? (plint) icfg.diagInterval
                                                               : (ade_VTI_iTer > 0 ? ade_VTI_iTer : 0);
                if (upsEvery > 0 && iT % upsEvery == 0)
                    complab_upscale_sample(iT, icfg, ade_dt, dx, nx, ny, nz,
                                           num_of_substrates, num_of_microbes,
                                           vec_substr_lattices, dC, maskLattice,
                                           vec_bFilm_lattices, vec_bFree_lattices,
                                           bmass_type, loctrack, pore_dynamics,
                                           vec_solute_bFilmD);
            }

            //   <upscaling><freeze_biomass>: throw the biomass increments away and keep the
            //   substrate ones. The reaction still runs, still consumes and still produces;
            //   only the catalyst is held while the concentration profile relaxes. Done here,
            //   after the upscaling sample above, so the sample sees the rate the organisms
            //   actually computed rather than a zeroed one.
            if (icfg.upsEnabled && icfg.upsFreezeBiomass) { dBp = dBp0; dBf = dBf0; }

            if (track_performance == 1) global::timer("rxn").restart();
            applyProcessingFunctional(new update_rxnLattices<T,RXNDES>(nx, num_of_substrates, num_of_microbes, no_dynamics, bounce_back),
                                      reactionBox, ptr_update_rxnLattices);
            accumulateRates();   /* [v1.3.1] the biotic share, before dC[] is reset below */
            if (track_performance == 1) { T rxntime=global::timer("rxn").getTime(); global::timer("rxn").stop(); if (kns_count>0) knstime+=rxntime; }
        }

        // Abiotic kinetics (substrate-only reactions without microbes).
        //   A learned abiotic law -- <symbolic><abiotic_file> or <graphnet><abiotic_file> --
        //   belongs in this block and not the biotic one, because it fires in every fluid
        //   voxel whether anything is alive there or not.  It therefore also has to be able
        //   to open the block on its own: a run with an abiotic .sym file and no
        //   <enable_abiotic_kinetics> is a perfectly reasonable configuration.
        const bool learnedAbiotic = complab_sym::haveAbiotic() || complab_gnn::haveAbiotic();
        if ((enable_abiotic_kinetics || learnedAbiotic) && num_of_substrates > 0) {
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
            if (enable_abiotic_kinetics) {
                if (precip_enabled && precip_surfaceOnly) {
                    applyProcessingFunctional(new surfaceAbioticKinetics3D<T,RXNDES>(nx, ny, nz, num_of_substrates, ade_dt, no_dynamics, bounce_back, 1),
                                              reactionBox, ptr_abiotic_kns_lattices);
                } else {
                    applyProcessingFunctional(new run_abiotic_kinetics<T,RXNDES>(nx, num_of_substrates, ade_dt, no_dynamics, bounce_back),
                                              reactionBox, ptr_abiotic_kns_lattices);
                }
            }
            //   The learned abiotic laws, into the SAME dC lattices, so a substrate touched
            //   by both a hand-written abiotic reaction and a fitted one gets one combined
            //   increment applied once.  Each has its own shared-budget clamp inside, so
            //   neither can drive a concentration negative on its own; the combination is
            //   bounded by update_abiotic_rxnLattices below.
            if (complab_sym::haveAbiotic()) {
                applyProcessingFunctional(new run_symbolic_abiotic3D<T,RXNDES>(
                                              nx, num_of_substrates, ade_dt, &mmcfg,
                                              no_dynamics, bounce_back),
                                          reactionBox, ptr_abiotic_kns_lattices);
            }
            if (complab_gnn::haveAbiotic()) {
                applyProcessingFunctional(new run_graphnet_abiotic3D<T,RXNDES>(
                                              nx, num_of_substrates, ade_dt, &mmcfg,
                                              no_dynamics, bounce_back),
                                          reactionBox, ptr_abiotic_kns_lattices);
            }
            // [DISSOL-VOP] Mineral dissolution, into the SAME dC lattices, before they
            //   are applied.  It runs here rather than in its own block so that a
            //   substrate both produced by dissolution and consumed by an abiotic
            //   reaction gets one combined increment, applied once.
            if (dissolCfg.enabled) {
                //   TWO passes, and the order matters. The first parks each product share in the
                //   mineral voxel's own increment slot; the second lets the water collect them.
                //   Between the two, Palabos refreshes the increment lattices' envelopes, which
                //   is what makes a share parked on one rank visible to the next -- and is why
                //   the result no longer depends on the processor count.
                applyProcessingFunctional(new surfaceDissolutionKinetics3D<T,RXNDES>(
                                              nx, ny, nz, num_of_substrates, ade_dt,
                                              no_dynamics, bounce_back, pore_dynamics, &dissolCfg),
                                          reactionBox, ptr_dissol);
                applyProcessingFunctional(new dissolutionGather3D<T,RXNDES>(
                                              nx, ny, nz, num_of_substrates,
                                              no_dynamics, bounce_back, &dissolCfg),
                                          reactionBox, ptr_dissol);
            }

            // Apply concentration changes
            applyProcessingFunctional(new update_abiotic_rxnLattices<T,RXNDES>(nx, num_of_substrates, no_dynamics, bounce_back, vec_immobile),
                                      reactionBox, ptr_abiotic_kns_lattices);
            accumulateRates();   /* [v1.3.1] the abiotic share, added on top of the biotic one */
            if (track_performance == 1) { knstime += global::timer("abiotic_kns").getTime(); global::timer("abiotic_kns").stop(); }
        }

        // [DISSOL-VOP] If the abiotic block above is switched off, dissolution still
        //   has to run somewhere.  Give it its own reset/apply pair so it works
        //   independently of <enable_abiotic_kinetics>.
        if (dissolCfg.enabled && !((enable_abiotic_kinetics || learnedAbiotic) && num_of_substrates > 0)) {
            dC = dC0;
            applyProcessingFunctional(new surfaceDissolutionKinetics3D<T,RXNDES>(
                                          nx, ny, nz, num_of_substrates, ade_dt,
                                          no_dynamics, bounce_back, pore_dynamics, &dissolCfg),
                                      reactionBox, ptr_dissol);
            applyProcessingFunctional(new dissolutionGather3D<T,RXNDES>(
                                          nx, ny, nz, num_of_substrates,
                                          no_dynamics, bounce_back, &dissolCfg),
                                      reactionBox, ptr_dissol);
            applyProcessingFunctional(new update_abiotic_rxnLattices<T,RXNDES>(nx, num_of_substrates, no_dynamics, bounce_back, vec_immobile),
                                      reactionBox, ptr_abiotic_kns_lattices);
            accumulateRates();   /* [v1.3.1] dissolution on its own, when the abiotic block is off */
        }

        /* [v1.3.1] The rate snapshot.
         *
         *  Written HERE, at the end of the reaction section, rather than up in the VTI
         *  block with the concentrations. The VTI block runs BEFORE the reaction, so a
         *  rate written there would be the PREVIOUS step's -- one interval stale, and
         *  identically zero in the very first file. Written here it carries the rate of
         *  the step whose number is in its filename, and it lands beside the
         *  concentrations of that same iteration because both use `iT`.
         *
         *  Equilibrium speciation is deliberately not counted. It redistributes a total
         *  between complexes rather than creating or destroying it, so folding it in
         *  would put a large number in a field labelled "reaction rate" for something
         *  that is not a reaction. */
        if (rateFieldsOn && ade_VTI_iTer > 0 && iT % ade_VTI_iTer == 0) {
            for (plint iS = 0; iS < num_of_substrates; ++iS)
                writeRateVTI(rateField[iS], iT, "rate_" + vec_subs_names[iS] + "_", (T)1./ade_dt);
        }

        // Equilibrium chemistry (runs regardless of enable_kinetics - controlled separately)
        if (useEquilibrium) {
            if (track_performance == 1) global::timer("eq").restart();
            applyProcessingFunctional(new run_equilibrium_biotic<T, RXNDES>(nx, num_of_substrates, eqSolver, no_dynamics, bounce_back),
                                      reactionBox, ptr_eq_lattices);
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
                /* computeDensity() cannot see mass held at a bounce-back wall or inside a
                 * grain, so a field that spreads out reports a FALLING total while it is in
                 * fact growing.  Example 07 is the clearest case: with the biomass on a
                 * lattice-Boltzmann solver the console line sat at 1.0800e+02 for the whole
                 * run while the biomass actually rose from 126.000 to 126.208, because the
                 * spreading patch moved 6 units of itself into wall-adjacent transit.  The
                 * conserved quantity is what the lattice holds, so that is what is printed. */
                T seen = T();
                const T totalB = complab_total_biomass(vec_bFilm_lattices, vec_bFree_lattices, &seen);
                const T held = totalB - seen;
                pcout << "│   Total biomass: " << std::scientific << totalB << std::fixed;
                if (std::fabs(held) > 1e-12 * std::fabs(totalB))
                    pcout << "   (" << std::scientific << seen << std::fixed
                          << " in open voxels, the rest held at walls)";
                pcout << "\n";
            }

            pcout << "└─────────────────────────────────────────────────────────────────────────┘\n";
        }

        // ════════════════════════════════════════════════════════════════════════
        // [NEW] <diagnostics>: one row of the summary CSV, and the conservation checks.
        //
        //   computeSum/Min/Max are Palabos reductions and are already global, which is
        //   what complab_diag::Diagnostics expects -- it knows nothing about MPI, which
        //   is what lets it be tested against synthetic fields.
        //
        //   <interval> overrides the VTI interval; 0 means "follow the VTI interval",
        //   because a scalar row is cheap and there is rarely a reason to want fewer of
        //   them than there are volumes.
        // ════════════════════════════════════════════════════════════════════════
        if (diag.active()) {
            const plint diagEvery = (icfg.diagInterval > 0) ? (plint) icfg.diagInterval
                                                            : (ade_VTI_iTer > 0 ? ade_VTI_iTer : 0);
            if (diagEvery > 0 && iT % diagEvery == 0) {
                complab_diag::Row row;
                row.iteration = (long) iT;

                //   EVERYTHING IS MEASURED OVER x = 1 .. nx-2, NOT THE WHOLE FIELD.
                //   readGeometry() duplicates the first and last slice of the file into the two
                //   ghost columns 0 and nx-1, so a reduction over the full bounding box counts
                //   the inlet and outlet slices twice.  For porosity that is a cosmetic error;
                //   for a conserved sum it is not -- the total would jump whenever the inlet
                //   concentration changed, and the mass-balance check would blame the chemistry.
                //   This is the same box saveGeometry() and the permeability calculation use.
                const Box3D physicalDomain(1, nx-2, 0, ny-1, 0, nz-1);

                //   "Open" is every material a solute can occupy: the pore materials and the
                //   biofilm materials.  Counting only the pore materials would make porosity
                //   fall as biofilm grows and would divide the substrate totals by the wrong
                //   volume, so the reported mean concentration would drift for a reason that
                //   has nothing to do with chemistry.
                //   The two lists can name the same material, so collect the distinct numbers
                //   first -- counting the same material twice would report a porosity above 1.
                std::vector<plint> openMat(pore_dynamics);
                for (size_t iB = 0; iB < bio_dynamics.size(); ++iB)
                    openMat.insert(openMat.end(), bio_dynamics[iB].begin(), bio_dynamics[iB].end());
                std::sort(openMat.begin(), openMat.end());
                openMat.erase(std::unique(openMat.begin(), openMat.end()), openMat.end());

                //   COUNTED ON THE LIVE MASK, NOT ON THE GEOMETRY FILE.
                //
                //   `geometry` is the material map as it was READ. Precipitation and dissolution
                //   change the pore space by writing the maskLattice, and neither touches
                //   `geometry` -- so counting there reports the starting geometry for the whole
                //   run, however much the pore space moves.
                //
                //   That made the headline number of the two cases the feature exists for a
                //   constant. Example 13 is documented as clogging, and its summary reported
                //   porosity 0.916667 at every interval from the first to the last. It was
                //   sealing voxels the whole time -- the conversion fires, the mask changes, the
                //   flow is re-solved -- and the CSV said nothing had happened.
                std::unique_ptr<MultiScalarField3D<T> > liveMask = computeDensity(maskLattice);
                plint openCount = 0;
                for (size_t k = 0; k < openMat.size(); ++k)
                    openCount += MaskedScalarCounts3D(physicalDomain, *liveMask, openMat[k]);
                row.openVoxels = (long) openCount;

                const double totalVoxels = (double) (nx - 2) * (double) ny * (double) nz;
                row.porosity = (totalVoxels > 0) ? (double) row.openVoxels / totalVoxels : 0.0;

                //   THE MASS THE DYNAMICS WILL NOT ADMIT TO.
                //
                //   computeDensity() asks each cell's dynamics.  BounceBack and NoDynamics both
                //   answer from a stored number and ignore the populations they are holding, so
                //   whatever is in flight at a wall, or resting inside a grain, is missing from
                //   the sum above.  It is not lost -- it streams back out next step -- but a
                //   conservation check written against `total` alone reads low by that share and
                //   sends the reader looking for a leak in the chemistry.  On a spread-out biomass
                //   patch it is about 9% of the field.
                //
                //   So it is measured, from the populations, and reported in its own column.
                //   Wall and grain voxels are parked at the substrate's background concentration
                //   at start-up (so a wall does not stream a spurious gradient into the water),
                //   and that parked value is subtracted here: what is left is the excess, which is
                //   the part that is genuinely in transit.
                for (plint iS = 0; iS < num_of_substrates; ++iS) {
                    complab_diag::FieldStat st;
                    st.name  = vec_subs_names[iS];
                    st.total = (double) computeSum(*computeDensity(vec_substr_lattices[iS]), physicalDomain);
                    st.minv  = (double) computeMin(*computeDensity(vec_substr_lattices[iS]), physicalDomain);
                    st.maxv  = (double) computeMax(*computeDensity(vec_substr_lattices[iS]), physicalDomain);
                    st.mean  = (row.openVoxels > 0) ? st.total / (double) row.openVoxels : 0.0;
                    //   Everything the lattice holds, read from the populations over the WHOLE
                    //   block, minus what the reported total was able to see.  Taking the
                    //   difference rather than masking a list of materials is what makes this
                    //   complete: a cell counts here whatever dynamics is attached to it, so the
                    //   correction cannot miss a category.  It picks up three at once --
                    //
                    //     mass in flight at a bounce-back wall, which BounceBack::computeDensity
                    //     will not report;
                    //     mass resting inside a grain, for the same reason;
                    //     mass sitting in the x=0 and x=nx-1 planes, which every reduction in this
                    //     program excludes because they are the boundary condition, and which for
                    //     a CLOSED boundary are part of the system rather than outside it.
                    //
                    //   The third was the largest.  With every boundary closed and every reaction
                    //   switched off -- a box in which nothing at all may change -- the reported
                    //   proton total fell by 7%, all of it sitting in those two planes.
                    //
                    //   An immobile species never streams, so nothing of it is ever in flight; the
                    //   difference is zero for it and the call is skipped.
                    st.held  = vec_immobile[iS] ? 0.0
                             : (PopulationSum3D(vec_substr_lattices[iS].getBoundingBox(),
                                                vec_substr_lattices[iS]) - st.total);
                    row.fields.push_back(st);
                }
                //   BIOMASS IS A CONSERVED FIELD TOO.
                //
                //   Until now the scalar record carried the substrates and nothing else, so the
                //   one field with its own transport solver -- and the one that showed the wall
                //   accounting most clearly, at about 9% of its total -- could not be checked
                //   from the CSV at all. Every population is added here, under its own name, so
                //   <conserve> can name a microbe beside a substrate.
                for (plint iM = 0; iM < num_of_microbes; ++iM) {
                    MultiBlockLattice3D<T,RXNDES> &bl = (bmass_type[iM] == 1)
                        ? vec_bFilm_lattices[loctrack[iM]] : vec_bFree_lattices[loctrack[iM]];
                    complab_diag::FieldStat st;
                    st.name  = vec_microbes_names[iM];
                    st.total = (double) computeSum(*computeDensity(bl), physicalDomain);
                    st.minv  = (double) computeMin(*computeDensity(bl), physicalDomain);
                    st.maxv  = (double) computeMax(*computeDensity(bl), physicalDomain);
                    st.mean  = (row.openVoxels > 0) ? st.total / (double) row.openVoxels : 0.0;
                    st.held  = PopulationSum3D(bl.getBoundingBox(), bl) - st.total;
                    row.fields.push_back(st);
                }

                const std::string dmsg = diag.record(row);
                if (!dmsg.empty()) pcout << dmsg;
            }
        }

        // CA biomass expansion
        if (ca_count > 0) {
            applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, caLlen, bounce_back, no_dynamics, bio_ca, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
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
                    applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, caLlen, bounce_back, no_dynamics, bio_ca, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_ca_lattices);
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
            applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, fdMaskLen, bounce_back, no_dynamics, bio_fd, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_fd_mask);
            for (plint iM=0; iM<bfilm_count; ++iM) vec_bFcopy_lattices[iM]=vec_bFilm_lattices[iM];
            for (plint iP=0; iP<bfree_count; ++iP) vec_bPcopy_lattices[iP]=vec_bFree_lattices[iP];
            applyProcessingFunctional(new fdDiffusion3D<T,RXNDES>(nx, ny, nz, fdLlen, 1, bioNUinPore[0]), vec_bFilm_lattices[0].getBoundingBox(), ptr_fd_lattices);
            applyProcessingFunctional(new updateLocalMaskNtotalLattices3D<T,RXNDES>(nx, ny, nz, fdMaskLen, bounce_back, no_dynamics, bio_fd, pore_dynamics, thrd_bFilmFrac, max_bMassRho), vec_bFilm_lattices[0].getBoundingBox(), ptr_fd_mask);
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
                        if (outletvel > thrd) {
                            /* [v1.3] The loop ran out of iterations and the outlet is still
                             * flowing, so the field coupled into every mobile solute below has
                             * NOT reached steady state. This used to clear the flag and carry on
                             * in silence -- the same failure the [FIX] at the initial NS solve was
                             * written to eliminate, still present on the two in-loop paths. */
                            ns_saturate = 0;
                            if (!ns_warned_unconverged) {
                                ns_warned_unconverged = true;
                                pcout << "\n  [NS] WARNING: the flow re-solve hit <ns_max_iT2> ("
                                      << ns_maxiTer_2 << ") without converging at iter=" << iT
                                      << ",\n         and the outlet is still flowing. The velocity "
                                      << "field advecting every solute from\n         here on is not "
                                      << "a steady state. Raise <ns_max_iT2>, or loosen\n         "
                                      << "<ns_convergence_iT2>. Reported once per run.\n";
                            }
                        }
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
                                      Box3D(1, nx-2, 0, ny-1, 0, nz-1), ptr_precip_conv);
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
                    if (outletvel > thrd) {
                        /* [v1.3] The loop ran out of iterations and the outlet is still
                         * flowing, so the field coupled into every mobile solute below has
                         * NOT reached steady state. This used to clear the flag and carry on
                         * in silence -- the same failure the [FIX] at the initial NS solve was
                         * written to eliminate, still present on the two in-loop paths. */
                        ns_saturate = 0;
                        if (!ns_warned_unconverged) {
                            ns_warned_unconverged = true;
                            pcout << "\n  [NS] WARNING: the flow re-solve hit <ns_max_iT2> ("
                                  << ns_maxiTer_2 << ") without converging at iter=" << iT
                                  << ",\n         and the outlet is still flowing. The velocity "
                                  << "field advecting every solute from\n         here on is not "
                                  << "a steady state. Raise <ns_max_iT2>, or loosen\n         "
                                  << "<ns_convergence_iT2>. Reported once per run.\n";
                        }
                    }
                    else { pcout << "\n  [PRECIP-VOP] Percolation limit reached (clogged) at iter=" << iT << "\n"; percolationFlag = 1; }
                }
                for (plint iS = 0; iS < num_of_substrates; ++iS) if (!vec_immobile[iS]) latticeToPassiveAdvDiff(nsLattice, vec_substr_lattices[iS], vec_substr_lattices[iS].getBoundingBox());
            }
            /* [FIX] A sealed voxel has to stop conducting solute, not just flow.  The substrate
             * lattices took their dynamics at start-up from the static `geometry` field, which
             * precipitation never touches, so without this line a fully clogged throat went on
             * diffusing at the full pore diffusivity and the case measured a breakthrough the
             * geometry no longer permits.  Cheap: it only walks the domain when the mask moved. */
            if (std::abs(pre_totMask - post_totMask) > 0) {
                applyProcessingFunctional(new updateSoluteSolidDynamics3D<T,RXNDES>(
                                              num_of_substrates, bounce_back, no_dynamics, pore_dynamics,
                                              substrOMEGAinbFilm, substrOMEGAinPore, vec_immobile),
                                          vec_substr_lattices[0].getBoundingBox(), substrate_lattices);
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
                                      reactionBox, ptr_dissol_conv);
            const plint post_totMask = util::roundToInt(computeAverage(*computeDensity(maskLattice))*nx*ny*nz);

            if (std::abs(pre_totMask - post_totMask) > 0 && Pe > thrd) {
                pcout << "  [DISSOL-VOP] pore space reopened at iter=" << iT << "; re-solving the flow\n";
                // A reopened voxel is ordinary pore again, so it takes the normal pore
                // relaxation, not the reduced precipitate one.
                applyProcessingFunctional(new updateNsLatticesDynamics3D<T,NSDES,T,RXNDES>(nsLatticeOmega, vec_permRatio.empty() ? (T)1 : vec_permRatio[0], pore_dynamics, no_dynamics, bounce_back),
                                          nsLattice.getBoundingBox(), nsLattice, maskLattice);
                ns_saturate = 0;                       // flow may percolate again
                plint dissol_ns_converged = 0;
                for (plint iT2 = 0; iT2 < ns_maxiTer_2; ++iT2) {
                    nsLattice.collideAndStream();
                    ns_convg2.takeValue(getStoredAverageEnergy(nsLattice),false);
                    if (ns_convg2.hasConverged()) { dissol_ns_converged = 1; break; }
                }
                /* [FIX] The precipitation block above reports a flow solve that ran out of
                 * iterations; this one used to exit silently and couple an unconverged velocity
                 * field into every mobile solute.  Say so instead. */
                if (dissol_ns_converged == 0)
                    pcout << "  [DISSOL-VOP] WARNING: the reopened flow field did not converge in "
                          << ns_maxiTer_2 << " iterations; the velocity coupled into the solutes "
                          << "below is not a steady field.\n";
                for (plint iS = 0; iS < num_of_substrates; ++iS)
                    if (!vec_immobile[iS]) latticeToPassiveAdvDiff(nsLattice, vec_substr_lattices[iS], vec_substr_lattices[iS].getBoundingBox());
                percolationFlag = 0;                   // reopening can undo a clog
            }
            /* [FIX] The mirror of the precipitation case: a voxel the mineral has vacated is
             * pore again, and a solute has to be able to enter it.  Without this the reopened
             * voxel kept the BounceBack it was given at start-up, acid could never reach the
             * fresh surface, and computeDensity there reported BounceBack's stored density
             * rather than what the voxel held. */
            if (std::abs(pre_totMask - post_totMask) > 0) {
                applyProcessingFunctional(new updateSoluteSolidDynamics3D<T,RXNDES>(
                                              num_of_substrates, bounce_back, no_dynamics, pore_dynamics,
                                              substrOMEGAinbFilm, substrOMEGAinPore, vec_immobile),
                                          vec_substr_lattices[0].getBoundingBox(), substrate_lattices);
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
        /* [v1.3.1] The closing rate snapshot, so every concentration file has a rate
         * file beside it at the same iteration. This one carries the rate of the LAST
         * step the loop ran, which is the most recent one there is -- the loop has
         * already exited, so there is no step numbered iT to take it from. */
        if (rateFieldsOn)
            for (plint iS = 0; iS < num_of_substrates; ++iS)
                writeRateVTI(rateField[iS], iT, "rate_" + vec_subs_names[iS] + "_", (T)1./ade_dt);
        /* [v1.3.2] and the closing energy snapshot, from the concentrations just written. */
        writeEnergyFields(iT);
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
        /* THE PEAK IS NOT THE TOTAL, AND IT IS NOT GROWTH.
         *
         * This line used to report (final peak - initial peak) / initial peak and call it
         * "Growth". For a solver that does not move biomass, the peak is a fair proxy. For one
         * that does, it is the opposite of the answer: example 06 spreads its patch from 108
         * voxels to 876, so the PEAK falls 4.05% while the total RISES 0.23%, and the run
         * reported -4.05% growth for a population that grew. Both are now printed, labelled
         * for what they are. */
        const T finalTotal = complab_total_biomass(vec_bFilm_lattices, vec_bFree_lattices);
        const T peakChange = (diag_initial_biomass > 0)
                           ? ((finalBmax - diag_initial_biomass) / diag_initial_biomass * 100.0) : 0.0;
        const T totalChange = (diag_initial_total_biomass > 0)
                           ? ((finalTotal - diag_initial_total_biomass) / diag_initial_total_biomass * 100.0) : 0.0;
        pcout << "╠══════════════════════════════════════════════════════════════════════════╣\n";
        pcout << "║ BIOMASS RESULTS:                                                         ║\n";
        pcout << "║   Peak density:     " << std::scientific << diag_initial_biomass
              << " -> " << finalBmax << " kg/m³  (" << std::fixed << peakChange
              << "%, falls when biomass spreads)\n" << std::scientific;
        pcout << "║   Total biomass:    " << diag_initial_total_biomass << " -> " << finalTotal
              << "  (" << std::fixed << totalChange << "%, this is the growth)\n";
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
        if (mmcfg.sym_count > 0) pcout << "│   Symbolic law:      " << symtime << " s\n";
        if (mmcfg.gnn_count > 0) pcout << "│   Graph network:     " << gnntime << " s\n";
        pcout << "└────────────────────────────────────────────────────────────────────────┘\n";
    }

    if (useEquilibrium) eqSolver.printStatistics();

    // [NEW] The scalar record's closing report: where summary.csv is, and whether any of the
    //   conserved sums drifted.  Then the surrogate's own account of how often it had to clamp
    //   to its training box -- both print nothing at all when the feature was off.
    if (diag.active()) pcout << diag.finalReport();
    //   The two learned paths report the same way: how many evaluations, and how many of them
    //   were outside the range the law or the network was fitted over.  A run that clamped a
    //   lot is a run whose results have to be looked at again, so it says so out loud.
    /* [v1.3] These four reports are built from counters incremented inside data processors, so
     * every one of them is PER RANK. pcout prints rank 0's, and a decomposition where rank 0
     * happens to hold little biomass understates the clamp fractions or -- since each report
     * returns the empty string at zero evaluations -- makes the report vanish altogether. Worse,
     * the thermodynamic gate's two verdicts ("closed everywhere, every time" / "never closed
     * anywhere") are global claims drawn from one rank's min and max. The dissolution counters
     * below have been reduced since v1.2 for exactly this reason; these were missed. */
    {
        double se = (double) complab_sym::runtime().evaluations, sc = (double) complab_sym::runtime().clamped;
        double ge = (double) complab_gnn::runtime().evaluations, gc = (double) complab_gnn::runtime().clamped;
        double re = (double) complab_srg::runtime().evaluations, rc = (double) complab_srg::runtime().clamped;
        double ce = (double) complab_srg::runtime().compiledEvaluations, cc = (double) complab_srg::runtime().compiledClamped;
        double te = (double) complab_thermo::runtime().evaluations, tb = (double) complab_thermo::runtime().blocked;
        double ts = complab_thermo::runtime().sumF;
        double tmin = complab_thermo::runtime().minF, tmax = complab_thermo::runtime().maxF;
        global::mpi().reduceAndBcast(se, MPI_SUM);  global::mpi().reduceAndBcast(sc, MPI_SUM);
        global::mpi().reduceAndBcast(ge, MPI_SUM);  global::mpi().reduceAndBcast(gc, MPI_SUM);
        global::mpi().reduceAndBcast(re, MPI_SUM);  global::mpi().reduceAndBcast(rc, MPI_SUM);
        global::mpi().reduceAndBcast(ce, MPI_SUM);  global::mpi().reduceAndBcast(cc, MPI_SUM);
        global::mpi().reduceAndBcast(te, MPI_SUM);  global::mpi().reduceAndBcast(tb, MPI_SUM);
        global::mpi().reduceAndBcast(ts, MPI_SUM);
        global::mpi().reduceAndBcast(tmin, MPI_MIN); global::mpi().reduceAndBcast(tmax, MPI_MAX);
        complab_sym::runtime().evaluations = (long) se;  complab_sym::runtime().clamped = (long) sc;
        complab_gnn::runtime().evaluations = (long) ge;  complab_gnn::runtime().clamped = (long) gc;
        complab_srg::runtime().evaluations = (long) re;  complab_srg::runtime().clamped = (long) rc;
        complab_srg::runtime().compiledEvaluations = (long) ce;
        complab_srg::runtime().compiledClamped     = (long) cc;
        complab_thermo::runtime().evaluations = (long) te;
        complab_thermo::runtime().blocked     = (long) tb;
        complab_thermo::runtime().sumF        = ts;
        /* A rank that evaluated the gate nowhere leaves minF at its 1.0 sentinel and maxF at 0.0,
         * which would otherwise widen the reduced range to [0,1] on every parallel run. */
        if (te > 0) { complab_thermo::runtime().minF = tmin; complab_thermo::runtime().maxF = tmax; }
    }
    pcout << complab_sym::runtimeReport();
    pcout << complab_gnn::runtimeReport();
    pcout << complab_srg::runtimeReport();
    pcout << complab_thermo::runtimeReport();

    /* The dissolution counters are PER RANK, and the loss this report exists to expose happens at
     * block interfaces -- so it is concentrated on exactly the ranks that are not rank 0, and a
     * report printed from rank 0 alone would understate it or miss it entirely. Reduce first. */
    {
        double rel = dissolutionRuntime().released, gat = dissolutionRuntime().gathered;
        double np  = (double) dissolutionRuntime().parked;
        double nc  = (double) dissolutionRuntime().collected;
        double ncl = (double) dissolutionRuntime().clamped;
        global::mpi().reduceAndBcast(rel, MPI_SUM);
        global::mpi().reduceAndBcast(gat, MPI_SUM);
        global::mpi().reduceAndBcast(np,  MPI_SUM);
        global::mpi().reduceAndBcast(nc,  MPI_SUM);
        global::mpi().reduceAndBcast(ncl, MPI_SUM);
        dissolutionRuntime().released  = rel;
        dissolutionRuntime().gathered  = gat;
        dissolutionRuntime().parked    = (long) (np  + 0.5);
        dissolutionRuntime().collected = (long) (nc  + 0.5);
        dissolutionRuntime().clamped   = (long) (ncl + 0.5);

        /* The per-substrate totals and the mineral removed need the same treatment, and the
         * per-substrate vector must be the same length on every rank before it is summed --
         * a rank whose blocks hold no mineral at all never calls note() and would otherwise
         * reduce a shorter vector. */
        double mrem = dissolutionRuntime().mineralRemoved;
        global::mpi().reduceAndBcast(mrem, MPI_SUM);
        dissolutionRuntime().mineralRemoved = mrem;

        double nSp = (double) dissolutionRuntime().perSpecies.size();
        global::mpi().reduceAndBcast(nSp, MPI_MAX);
        dissolutionRuntime().perSpecies.resize((size_t) (nSp + 0.5), 0.0);
        for (size_t k = 0; k < dissolutionRuntime().perSpecies.size(); ++k) {
            double v = dissolutionRuntime().perSpecies[k];
            global::mpi().reduceAndBcast(v, MPI_SUM);
            dissolutionRuntime().perSpecies[k] = v;
        }
    }
    pcout << dissolutionRuntimeReport();

    /* The upscaling record, if this run asked for one. Written from rank 0; every number in it is
     * already a Palabos reduction and therefore global. */
    if (icfg.upsEnabled) {
        const std::string upath = str_outputDir + icfg.upsCsv;
        complab_upscale::writeCsv(upath, global::mpi().isMainProcessor());
        pcout << complab_upscale::report(upath);
    }

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
