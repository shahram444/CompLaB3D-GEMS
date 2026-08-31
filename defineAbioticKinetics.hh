/* =============================================================================
 * defineAbioticKinetics.hh   --   ABIOTIC RATE LAWS.  This file is yours to edit.
 *
 * Two hooks live here, and they are called by different sweeps. Both are shown
 * live, with generic placeholder values. Nothing is commented out.
 *
 *   defineAbioticRxnKinetics   chemistry in the water. Called in every fluid
 *                              voxel, whether or not any organism is present.
 *                              Switched on by <enable_abiotic_kinetics>.
 *
 *   defineDissolutionRate      a mineral going back into solution. Called only
 *                              on voxels that carry a phase id, and only when
 *                              <dissolution><enabled> is true.
 *
 * The substrate ORDER is the order of <substrate0>, <substrate1> ... in
 * CompLaB.xml.
 *
 * UNITS
 *     C[i]       mol/L        concentration
 *     subsR[i]   mol/L/s      rate of change of substrate i
 *     mineral    mol/L        how much mineral this voxel still holds
 *     mineralR   mol/L/s      rate of change of the mineral. MUST be <= 0 here:
 *                             a positive value is precipitation, which is the
 *                             <precipitation> block's job, not this one's.
 *     rxnR[k]    mol/L/s      optional, one entry per reaction
 * ============================================================================= */
#ifndef DEFINE_ABIOTIC_KINETICS_HH
#define DEFINE_ABIOTIC_KINETICS_HH

#include <vector>
#include <cmath>
#include <algorithm>

/* -----------------------------------------------------------------------------
 * REACTION RATE OUTPUT  -  the same contract, opted into separately.
 * The two files are separate because they are swapped separately.
 * --------------------------------------------------------------------------- */
#define COMPLAB_HAS_ABIO_RXN_RATES 1
#include <string>
namespace AbioRxnRates {
    inline const std::vector<std::string>& names() {
        static const std::vector<std::string> n = {
            "A_plus_B"        /* rxnR[0] : the single abiotic reaction below */
        };
        return n;
    }
}

namespace AbioticParams {
    /* -------- generic placeholders: replace with your own ------------------ */
    const double k_second  = 1.0;       /* L/mol/s  second order rate constant          */
    const double k_dissol  = 5.0e-2;    /* L/mol/s  dissolution, per unit acid          */
    const double s_B_per_A = 1.0;       /* mol B consumed per mol A                     */
    const double s_P_per_A = 1.0;       /* mol product formed per mol A                 */
    const double MAX_FRAC  = 0.5;       /* never consume more than half a pool per step */
    const double dt_kin    = 1.0e-2;    /* s   the reaction step; match <ade_dt>        */
    const double eps       = 1.0e-30;
}

namespace AbioticIdx {
    const int A    = 0;   /* C[0]  first reactant,  e.g. Fe2+          */
    const int Bx   = 1;   /* C[1]  second reactant, e.g. HS-           */
    const int PROD = 2;   /* C[2]  product, often the immobile mineral */
}

/* -----------------------------------------------------------------------------
 * Chemistry in the water:      A + s_B_per_A * B  ->  s_P_per_A * P
 * Second order, so both concentrations matter and a mis-ordered index changes
 * the answer rather than leaving it alone.
 * --------------------------------------------------------------------------- */
void defineAbioticRxnKinetics(
    std::vector<double> C,
    std::vector<double>& subsR,
    plb::plint mask,
    std::vector<double>* rxnR = 0
) {
    using namespace AbioticParams;
    using namespace AbioticIdx;

    for (size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    if (rxnR) for (size_t i = 0; i < rxnR->size(); ++i) (*rxnR)[i] = 0.0;

    (void) mask;   /* available if you want surface chemistry to differ from bulk */

    if ((int) C.size() <= PROD) return;

    const double a = std::max(C[A],  0.0);
    const double b = std::max(C[Bx], 0.0);

    double R = k_second * a * b;                       /* mol/L/s */

    const double capA = a * MAX_FRAC / dt_kin;
    const double capB = b * MAX_FRAC / (dt_kin * s_B_per_A);
    R = std::min(R, std::min(capA, capB));
    if (R < 0.0) R = 0.0;

    subsR[A]    = -R;
    subsR[Bx]   = -R * s_B_per_A;
    subsR[PROD] = +R * s_P_per_A;

    if (rxnR && rxnR->size() > 0) (*rxnR)[0] = R;
}

/* -----------------------------------------------------------------------------
 * A mineral going back into solution:   P + acid  ->  A + B
 *
 * C here is the water in the OPEN NEIGHBOURS of the mineral voxel, averaged, not
 * the mineral voxel itself: a solid voxel holds no solute to react with. The
 * products written into subsR are released into those same neighbours, and the
 * solver scales them by exactly the factor it applies to mineralR, so whatever
 * limiter fires, the element balance survives it.
 * --------------------------------------------------------------------------- */
inline void defineDissolutionRate(plb::plint phaseId,
                                  const std::vector<double>& C,
                                  double mineral,
                                  std::vector<double>& subsR,
                                  double& mineralR,
                                  plb::plint mask)
{
    using namespace AbioticParams;
    using namespace AbioticIdx;

    for (size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    mineralR = 0.0;

    (void) mask;   /* available if you want a different law at a wall */

    /* One phase here. With several, switch on phaseId: it is the 1-based
     * position of the <phase0>, <phase1> ... block that declared it. */
    if (phaseId != 1) return;
    if ((int) C.size() <= PROD || mineral <= 0.0) return;

    const double acid = std::max(C[Bx], 0.0);          /* whatever attacks the mineral */

    const double R = k_dissol * acid * mineral;        /* mol/L/s, positive */

    mineralR    = -R;                                  /* must be negative */
    subsR[A]    = +R * s_P_per_A;
    subsR[Bx]   = -R * s_B_per_A;
}

#endif
