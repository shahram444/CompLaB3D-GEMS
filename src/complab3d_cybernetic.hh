/* This file is a part of the CompLaB program.  AGPL-3.0-or-later.
 * Meile Lab, University of Georgia.  shahram.asgari@uga.edu
*/

/* ================================================================================================
 * complab3d_cybernetic.hh  --  METABOLIC SWITCHING BETWEEN ALTERNATIVE CARBON SOURCES
 * ================================================================================================
 *
 *  WHAT PROBLEM THIS SOLVES
 *
 *  An organism offered lactate, pyruvate and acetate does not eat all three at once and does not
 *  eat them in a fixed ratio. It eats the best one available, and switches when that one runs out.
 *  Shewanella oneidensis grows on lactate while excreting pyruvate and acetate, then turns round
 *  and consumes the pyruvate it made, then the acetate.
 *
 *  Flux balance analysis alone cannot produce that behaviour, and neither can a Monod term. The
 *  reason is the same in both cases: the same compound is a product in one regime and a substrate
 *  in the next, and nothing in a concentration alone says which. A single linear program handed
 *  all three carbon sources at once will simply pick whichever gives the highest growth and eat it
 *  together with the others, which is not what the organism does.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE METHOD
 *
 *  Treat the organism as a competition between several separate growth options, one per carbon
 *  source, and let the outcome of that competition decide how much of each is realised.
 *
 *    1.  Solve the metabolic model ONCE PER CARBON SOURCE, each time with only that carbon source
 *        open and the others closed. Three carbon sources means three flux vectors -- three
 *        different ways the organism could be making a living right now.
 *
 *    2.  Score each option by how much carbon it would bring in, using an ordinary Monod rate on
 *        that source alone, weighted by the number of carbon atoms the molecule carries:
 *
 *              r_k   =  k_k * C_k / (K_k + C_k)                      the unregulated rate
 *              u_k   =  n_k r_k  /  SUM_j ( n_j r_j )                the cybernetic variable
 *
 *        n_k is the carbon count: 3 for lactate, 3 for pyruvate, 2 for acetate. The u_k are
 *        non-negative and sum to one, so they are a partition of the organism's effort.
 *
 *    3.  Blend the flux vectors with those weights:
 *
 *              v  =  SUM_k  u_k * v_k                                growth and every exchange
 *
 *  When lactate is plentiful u_Lac is near one and the organism behaves as the lactate solution
 *  says -- growing, and EXCRETING pyruvate and acetate, because that is what that flux vector
 *  does. As lactate depletes u_Lac falls, u_Pyr rises, and the pyruvate solution -- which consumes
 *  pyruvate -- takes over. The switch is continuous, it needs no threshold, and it introduces no
 *  parameter that was not already in the Monod rates.
 *
 *  This is the cybernetic approach of Ramkrishna and Song (AIChE J 58:986, 2012; Cybernetic
 *  Modeling for Bioreaction Engineering, CUP 2018), in the form Song et al. (2025), Sci Rep
 *  15:6042 used to couple it to reactive transport. The specific choice of carbon uptake rate as
 *  the objective of the competition is theirs.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT IT COSTS, AND WHY IT IS NOT AS BAD AS IT LOOKS
 *
 *  m carbon sources is m linear programs per organism per voxel per step instead of one. That is
 *  the honest cost and there is no way round it: the whole method rests on having m separate flux
 *  vectors to blend.
 *
 *  Two things make it cheaper than the factor of m suggests. The m programs differ only in which
 *  exchange bounds are open, so every one of them warm-starts from the basis the previous left.
 *  And a source whose weight u_k is below `weight_floor` is not solved at all -- its contribution
 *  would be multiplied by very nearly zero. In a run where one substrate dominates most of the
 *  domain most of the time, that skip is most of the cost.
 *
 *  The surrogate path removes the cost entirely: this layer calls a solve callback, not GLPK, so
 *  the same blending works over a trained network. That is the combination Song et al. actually
 *  ran, and it is the reason this file is written against a callback rather than against a solver.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT IT DOES NOT DO
 *
 *  The weights depend on concentration only. Nothing here carries enzyme state between time steps,
 *  so an organism switches as fast as the concentrations move and there is no lag or diauxic
 *  plateau of the kind a full cybernetic model with inducible enzyme levels would produce. That is
 *  the same simplification Song et al. made, and it is adequate exactly when the switch is fast
 *  compared with transport, which is the pore-scale case.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHEN IT IS OFF
 *
 *  Fewer than two sources means the feature is inactive and every caller falls through to the
 *  ordinary single-solve path. Nothing in this file executes unless <cybernetic> appears in the
 *  input.
 * ================================================================================================
 */
#ifndef COMPLAB3D_CYBERNETIC_HH
#define COMPLAB3D_CYBERNETIC_HH

#include <string>
#include <vector>
#include <cstdio>
#include <algorithm>

namespace complab_cyb {

/* ------------------------------------------------------------------------------------------------
 *  One growth option: one carbon source the organism can live on.
 *
 *    subsIndex  index into the simulation's substrate list -- the CompLaB substrate, not the model
 *               reaction. This is what ties the option to a concentration on the lattice.
 *    carbon     number of carbon atoms in the molecule. Lactate 3, pyruvate 3, acetate 2. This is
 *               the n_k that makes the competition about carbon rather than about molecules.
 *    kmax, Ks   Monod parameters for the UNREGULATED uptake rate of this source. They set the
 *               weights only; the flux that is actually realised comes from the metabolic model.
 *    exclusive  substrate indices to close while this option is being solved -- normally the other
 *               carbon sources. Leaving them open would let the model eat two at once, which is
 *               the behaviour this whole file exists to prevent.
 * ---------------------------------------------------------------------------------------------- */
struct CyberneticSource {
    std::string        name;
    int                subsIndex;
    double             carbon;
    double             kmax;
    double             Ks;
    std::vector<int>   exclusive;

    CyberneticSource() : subsIndex(-1), carbon(1.0), kmax(0.0), Ks(1.0) {}
};

/* ------------------------------------------------------------------------------------------------
 *  The plan for one microbe.
 *
 *    weightFloor  a source whose weight falls below this is not solved. 0 solves every source
 *                 every time, which is the reference behaviour; 1e-3 is a good production value
 *                 and changes the answer in the fourth decimal at most.
 * ---------------------------------------------------------------------------------------------- */
struct CyberneticPlan {
    std::vector<CyberneticSource> sources;
    double                        weightFloor;

    CyberneticPlan() : weightFloor(1e-3) {}

    bool   active() const { return sources.size() > 1; }
    size_t size()   const { return sources.size(); }

    std::string describe() const {
        if (!active()) return "off (single growth option)";
        std::string s;
        for (size_t i = 0; i < sources.size(); ++i) {
            if (i) s += " / ";
            char buf[96];
            std::snprintf(buf, sizeof buf, "%s(C%.0f, k=%.3g, K=%.3g)",
                          sources[i].name.c_str(), sources[i].carbon,
                          sources[i].kmax, sources[i].Ks);
            s += buf;
        }
        return s;
    }
};

/* ================================================================================================
 *  weights
 *  ------------------------------------------------------------------------------------------------
 *  The cybernetic variables u_k for one voxel, from the local concentrations.
 *
 *      u_k  =  n_k r_k / SUM_j n_j r_j ,    r_k = kmax_k C_k / (K_k + C_k)
 *
 *  `avail` is indexed by CompLaB substrate index and is whatever the caller decided the organism
 *  can see -- free concentration, or the equilibrium total, per <substrate_basis>. Using the same
 *  quantity the uptake bound uses is what keeps the weights consistent with the fluxes.
 *
 *  Returns false when every source scores zero -- no carbon anywhere in this voxel. The weights
 *  are then all zero and the caller should skip the organism entirely rather than blend nothing.
 * ============================================================================================== */
template<typename T>
inline bool weights (const CyberneticPlan &plan, const std::vector<T> &avail,
                     std::vector<double> &u)
{
    const size_t m = plan.sources.size();
    u.assign(m, 0.0);

    double total = 0.0;
    for (size_t k = 0; k < m; ++k) {
        const CyberneticSource &S = plan.sources[k];
        if (S.subsIndex < 0 || S.subsIndex >= (int) avail.size()) continue;

        const double C = (double) avail[S.subsIndex];
        if (C <= 0.0) continue;

        const double denom = S.Ks + C;
        if (denom <= 0.0) continue;

        const double r = S.kmax * C / denom;      /* the unregulated (kinetic) rate */
        u[k]  = S.carbon * r;                     /* scored in carbon, not in molecules */
        total += u[k];
    }

    if (total <= 0.0) { u.assign(m, 0.0); return false; }

    for (size_t k = 0; k < m; ++k) u[k] /= total;
    return true;
}

/* ================================================================================================
 *  blend
 *  ------------------------------------------------------------------------------------------------
 *  Combine the per-source answers into the one answer the transport step will see.
 *
 *      mu       = SUM_k u_k mu_k
 *      flux[s]  = SUM_k u_k flux_k[s]
 *
 *  `fluxPer[k][s]` is the flux of substrate s in the flux vector obtained with source k open, and
 *  `gratePer[k]` the growth rate of that same solution. Sources that were skipped for low weight
 *  must have been left at zero by the caller; they contribute nothing, which is the point.
 *
 *  A NOTE ON WHAT SURVIVES THE BLEND. Each v_k satisfied S v_k = 0, and a convex combination of
 *  vectors in the null space of S is itself in the null space of S. So the blended vector is
 *  stoichiometrically consistent exactly, not approximately -- no mass is created or destroyed by
 *  the averaging. That is the property that makes this blend legitimate and it is worth stating,
 *  because averaging model outputs is not usually safe.
 * ============================================================================================== */
template<typename T>
inline void blend (const std::vector<double> &u,
                   const std::vector<T>                 &gratePer,
                   const std::vector< std::vector<T> >  &fluxPer,
                   T &grate, std::vector<T> &flux)
{
    grate = T();
    std::fill(flux.begin(), flux.end(), T());

    for (size_t k = 0; k < u.size() && k < fluxPer.size(); ++k) {
        if (u[k] <= 0.0) continue;
        const T w = (T) u[k];
        grate += w * gratePer[k];
        for (size_t s = 0; s < flux.size() && s < fluxPer[k].size(); ++s)
            flux[s] += w * fluxPer[k][s];
    }
}

}  // namespace complab_cyb

#endif  /* COMPLAB3D_CYBERNETIC_HH */
