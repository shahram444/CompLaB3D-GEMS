/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =================================================================================================
 * complab3d_upscale.hh  --  from one aggregate to a number a continuum model can use
 *
 * THE QUESTION THIS ANSWERS
 *
 *   A pore-scale run resolves an aggregate: substrate arrives at the rim, is consumed on the way
 *   in, and the middle sees less than the water outside. A reservoir-scale or column-scale model
 *   cannot afford that. It carries one concentration per grid block and needs one rate.
 *
 *   The wrong way to get it is to take the bulk concentration and put it through the same rate
 *   law, because the average of a rate is not the rate of the average. The right way is the
 *   effectiveness factor: the ratio of what the aggregate actually does to what it would do if
 *   every point inside it saw the bulk.
 *
 *       eta = < r >_aggregate / r(C_bulk)
 *
 *   eta is 1 for an aggregate small or slow enough that nothing is depleted inside it, and falls
 *   towards zero as the reaction outruns diffusion. Which regime you are in is set by the Thiele
 *   modulus, the ratio of the reaction rate to the diffusion rate over the aggregate's own size:
 *
 *       phi = R * sqrt(k / D)          k = r(C_bulk) / C_bulk, a first-order rate constant, 1/s
 *
 *   For a sphere with first-order kinetics the classical result, written for THAT modulus, is
 *
 *       eta = (3/phi^2) * ( phi/tanh(phi) - 1 )
 *
 *   which is printed alongside the measured value. They agree where the assumptions hold, and
 *   the interesting cases are the ones where they do not. The modulus and the formula have to be
 *   written for the same characteristic length or the curve is evaluated at the wrong argument;
 *   complab3d_thiele.hh says why, and what went wrong when they were not.
 *
 * WHY THIS NEEDED A THERMODYNAMIC GATE TO BE WORTH BUILDING
 *
 *   With Monod kinetics alone, eta falls smoothly and the classical curve is a fair guide. With
 *   the gate on, the core of the aggregate does not merely run slowly -- it stops, because its own
 *   products have raised the free energy of the reaction to the point where the organism cannot
 *   make ATP. That is a moving internal boundary, not a gradient, and no first-order Thiele
 *   analysis anticipates it.
 *
 *   So a second group is reported beside phi: the gate evaluated at the bulk composition,
 *   F_T(C_bulk). An aggregate whose bulk sits near the threshold behaves very differently from one
 *   whose bulk is far above it, at the same phi. Reporting eta against BOTH is the whole point --
 *   and the reason evaluating F_T at the bulk concentration is not itself a shortcut, which this
 *   diagnostic demonstrates rather than asserts: r(C_bulk) includes F_T(C_bulk), so any departure
 *   of eta from 1 that survives is exactly the error you would make by upscaling that way.
 *
 * WHAT IS MEASURED, AND WHERE FROM
 *
 *   The per-voxel reaction rate is not recomputed here. Every rate path already writes its result
 *   into the dC increment lattices, and this samples those, once per diagnostic interval, after
 *   the rate processors have run and before update_rxnLattices applies them. r = dC / dt, exactly
 *   the number the solver is about to use. Averaging anything else would be measuring a different
 *   model from the one that ran.
 *
 *   The aggregate is named by material number -- the biofilm materials, normally -- and the bulk
 *   is the open pore OUTSIDE it. Both are counted on the live mask, so an aggregate that grows or
 *   dissolves is followed rather than assumed.
 *
 * WHAT THIS IS NOT
 *
 *   One run gives one point. A rate law for a continuum model needs the curve, which means a sweep
 *   over aggregate radius and bulk composition; examples/20_upscaling/offline/upscale.py
 *   drives that and fits the result. And eta is a steady-state quantity: measured during a
 *   transient it is still moving, so the run has to be long enough that it stops. The report says
 *   which of the last intervals agreed, so that is visible rather than assumed.
 * ================================================================================================= */
#ifndef COMPLAB3D_UPSCALE_HH
#define COMPLAB3D_UPSCALE_HH

#include <cmath>
#include "complab3d_thiele.hh"
#include <cstdio>
#include <string>
#include <vector>

namespace complab_upscale {

struct Config {
    bool enabled;
    std::vector<plint> aggregateMat;   /* which materials are the aggregate */
    plint species;                     /* which substrate's rate defines the reaction */
    double radius;                     /* aggregate radius, micrometres; 0 = derive from volume */
    double diffusivity;                /* m2/s inside the aggregate; 0 = take the substrate's */
    std::string csv;                   /* where the per-interval record goes */
    Config() : enabled(false), species(0), radius(0), diffusivity(0), csv("upscaling.csv") {}
};

struct Point {
    long   iteration;
    double aggregateVoxels, bulkVoxels;
    double rateMean;        /* < r > over the aggregate, mol/L/s, positive = consumption */
    double bulkConc;        /* mean of the reacting substrate outside the aggregate, mol/L */
    double rateAtBulk;      /* r(C_bulk): the rate a continuum model would have used */
    double eta;             /* rateMean / rateAtBulk */
    double thiele;
    double etaClassical;
    double gateAtBulk;      /* F_T(C_bulk), 1 when no gate is loaded */
    Point() : iteration(0), aggregateVoxels(0), bulkVoxels(0), rateMean(0), bulkConc(0),
              rateAtBulk(0), eta(0), thiele(0), etaClassical(0), gateAtBulk(1) {}
};

/* The classical effectiveness factor lives in its own header so it can be tested without
 * Palabos.  It is the RADIUS convention, phi = R sqrt(k/D), matching what complab.cpp
 * computes and what every README prints; complab3d_thiele.hh explains why that matters. */

/* The volume-averaged reaction rate inside the named materials, read from the increment lattice.
 * `dt` converts an increment back into a rate. The sign is flipped so that consumption reports
 * positive, which is the convention every effectiveness-factor treatment uses. */
template <typename T1, template <typename U1> class Descriptor1, typename T2>
class RegionRateFunctional3D : public ReductiveBoxProcessingFunctional3D_LS<T1,Descriptor1,T2>
{
public:
    RegionRateFunctional3D(std::vector<plint> mats_, T1 dt_)
        : sumId(this->getStatistics().subscribeSum()),
          cntId(this->getStatistics().subscribeSum()),
          mats(mats_), dt(dt_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor1> &inc, ScalarField3D<T2> &mask)
    {
        BlockStatistics &st = this->getStatistics();
        Dot3D ofs = computeRelativeDisplacement(inc, mask);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX)
            for (plint iY=domain.y0; iY<=domain.y1; ++iY)
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    const plint m = util::roundToInt(mask.get(iX+ofs.x, iY+ofs.y, iZ+ofs.z));
                    bool in = false;
                    for (size_t k=0; k<mats.size(); ++k) if (m == mats[k]) { in = true; break; }
                    if (!in) continue;
                    /* An increment lattice is read by update_rxnLattices through
                     * computeDensity(), which for these descriptors is sum(f) + 1: a zeroed
                     * increment has sum(f) = -1, not 0. Reading sum(f) alone reports every
                     * quiet voxel as an increment of -1, which on this case came out as a
                     * mean rate of 200 mol/L/s against a bulk rate of 0.0018. Read it the
                     * same way the solver does. */
                    T1 d = T1();
                    for (plint iPop=0; iPop<Descriptor1<T1>::q; ++iPop) d += inc.get(iX,iY,iZ)[iPop];
                    d += (T1) 1;
                    st.gatherSum(sumId, (double) (-d / dt));
                    st.gatherSum(cntId, 1.0);
                }
    }
    virtual RegionRateFunctional3D<T1,Descriptor1,T2>* clone() const {
        return new RegionRateFunctional3D<T1,Descriptor1,T2>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    double getSum()   const { return this->getStatistics().getSum(sumId); }
    double getCount() const { return this->getStatistics().getSum(cntId); }
private:
    plint sumId, cntId;
    std::vector<plint> mats;
    T1 dt;
};

/* The mean concentration over the named materials. Same shape, on the concentration lattice. */
template <typename T1, template <typename U1> class Descriptor1, typename T2>
class RegionMeanFunctional3D : public ReductiveBoxProcessingFunctional3D_LS<T1,Descriptor1,T2>
{
public:
    RegionMeanFunctional3D(std::vector<plint> mats_)
        : sumId(this->getStatistics().subscribeSum()),
          cntId(this->getStatistics().subscribeSum()), mats(mats_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor1> &lat, ScalarField3D<T2> &mask)
    {
        BlockStatistics &st = this->getStatistics();
        Dot3D ofs = computeRelativeDisplacement(lat, mask);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX)
            for (plint iY=domain.y0; iY<=domain.y1; ++iY)
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    const plint m = util::roundToInt(mask.get(iX+ofs.x, iY+ofs.y, iZ+ofs.z));
                    bool in = false;
                    for (size_t k=0; k<mats.size(); ++k) if (m == mats[k]) { in = true; break; }
                    if (!in) continue;
                    st.gatherSum(sumId, (double) lat.get(iX,iY,iZ).computeDensity());
                    st.gatherSum(cntId, 1.0);
                }
    }
    virtual RegionMeanFunctional3D<T1,Descriptor1,T2>* clone() const {
        return new RegionMeanFunctional3D<T1,Descriptor1,T2>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    double getSum()   const { return this->getStatistics().getSum(sumId); }
    double getCount() const { return this->getStatistics().getSum(cntId); }
private:
    plint sumId, cntId;
    std::vector<plint> mats;
};

/* The record. Held here rather than written straight out so the closing report can say whether the
 * last few intervals agreed -- an effectiveness factor still moving is a transient, not a result. */
inline std::vector<Point> &record()
{
    static std::vector<Point> v;
    return v;
}

inline void writeCsv(const std::string &path, bool master)
{
    if (!master || record().empty()) return;
    std::FILE *f = std::fopen(path.c_str(), "w");
    if (!f) return;
    std::fprintf(f, "# CompLB3D upscaling record: one row per diagnostic interval\n");
    std::fprintf(f, "# eta = <r>_aggregate / r(C_bulk).  phi = R sqrt(k/D), k = r(C_bulk)/C_bulk.\n");
    std::fprintf(f, "# eta_classical is the sphere with first-order kinetics, for comparison only.\n");
    std::fprintf(f, "# F_T_bulk is the thermodynamic gate at the bulk composition, 1 with no gate.\n");
    std::fprintf(f, "iteration,aggregate_voxels,bulk_voxels,rate_mean,bulk_conc,"
                    "rate_at_bulk,eta,thiele,eta_classical,F_T_bulk\n");
    for (size_t i = 0; i < record().size(); ++i) {
        const Point &p = record()[i];
        std::fprintf(f, "%ld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                     p.iteration, p.aggregateVoxels, p.bulkVoxels, p.rateMean, p.bulkConc,
                     p.rateAtBulk, p.eta, p.thiele, p.etaClassical, p.gateAtBulk);
    }
    std::fclose(f);
}

inline std::string report(const std::string &path)
{
    if (record().empty()) return std::string();
    const Point &p = record().back();
    char b[900];
    std::string s = "\n  [UPSCALE] one aggregate, reduced to the two numbers a continuum model needs\n";
    std::sprintf(b, "  [UPSCALE]   aggregate %.0f voxel(s), bulk %.0f voxel(s)\n",
                 p.aggregateVoxels, p.bulkVoxels);
    s += b;
    std::sprintf(b, "  [UPSCALE]   bulk concentration      %.6g mol/L\n", p.bulkConc);   s += b;
    std::sprintf(b, "  [UPSCALE]   rate at bulk            %.6g mol/L/s\n", p.rateAtBulk); s += b;
    std::sprintf(b, "  [UPSCALE]   measured mean rate      %.6g mol/L/s\n", p.rateMean);  s += b;
    std::sprintf(b, "  [UPSCALE]   effectiveness factor    %.6g\n", p.eta);               s += b;
    std::sprintf(b, "  [UPSCALE]   Thiele modulus          %.6g   (classical eta %.6g)\n",
                 p.thiele, p.etaClassical); s += b;
    std::sprintf(b, "  [UPSCALE]   gate at bulk            %.6g\n", p.gateAtBulk);        s += b;

    /* STEADY STATE, CHECKED RATHER THAN ASSUMED.  An effectiveness factor read off a transient is
     * a number that is still moving.  The last three intervals are compared; if they have not
     * settled, the run was too short and the figure above should not be quoted. */
    if (record().size() >= 3) {
        const size_t n = record().size();
        double lo = record()[n-3].eta, hi = lo;
        for (size_t i = n-3; i < n; ++i) {
            if (record()[i].eta < lo) lo = record()[i].eta;
            if (record()[i].eta > hi) hi = record()[i].eta;
        }
        const double spread = (std::fabs(hi) > 0) ? (hi - lo) / std::fabs(hi) : 0.0;
        if (spread < 0.01) {
            std::sprintf(b, "  [UPSCALE]   steady: the last three intervals agree to %.2g%%\n",
                         100.0 * spread);
            s += b;
        } else {
            std::sprintf(b, "  [UPSCALE]   NOT STEADY: the last three intervals span %.2g%%. This is a\n"
                            "  [UPSCALE]   transient, and eta is still moving. Run longer before quoting it.\n",
                         100.0 * spread);
            s += b;
        }
    }
    s += "  [UPSCALE]   one run is one point. examples/20_upscaling/offline/upscale.py\n"
         "  [UPSCALE]   sweeps radius and bulk composition and fits the curve.\n";
    if (!path.empty()) { std::sprintf(b, "  [UPSCALE]   record written to %s\n", path.c_str()); s += b; }
    return s;
}

}  // namespace complab_upscale

#endif
