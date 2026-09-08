/* ================================================================================================
 *  complab3d_energyfield.hh  --  THE ENERGY SNAPSHOT
 *
 *  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
 *  Meile Lab, University of Georgia.  shahram.asgari@uga.edu
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT THIS IS FOR
 *
 *  complab3d_thermo.hh computes, in every fluid voxel and at every step, how much energy the
 *  reaction actually releases at the local composition (dG, kJ/mol) and what fraction of the
 *  unrestricted rate that leaves (F_T, dimensionless, 0 to 1).  Until now neither number left the
 *  inner loop.  The whole run condensed to one line:
 *
 *      [THM] 92160024 gate evaluations, mean F_T 0.7933, range 0.0058..0.9993
 *
 *  which says the gate did something, and says nothing about WHERE.  The question that line cannot
 *  answer is the pore-scale question: the gate closes because products build up, products build up
 *  where transport is slow, so the throttling is a spatial pattern and a mean over 92 million
 *  evaluations is exactly the wrong summary of it.  A user looking at a concentration field and a
 *  rate field can see that the rate is low in the middle of an aggregate, but not whether that is
 *  because the substrate ran out (a transport limit) or because the reaction stopped paying (an
 *  energy limit).  Those two have different consequences and the same appearance.
 *
 *  This header writes both numbers as fields, so the answer is a picture:
 *
 *      dG_<microbe>_<iter>.vti     kJ/mol, negative where the reaction releases energy
 *      FT_<microbe>_<iter>.vti     0 to 1, the factor the rate was multiplied by
 *
 *  Read together with rate_<species>_<iter>.vti and <species>_<iter>.vti from the same iteration:
 *  substrate present and F_T near 1 and rate low means the organism is simply slow there;
 *  substrate present and F_T near 0 means the reaction is energetically shut, which no Monod term
 *  in the model can express.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHY IT IS A SNAPSHOT PASS AND NOT AN ACCUMULATOR
 *
 *  The rate fields in complab.cpp accumulate, because a rate is what dC[] holds for one step and
 *  dC[] is cleared twice per step; the field has to catch each contribution as it is made.
 *
 *  dG and F_T are not like that.  They are functions of the concentrations standing in the voxel at
 *  one instant, and those concentrations are still there when the output interval comes round.  So
 *  the field can be computed from scratch at write time, from the same lattices the gate reads,
 *  and the run pays nothing at all on the steps in between.  For a case that writes every 100
 *  iterations that is a hundredfold difference in cost, for an identical file.
 *
 *  One consequence worth stating plainly: the fields are an instantaneous picture at iteration iT,
 *  not an average over the interval that ended there.  A rate field beside it IS such an average.
 *  Comparing them voxel by voxel is still meaningful, because F_T moves slowly once the profile
 *  has developed, but during a fast transient the two carry different time meanings.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHICH CONCENTRATIONS
 *
 *  The gate is fed exactly what this pass feeds it: the lattice density of each substrate, scaled
 *  by <concentration_scale>, through complab_thermo::fillConc().  For the kinetics, symbolic and
 *  graph-network paths that is the same vector those paths gated with, so the field reproduces the
 *  gate the run applied.  For an FBA path running on a TOTAL rather than a FREE basis
 *  (<fba_concentration_basis>), the solver gated on totals from the speciation solver while this
 *  pass sees the lattice values; the two agree wherever complexation is weak and differ where it
 *  is not.  The distinction is noted in the log line this header prints at startup, so a user
 *  reading a dG field knows which basis produced it.
 * ================================================================================================ */

#ifndef COMPLAB3D_ENERGYFIELD_HH
#define COMPLAB3D_ENERGYFIELD_HH

#include "complab3d_thermo.hh"
#include <string>
#include <vector>

namespace complab_efield {

/* ------------------------------------------------------------------------------------------------
 *  THE FUNCTIONAL
 *
 *  Block order, which the caller in complab.cpp must match:
 *
 *      [0 .. subsNum-1]   the substrate lattices, in model order
 *      [subsNum]          the mask lattice
 *      [subsNum+1]        the dG output field
 *      [subsNum+2]        the F_T output field
 *
 *  Solid, bounce-back and boundary-column voxels are written as zero rather than left alone.  Zero
 *  is not a value dG takes in a live voxel (a reaction at exactly zero free energy is a measure-zero
 *  coincidence), so a zero in the file reads unambiguously as "no chemistry here", and ParaView's
 *  threshold filter separates the two with one click.
 * ---------------------------------------------------------------------------------------------- */
template<typename T, template<typename U> class Descriptor>
class computeEnergyField3D : public BoxProcessingFunctional3D
{
public:
    computeEnergyField3D(plint nx_, plint subsNum_, int globalMicrobe_, plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), globalMicrobe(globalMicrobe_), solid(solid_), bb(bb_)
    {}

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
    {
        const size_t nLat = (size_t) subsNum + 1;              /* substrates + mask */
        if (blocks.size() < nLat + 2) return;

        std::vector<BlockLattice3D<T,Descriptor>*> lat(nLat, (BlockLattice3D<T,Descriptor>*) 0);
        for (size_t i = 0; i < nLat; ++i)
            lat[i] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[i]);
        ScalarField3D<T>* fdG = dynamic_cast<ScalarField3D<T>*>(blocks[nLat]);
        ScalarField3D<T>* fFT = dynamic_cast<ScalarField3D<T>*>(blocks[nLat + 1]);
        if (!lat[0] || !fdG || !fFT) return;

        const Dot3D absoluteOffset = lat[0]->getLocation();
        std::vector<Dot3D> off(nLat);
        for (size_t i = 0; i < nLat; ++i) off[i] = computeRelativeDisplacement(*lat[0], *lat[i]);
        const Dot3D offG = computeRelativeDisplacement(*lat[0], *fdG);
        const Dot3D offF = computeRelativeDisplacement(*lat[0], *fFT);

        std::vector<T> conc((size_t) (subsNum > 0 ? subsNum : 0), (T) 0);
        std::vector<double> tconc;
        double dG = 0.0, ft = 0.0;

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            const plint absX = iX + absoluteOffset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                    T &outG = fdG->get(iX + offG.x, iY + offG.y, iZ + offG.z);
                    T &outF = fFT->get(iX + offF.x, iY + offF.y, iZ + offF.z);
                    outG = (T) 0;
                    outF = (T) 0;

                    if (absX <= 0 || absX >= nx - 1) continue;

                    const plint m = util::roundToInt(
                        lat[subsNum]->get(iX + off[subsNum].x,
                                          iY + off[subsNum].y,
                                          iZ + off[subsNum].z).computeDensity());
                    if (m == solid || m == bb) continue;

                    for (plint iS = 0; iS < subsNum; ++iS) {
                        const T c = lat[iS]->get(iX + off[iS].x,
                                                 iY + off[iS].y,
                                                 iZ + off[iS].z).computeDensity();
                        conc[(size_t) iS] = (c > (T) 0) ? c : (T) 0;
                    }

                    complab_thermo::fillConc(conc, (int) subsNum, tconc);
                    if (complab_thermo::gateProbe(globalMicrobe, tconc, dG, ft)) {
                        outG = (T) dG;
                        outF = (T) ft;
                    }
                }
            }
        }
    }

    virtual computeEnergyField3D<T,Descriptor>* clone() const
    {
        return new computeEnergyField3D<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        for (size_t i = 0; i < modified.size(); ++i) modified[i] = modif::nothing;
        if (modified.size() >= 2) {
            modified[modified.size() - 2] = modif::staticVariables;
            modified[modified.size() - 1] = modif::staticVariables;
        }
    }

    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }

private:
    plint nx;
    plint subsNum;
    int   globalMicrobe;
    plint solid;
    plint bb;
};

/* ------------------------------------------------------------------------------------------------
 *  WHICH ORGANISMS GET FIELDS
 *
 *  One pair per organism that actually has a block in the .thm file.  An ungated organism in a
 *  partly gated run would write a pair of files full of zeros, which is worse than no file: it
 *  looks like a shut reaction rather than an absent one.
 *
 *  The kinetics path is the exception the loader already knows about.  defineKinetics.hh returns
 *  one combined rate vector for the whole voxel, so integ::prepareThermo() refuses a .thm file with
 *  more than one block whenever that path is in use and the run gates through gateFor(0).  That
 *  single block is reported under the name of the first gated organism, and the log says so.
 * ---------------------------------------------------------------------------------------------- */
inline std::vector<int> gatedMicrobes(int nMicrobes)
{
    std::vector<int> out;
    if (!complab_thermo::enabled()) return out;
    for (int iM = 0; iM < nMicrobes; ++iM)
        if (complab_thermo::specOf(iM) >= 0) out.push_back(iM);
    return out;
}

}  /* namespace complab_efield */

#endif  /* COMPLAB3D_ENERGYFIELD_HH */
