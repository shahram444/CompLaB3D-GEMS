/* This file is a part of the CompLaB program.
 *
 * The CompLaB3D software (3-D pore-scale extension) is developed since 2024
 * by the University of Georgia (United States, Meile Lab, Department of
 * Marine Sciences). The original 2-D CompLaB v1.0 was a collaboration of
 * the University of Georgia and Chungnam National University (South Korea).
 * 
 * The CompLaB3D extension (3-D pore-scale model with biotic/abiotic
 * kinetics, Anderson-Accelerated Newton-Raphson PCF equilibrium chemistry,
 * and the CompLaB Studio GUI) is developed by Shahram Asgari and
 * Christof Meile (Meile Lab) at the University of Georgia. The 2-D
 * predecessor (CompLaB v1.0) was developed by Heewon Jung et al.
 *
 * Contact:
 * Shahram Asgari
 * Department of Marine Sciences (Meile Lab)
 * University of Georgia
 * Athens, GA 30602, USA
 * shahram.asgari@uga.edu
 *
 * Christof Meile
 * Department of Marine Sciences
 * University of Georgia
 * Athens, GA 30602, USA
 *
 * The most recent release of CompLaB can be downloaded at 
 * https://bitbucket.org/MeileLab/complab/downloads/
 *
 * CompLaB is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef COMPLAB3D_PROCESSORS_PART3_HH
#define COMPLAB3D_PROCESSORS_PART3_HH

using namespace plb;
typedef double T;

/* ===============================================================================================================
   ================================= COPY AND INITIALIZATION DATA PROCESSORS =====================================
   =============================================================================================================== */

// Link geometry scalar numbers and maskLattice
template<typename T1, template<typename U> class Descriptor, typename T2>
class CopyGeometryScalar2maskLattice3D : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2>
{
public:
    CopyGeometryScalar2maskLattice3D(std::vector< std::vector<plint> > mask0_) : mask0(mask0_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice, ScalarField3D<T2> &field) {
        Dot3D offset = computeRelativeDisplacement(lattice, field);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iX1 = iX+offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iY1 = iY+offset.y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iZ1 = iZ+offset.z;
                    bool flag = 0; T mask1 = field.get(iX1,iY1,iZ1); T mask2 = -1;
                    for (size_t iM=0; iM<mask0.size(); ++iM) {
                        for (size_t iN=0; iN<mask0[iM].size(); ++iN) {
                            if ( mask1 == mask0[iM][iN] ) {
                                flag = 1;
                                mask2 = mask0[iM][0];
                                break;
                            }
                        }
                        if (flag == 1) {
                            break;
                        }
                    }
                    Array<T,7> g;
                    if (flag == 0) {
                        g[0]=(T) (mask1-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (mask1-1)/8;
                    }
                    else {
                        g[0]=(T) (mask2-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (mask2-1)/8;
                    }
                    lattice.get(iX,iY,iZ).setPopulations(g); // allocate the mask number
                }
            }
        }
    }
    virtual CopyGeometryScalar2maskLattice3D<T1,Descriptor,T2>* clone() const {
        return new CopyGeometryScalar2maskLattice3D<T1,Descriptor,T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
private:
    std::vector<std::vector<plint>> mask0;
};

// Link geometry scalar numbers and ageLattice
template<typename T1, template<typename U> class Descriptor, typename T2>
class CopyGeometryScalar2ageLattice3D : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2>
{
public:
    CopyGeometryScalar2ageLattice3D()
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice, ScalarField3D<T2> &field) {
        Dot3D offset = computeRelativeDisplacement(lattice, field);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iX1 = iX+offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iY1 = iY+offset.y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iZ1 = iZ+offset.z;
                    plint mask = field.get(iX1,iY1,iZ1);
                    if (mask <0) {mask = -1;}
                    Array<T,7> g;
                    g[0]=(T) (mask-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (mask-1)/8;
                    lattice.get(iX,iY,iZ).setPopulations(g); // allocate the mask number
                }
            }
        }
    }
    virtual CopyGeometryScalar2ageLattice3D<T1,Descriptor,T2>* clone() const {
        return new CopyGeometryScalar2ageLattice3D<T1,Descriptor,T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
private:
};

// Link geometry scalar numbers and distLattice
template<typename T1, template<typename U> class Descriptor, typename T2>
class CopyGeometryScalar2distLattice3D : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2>
{
public:
    CopyGeometryScalar2distLattice3D()
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice, ScalarField3D<T2> &field) {
        Dot3D offset = computeRelativeDisplacement(lattice, field);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iX1 = iX+offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iY1 = iY+offset.y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iZ1 = iZ+offset.z;
                    plint dist = field.get(iX1,iY1,iZ1);
                    Array<T,7> g;
                    g[0]=(T) (dist-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (dist-1)/8;
                    lattice.get(iX,iY,iZ).setPopulations(g);
                }
            }
        }
    }
    virtual CopyGeometryScalar2distLattice3D<T1,Descriptor,T2>* clone() const {
        return new CopyGeometryScalar2distLattice3D<T1,Descriptor,T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
private:
};

// Copy maskLattice to geometry field
template<typename T1, template<typename U> class Descriptor, typename T2>
class CopyLattice2ScalarField3D : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2>
{
public:
    CopyLattice2ScalarField3D()
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice, ScalarField3D<T2> &field) {
        Dot3D offset = computeRelativeDisplacement(lattice, field);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iX1 = iX+offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iY1 = iY+offset.y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iZ1 = iZ+offset.z;
                    field.get(iX1,iY1,iZ1)=util::roundToInt( lattice.get(iX,iY,iZ).computeDensity() );
                }
            }
        }
    }
    virtual CopyLattice2ScalarField3D<T1,Descriptor,T2>* clone() const {
        return new CopyLattice2ScalarField3D<T1,Descriptor,T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
private:
    std::vector<std::vector<plint>> mask0;
};

// initialize scalar biomass lattice
template<typename T1, template<typename U> class Descriptor, typename T2>
class initializeScalarLattice3D : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2>
{
public:
    initializeScalarLattice3D(std::vector<T> b0_, std::vector<plint> mask0_): b0(b0_), mask0(mask0_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice, ScalarField3D<T2> &field) {
        Dot3D offset = computeRelativeDisplacement(lattice, field);
        if (b0.size() != mask0.size()) {
            std::cout << "ERROR: the size of vectors b0 and mask0 should be the same.\n";
            exit(EXIT_FAILURE);
        }
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iX1 = iX + offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iY1 = iY + offset.y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iZ1 = iZ + offset.z;
                    T2 mask1 = field.get(iX1,iY1,iZ1);
                    for (size_t iN=0; iN<b0.size(); ++iN) {
                        if (mask1 == mask0[iN]) {
                            Array<T,7> g;
                            g[0]=(T) (b0[iN]-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (b0[iN]-1)/8;
                            lattice.get(iX,iY,iZ).setPopulations(g);
                        }
                    }
                }
            }
        }
    }
    virtual initializeScalarLattice3D<T1,Descriptor,T2>* clone() const {
        return new initializeScalarLattice3D<T1,Descriptor,T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
private:
    std::vector<T> b0;
    std::vector<plint> mask0;
};

// stabilize ADE lattice
template<typename T1, template<typename U> class Descriptor, typename T2>
class stabilizeADElattice3D : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2>
{
public:
    stabilizeADElattice3D(T c0_, std::vector<plint> pore_, std::vector< std::vector<plint> > bio_): c0(c0_), pore(pore_), bio(bio_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice, ScalarField3D<T2> &field) {
        Dot3D offset = computeRelativeDisplacement(lattice, field);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iX1 = iX + offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iY1 = iY + offset.y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iZ1 = iZ + offset.z;
                    T2 mask = field.get(iX1,iY1,iZ1);
                    bool chk = 0;
                    for (size_t iP=0; iP<pore.size(); ++iP) {
                        if ( mask==pore[iP] ) { chk = 1; break;}
                    }
                    for (size_t iB0=0; iB0<bio.size(); ++iB0) {
                        if (chk == 1) {break;}
                        for (size_t iB1=0; iB1<bio[iB0].size(); ++iB1) {
                            if ( mask==bio[iB0][iB1] ) { chk = 1; break;}
                        }
                    }
                    if (chk == 1) {
                        if (c0<thrd &&c0>-thrd) {c0 = 0;}
                        Array<T,7> g;
                        g[0]=(T) (c0-1)/4; g[1]=(T) (c0-1)/8; g[2]=(T) (c0-1)/8; g[3]=(T) (c0-1)/8; g[4]=(T) (c0-1)/8; g[5]=(T) (c0-1)/8; g[6]=(T) (c0-1)/8;
                        lattice.get(iX,iY,iZ).setPopulations(g);
                    }
                }
            }
        }
    }
    virtual stabilizeADElattice3D<T1,Descriptor,T2>* clone() const {
        return new stabilizeADElattice3D<T1,Descriptor,T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
private:
    T c0;
    std::vector<plint> pore;
    std::vector< std::vector<plint> > bio;
};

// create a domain distance scalarfield3d
template<typename T1>
class createDistanceDomain3D : public BoxProcessingFunctional3D_S<T1>
{
public:
    createDistanceDomain3D(std::vector<std::vector<std::vector<plint>>> distVec_): distVec(distVec_)
    {}
    virtual void process(Box3D domain, ScalarField3D<T1> &field) {
        Dot3D absoluteOffset = field.getLocation();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint absY = iY + absoluteOffset.y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint absZ = iZ + absoluteOffset.z;
                    field.get(iX,iY,iZ) = distVec[absX][absY][absZ];
                }
            }
        }
    }
    virtual createDistanceDomain3D<T1>* clone() const {
        return new createDistanceDomain3D<T1>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
private:
    std::vector<std::vector<std::vector<plint>>> distVec;
};

// create a domain age scalarfield3d
template<typename T1>
class createAgeDomain3D : public BoxProcessingFunctional3D_S<T1>
{
public:
    createAgeDomain3D(std::vector<plint> pore_, plint bb_, plint solid_): pore(pore_), bb(bb_), solid(solid_)
    {}
    virtual void process(Box3D domain, ScalarField3D<T1> &field) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint mask=field.get(iX,iY,iZ);
                    if (mask == solid || mask == bb) { field.get(iX,iY,iZ)=-1; }
                    else {
                        bool poreflag = 0;
                        for (size_t iP=0; iP<pore.size(); ++iP) { if (mask == pore[iP]) {poreflag = 1; break;} }
                        if (poreflag == 1) { field.get(iX,iY,iZ)=0; }
                        else { field.get(iX,iY,iZ)=1; }
                    }
                }
            }
        }
    }
    virtual createAgeDomain3D<T1>* clone() const {
        return new createAgeDomain3D<T1>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
private:
    std::vector<plint> pore;
    plint bb, solid;
};

/* ===============================================================================================================
   ========================================== REDUCTIVE DATA PROCESSORS ==========================================
   =============================================================================================================== */

template<typename T1>
class MaskedBoxScalarCountFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T1>
{
public:
    MaskedBoxScalarCountFunctional3D(plint mask_) : countId(this->getStatistics().subscribeSum()), mask(mask_)
    {}
    virtual void process(Box3D domain, ScalarField3D<T1>& scalar) {
        BlockStatistics& statistics = this->getStatistics();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint tmpMask = util::roundToInt( scalar.get(iX,iY,iZ) );
                    if (tmpMask == mask) {
                        statistics.gatherSum(countId, (int) 1);
                    }
                }
            }
        }
    }
    virtual MaskedBoxScalarCountFunctional3D<T1>* clone() const {
        return new MaskedBoxScalarCountFunctional3D<T1>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    plint getCount() const {
        // The sum is internally computed on floating-point values.
        // If T is integer, the value must be rounded at the end.
        plint doubleSum = this->getStatistics().getSum(countId);
        return (T) doubleSum;
    }
private:
    plint countId;
    plint mask;
};

template<typename T1>
plint MaskedScalarCounts3D(Box3D domain, MultiScalarField3D<T1>& field, plint mask) {
    MaskedBoxScalarCountFunctional3D<T1> functional = MaskedBoxScalarCountFunctional3D<T1> ( mask );
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

// calculate RMSE for convergence checking
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
class BoxLatticeRMSEFunctional3D : public ReductiveBoxProcessingFunctional3D_LL<T1,Descriptor1,T2,Descriptor2>
{
public:
    BoxLatticeRMSEFunctional3D() : sumId(this->getStatistics().subscribeSum())
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor1> &lattice0, BlockLattice3D<T2,Descriptor2> &lattice1) {
        BlockStatistics& statistics = this->getStatistics();
        Dot3D offset_01 = computeRelativeDisplacement(lattice0,lattice1);
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    plint iX1 = iX0 + offset_01.x; plint iY1 = iY0 + offset_01.y; plint iZ1 = iZ0 + offset_01.z;
                    T deltaC = lattice0.get(iX0,iY0,iZ0).computeDensity() - lattice1.get(iX1,iY1,iZ1).computeDensity();
                    T RMSE = deltaC * deltaC;
                    statistics.gatherSum(sumId, RMSE);
                }
            }
        }
    }
    virtual BoxLatticeRMSEFunctional3D<T1,Descriptor1,T2,Descriptor2>* clone() const {
        return new BoxLatticeRMSEFunctional3D<T1,Descriptor1,T2,Descriptor2>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getCount() const {
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        double doubleSum = this->getStatistics().getSum(sumId);
        if (std::numeric_limits<T>::is_integer) {
            return (T) util::roundToInt(doubleSum);
        }
        return (T) doubleSum;
    }
private:
    plint sumId;
};

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
T computeRMSE3D(Box3D domain, MultiBlockLattice3D<T1,Descriptor1>& lattice0, MultiBlockLattice3D<T2,Descriptor2>& lattice1, T poreLen) {
    BoxLatticeRMSEFunctional3D<T1,Descriptor1,T2,Descriptor2> functional;
    applyProcessingFunctional(functional, domain, lattice0, lattice1);
    return std::sqrt(functional.getCount()/poreLen);
}

#endif // COMPLAB3D_PROCESSORS_PART3_HH
