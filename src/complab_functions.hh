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


#ifndef COMPLAB_FUNCTIONS_HH
#define COMPLAB_FUNCTIONS_HH

// ADD THESE LINES:
#include "palabos3D.h"
#include "palabos3D.hh"

#include <string>
#include <vector>


using namespace plb;
typedef double T;

// [UPD_RXN] debug-print toggle, defined in complab3d_processors_part1.hh.
//   Set here from CompLaB.xml <IO><debug_updRxn> during initialize_complab().
extern bool g_updRxn_debug;

#define NSDES descriptors::D3Q19Descriptor // Cs2 = 1/3
#define RXNDES descriptors::AdvectionDiffusionD3Q7Descriptor // Cs2 = 1/3
#define COMPLAB_THRD 1e-14
// #define thrd DBL_EPSILON

// Optional metabolic layer: needed here only for the rxntype:: vocabulary that
// the <reaction_type> parser below uses.  The header is self-contained and
// compiles to almost nothing when GLPK and Python are not built in.
#include "complab3d_metabolic.hh"

class PressureGradient {
public:
    PressureGradient(T deltaP_, plint nx_) : deltaP(deltaP_), nx(nx_)
    { }
    void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const
    {
        velocity.resetToZero();
        density = (T)1 - deltaP*NSDES<T>::invCs2 / (T)(nx-1) * (T)iX;
    }
private:
    T deltaP;
    plint nx;
};

void writeNsVTI(MultiBlockLattice3D<T,NSDES>& lattice, plint iter, std::string nameid)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();
    VtkImageOutput3D<T> vtkOut(createFileName(nameid, iter, 7), 1.);

    vtkOut.writeData<float>(*computeVelocityNorm(lattice, Box3D(1,nx-2,0,ny-1,0,nz-1)), "velocityNorm", 1.);
    vtkOut.writeData<3,float>(*computeVelocity(lattice, Box3D(1,nx-2,0,ny-1,0,nz-1)), "velocity", 1.);
}

void writePorousMediumVTI(MultiScalarField3D<int>& geometry, plint iter, std::string nameid)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    const plint nz = geometry.getNz();
    VtkImageOutput3D<T> vtkOut(createFileName(nameid, iter, 7), 1.);
    vtkOut.writeData<float>(*copyConvert<int,T>(geometry, Box3D(1,nx-2,0,ny-1,0,nz-1)), "tag", 1.0);
}

void writeAdvVTI(MultiBlockLattice3D<T,RXNDES>& lattice, plint iter, std::string nameid)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    VtkImageOutput3D<T> vtkOut(createFileName(nameid, iter, 7), 1.);
    // vtkOut.writeData<T>(*computeDensity(lattice, Box3D(0,nx-1,0,ny-1,0,nz-1)), "Density", 1.);
    vtkOut.writeData<T>(*computeDensity(lattice, Box3D(1,nx-2,0,ny-1,0,nz-1)), "Density", 1.);
}

void writeScalarVTI(MultiScalarField3D<int>& field)
{
    const plint nx = field.getNx();
    const plint ny = field.getNy();
    const plint nz = field.getNz();

    VtkImageOutput3D<T> vtkOut("distanceDomain", 1.);
    vtkOut.writeData<float>(*copyConvert<int,T>(field, Box3D (1,nx-2, 0,ny-1, 0,nz-1)), "tag", 1.0);
}

void readGeometry(std::string fNameIn, MultiScalarField3D<int>& geometry)
{
    pcout << "Reading the geometry file (" << fNameIn << ").\n";
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    const plint nz = geometry.getNz();

    Box3D sliceBox(0,0, 0,ny-1, 0,nz-1);
    std::unique_ptr< MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);

    plb_ifstream geometryFile(fNameIn.c_str());
    if (!geometryFile.is_open()) {
        pcout << "Error: could not open geometry file " << fNameIn << std::endl;
        exit(EXIT_FAILURE);
    }
    for (plint iX=1; iX<nx-1; ++iX) {
        geometryFile >> *slice;
        if (iX==1) {
            copy(*slice, slice->getBoundingBox(), geometry, Box3D(0,0, 0,ny-1, 0,nz-1));
            copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,ny-1, 0,nz-1));
        }
        else if (iX==nx-2) {
            copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,ny-1, 0,nz-1));
            copy(*slice, slice->getBoundingBox(), geometry, Box3D(nx-1,nx-1, 0,ny-1, 0,nz-1));
        }
        else {
            copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,ny-1, 0,nz-1));
        }
    }
}

void saveGeometry(std::string fNameIn, MultiScalarField3D<int>& geometry)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    const plint nz = geometry.getNz();

    pcout << "Saving geometry vti file ("<< fNameIn <<").\n";
    VtkImageOutput3D<T> vtkOut(fNameIn, 1.0);
    vtkOut.writeData<float>(*copyConvert<int,T>(geometry, Box3D (1,nx-2, 0,ny-1, 0,nz-1)), "tag", 1.0);
}

void calculateDistanceFromSolid(MultiScalarField3D<int> distance, plint nodymcs, plint bb, std::vector< std::vector< std::vector<plint> > > &distVec)
{
    const plint nx = distance.getNx();
    const plint ny = distance.getNy();
    const plint nz = distance.getNz();

    for (plint iX=0; iX<nx-1; ++iX) {
        for (plint iY=0; iY<ny-1; ++iY) {
            for (plint iZ=0; iZ<nz-1; ++iZ) {
                plint mask = distance.get(iX,iY,iZ);
                if (mask == nodymcs) { distVec[iX][iY][iZ] = -1; }
                else if (mask == bb) { distVec[iX][iY][iZ] = 0; }
                else { distVec[iX][iY][iZ] = 1; }
            }
        }
    }
    for (plint iX=0; iX<nx-1; ++iX) {
        for (plint iY=0; iY<ny-1; ++iY) {
            for (plint iZ=0; iZ<nz-1; ++iZ) {
                if ( distVec[iX][iY][iZ]==1 ) {
                    plint lp=1, r=0, dist=0;
                    while ( lp==1 ) {
                        ++r;
                        std::vector<plint> vx(r+1), vy(r+1), vz(r+1);
                        for (plint tmp=0; tmp<r+1; ++tmp) { vx[tmp]=tmp; vy[tmp]=r-tmp; vz[tmp]=0; }
                        for (plint it=0; it<r+1; ++it) {
                            plint xp = iX+vx[it], yp = iY+vy[it], zp = iZ+vz[it], xn = iX-vx[it], yn = iY-vy[it], zn = iZ-vz[it];
                            if (xp>=0 && yp>=0 && zp>=0 && xp<nx && yp<ny && zp<nz) {
                                if ( distVec[xp][yp][zp]==0 ) {
                                    dist = r; lp = 0; break;
                                }
                            }
                            if (xp>=0 && yn>=0 && zp>=0 && xp<nx && yn<ny && zp<nz) {
                                if ( distVec[xp][yn][zp]==0 ) {
                                    dist = r; lp = 0; break;
                                }
                            }
                            if (xn>=0 && yp>=0 && zp>=0 && xn<nx && yp<ny && zp<nz) {
                                if ( distVec[xn][yp][zp]==0 ) {
                                    dist = r; lp = 0; break;
                                }
                            }
                            if (xn>=0 && yn>=0 && zp>=0 && xn<nx && yn<ny && zp<nz) {
                                if ( distVec[xn][yn][zp]==0 ) {
                                    dist = r; lp=0; break;
                                }
                            }
                            if (xp>=0 && yp>=0 && zn>=0 && xp<nx && yp<ny && zn<nz) {
                                if ( distVec[xp][yp][zn]==0 ) {
                                    dist = r; lp = 0; break;
                                }
                            }
                            if (xp>=0 && yn>=0 && zn>=0 && xp<nx && yn<ny && zn<nz) {
                                if ( distVec[xp][yn][zn]==0 ) {
                                    dist = r; lp = 0; break;
                                }
                            }
                            if (xn>=0 && yp>=0 && zn>=0 && xn<nx && yp<ny && zn<nz) {
                                if ( distVec[xn][yp][zn]==0 ) {
                                    dist = r; lp = 0; break;
                                }
                            }
                            if (xn>=0 && yn>=0 && zn>=0 && xn<nx && yn<ny && zn<nz) {
                                if ( distVec[xn][yn][zn]==0 ) {
                                    dist = r; lp=0; break;
                                }
                            }
                        }
                    }
                    if (lp==0) { distVec[iX][iY][iZ] = dist; }
                }
            }
        }
    }
}

void NSdomainSetup(MultiBlockLattice3D<T,NSDES>& lattice, OnLatticeBoundaryCondition3D<T,NSDES>* boundaryCondition, MultiScalarField3D<int>& geometry, T deltaP,
                    T fluidOmega, std::vector<plint> pore, plint bounceback, plint nodymcs, std::vector<std::vector<plint>> bio_dynamics, std::vector<T> permRatio)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    Box3D west   (0,0, 0,ny-1, 0,nz-1);
    Box3D east   (nx-1,nx-1, 0,ny-1, 0,nz-1);

    // default. initialize the entire domain. may be redundant
    defineDynamics(lattice, lattice.getBoundingBox(), new IncBGKdynamics<T,NSDES>(fluidOmega));

    // pore space
    for (size_t iP=0; iP<pore.size(); ++iP) {
        if (pore[iP] > 0) { defineDynamics(lattice, geometry, new IncBGKdynamics<T,NSDES>(fluidOmega), pore[iP]); }
    }
    // bounce-back boundary
    if (bounceback > 0) {
        defineDynamics(lattice, geometry, new BounceBack<T,NSDES>(), bounceback);
    }
    // no dynamics
    if (nodymcs >= 0) {
        defineDynamics(lattice, geometry, new NoDynamics<T,NSDES>(),nodymcs);
    }
    // microbial material number
    for (size_t iM = 0; iM < bio_dynamics.size(); ++iM) {
        T bioOmega = 1/(permRatio[iM]*(1/fluidOmega-.5)+.5);
        for (size_t iB = 0; iB < bio_dynamics[iM].size(); ++iB) {
            if (bio_dynamics[iM][iB] > 0) {
                if (permRatio[iM] > COMPLAB_THRD) { // permeable biofilm
                    defineDynamics(lattice, geometry, new IncBGKdynamics<T,NSDES>(bioOmega), bio_dynamics[iM][iB]);
                }
                else { // impermeable biofilm
                    defineDynamics(lattice, geometry, new BounceBack<T,NSDES>(), bio_dynamics[iM][iB]);
                }
            }
        }
    }

    boundaryCondition->addPressureBoundary0N(west, lattice);
    setBoundaryDensity(lattice, west, (T) 1.);
    boundaryCondition->addPressureBoundary0P(east, lattice);
    setBoundaryDensity(lattice, east, (T) 1. - deltaP*NSDES<T>::invCs2);

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), PressureGradient(deltaP, nx));

    lattice.initialize();
    delete boundaryCondition;
}

void soluteDomainSetup( MultiBlockLattice3D<T,RXNDES> &lattice, OnLatticeAdvectionDiffusionBoundaryCondition3D<T, RXNDES>* boundaryCondition, MultiScalarField3D<int>& geometry,
                        T substr_bMassOmega, T substrOmega, std::vector<plint> pore, plint bounceback, plint nodymcs, std::vector<std::vector<plint>> bio_dynamics,
                        T rho0, bool left_btype, bool right_btype, T left_BC, T right_BC,
                        bool immobile = false )
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    Box3D west   (0,0, 0,ny-1, 0,nz-1);
    Box3D east   (nx-1,nx-1, 0,ny-1, 0,nz-1);
    plint processorLevelBC = 1;

    // default. initialize the entire domain. may be redundant
    defineDynamics(lattice, lattice.getBoundingBox(), new AdvectionDiffusionBGKdynamics<T,RXNDES>(substrOmega));

    // pore space
    for (size_t iP=0; iP<pore.size(); ++iP) {
        if (pore[iP] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(substrOmega), pore[iP]); }
    }
    // bounceback boundary
    if (bounceback > 0) { defineDynamics(lattice, geometry, new BounceBack<T,RXNDES>(), bounceback); }
    // no dynamics
    //   [FIX] ... except for an IMMOBILE substrate.  Palabos NoDynamics::computeDensity()
    //   returns its stored rho, which defaults to zero, so a solid voxel of a NoDynamics lattice
    //   reads back as 0 whatever was written into its populations.  For a transported solute that
    //   is correct -- there is no solute inside a grain.  For a mineral inventory it is fatal:
    //   seeding a grain phase appeared to work, every later read returned 0, the phase could never
    //   dissolve, and the reopen test (mineral < reopen_fraction * full_density) was satisfied at
    //   iteration 0 for every voxel of the phase, so the whole declared grain pack turned to pore
    //   on the first step.
    //
    //   An immobile substrate is never collided and never streamed, so leaving BGK dynamics on its
    //   solid voxels costs nothing and is what makes the stored value readable.
    if (nodymcs >= 0 && !immobile) { defineDynamics(lattice, geometry, new NoDynamics<T,RXNDES>(),nodymcs); }
    // microbial material number
    for (size_t iM = 0; iM < bio_dynamics.size(); ++iM) {
        for (size_t iB = 0; iB < bio_dynamics[iM].size(); ++iB) {
            if (bio_dynamics[iM][iB] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(substr_bMassOmega), bio_dynamics[iM][iB]); }
        }
    }

    // Set the boundary-conditions
    boundaryCondition->addTemperatureBoundary0N(west, lattice);
    if (left_btype == 0) { setBoundaryDensity(lattice, west, left_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional3D<T,RXNDES, 0, -1>, west, lattice, processorLevelBC); }
    boundaryCondition->addTemperatureBoundary0P(east, lattice);
    if (right_btype == 0) { setBoundaryDensity(lattice, east, right_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional3D<T,RXNDES, 0, +1>, east, lattice, processorLevelBC); }

    // Init lattice
    Array<T,3> u0(0., 0., 0.);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);

    lattice.initialize();
    delete boundaryCondition;
}

void soluteDeltaSetup( MultiBlockLattice3D<T,RXNDES> &lattice, OnLatticeAdvectionDiffusionBoundaryCondition3D<T, RXNDES>* boundaryCondition, MultiScalarField3D<int>& geometry,
                        T substr_bMassOmega, T substrOmega, std::vector<plint> pore, plint bounceback, plint nodymcs, std::vector<std::vector<plint>> bio_dynamics,
                        T rho0, bool left_btype, bool right_btype, T left_BC, T right_BC )
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    Box3D west   (0,0, 0,ny-1, 0,nz-1);
    Box3D east   (nx-1,nx-1, 0,ny-1, 0,nz-1);
    plint processorLevelBC = 1;

    // default. initialize the entire domain. may be redundant
    defineDynamics(lattice, lattice.getBoundingBox(), new AdvectionDiffusionBGKdynamics<T,RXNDES>(substrOmega));

    // pore space
    for (size_t iP=0; iP<pore.size(); ++iP) {
        if (pore[iP] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(substrOmega), pore[iP]); }
    }
    // bounceback boundary
    if (bounceback > 0) { defineDynamics(lattice, geometry, new BounceBack<T,RXNDES>(), bounceback); }
    // no dynamics
    if (nodymcs >= 0) { defineDynamics(lattice, geometry, new NoDynamics<T,RXNDES>(),nodymcs); }
    // microbial material number
    for (size_t iM = 0; iM < bio_dynamics.size(); ++iM) {
        for (size_t iB = 0; iB < bio_dynamics[iM].size(); ++iB) {
            if (bio_dynamics[iM][iB] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(substr_bMassOmega), bio_dynamics[iM][iB]); }
        }
    }

    // Set the boundary-conditions
    boundaryCondition->addTemperatureBoundary0N(west, lattice);
    if (left_btype == 0) { setBoundaryDensity(lattice, west, left_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional3D<T,RXNDES, 0, -1>, west, lattice, processorLevelBC); }
    boundaryCondition->addTemperatureBoundary0P(east, lattice);
    if (right_btype == 0) { setBoundaryDensity(lattice, east, right_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional3D<T,RXNDES, 0, +1>, east, lattice, processorLevelBC); }

    // Init lattice
    Array<T,3> u0(0., 0., 0.);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);

    lattice.initialize();
    delete boundaryCondition;
}

void bmassDomainSetup(MultiBlockLattice3D<T,RXNDES> &lattice, OnLatticeAdvectionDiffusionBoundaryCondition3D<T, RXNDES>* boundaryCondition,
                    MultiScalarField3D<int>& geometry, T bioOmegaPore, T bioOmegaFilm, std::vector<plint> pore, plint bounceback, plint nodymcs, std::vector<std::vector<plint>> bio_dynamics,
                    bool left_btype, bool right_btype, T left_BC, T right_BC )
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();
    plint processorLevelBC = 1;

    Box3D west   (0,0, 0,ny-1, 0,nz-1);
    Box3D east   (nx-1,nx-1, 0,ny-1, 0,nz-1);

    // default. initialize the entire domain. may be redundant
    defineDynamics(lattice, lattice.getBoundingBox(), new AdvectionDiffusionBGKdynamics<T,RXNDES>(bioOmegaPore));

    // pore space
    for (size_t iP=0; iP<pore.size(); ++iP) {
        if (pore[iP] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(bioOmegaPore), pore[iP]); }
    }
    // bounceback boundary
    if (bounceback > 0) { defineDynamics(lattice, geometry, new BounceBack<T,RXNDES>(), bounceback); }
    // no dynamics
    if (nodymcs >= 0) { defineDynamics(lattice, geometry, new NoDynamics<T,RXNDES>(),nodymcs); }
    // microbial material number
    for (size_t iM = 0; iM < bio_dynamics.size(); ++iM) {
        for (size_t iB = 0; iB < bio_dynamics[iM].size(); ++iB) {
            if (bio_dynamics[iM][iB] > 0) { defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(bioOmegaFilm), bio_dynamics[iM][iB]); }
        }
    }

    // Set the boundary-conditions
    boundaryCondition->addTemperatureBoundary0N(west, lattice);
    if (left_btype == 0) { setBoundaryDensity(lattice, west, left_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional3D<T,RXNDES, 0, -1>, west, lattice, processorLevelBC); }
    boundaryCondition->addTemperatureBoundary0P(east, lattice);
    if (right_btype == 0) { setBoundaryDensity(lattice, east, right_BC); }
    else { integrateProcessingFunctional(new FlatAdiabaticBoundaryFunctional3D<T,RXNDES, 0, +1>, east, lattice, processorLevelBC); }

    // Init lattice
    Array<T,3> u0(0., 0., 0.);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 0., u0);

    lattice.initialize();
    delete boundaryCondition;
}

void scalarDomainDynamicsSetupFromVectors(MultiBlockLattice3D<T,RXNDES> &lattice, MultiScalarField3D<int>& geometry, std::vector<plint> mtrvec, std::vector<T> omegavec)
{
    // default. initialize the entire domain. may be redundant
    defineDynamics(lattice, lattice.getBoundingBox(), new AdvectionDiffusionBGKdynamics<T,RXNDES>(0.));

    if ( mtrvec.size() != omegavec.size() ) {
        pcout << "ERROR: the length of input vectors (mtrvec and omegavec) must be the same.\n";
        exit(EXIT_FAILURE);
    }
    // assign lattice omegas (dynamics) for each mask number
    for (size_t iT=0; iT<mtrvec.size(); ++iT) {
        defineDynamics(lattice, geometry, new AdvectionDiffusionBGKdynamics<T,RXNDES>(omegavec[iT]), mtrvec[iT]);
    }
    // Init lattice
    Array<T,3> jEq(0., 0., 0.);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 0., jEq);

    lattice.initialize();
}

void scalarDomainDynamicsSetupFromGeometry(MultiBlockLattice3D<T,RXNDES> &lattice, MultiScalarField3D<int>& geometry, plint nx, plint ny, plint nz)
{
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                plint geom = geometry.get(iX,iY,iZ);
                defineDynamics( lattice, iX, iY, iZ, new AdvectionDiffusionBGKdynamics<T,RXNDES>((T) geom) );
            }
        }
    }
    // Init lattice
    Array<T,3> jEq(0., 0., 0.);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 0., jEq);

    lattice.initialize();
}

void gridSearch(MultiScalarField3D<int> geometry, std::vector< std::vector< std::vector<plint> > > &distVec, plint bb, plint solid, std::vector< std::vector<plint> > bio, std::vector<plint> pore)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    const plint nz = geometry.getNz();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                bool flag0 = 0;
                plint geom = geometry.get(iX,iY,iZ);
                for (size_t iB0=0; iB0<bio.size(); ++iB0) {
                    for (size_t iB1=0; iB1<bio[iB0].size(); ++iB1) {
                        if (geom == bio[iB0][iB1]) { flag0 = 1; }
                    }
                }
                if ( flag0 == 1 ) {
                    plint iR = 0; bool flag1 = 0;
                    while (flag1 == 0) {
                        ++iR;
                        for (plint rx=0; rx<iR; ++rx) {
                            plint ry=iR-rx;
                            for (plint rz=0; rz<=ry; ++rz) {
                                plint ry2=ry-rz;
                                if (iX+rx<nx && iY+ry2<ny && iZ+rz<nz) {
                                    plint mask = geometry.get(iX+rx,iY+ry2,iZ+rz);
                                    for (size_t iP=0; iP<pore.size(); ++iP) {
                                        if (mask == pore[iP]) {
                                            flag1 = 1; distVec[iX][iY][iZ]=iR; break;
                                        }
                                    }
                                    if (flag1==1) { break; }
                                }
                                if (iX+rx<nx && iY-ry2>0 && iZ+rz<nz) {
                                    plint mask = geometry.get(iX+rx,iY-ry2,iZ+rz);
                                    for (size_t iP=0; iP<pore.size(); ++iP) {
                                        if (mask == pore[iP]) {
                                            flag1 = 1; distVec[iX][iY][iZ]=iR; break;
                                        }
                                    }
                                    if (flag1==1) { break; }
                                }
                                if (iX-rx>0 && iY+ry2<ny && iZ+rz<nz) {
                                    plint mask = geometry.get(iX-rx,iY+ry2,iZ+rz);
                                    for (size_t iP=0; iP<pore.size(); ++iP) {
                                        if (mask == pore[iP]) {
                                            flag1 = 1; distVec[iX][iY][iZ]=iR; break;
                                        }
                                    }
                                    if (flag1==1) { break; }
                                }
                                if (iX-rx<0 && iY-ry2>0 && iZ+rz<nz) {
                                    plint mask = geometry.get(iX-rx,iY-ry2,iZ+rz);
                                    for (size_t iP=0; iP<pore.size(); ++iP) {
                                        if (mask == pore[iP]) {
                                            flag1 = 1; distVec[iX][iY][iZ]=iR; break;
                                        }
                                    }
                                    if (flag1==1) { break; }
                                }
                                if (iX+rx<nx && iY+ry2<ny && iZ-rz>0) {
                                    plint mask = geometry.get(iX+rx,iY+ry2,iZ-rz);
                                    for (size_t iP=0; iP<pore.size(); ++iP) {
                                        if (mask == pore[iP]) {
                                            flag1 = 1; distVec[iX][iY][iZ]=iR; break;
                                        }
                                    }
                                    if (flag1==1) { break; }
                                }
                                if (iX+rx<nx && iY-ry2>0 && iZ-rz>0) {
                                    plint mask = geometry.get(iX+rx,iY-ry2,iZ-rz);
                                    for (size_t iP=0; iP<pore.size(); ++iP) {
                                        if (mask == pore[iP]) {
                                            flag1 = 1; distVec[iX][iY][iZ]=iR; break;
                                        }
                                    }
                                    if (flag1==1) { break; }
                                }
                                if (iX-rx>0 && iY+ry2<ny && iZ-rz>0) {
                                    plint mask = geometry.get(iX-rx,iY+ry2,iZ-rz);
                                    for (size_t iP=0; iP<pore.size(); ++iP) {
                                        if (mask == pore[iP]) {
                                            flag1 = 1; distVec[iX][iY][iZ]=iR; break;
                                        }
                                    }
                                    if (flag1==1) { break; }
                                }
                                if (iX-rx<0 && iY-ry2>0 && iZ-rz>0) {
                                    plint mask = geometry.get(iX-rx,iY-ry2,iZ-rz);
                                    for (size_t iP=0; iP<pore.size(); ++iP) {
                                        if (mask == pore[iP]) {
                                            flag1 = 1; distVec[iX][iY][iZ]=iR; break;
                                        }
                                    }
                                    if (flag1==1) { break; }
                                }
                            }
                            if (flag1==1) { break; }
                        }
                        // if (iR == 2) {
                        //     pcout << "iX = " << iX << ", iY = " << iY << ", iZ = " << iZ << ", iX+1 = " << geometry.get(iX+1,iY,iZ) << ", iX-1 = " << geometry.get(iX-1,iY,iZ) << ", iY+1 = " << geometry.get(iX,iY+1,iZ) << ", iY-1 = " << geometry.get(iX,iY-1,iZ) << std::endl;
                        // }
                    }
                }
                else if (geom == bb || geom == solid) { distVec[iX][iY][iZ]=-1; }
                else { distVec[iX][iY][iZ]=0; }
            }
        }
    }
}

void initTotalbFilmLatticeDensity(MultiBlockLattice3D<T,RXNDES> &lattice1, MultiBlockLattice3D<T,RXNDES> &lattice2)
{
    const plint nx = lattice1.getNx();
    const plint ny = lattice1.getNy();
    const plint nz = lattice1.getNz();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                T bmass = lattice1.get(iX,iY,iZ).computeDensity();
                Array<T,7> g;
                lattice2.get(iX,iY,iZ).getPopulations(g);
                g[0]+=(T) (bmass)/4; g[1]+=(T) (bmass)/8; g[2]+=(T) (bmass)/8; g[3]+=(T) (bmass)/8; g[4]+=(T) (bmass)/8; g[5]+=(T) (bmass)/8; g[6]+=(T) (bmass)/8;
                lattice2.get(iX,iY,iZ).setPopulations(g);
            }
        }
    }
}

void defineMaskLatticeDynamics(MultiBlockLattice3D<T,RXNDES> &lattice1, MultiBlockLattice3D<T,RXNDES> &lattice2, T fbM)
{
    const plint nx = lattice1.getNx();
    const plint ny = lattice1.getNy();
    const plint nz = lattice1.getNz();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                T bmass = lattice1.get(iX,iY,iZ).computeDensity();
                T omega = 0.;
                if ( bmass > fbM ) { omega = 1.; }
                defineDynamics( lattice2, iX, iY, iZ, new AdvectionDiffusionBGKdynamics<T,RXNDES>(omega) );
            }
        }
    }
    // Init lattice
    Array<T,3> jEq(0., 0., 0.);
    initializeAtEquilibrium(lattice2, lattice2.getBoundingBox(), 0., jEq);

    lattice2.initialize();
}

int initialize_complab( char *&main_path, char *&src_path, char *&input_path, char *&output_path, char *&ns_filename, std::string &ade_filename, std::string &bio_filename, std::string &geom_filename,
    std::string &mask_filename, bool &read_NS_file, plint &ns_rerun_iT0, T &ns_converge_iT1, T &ns_converge_iT2, plint &ns_maxiTer_1, plint &ns_maxiTer_2, plint &ns_update_interval, plint &ade_update_interval,
    bool &read_ADE_file, plint &ade_rerun_iT0, plint &ade_VTK_iTer, plint &ade_CHK_iTer, T &ade_converge_iT, plint &ade_maxiTer, plint &nx, plint &ny, plint &nz, T &dx, T &dy, T &dz, T &delta_P, T &tau,
    T &Pe, T &charcs_length, std::vector<T> &solute_poreD, std::vector<T> &solute_bFilmD, std::vector<T> &bMass_poreD, std::vector<T> &bMass_bFilmD, bool &soluteDindex, bool &bmassDindex,
    T &thrd_bFilmFrac, std::vector<T> &vec_permRatio, T &max_bMassRho, std::vector<plint> &pore_dynamics, plint &bounce_back, plint &no_dynamics, std::vector< std::vector<plint> > &bio_dynamics,
    plint &num_of_microbes, plint &num_of_substrates, std::vector<std::string> &vec_subs_names, std::vector<std::string>& vec_microbes_names, std::vector<plint> &solver_type, plint &fd_count,
    plint &lb_count, plint &ca_count, plint &bfilm_count, plint &bfree_count, plint &kns_count, std::vector<plint> &reaction_type,  // ADDED
    std::vector<T> &vec_c0, std::vector<bool> &left_btype, std::vector<bool> &right_btype, std::vector<T> &vec_leftBC, std::vector<T> &vec_rightBC, std::vector< std::vector<T> > &vec_b0_all,
    std::vector<bool> &bio_left_btype, std::vector<bool> &bio_right_btype, std::vector<T> &bio_leftBC, std::vector<T> &bio_rightBC,
    std::vector< std::vector<T> > &vec_Kc, std::vector< std::vector<T> > &vec_Kc_kns,
    std::vector<T> &vec_mu, std::vector<T> &vec_mu_kns,
    std::vector<bool> &bmass_type, std::vector<T> &vec_b0_free, std::vector< std::vector<T> > &vec_b0_film,
    std::vector<std::vector<T>> &vec_Vmax, std::vector<std::vector<T>> &vec_Vmax_kns,
    bool &track_performance, bool &halfflag,
    // Equilibrium chemistry parameters
    bool &useEquilibrium, std::vector<std::string> &eq_component_names,
    std::vector<T> &eq_logK_values, std::vector<std::vector<T>> &eq_stoich_matrix,
    // NEW: Biotic/Abiotic and Kinetics control
    bool &biotic_mode, bool &enable_kinetics,
    // NEW: Abiotic kinetics (chemical reactions without microbes)
    bool &enable_abiotic_kinetics,
    // NEW: Validation diagnostics option
    bool &enable_validation_diagnostics)

{


    try {
        std::string fin("CompLaB.xml");
        XMLreader doc(fin);

        // ════════════════════════════════════════════════════════════════════════════
        // BIOTIC/ABIOTIC MODE AND KINETICS CONTROL
        // ════════════════════════════════════════════════════════════════════════════
        // Default values
        biotic_mode = true;      // Default: biotic simulation with microbes
        enable_kinetics = true;  // Default: kinetics enabled

        try {
            std::string tmp;
            doc["parameters"]["simulation_mode"]["biotic_mode"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0 || tmp.compare("biotic")==0) {
                biotic_mode = true;
                pcout << "\n╔══════════════════════════════════════════════════════════════════════╗\n";
                pcout << "║ SIMULATION MODE: BIOTIC (with microbes/biomass)                      ║\n";
                pcout << "╚══════════════════════════════════════════════════════════════════════╝\n";
            }
            else if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0 || tmp.compare("abiotic")==0) {
                biotic_mode = false;
                pcout << "\n╔══════════════════════════════════════════════════════════════════════╗\n";
                pcout << "║ SIMULATION MODE: ABIOTIC (no microbes - transport only)              ║\n";
                pcout << "╚══════════════════════════════════════════════════════════════════════╝\n";
            }
            else {
                pcout << "biotic_mode (" << tmp << ") should be true/false or biotic/abiotic. Defaulting to biotic.\n";
                biotic_mode = true;
            }
        }
        catch (PlbIOException& exception) {
            biotic_mode = true;  // Default to biotic mode
        }

        try {
            std::string tmp;
            doc["parameters"]["simulation_mode"]["enable_kinetics"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0) {
                enable_kinetics = true;
                pcout << "Kinetics reactions: ENABLED\n";
            }
            else if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0) {
                enable_kinetics = false;
                pcout << "Kinetics reactions: DISABLED (equilibrium solver only)\n";
            }
            else {
                pcout << "enable_kinetics (" << tmp << ") should be true or false. Defaulting to true.\n";
                enable_kinetics = true;
            }
        }
        catch (PlbIOException& exception) {
            enable_kinetics = true;  // Default to kinetics enabled
        }

        // ════════════════════════════════════════════════════════════════════════════
        // ABIOTIC KINETICS OPTION
        // ════════════════════════════════════════════════════════════════════════════
        // Abiotic kinetics = chemical reactions between substrates WITHOUT microbes
        // Examples: first-order decay, bimolecular reactions, radioactive decay
        enable_abiotic_kinetics = false;  // Default: off

        try {
            std::string tmp;
            doc["parameters"]["simulation_mode"]["enable_abiotic_kinetics"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0) {
                enable_abiotic_kinetics = true;
                pcout << "╔══════════════════════════════════════════════════════════════════════╗\n";
                pcout << "║ ABIOTIC KINETICS: ENABLED                                            ║\n";
                pcout << "║ Chemical reactions between substrates (no microbes)                  ║\n";
                pcout << "╚══════════════════════════════════════════════════════════════════════╝\n";
            }
        }
        catch (PlbIOException& exception) {
            enable_abiotic_kinetics = false;
        }

        // If abiotic mode, disable biotic kinetics but allow abiotic kinetics
        if (!biotic_mode) {
            enable_kinetics = false;  // Disable biotic (Monod) kinetics
            if (!enable_abiotic_kinetics) {
                pcout << "Note: Biotic kinetics disabled (abiotic mode)\n";
                pcout << "      Set enable_abiotic_kinetics=true for substrate reactions\n\n";
            }
        }

        // NOTE: Validation that abiotic kinetics requires substrates is done
        // after chemistry section parsing (see below around line 1023)

        // ════════════════════════════════════════════════════════════════════════════
        // VALIDATION DIAGNOSTICS OPTION
        // ════════════════════════════════════════════════════════════════════════════
        // When enabled, prints detailed data flow and calculation verification at each iteration
        enable_validation_diagnostics = false;  // Default: off (performance mode)

        try {
            std::string tmp;
            doc["parameters"]["simulation_mode"]["enable_validation_diagnostics"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0) {
                enable_validation_diagnostics = true;
                pcout << "╔══════════════════════════════════════════════════════════════════════╗\n";
                pcout << "║  VALIDATION DIAGNOSTICS: ENABLED                                     ║\n";
                pcout << "║  Detailed per-iteration output for data flow verification            ║\n";
                pcout << "║  WARNING: This adds overhead - use for debugging only!               ║\n";
                pcout << "╚══════════════════════════════════════════════════════════════════════╝\n\n";
            }
            else if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0) {
                enable_validation_diagnostics = false;
            }
        }
        catch (PlbIOException& exception) {
            enable_validation_diagnostics = false;  // Default to off
        }

        // terminate the simulation if inputs are undefined.
        try {
            // LB_numerics
            doc["parameters"]["LB_numerics"]["domain"]["nx"].read(nx); nx+=2;
            doc["parameters"]["LB_numerics"]["domain"]["ny"].read(ny);
            doc["parameters"]["LB_numerics"]["domain"]["nz"].read(nz);
            doc["parameters"]["LB_numerics"]["domain"]["dx"].read(dx);
            doc["parameters"]["LB_numerics"]["domain"]["filename"].read(geom_filename);

            // chemistry
            doc["parameters"]["chemistry"]["number_of_substrates"].read(num_of_substrates);
            for (plint iT=0; iT<num_of_substrates; ++iT) {
                T bc0, bc1;
                std::string tmp0, tmp1;
                std::string chemname = "substrate"+std::to_string(iT);
                doc["parameters"]["chemistry"][chemname]["left_boundary_type"].read(tmp0);
                std::transform(tmp0.begin(), tmp0.end(), tmp0.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp0.compare("dirichlet")==0) { left_btype.push_back(0); }
                else if (tmp0.compare("neumann")==0) { left_btype.push_back(1); }
                else { pcout << "left_boundary_type (" << tmp0 << ") should be either Dirichlet or Neumann. Terminating the simulation.\n"; return -1; }
                doc["parameters"]["chemistry"][chemname]["right_boundary_type"].read(tmp1);
                std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp1.compare("dirichlet")==0) { right_btype.push_back(0); }
                else if (tmp1.compare("neumann")==0) { right_btype.push_back(1); }
                else { pcout << "right_boundary_type (" << tmp1 << ") should be either Dirichlet or Neumann. Terminating the simulation.\n"; return -1; }
                doc["parameters"]["chemistry"][chemname]["left_boundary_condition"].read(bc0); vec_leftBC.push_back(bc0);
                doc["parameters"]["chemistry"][chemname]["right_boundary_condition"].read(bc1); vec_rightBC.push_back(bc1);
            }
            if (left_btype.size()!=(unsigned)num_of_substrates) { pcout << "The length of left_boundary_type vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
            if (right_btype.size()!=(unsigned)num_of_substrates) { pcout << "The length of right_boundary_type vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
            if (vec_leftBC.size()!=(unsigned)num_of_substrates) { pcout << "The length of left_boundary_condition vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
            if (vec_rightBC.size()!=(unsigned)num_of_substrates) { pcout << "The length of right_boundary_condition vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}

            // microbiology
            // In abiotic mode, skip microbe parsing entirely
            if (!biotic_mode) {
                // Check if user has microbes defined in XML but biotic_mode=false
                try {
                    plint xml_microbes = 0;
                    doc["parameters"]["microbiology"]["number_of_microbes"].read(xml_microbes);
                    if (xml_microbes > 0) {
                        pcout << "╔══════════════════════════════════════════════════════════════════════╗\n";
                        pcout << "║ WARNING: " << xml_microbes << " microbe(s) defined in XML but biotic_mode=false    ║\n";
                        pcout << "║ Microbes will be IGNORED. Set biotic_mode=true to use them.         ║\n";
                        pcout << "╚══════════════════════════════════════════════════════════════════════╝\n";
                    }
                } catch (...) { }  // No microbes section defined, which is expected in abiotic mode

                num_of_microbes = 0;
                kns_count = 0; fd_count = 0; ca_count = 0; lb_count = 0;
                pcout << "Abiotic mode: Skipping microbiology section\n";
            }
            else {
                doc["parameters"]["microbiology"]["number_of_microbes"].read(num_of_microbes);
            }
            kns_count = 0; fd_count = 0; ca_count = 0; lb_count = 0;  // Initialize all counters

            for (plint iT=0; iT<num_of_microbes; ++iT) {
                std::string tmp1, tmp2;
                std::vector<T> b0;
                std::string bioname = "microbe"+std::to_string(iT);
                
                // Parse solver_type (CA, LBM, FD)
                doc["parameters"]["microbiology"][bioname]["solver_type"].read(tmp1);
                std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp1.compare("fd")==0 || tmp1.compare("finite difference")==0 || tmp1.compare("finite_difference")==0) { 
                    solver_type.push_back(1); ++fd_count; 
                }
                else if (tmp1.compare("ca")==0 || tmp1.compare("cellular automata")==0 || tmp1.compare("cellular_automata")==0) { 
                    solver_type.push_back(2); ++ca_count; 
                }
                else if (tmp1.compare("lbm")==0 || tmp1.compare("lb")==0 || tmp1.compare("lattice boltzmann")==0) { 
                    solver_type.push_back(3); ++lb_count; 
                }
                else {
                    pcout << "Palabos IO exception: Element solver_type " << tmp1 << " is not defined. Use FD, CA, or LBM. Terminating the simulation.\n"; 
                    return -1;
                }
                
                // Parse reaction_type.
                //
                // The vocabulary was extended when the optional metabolic layer was added
                // back (see complab3d_metabolic.hh).  It now accepts:
                //
                //     none                     0   no reaction at all
                //     kinetics                 1   defineKinetics.hh          <-- original
                //     glpk                     2   flux balance via GLPK
                //     surrogate                3   surrogateModel.hh
                //     cobrapy                  4   flux balance via COBRApy
                //     glpk_and_kinetics        5   both
                //     surrogate_and_kinetics   6   both
                //     cobrapy_and_kinetics     7   both
                //
                // 0 and 1 keep the meaning they always had in CompLB3D, so existing input
                // files are unaffected.  rxntype::parse() owns the spelling table; do not
                // duplicate it here.
                try {
                    doc["parameters"]["microbiology"][bioname]["reaction_type"].read(tmp2);
                    const plint rt = rxntype::parse(tmp2);
                    if (rt < 0) {
                        pcout << "Palabos IO exception: Element reaction_type " << tmp2 << " is not defined.\n"
                              << "  Use one of: none, kinetics, glpk, surrogate, cobrapy,\n"
                              << "              glpk_and_kinetics, surrogate_and_kinetics, cobrapy_and_kinetics.\n"
                              << "  Terminating the simulation.\n";
                        return -1;
                    }
                    reaction_type.push_back(rt);
                    if (rxntype::usesKinetics(rt)) ++kns_count;
                }
                catch (PlbIOException& exception) {
                    // Default: assume kinetics if not specified
                    reaction_type.push_back(rxntype::KINETICS); ++kns_count;
                    pcout << "WARNING: reaction_type not specified for " << bioname << ". Defaulting to 'kinetics'.\n";
                }
                
                doc["parameters"]["microbiology"][bioname]["initial_densities"].read(b0); 
                vec_b0_all.push_back(b0);
            }

            // Validation
            if (solver_type.size() != (unsigned)num_of_microbes) { 
                pcout << "The length of solver_type vector does not match the num_of_microbes. Terminating the simulation.\n"; 
                return -1; 
            }
            if (reaction_type.size() != (unsigned)num_of_microbes) { 
                pcout << "The length of reaction_type vector does not match the num_of_microbes. Terminating the simulation.\n"; 
                return -1; 
            }
           if (vec_b0_all.size() != (unsigned)num_of_microbes) { 
                pcout << "The length of initial_densities vector does not match the num_of_microbes. Terminating the simulation.\n"; 
                return -1; 
            }
        }
        catch (PlbIOException& exception) {
            pcout << exception.what() << " Required parameter missing. Terminating the simulation.\n" << std::endl;
            return -1;
        }

        // parameters with default values
        // define paths
        try {
            std::string item;
            doc["parameters"]["path"]["src_path"].read(item);
            src_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { src_path[i] = item[i]; }
            src_path[item.size()+1] = '\0';
        }
        catch (PlbIOException& exception) {
            std::string item = "src";
            src_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { src_path[i] = item[i]; }
            src_path[item.size()+1] = '\0';
        }
        try {
            std::string item;
            doc["parameters"]["path"]["input_path"].read(item);
            input_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { input_path[i] = item[i]; }
            input_path[item.size()+1] = '\0';
        }
        catch (PlbIOException& exception) {
            std::string item = "input";
            input_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { input_path[i] = item[i]; }
            input_path[item.size()+1] = '\0';
        }
        try {
            std::string item;
            doc["parameters"]["path"]["output_path"].read(item);
            output_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { output_path[i] = item[i]; }
            output_path[item.size()+1] = '\0';
        }
        catch (PlbIOException& exception) {
            std::string item = "output";
            output_path = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { output_path[i] = item[i]; }
            output_path[item.size()+1] = '\0';
        }

        // LB_numerics
        try { doc["parameters"]["LB_numerics"]["delta_P"].read(delta_P); }
        catch (PlbIOException& exception) { delta_P=0; }
        try {
            std::string tmp;
            doc["parameters"]["LB_numerics"]["track_performance"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0 ) { track_performance = 0; }
            else if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0 ) { track_performance = 1; }
            else { pcout << "track_performance (" << tmp << ") should be either true or false. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { track_performance=0; }
        try { doc["parameters"]["LB_numerics"]["Peclet"].read(Pe); }
        catch (PlbIOException& exception) { Pe=0; }
        if (delta_P < COMPLAB_THRD) {Pe=0;}
        try { doc["parameters"]["LB_numerics"]["tau"].read(tau); }
        catch (PlbIOException& exception) { tau=0.8; }
        try { doc["parameters"]["LB_numerics"]["domain"]["dy"].read(dy); }
        catch (PlbIOException& exception) { dy = dx; }
        try { doc["parameters"]["LB_numerics"]["domain"]["dz"].read(dz); }
        catch (PlbIOException& exception) { dz = dx; }
        try { doc["parameters"]["LB_numerics"]["domain"]["characteristic_length"].read(charcs_length); }
        catch (PlbIOException& exception) {
            charcs_length = 0;
            if (Pe>COMPLAB_THRD) {
                pcout << "charcs_length must be defined when for transport simulations (Pe > 0). Terminating the simulation.\n"; return -1;
            }
        }
        try {
            std::string unit;
            doc["parameters"]["LB_numerics"]["domain"]["unit"].read(unit);
            if (unit == "m") {  charcs_length /= dx; /* do nothing */ }
            else if (unit == "mm") { charcs_length /= dx; dx *= 1e-3; }
            else if (unit == "um") { charcs_length /= dx; dx *= 1e-6; }
            else { pcout << "unit ("<< unit << ") must be either m, mm, or um. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { charcs_length /= dx; dx *= 1e-6; }
        try { doc["parameters"]["LB_numerics"]["domain"]["material_numbers"]["pore"].read(pore_dynamics); }
        catch (PlbIOException& exception) { pore_dynamics.push_back(2); }
        try { doc["parameters"]["LB_numerics"]["domain"]["material_numbers"]["solid"].read(no_dynamics); }
        catch (PlbIOException& exception) { no_dynamics=0; }
        try { doc["parameters"]["LB_numerics"]["domain"]["material_numbers"]["bounce_back"].read(bounce_back); }
        catch (PlbIOException& exception) { bounce_back=1; }
        bfilm_count=0; bfree_count=0;
        for (plint iT=0; iT<num_of_microbes; ++iT) {
            std::string name = "microbe" + std::to_string(iT);
            try {
                std::vector<plint> tmp;
                doc["parameters"]["LB_numerics"]["domain"]["material_numbers"][name].read(tmp);
                bio_dynamics.push_back(tmp);
                bmass_type.push_back(1);
                vec_b0_film.push_back(vec_b0_all[iT]);
                if (tmp.size()!=vec_b0_all[iT].size()) {pcout << "The length of microbial material_number vector is not consistent with the length of initial_densities vector. Terminating the simulation.\n"; return -1;}
                ++bfilm_count;
            }
            catch (PlbIOException& exception) { bmass_type.push_back(0); vec_b0_free.push_back(vec_b0_all[iT][0]); ++bfree_count; }
        }
        try {
            std::string tmp;
            doc["parameters"]["IO"]["read_NS_file"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0 ) { read_NS_file = 0; }
            else if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0 ) { read_NS_file = 1; }
            else { pcout << "read_NS_file (" << tmp << ") should be either true or false. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { read_NS_file=0; }
        try {
            doc["parameters"]["LB_numerics"]["iteration"]["ns_rerun_iT0"].read(ns_rerun_iT0);
            if (ns_rerun_iT0 < 0) {
                pcout << "ns_rerun_iT0 ("<< ns_rerun_iT0 << ") must be a positive number. Terminating the simulation.\n";
                return -1;
            }
        }
        catch (PlbIOException& exception) { if (read_NS_file == 1) { pcout << "WARNING: NS checkpoint file is loaded but ns_rerun_iT0 is not provided. Assume no further flow simulation.\n"; ns_rerun_iT0=0; } }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_update_interval"].read(ns_update_interval); }
        catch (PlbIOException& exception) { ns_update_interval=1; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ade_update_interval"].read(ade_update_interval); }
        catch (PlbIOException& exception) { ade_update_interval=1; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_max_iT1"].read(ns_maxiTer_1); }
        catch (PlbIOException& exception) { ns_maxiTer_1=100000; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_max_iT2"].read(ns_maxiTer_2); }
        catch (PlbIOException& exception) { ns_maxiTer_2=100000; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_converge_iT1"].read(ns_converge_iT1); }
        catch (PlbIOException& exception) { ns_converge_iT1=1e-8; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ns_converge_iT2"].read(ns_converge_iT2); }
        catch (PlbIOException& exception) { ns_converge_iT2=1e-6; }
        try {
            doc["parameters"]["LB_numerics"]["iteration"]["ade_rerun_iT0"].read(ade_rerun_iT0);
            if (ade_rerun_iT0 < 0) {
                pcout << "ade_rerun_iT0 ("<< ade_rerun_iT0 << ") must be a positive number. Terminating the simulation.\n";
                return -1;
            }
        }
        catch (PlbIOException& exception) { if (read_ADE_file == 1) { pcout << "WARNING: ADE checkpoint file is loaded but ade_rerun_iT0 is not provided. Assume no further flow simulation.\n"; ade_rerun_iT0=0; } }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ade_max_iT"].read(ade_maxiTer); }
        catch (PlbIOException& exception) { ade_maxiTer=10000000; }
        try { doc["parameters"]["LB_numerics"]["iteration"]["ade_converge_iT"].read(ade_converge_iT); }
        catch (PlbIOException& exception) { ade_converge_iT=1e-8; }

        // chemistry
        soluteDindex = 0;
        for (plint iT=0; iT<num_of_substrates; ++iT) {
            T D0, D1, c0;
            std::string chemname = "substrate"+std::to_string(iT);
            try { doc["parameters"]["chemistry"][chemname]["name_of_substrates"].read(vec_subs_names); }
            catch (PlbIOException& exception) { vec_subs_names.push_back("substrate_"+std::to_string(iT)); }
            try { doc["parameters"]["chemistry"][chemname]["substrate_diffusion_coefficients"]["in_pore"].read(D0); solute_poreD.push_back(D0);}
            catch (PlbIOException& exception) { solute_poreD.push_back(1e-9); }
            try { doc["parameters"]["chemistry"][chemname]["substrate_diffusion_coefficients"]["in_biofilm"].read(D1); solute_bFilmD.push_back(D1);}
            catch (PlbIOException& exception) { solute_bFilmD.push_back(2e-10); }
            if ( std::abs(D1-D0)> COMPLAB_THRD ) { soluteDindex = 1; }
            try { doc["parameters"]["chemistry"][chemname]["initial_concentration"].read(c0); vec_c0.push_back(c0); }
            catch (PlbIOException& exception) { vec_c0.push_back(0.0); }
        }
        if (vec_c0.size()!=(unsigned)num_of_substrates) { pcout << "The length of initial_concentration vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1;}
        if (solute_poreD.size()!=(unsigned)num_of_substrates) { pcout << "The length of substrate_diffusion_coefficients in_pore vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1; }
        if (solute_bFilmD.size()!=(unsigned)num_of_substrates) { pcout << "The length of substrate_diffusion_coefficients in_biofilm vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1; }
        if (vec_subs_names.size()!=(unsigned)num_of_substrates) { pcout << "The length of name_of_substrates vector does not match the num_of_substrates. Terminating the simulation.\n"; return -1; }

        // Validate abiotic kinetics requires at least one substrate
        if (enable_abiotic_kinetics && num_of_substrates < 1) {
            pcout << "╔══════════════════════════════════════════════════════════════════════╗\n";
            pcout << "║ ERROR: enable_abiotic_kinetics=true but no substrates defined!       ║\n";
            pcout << "║ Abiotic kinetics requires at least one substrate to react.           ║\n";
            pcout << "║ Add substrates in <chemistry> section or set enable_abiotic_kinetics=false.║\n";
            pcout << "╚══════════════════════════════════════════════════════════════════════╝\n";
            return -1;
        }

        // microbiology
        bmassDindex = 0;
        // Skip microbiology parsing in abiotic mode
        if (biotic_mode) {
        for (plint iT=0; iT<num_of_microbes; ++iT) {
            std::string bioname = "microbe"+std::to_string(iT);
            try {
                std::string tmp;
                doc["parameters"]["microbiology"][bioname]["name_of_microbes"].read(tmp);
                vec_microbes_names.push_back(tmp);
            }
            catch (PlbIOException& exception) { vec_microbes_names.push_back(bioname); }
            try {
                T mu;
                doc["parameters"]["microbiology"][bioname]["decay_coefficient"].read(mu);
                vec_mu.push_back(mu);
            }
            catch (PlbIOException& exception) { vec_mu.push_back(0.); }
            try {
                std::string tmp0, tmp1;
                doc["parameters"]["microbiology"][bioname]["left_boundary_type"].read(tmp0); tmp1=tmp0;
                std::transform(tmp0.begin(), tmp0.end(), tmp0.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp0.compare("dirichlet")==0) { bio_left_btype.push_back(0); }
                else if (tmp0.compare("neumann")==0) { bio_left_btype.push_back(1); }
                else { pcout << "left_boundary_type (" << tmp1 << ") should be either Dirichlet or Neumann. Terminating the simulation.\n"; return -1; }
            }
            catch (PlbIOException& exception) { bio_left_btype.push_back(1); }
            try {
                std::string tmp0, tmp1;
                doc["parameters"]["microbiology"][bioname]["right_boundary_type"].read(tmp0); tmp1=tmp0;
                std::transform(tmp0.begin(), tmp0.end(), tmp0.begin(), [](unsigned char c){ return std::tolower(c); });
                if (tmp0.compare("dirichlet")==0) { bio_right_btype.push_back(0); }
                else if (tmp0.compare("neumann")==0) { bio_right_btype.push_back(1); }
                else { pcout << "right_boundary_type (" << tmp1 << ") should be either Dirichlet or Neumann. Terminating the simulation.\n"; return -1; }
            }
            catch (PlbIOException& exception) { bio_right_btype.push_back(1); }
            try {
                T bc;
                doc["parameters"]["microbiology"][bioname]["left_boundary_condition"].read(bc);
                bio_leftBC.push_back(bc);
            }
            catch (PlbIOException& exception) { bio_leftBC.push_back(0.); }
            try {
                T bc;
                doc["parameters"]["microbiology"][bioname]["right_boundary_condition"].read(bc);
                bio_rightBC.push_back(bc);
            }
            catch (PlbIOException& exception) { bio_rightBC.push_back(0.); }
            try {
                T D0;
                doc["parameters"]["microbiology"][bioname]["biomass_diffusion_coefficients"]["in_pore"].read(D0);
                bMass_poreD.push_back(D0);
            }
            catch (PlbIOException& exception) {
                if ( solver_type[iT]==1 ) { pcout << exception.what() << " for microbe" << iT << ". It must be defined when solver_type is Finite Difference. Terminating the simulation.\n"; return -1; }
                else { bMass_poreD.push_back(-99); }
            }
            try {
                T D0;
                doc["parameters"]["microbiology"][bioname]["biomass_diffusion_coefficients"]["in_biofilm"].read(D0);
                bMass_bFilmD.push_back(D0);
            }
            catch (PlbIOException& exception) {
                if ( solver_type[iT]==1 ) { pcout << exception.what() << " for microbe" << iT << ". It must be defined when solver_type is Finite Difference. Terminating the simulation.\n"; return -1; }
                else { bMass_bFilmD.push_back(-99); }
            }
            if ( (bMass_poreD[iT]>0&&bMass_poreD[iT]>0) && (std::abs(bMass_poreD[iT]-bMass_bFilmD[iT])>COMPLAB_THRD) ) { bmassDindex = 1; }
            try {
                T nu0;
                doc["parameters"]["microbiology"][bioname]["viscosity_ratio_in_biofilm"].read(nu0);
                if (nu0 > COMPLAB_THRD) { vec_permRatio.push_back(1/nu0); }
                else { vec_permRatio.push_back(-99); }
            }
            catch (PlbIOException& exception) {
                if (solver_type[iT]==2 ) { pcout << exception.what() << " It must be defined when solver_type is Cellular Automata. Terminating the simulation.\n"; return -1; }
            }
            try {
                std::vector<T> tmp;
                doc["parameters"]["microbiology"][bioname]["half_saturation_constants"].read(tmp);
                if (tmp.size()!=(unsigned)num_of_substrates) { pcout << "The length of half_saturation_constants should be equal to the number_of_substates. Terminating the simulation.\n"; }
                else { vec_Kc.push_back(tmp); }
            }
            catch (PlbIOException& exception) {
                vec_Kc.push_back( std::vector<T> (1,-99) );
            }
            try {
                std::vector<T> vmax;
                doc["parameters"]["microbiology"][bioname]["maximum_uptake_flux"].read(vmax);
                vec_Vmax.push_back(vmax);
            }
            catch (PlbIOException& exception) { vec_Vmax.push_back( std::vector<T> (num_of_substrates,0) ); }
        }
        try { doc["parameters"]["microbiology"]["thrd_biofilm_fraction"].read(thrd_bFilmFrac); }
        catch (PlbIOException& exception) { if (ca_count>0 ) { pcout << exception.what() << " It must be defined when solver_type is Cellular Automata. Terminating the simulation.\n"; return -1; } }
        if (vec_microbes_names.size()!=(unsigned)num_of_microbes) { pcout << "The length of name_of_microbes vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
        if (vec_mu.size()!=(unsigned)num_of_microbes) { pcout << "The length of decay_coefficient vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
        if (bMass_poreD.size()!=(unsigned)num_of_microbes) { pcout << "The length of biomass_diffusion_coefficients in_pore vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
        if (bMass_bFilmD.size()!=(unsigned)num_of_microbes) { pcout << "The length of biomass_diffusion_coefficients in_biofilm vector does not match the num_of_microbes. Terminating the simulation.\n"; return -1; }
        if (vec_permRatio.size()!=(unsigned)bfilm_count) { pcout << "The length of biofilm_permeability_ratio vector does not match the length of the material number vector for microbes. Terminating the simulation.\n"; return -1; }

        try { doc["parameters"]["microbiology"]["maximum_biomass_density"].read(max_bMassRho); }
        catch (PlbIOException& exception) {
            if (ca_count==1 ) { pcout << exception.what() << " It must be defined when solver_type is Cellular Automata. Terminating the simulation.\n"; return -1; }
            else { max_bMassRho = 999999999.; }
        }
        try {
            std::string tmp;
            doc["parameters"]["microbiology"]["CA_method"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("fraction")==0 || tmp.compare("0")==0 || tmp.compare("no")==0 ) { halfflag=0; }
            else if (tmp.compare("half")==0 || tmp.compare("1")==0 || tmp.compare("yes")==0 ) { halfflag=1; }
            else { pcout << "halfflag (" << tmp << ") should be either half or fraction. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { halfflag=0; }
        } // END if (biotic_mode) - microbiology section

        // IO
        try {
            std::string tmp;
            doc["parameters"]["IO"]["read_ADE_file"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0 ) { read_ADE_file = 0; }
            else if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0 ) { read_ADE_file = 1; }
            else { pcout << "read_ADE_file (" << tmp << ") should be either true or false. Terminating the simulation.\n"; return -1; }
        }
        catch (PlbIOException& exception) { read_ADE_file=0; }
        try {
            std::string item;
            doc["parameters"]["IO"]["ns_filename"].read(item);
            ns_filename = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { ns_filename[i] = item[i]; }
            ns_filename[item.size()+1] = '\0';
        }
        catch (PlbIOException& exception) {
            std::string item = "nsLattice";
            ns_filename = (char *) calloc(item.size()+1,sizeof(char));
            for (size_t i = 0; i < item.size(); ++i) { ns_filename[i] = item[i]; }
            ns_filename[item.size()+1] = '\0';
        }
        try { doc["parameters"]["IO"]["mask_filename"].read(mask_filename); }
        catch (PlbIOException& exception) { mask_filename="maskLattice"; }
        try { doc["parameters"]["IO"]["subs_filename"].read(ade_filename); }
        catch (PlbIOException& exception) { ade_filename="subsLattice"; }
        try { doc["parameters"]["IO"]["bio_filename"].read(bio_filename); }
        catch (PlbIOException& exception) { bio_filename="bioLattice"; }
        try { doc["parameters"]["IO"]["save_VTK_interval"].read(ade_VTK_iTer); }
        catch (PlbIOException& exception) { ade_VTK_iTer=1000; }
        try { doc["parameters"]["IO"]["save_CHK_interval"].read(ade_CHK_iTer); }
        catch (PlbIOException& exception) { ade_CHK_iTer=1000000; }
        try { std::string dbg; doc["parameters"]["IO"]["debug_updRxn"].read(dbg);
              std::transform(dbg.begin(), dbg.end(), dbg.begin(), [](unsigned char c){ return std::tolower(c); });
              g_updRxn_debug = (dbg=="1"||dbg=="true"||dbg=="yes"||dbg=="on"); }
        catch (PlbIOException& exception) { g_updRxn_debug=false; }

        // ════════════════════════════════════════════════════════════════════════════
        // EQUILIBRIUM CHEMISTRY PARSING
        // ════════════════════════════════════════════════════════════════════════════
        // XML structure:
        // <equilibrium>
        //     <enabled>true</enabled>
        //     <components>HCO3- H+</components>
        //     <species>
        //         <species0 logK="0.0">0.0 0.0</species0>        <!-- O2: not in equilibrium -->
        //         <species1 logK="0.0">0.0 0.0</species1>        <!-- CH4: not in equilibrium -->
        //         <species2 logK="6.35">1.0 1.0</species2>       <!-- H2CO3 -->
        //         <species3 logK="0.0">1.0 0.0</species3>        <!-- HCO3- (component) -->
        //         <species4 logK="-10.33">1.0 -1.0</species4>    <!-- CO3-- -->
        //         <species5 logK="0.0">0.0 1.0</species5>        <!-- H+ (component) -->
        //     </species>
        // </equilibrium>
        useEquilibrium = false;
        try {
            std::string tmp;
            doc["parameters"]["equilibrium"]["enabled"].read(tmp);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){ return std::tolower(c); });
            if (tmp.compare("yes")==0 || tmp.compare("true")==0 || tmp.compare("1")==0 ) { 
                useEquilibrium = true;
                pcout << "\nEquilibrium chemistry: ENABLED\n";
            }
            else if (tmp.compare("no")==0 || tmp.compare("false")==0 || tmp.compare("0")==0 ) { 
                useEquilibrium = false;
                pcout << "\nEquilibrium chemistry: DISABLED\n";
            }
            else { 
                pcout << "equilibrium/enabled (" << tmp << ") should be either true or false. Defaulting to false.\n"; 
                useEquilibrium = false;
            }
        }
        catch (PlbIOException& exception) { 
            useEquilibrium = false; 
        }

        if (useEquilibrium) {
            // Parse component names
            try {
                std::vector<std::string> components;
                doc["parameters"]["equilibrium"]["components"].read(components);
                eq_component_names = components;
                pcout << "Equilibrium components (" << eq_component_names.size() << "): ";
                for (const auto& name : eq_component_names) { pcout << name << " "; }
                pcout << std::endl;
            }
            catch (PlbIOException& exception) {
                pcout << "WARNING: equilibrium/components not specified. Equilibrium solver may not work correctly.\n";
            }

            // Parse stoichiometry matrix and logK values for each species
            eq_stoich_matrix.resize(num_of_substrates);
            eq_logK_values.resize(num_of_substrates, 0.0);
            
            plint num_components = eq_component_names.size();
            for (plint iS = 0; iS < num_of_substrates; ++iS) {
                eq_stoich_matrix[iS].resize(num_components, 0.0);
            }

            for (plint iS = 0; iS < num_of_substrates; ++iS) {
                std::string speciesTag = "species" + std::to_string(iS);
                
                // Read stoichiometry row
                try {
                    std::vector<T> stoich_row;
                    doc["parameters"]["equilibrium"]["stoichiometry"][speciesTag].read(stoich_row);
                    if (stoich_row.size() == (size_t)num_components) {
                        eq_stoich_matrix[iS] = stoich_row;
                    }
                    else {
                        pcout << "WARNING: Stoichiometry for " << speciesTag << " has wrong size (" 
                              << stoich_row.size() << " vs " << num_components << " components). Using zeros.\n";
                    }
                }
                catch (PlbIOException& exception) {
                    // Species not defined in equilibrium - keep zeros (not participating)
                }

                // Read logK value
                try {
                    T logK_val;
                    std::string logKTag = "species" + std::to_string(iS);
                    doc["parameters"]["equilibrium"]["logK"][logKTag].read(logK_val);
                    eq_logK_values[iS] = logK_val;
                }
                catch (PlbIOException& exception) {
                    eq_logK_values[iS] = 0.0;
                }
            }

            // Print parsed equilibrium data
            pcout << "Equilibrium stoichiometry matrix and logK values:\n";
            pcout << "  Species";
            for (const auto& comp : eq_component_names) { pcout << "\t" << comp; }
            pcout << "\tlogK\n";
            for (plint iS = 0; iS < num_of_substrates; ++iS) {
                pcout << "  " << vec_subs_names[iS];
                for (plint iC = 0; iC < num_components; ++iC) {
                    pcout << "\t" << eq_stoich_matrix[iS][iC];
                }
                pcout << "\t" << eq_logK_values[iS] << "\n";
            }
            pcout << std::endl;
        }
        // ════════════════════════════════════════════════════════════════════════════

        if (vec_Kc.size()>0) {
            for (plint iT=0; iT<num_of_microbes; ++iT) {
                std::vector<T> Vmax0=vec_Vmax[iT]; std::vector<T> Kc0=vec_Kc[iT]; T mu0=vec_mu[iT];
                vec_Kc_kns.push_back(Kc0); vec_mu_kns.push_back(mu0); vec_Vmax_kns.push_back(Vmax0);
            }
        }
    }
    catch (PlbIOException& exception) {
        pcout << exception.what() << " Terminating the simulation.\n" << std::endl;
        return -1;
    }
    return 0;
}

#endif
