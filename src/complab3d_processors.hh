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

#ifndef COMPLAB3D_PROCESSORS_HH
#define COMPLAB3D_PROCESSORS_HH

#include <cmath>
#include <limits>
#include <vector>
#include <random>
#include <iostream>

using namespace plb;
typedef double T;

/* ===============================================================================================================
   ================================================ CONSTANTS ====================================================
   =============================================================================================================== */

// Threshold value for numerical comparisons
const T thrd = 1e-12;

/* ===============================================================================================================
   ============================================ LATTICE DESCRIPTORS ==============================================
   =============================================================================================================== */

#define NSDES descriptors::D3Q19Descriptor      // Navier-Stokes: 
#define RXNDES descriptors::AdvectionDiffusionD3Q7Descriptor  // Reaction-Diffusion: 

/* ===============================================================================================================
   =========================================== INCLUDE PROCESSOR PARTS ===========================================
   =============================================================================================================== */

// Part 1: Kinetics and biomass redistribution processors
//   - run_kinetics
//   - update_rxnLattices
//   - pushExcessBiomass3D
//   - halfPushExcessBiomass3D
//   - pullExcessBiomass3D
#include "complab3d_processors_part1.hh"

// Part 2: Mask and dynamics update processors
//   - updateLocalMaskNtotalLattices3D
//   - fdDiffusion3D
//   - updateSoluteDynamics3D
//   - updateBiomassDynamics3D
//   - updateNsLatticesDynamics3D
//   - updateAgeDistance3D
#include "complab3d_processors_part2.hh"

// Part 3: Copy/initialize processors and reductive functionals
//   - CopyGeometryScalar2maskLattice3D
//   - CopyGeometryScalar2ageLattice3D
//   - CopyGeometryScalar2distLattice3D
//   - CopyLattice2ScalarField3D
//   - initializeScalarLattice3D
//   - stabilizeADElattice3D
//   - createDistanceDomain3D
//   - createAgeDomain3D
//   - MaskedBoxScalarCountFunctional3D
//   - BoxLatticeRMSEFunctional3D
//   - MaskedScalarCounts3D (helper function)
//   - computeRMSE3D (helper function)
#include "complab3d_processors_part3.hh"

// Part 4: Equilibrium chemistry solver and processors (Anderson Acceleration + PCF)
//   - EquilibriumChemistry<T> class
//   - run_equilibrium_biotic
//   - run_equilibrium_full  
//   - update_equilibriumLattices
//   - reset_deltaLattices
#include "complab3d_processors_part4_eqsolver.hh"

// Part 5: OPTIONAL metabolic layer -- flux balance analysis and surrogate models.
//   Nothing below is compiled unless the matching CMake option is on, and
//   nothing runs unless the matching switch in <simulation_mode> is on.  A run
//   that asks for none of them behaves exactly like CompLB3D before this layer
//   existed.  See complab3d_metabolic.hh for the units and the XML reference.
//
//   - MetabolicConfig / MetabolicTotals   configuration, model loading, unit conversion
//   - runFBA_glpk3D                       one LP per voxel via GLPK        (ENABLE_GLPK)
//   - runFBA_cobrapy3D                    one LP per voxel via COBRApy     (ENABLE_COBRAPY)
//   - run_surrogate3D                     pre-trained stand-in for FBA     (always built)
#include "complab3d_metabolic.hh"
#include "complab3d_processors_fba.hh"
#include "complab3d_processors_surrogate.hh"

/* ===============================================================================================================
   ============================================= PROCESSOR SUMMARY ================================================
   ===============================================================================================================
   
   DATA PROCESSORS (LatticeBoxProcessingFunctional3D):
   ---------------------------------------------------
   1.  run_kinetics                    - Execute user-defined reaction kinetics
   2.  update_rxnLattices              - Update reaction lattices to concentration/biomass
   3.  pushExcessBiomass3D             - Redistribute excess biomass (push to neighbors)
   4.  halfPushExcessBiomass3D         - Redistribute half of excess biomass
   5.  pullExcessBiomass3D             - Redistribute biomass (pull from neighbors, MPI)
   6.  updateLocalMaskNtotalLattices3D - Update mask numbers and total biomass
   7.  fdDiffusion3D                   - Finite difference diffusion for biomass
   8.  updateSoluteDynamics3D          - Update solute diffusivity at biomass voxels
   9.  updateBiomassDynamics3D         - Update planktonic biomass diffusivity
   10. updateAgeDistance3D             - Update age/distance lattice for biofilm tracking
   
   BOX PROCESSING FUNCTIONALS (Two Lattices):
   ------------------------------------------
   11. updateNsLatticesDynamics3D      - Update Navier-Stokes flow field dynamics
   
   BOX PROCESSING FUNCTIONALS (Lattice + Scalar):
   ----------------------------------------------
   12. CopyGeometryScalar2maskLattice3D - Copy geometry scalar to mask lattice
   13. CopyGeometryScalar2ageLattice3D  - Copy geometry scalar to age lattice
   14. CopyGeometryScalar2distLattice3D - Copy geometry scalar to distance lattice
   15. CopyLattice2ScalarField3D        - Copy lattice values to scalar field
   16. initializeScalarLattice3D        - Initialize scalar biomass lattice
   17. stabilizeADElattice3D            - Stabilize advection-diffusion lattice
   
   BOX PROCESSING FUNCTIONALS (Scalar Only):
   -----------------------------------------
   18. createDistanceDomain3D           - Create distance domain scalar field
   19. createAgeDomain3D                - Create age domain scalar field
   
   REDUCTIVE DATA PROCESSORS:
   --------------------------
   20. MaskedBoxScalarCountFunctional3D - Count scalar values matching a mask
   21. BoxLatticeRMSEFunctional3D       - Calculate RMSE between two lattices
   
   HELPER FUNCTIONS:
   -----------------
   22. MaskedScalarCounts3D             - Wrapper for MaskedBoxScalarCountFunctional3D
   23. computeRMSE3D                    - Wrapper for BoxLatticeRMSEFunctional3D

   EQUILIBRIUM CHEMISTRY (Part 4):
   -------------------------------
   24. EquilibriumChemistry<T>          - Anderson Acceleration + PCF solver class
   25. run_equilibrium_biotic           - Direct equilibrium update (RECOMMENDED)
   26. run_equilibrium_full             - Equilibrium with delta lattices
   27. update_equilibriumLattices       - Apply delta changes from equilibrium
   28. reset_deltaLattices              - Reset delta lattices to zero
   
   ===============================================================================================================
   ========================================== D3Q7 POPULATION WEIGHTS ============================================
   ===============================================================================================================
   
   For D3Q7 lattice (AdvectionDiffusionD3Q7Descriptor):
   - Center (i=0):      w_0 = 1/4
   - Face neighbors (i=1-6): w_i = 1/8
   
   Population distribution for density rho:
   g[0] = (rho - 1) * 1/4
   g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (rho - 1) * 1/8
   
   Velocity directions for D3Q7:
   i=0: ( 0,  0,  0)  center
   i=1: ( 1,  0,  0)  +x
   i=2: (-1,  0,  0)  -x
   i=3: ( 0,  1,  0)  +y
   i=4: ( 0, -1,  0)  -y
   i=5: ( 0,  0,  1)  +z
   i=6: ( 0,  0, -1)  -z
   
   =============================================================================================================== */

#endif // COMPLAB3D_PROCESSORS_HH
