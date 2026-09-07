/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =================================================================================================
 * complab3d_outerfaces.hh  --  what happens at the four sides of the box
 *
 * THE DEFECT, STATED PLAINLY
 *
 *   The solver gives the domain a boundary condition on two of its six faces. west (x = 0) and
 *   east (x = nx-1) get an inlet and an outlet, set per substrate in CompLaB.xml. The other four
 *   -- y = 0, y = ny-1, z = 0, z = nz-1 -- get nothing at all. Not a wall, not a symmetry plane,
 *   not periodicity. Nothing.
 *
 *   Where the geometry happens to put a wall voxel on one of those faces the omission does not
 *   matter, because bounce-back closes it. Where it does not, the advection-diffusion lattice
 *   streams populations off the edge of the block every step, and the populations that should
 *   arrive from outside are read out of an envelope nothing updates. Mass leaves and what comes
 *   back is stale.
 *
 *   Every example shipped with this program leaves at least one outer face open. The generator
 *   for examples 05 to 08 even describes its grain blocks as sitting "against the y = 0 wall" --
 *   a wall it never writes.
 *
 * HOW BIG IT IS
 *
 *   Measured on example 07 with every reaction switched off, so the only correct answer is that
 *   the patch spreads and the total never changes. Biomass total over 200 steps, at a relaxation
 *   time well away from 0.5 so the ringing described in complab.cpp is not in the way:
 *
 *       faces open      108.000 -> 104.000 , and flat from there
 *       faces walled     72.000 ->  65.039 , and flat from there
 *
 *   Both settle rather than draining. Most of both deficits turned out not to be lost mass at all
 *   but a measurement artefact -- see the `held` column in the scalar record, which is what
 *   computeDensity() cannot see. With that accounted for, and with the four faces actually
 *   closed, the same case conserves biomass to 1.3e-14 over 1000 steps.
 *
 *   So this is a defect of definition more than of magnitude: an unconditioned face is not a
 *   boundary condition, and a result that depends on what the envelope happened to hold is not
 *   reproducible, even when it is close.
 *
 * HOW IT IS FIXED, AND WHERE
 *
 *   Not here. Every example in this repository now draws its own confining wall: each
 *   preprocess.py calls pad_closed_faces(), which ADDS one layer of inert wall outside any face
 *   that is still open, and grows NY and NZ by two to hold it. Adding a layer rather than
 *   converting one is the whole point -- the pore space, the grains and the patches keep exactly
 *   the volumes the case declares, where walling the two z faces of a six-deep slab would have
 *   spent a third of the pore space on the boundary condition. The x faces need nothing: they
 *   already carry the inlet and the outlet, and a species that should not leave through them now
 *   asks for <left_boundary_type>closed</left_boundary_type>, which is bounce-back and exactly
 *   conservative.
 *
 *   What is left in this file is the part that cannot live in a generator: a run reads whatever
 *   geometry it is given, including one written years ago or by hand, so the solver counts the
 *   open voxels on those four faces at start-up and says what follows from them.
 *
 * WHAT THIS FILE DOES
 *
 *       <LB_numerics><domain><outer_faces>open</outer_faces></domain></LB_numerics>
 *
 *         open    (default) leave the geometry exactly as written. Every open voxel on the four
 *                 faces is counted and reported at start-up, with what it means.
 *         sealed  turn every voxel on those faces that is not already wall or solid into inert
 *                 wall, on the geometry field itself, before any lattice is built from it.
 *
 *   `open` is the default because `sealed` DELETES voxels the user declared, and the right answer
 *   is almost always to pad the geometry instead. `sealed` is for a geometry you cannot regenerate.
 *
 * WHAT IS NOT DONE, AND WHY
 *
 *   A zero-gradient closure on those four faces -- exact for a field that does not vary across
 *   the closed direction, and free of any volume cost -- was tried twice and is not here. Palabos
 *   has the machinery: addTemperatureBoundary1N/1P/2N/2P with boundary::neumann, the same calls
 *   the x faces use. The obstacle is that the boundary dynamics is a composite built on
 *   NoDynamics, which answers computeDensity() from a stored number and ignores the populations,
 *   so initializeAtEquilibrium cannot reach those cells. Left alone they hold zero populations,
 *   which in the deviation form the advection-diffusion descriptors use reads as density 1, and a
 *   face at density 1 is a source: example 07's biomass ran 108 -> 3310 with every voxel
 *   saturated. Setting the stored density instead pins the face, and the same case drained
 *   72 -> 0.02. Padding the geometry gives the same physics with none of that, which is why the
 *   generators do it and this file does not.
 *
 *   Periodicity is not offered either. Palabos would wrap the lattices happily, but the cellular
 *   automaton, the finite-difference biomass step and the dissolution gather all walk their
 *   neighbours with hand-written index arithmetic that stops at the edge and does not wrap. A
 *   periodic option would be true for the lattice-Boltzmann fields and false for everything else
 *   in the same run, which is worse than not offering it.
 * ================================================================================================= */
#ifndef COMPLAB3D_OUTERFACES_HH
#define COMPLAB3D_OUTERFACES_HH

#include <string>
#include <vector>

namespace complab_faces {

/* Turn every voxel in `domain` that is neither wall nor solid into wall. Applied only to the four
 * face slabs, so the cost is four thin boxes and not a sweep of the field. */
template <typename T1>
class sealFace3D : public BoxProcessingFunctional3D_S<T1>
{
public:
    sealFace3D(plint bb_, plint solid_) : bb(bb_), solid(solid_) {}

    virtual void process(Box3D domain, ScalarField3D<T1> &field)
    {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    const plint m = (plint) field.get(iX, iY, iZ);
                    if (m != bb && m != solid) field.get(iX, iY, iZ) = (T1) bb;
                }
    }
    virtual sealFace3D<T1>* clone() const { return new sealFace3D<T1>(*this); }
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; }
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
private:
    plint bb, solid;
};

/* The four faces, in one place, so the count and the seal can never disagree about which they are. */
inline void faceBoxes(plint nx, plint ny, plint nz, Box3D out[4])
{
    /* [v1.3] The y=0 and y=ny-1 boxes used to span the whole z range while the z=0 and z=nz-1
     * boxes spanned the whole y range, so all four met on the four x-parallel edges and every
     * voxel there was counted twice -- 4*nx of them. Sealing is idempotent so the seal itself was
     * fine, but both numbers in the "N of M voxels are open" line were inflated, and unequally:
     * the denominator always by 4*nx, the numerator only by whichever edge voxels were open. The
     * z boxes now stop one voxel short of each y face, so the four are disjoint and their volumes
     * sum to the true face-voxel count. */
    out[0] = Box3D(0, nx-1, 0,    0,    0,    nz-1);
    out[1] = Box3D(0, nx-1, ny-1, ny-1, 0,    nz-1);
    out[2] = Box3D(0, nx-1, 1,    ny-2, 0,    0   );
    out[3] = Box3D(0, nx-1, 1,    ny-2, nz-1, nz-1);
}

/* How many voxels on those faces are open. MaskedScalarCounts3D is the reduction the rest of the
 * program already uses for porosity, so this number and the porosity come from one counter and
 * are correct on any number of ranks. */
inline plint countOpenOnFaces(MultiScalarField3D<int> &geometry,
                              plint nx, plint ny, plint nz, plint bb, plint solid)
{
    Box3D f[4];
    faceBoxes(nx, ny, nz, f);
    plint open = 0;
    for (int i = 0; i < 4; ++i) {
        const plint total = (f[i].x1-f[i].x0+1) * (f[i].y1-f[i].y0+1) * (f[i].z1-f[i].z0+1);
        plint closed = 0;
        if (bb    >= 0)                 closed += MaskedScalarCounts3D(f[i], geometry, bb);
        if (solid >= 0 && solid != bb)  closed += MaskedScalarCounts3D(f[i], geometry, solid);
        open += total - closed;
    }
    return open;
}

/* mode is whatever <outer_faces> said. Returns false only for a word this program does not know,
 * so the caller stops rather than running something the user did not ask for. */
inline bool apply(MultiScalarField3D<int> &geometry,
                  plint nx, plint ny, plint nz, plint bb, plint solid,
                  const std::string &raw)
{
    std::string mode = raw;
    for (size_t i = 0; i < mode.size(); ++i)
        mode[i] = (char) std::tolower((unsigned char) mode[i]);
    if (mode != "open" && mode != "sealed") {
        pcout << "  <outer_faces> must be `open` or `sealed`, not `" << raw << "`. Terminating.\n";
        return false;
    }

    const plint open = countOpenOnFaces(geometry, nx, ny, nz, bb, solid);
    /* [v1.3] 2*nx*nz + 2*nx*ny counted the four x-parallel edges twice. faceBoxes() is now
     * disjoint, so the true count is the sum of the four box volumes: the two y faces whole, plus
     * the two z faces less the row they share with each y face. */
    const plint faceTotal = 2*nx*nz + 2*nx*(ny-2);

    if (open == 0) {
        pcout << "  [GEOM] outer faces y=0, y=" << (ny-1) << ", z=0, z=" << (nz-1)
              << ": closed by the geometry. Nothing to do.\n";
        return true;
    }

    if (mode == "sealed") {
        if (bb < 0) {
            pcout << "  <outer_faces> is `sealed` but no <bounce_back> material number is defined, "
                  << "so there is nothing\n  to seal the faces WITH. Give <material_numbers> a "
                  << "<bounce_back> entry, or leave <outer_faces> at `open`.\n  Terminating.\n";
            return false;
        }
        Box3D f[4];
        faceBoxes(nx, ny, nz, f);
        for (int i = 0; i < 4; ++i)
            applyProcessingFunctional(new sealFace3D<int>(bb, solid), f[i], geometry);
        pcout << "  [GEOM] outer faces sealed: " << open << " of " << faceTotal
              << " voxel(s) on y=0, y=" << (ny-1) << ", z=0 and z=" << (nz-1)
              << " turned into wall (material " << bb << ").\n"
              << "         Those voxels are gone from the pore space, so porosity and every total "
              << "below are for the\n         closed domain, not for the geometry file. x=0 and x="
              << (nx-1) << " keep the inlet and outlet.\n";
        return true;
    }

    pcout << "  [GEOM] WARNING: " << open << " of " << faceTotal
          << " voxel(s) on y=0, y=" << (ny-1) << ", z=0 and z=" << (nz-1)
          << " are open pore and have\n"
          << "         NO boundary condition. Only x=0 and x=" << (nx-1)
          << " are given one. The advection-diffusion\n"
          << "         lattices stream off the block at those voxels and read back an envelope "
          << "nothing updates,\n"
          << "         so what crosses them is undefined: mass is not conserved there and the "
          << "answer depends on\n"
          << "         the block decomposition. Measured on example 07 with no reaction, the "
          << "biomass total settles\n"
          << "         about 4% below where it started and stays there.\n"
          << "         Two ways to make the domain mean something:\n"
          << "           1. draw the confining wall in the geometry, which is almost always what "
          << "was intended;\n"
          << "           2. set <outer_faces>sealed</outer_faces> and let the solver wall them -- "
          << "but read the count\n"
          << "              above first, because on a thin slab those voxels are a large part of "
          << "the pore space.\n";
    return true;
}

}  // namespace complab_faces

#endif
