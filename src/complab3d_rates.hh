/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * complab3d_rates.hh  --  OPTIONAL REACTION RATE OUTPUT
 *
 * WHAT THIS IS FOR
 *   CompLaB writes concentrations and the flow field. It does not write how
 *   fast anything reacted. The reaction terms are computed inside
 *   defineKinetics.hh and defineAbioticKinetics.hh, used to move the
 *   concentrations, and then discarded. Anyone who wants a rate afterwards has
 *   to reconstruct it by differencing the saved snapshots, which gives
 *   transport plus reaction rather than reaction, at the output interval rather
 *   than at the time step, after a clamp that leaves no trace, and never
 *   separated per reaction. This file writes the rate directly instead.
 *
 * WHAT IT DOES NOT DO
 *   Nothing here is switched on by default. With <save_reaction_rates> absent
 *   or false, not one lattice is allocated and not one extra line of work is
 *   done, so every existing input file behaves exactly as it did before.
 *
 * HOW IT STAYS PROBLEM INDEPENDENT
 *   This file knows no chemical name, no reaction name and no rate law. The
 *   chemical and microbe names come from CompLaB.xml, the same place the
 *   solver already gets them. The reaction names, if there are any, come from
 *   the user's own kinetics header through the small opt-in contract below.
 *   Changing the chemistry means editing that header and nothing else.
 *
 * THE OPT-IN CONTRACT  (three lines in the user's defineKinetics.hh)
 *
 *     #define COMPLAB_HAS_RXN_RATES 1
 *     namespace RxnRates {
 *         inline const std::vector<std::string>& names() {
 *             static const std::vector<std::string> n = { "Geo_Fe", "Geo_U" };
 *             return n;                       // as many as you have, your order
 *         }
 *     }
 *     void defineRxnKinetics(std::vector<double> B, std::vector<double> C,
 *                            std::vector<double>& subsR, std::vector<double>& bioR,
 *                            plb::plint mask, std::vector<double>* rxnR);
 *
 *   and inside that function, once the rate laws have been evaluated:
 *
 *     if (rxnR) { (*rxnR)[0] = r_first; (*rxnR)[1] = r_second; }
 *
 *   defineAbioticKinetics.hh has the same contract under the names
 *   COMPLAB_HAS_ABIO_RXN_RATES and AbioRxnRates. The two files are opted in
 *   separately because they are swapped separately.
 *
 *   A header that declares none of this still gets the per chemical and per
 *   microbe rates, because those are sized from CompLaB.xml. It simply gets no
 *   per reaction channels. Nothing has to change in any existing header.
 *
 * ONE MORE RULE, AND IT MATTERS
 *   A call with rxnR set is a REPORTING call made by the sampler below, not a
 *   step of the simulation. The rates it returns are the same, but any running
 *   total the header keeps must not count it: those totals are printed as run
 *   statistics and a reporting call would inflate them. So a header that keeps
 *   statistics guards them:
 *
 *     if (!rxnR) KineticsStats::accumulate(...);
 *
 *   The sampler always passes a non null rxnR, even when the header declares no
 *   reactions, precisely so this test works everywhere.
 *
 * WHERE THE NUMBERS COME FROM
 *   The sampler calls the SAME functions the solver calls, on the SAME
 *   concentrations, at the point in the iteration where the solver is about to
 *   call them. It is not an estimate and not a re-run at a different state. It
 *   is the rate the solver used, one line before the solver used it.
 *
 * UNITS
 *   Whatever unit CompLaB.xml used for concentration, per second. The solver
 *   never converts and this file does not guess. <rate_unit> in the IO block
 *   is free text, recorded in the channel list beside the output so a reader a
 *   year later knows what the numbers are.
 *
 * COST
 *   One extra evaluation of the rate laws per voxel, on output iterations only,
 *   and one lattice per channel held for the length of the run. The lattice is
 *   the D3Q7 type the rest of the code uses, which is 7 doubles per voxel. The
 *   startup log prints the channel count and the memory before allocating, so a
 *   run that cannot afford it says so rather than failing later.
 * ============================================================================= */

#ifndef COMPLAB3D_RATES_HH
#define COMPLAB3D_RATES_HH

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

/* -----------------------------------------------------------------------------
 * How many reaction channels the two kinetics headers declare.
 * Zero when a header has not opted in, which is the default for every header
 * that existed before this feature.
 * ----------------------------------------------------------------------------- */
inline plb::plint complabBioRxnCount() {
#ifdef COMPLAB_HAS_RXN_RATES
    return (plb::plint) RxnRates::names().size();
#else
    return (plb::plint) 0;
#endif
}

inline plb::plint complabAbioRxnCount() {
#ifdef COMPLAB_HAS_ABIO_RXN_RATES
    return (plb::plint) AbioRxnRates::names().size();
#else
    return (plb::plint) 0;
#endif
}

inline std::string complabBioRxnName(plb::plint k) {
#ifdef COMPLAB_HAS_RXN_RATES
    return RxnRates::names()[(size_t) k];
#else
    (void) k; return std::string("");
#endif
}

inline std::string complabAbioRxnName(plb::plint k) {
#ifdef COMPLAB_HAS_ABIO_RXN_RATES
    return AbioRxnRates::names()[(size_t) k];
#else
    (void) k; return std::string("");
#endif
}

/* =============================================================================
 * RateOutput
 *
 *   Everything the rest of the program needs to know about the rate output:
 *   whether it is on, how often, what the channels are called and in which
 *   order they sit in the lattice vector. Built once at startup by
 *   buildRateOutput() and then read only.
 *
 *   CHANNEL ORDER, fixed, and the sampler and the writer both rely on it:
 *
 *       [0            .. S-1        ]  bio_<chemical>     biotic, per chemical
 *       [S            .. S+M-1      ]  bio_<microbe>      biotic, per microbe
 *       [S+M          .. S+M+nrb-1  ]  R_<reaction>       biotic, per reaction
 *       [S+M+nrb      .. +S-1       ]  abio_<chemical>    abiotic, per chemical
 *       [                .. +nra-1  ]  Rabio_<reaction>   abiotic, per reaction
 *
 *   A block that is switched off contributes no channels at all, so the
 *   offsets below are the only correct way to find a channel.
 * ============================================================================= */
struct RateOutput {
    bool        enabled;
    bool        doBio;          // biotic kinetics contributes channels
    bool        doAbio;         // abiotic kinetics contributes channels
    int         surfaceOnly;    // mirror of <precipitation><surface_only>
    plb::plint  interval;       // iterations between rate snapshots
    std::string filename;       // base name of the .vti files
    std::string unit;           // free text, recorded only

    plb::plint  S;              // chemicals
    plb::plint  M;              // microbes
    plb::plint  nrb;            // biotic reactions declared by the header
    plb::plint  nra;            // abiotic reactions declared by the header

    std::vector<std::string> channelName;
    std::vector<std::string> channelSource;   // "biotic" or "abiotic"
    std::vector<std::string> channelKind;     // "chemical", "microbe" or "reaction"

    RateOutput()
      : enabled(false), doBio(false), doAbio(false), surfaceOnly(0),
        interval(0), filename("rateLattice"), unit(""),
        S(0), M(0), nrb(0), nra(0) {}

    plb::plint nch()      const { return (plb::plint) channelName.size(); }
    plb::plint offBioC()  const { return 0; }
    plb::plint offBioM()  const { return S; }
    plb::plint offBioR()  const { return S + M; }
    plb::plint offAbioC() const { return doBio ? (S + M + nrb) : 0; }
    plb::plint offAbioR() const { return offAbioC() + S; }
};

/* -----------------------------------------------------------------------------
 * readRateSettings
 *
 *   Reads the four new entries of the IO block from CompLaB.xml. Every one of
 *   them is optional and every one has a default that reproduces the old
 *   behaviour exactly, so an input file written before this feature existed
 *   reads without a word of complaint and switches nothing on.
 *
 *       <IO>
 *           <save_reaction_rates>true</save_reaction_rates>   default false
 *           <rate_filename>rateLattice</rate_filename>        default rateLattice
 *           <rate_interval>0</rate_interval>   0 means follow save_VTK_interval
 *           <rate_unit>mol/L/s</rate_unit>     free text, recorded only
 *       </IO>
 *
 *   Read here from its own XMLreader rather than threaded through
 *   initialize_complab, which is how the immobile, precipitation and
 *   dissolution blocks are already read, and which keeps a forty argument
 *   signature from growing four more.
 * ----------------------------------------------------------------------------- */
inline void readRateSettings(bool& enabled, plb::plint& interval,
                             std::string& filename, std::string& unit)
{
    enabled  = false;
    interval = 0;
    filename = "rateLattice";
    unit     = "";

    try {
        plb::XMLreader rdoc("CompLaB.xml");

        try {
            std::string s;
            rdoc["parameters"]["IO"]["save_reaction_rates"].read(s);
            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            enabled = (s == "1" || s == "true" || s == "yes" || s == "on");
        } catch (plb::PlbIOException&) { enabled = false; }

        try { rdoc["parameters"]["IO"]["rate_filename"].read(filename); }
        catch (plb::PlbIOException&) { filename = "rateLattice"; }

        try { rdoc["parameters"]["IO"]["rate_interval"].read(interval); }
        catch (plb::PlbIOException&) { interval = 0; }

        try { rdoc["parameters"]["IO"]["rate_unit"].read(unit); }
        catch (plb::PlbIOException&) { unit = ""; }
    }
    catch (plb::PlbIOException&) {
        enabled = false;
    }
}

/* -----------------------------------------------------------------------------
 * buildRateOutput
 *
 *   Assembles the channel list. Called once, after CompLaB.xml has been read,
 *   so every name in it comes either from the input file or from the kinetics
 *   header. Nothing is hard coded.
 * ----------------------------------------------------------------------------- */
inline RateOutput buildRateOutput(
        bool enabled, plb::plint interval, std::string const& filename, std::string const& unit,
        bool bioActive, bool abioActive, int surfaceOnly,
        std::vector<std::string> const& subsNames,
        std::vector<std::string> const& microbeNames)
{
    RateOutput ro;
    ro.enabled     = enabled;
    ro.interval    = interval;
    ro.filename    = filename;
    ro.unit        = unit;
    ro.surfaceOnly = surfaceOnly;
    ro.doBio       = enabled && bioActive;
    ro.doAbio      = enabled && abioActive;
    ro.S           = (plb::plint) subsNames.size();
    ro.M           = (plb::plint) microbeNames.size();
    ro.nrb         = ro.doBio  ? complabBioRxnCount()  : 0;
    ro.nra         = ro.doAbio ? complabAbioRxnCount() : 0;

    if (!ro.enabled) return ro;

    if (ro.doBio) {
        for (plb::plint i = 0; i < ro.S; ++i) {
            ro.channelName.push_back("bio_" + subsNames[(size_t) i]);
            ro.channelSource.push_back("biotic");
            ro.channelKind.push_back("chemical");
        }
        for (plb::plint i = 0; i < ro.M; ++i) {
            ro.channelName.push_back("bio_" + microbeNames[(size_t) i]);
            ro.channelSource.push_back("biotic");
            ro.channelKind.push_back("microbe");
        }
        for (plb::plint i = 0; i < ro.nrb; ++i) {
            ro.channelName.push_back("R_" + complabBioRxnName(i));
            ro.channelSource.push_back("biotic");
            ro.channelKind.push_back("reaction");
        }
    }
    if (ro.doAbio) {
        for (plb::plint i = 0; i < ro.S; ++i) {
            ro.channelName.push_back("abio_" + subsNames[(size_t) i]);
            ro.channelSource.push_back("abiotic");
            ro.channelKind.push_back("chemical");
        }
        for (plb::plint i = 0; i < ro.nra; ++i) {
            ro.channelName.push_back("Rabio_" + complabAbioRxnName(i));
            ro.channelSource.push_back("abiotic");
            ro.channelKind.push_back("reaction");
        }
    }
    return ro;
}

/* -----------------------------------------------------------------------------
 * writeRateChannelList
 *
 *   One small text file written once beside the results, saying what each
 *   channel is. The .vti files carry the names too, but a reader that wants to
 *   know how many reactions a campaign had, and in what unit, should not have
 *   to parse a binary file to find out.
 * ----------------------------------------------------------------------------- */
inline void writeRateChannelList(RateOutput const& ro, std::string const& outputDir)
{
    if (!ro.enabled || ro.nch() == 0) return;
    if (!plb::global::mpi().isMainProcessor()) return;

    std::string path = outputDir + ro.filename + "_channels.txt";
    std::ofstream f(path.c_str());
    if (!f.is_open()) return;

    f << "# CompLaB reaction rate channels\n";
    f << "# one line per array in " << ro.filename << "_*.vti\n";
    f << "# unit of every channel: "
      << (ro.unit.empty() ? std::string("unstated, concentration unit of CompLaB.xml per second")
                          : ro.unit) << "\n";
    f << "# columns: index  name  source  kind\n";
    for (plb::plint k = 0; k < ro.nch(); ++k) {
        f << k << "\t" << ro.channelName[(size_t) k] << "\t"
          << ro.channelSource[(size_t) k] << "\t" << ro.channelKind[(size_t) k] << "\n";
    }
    f.close();
}

/* -----------------------------------------------------------------------------
 * writeRateVTI
 *
 *   One file per snapshot holding every channel, rather than one file per
 *   channel. VtkImageOutput3D takes as many arrays as you give it, which is
 *   what writeNsVTI already relies on for velocityNorm and velocity.
 *
 *   nx, ny and nz are passed in rather than read back from the lattice,
 *   deliberately: the four older writers in complab_functions.hh read ny from
 *   getNx(), which is only right on a cubic domain.
 * ----------------------------------------------------------------------------- */
template<typename T, template<typename U> class Descriptor>
void writeRateVTI(std::vector< MultiBlockLattice3D<T,Descriptor> >& rateLat,
                  RateOutput const& ro, plint iter, plint nx, plint ny, plint nz)
{
    if (!ro.enabled || ro.nch() == 0) return;

    VtkImageOutput3D<T> vtkOut(createFileName(ro.filename + "_", iter, 7), 1.);
    for (plint k = 0; k < ro.nch(); ++k) {
        vtkOut.template writeData<float>(
            *computeDensity(rateLat[(size_t) k], Box3D(1,nx-2, 0,ny-1, 0,nz-1)),
            ro.channelName[(size_t) k], 1.0);
    }
}

/* =============================================================================
 * sample_reaction_rates
 *
 *   Reads the concentrations and the biomass, calls the same two kinetics
 *   functions the solver calls, and stores what they returned in one lattice
 *   per channel. It changes nothing else: the concentration lattices and the
 *   dC lattices are read only here.
 *
 *   LATTICE ORDER, and complab.cpp builds the vector in exactly this order:
 *
 *       [0        .. S-1      ]  substrate lattices        (read)
 *       [S        .. S+M-1    ]  biomass lattices          (read)
 *       [S+M      .. S+M+NCH-1]  rate lattices             (written)
 *       [S+M+NCH             ]   mask lattice              (read)
 *
 *   The rate lattices are cleared here rather than by a copy of a zero lattice,
 *   so no second set has to be held in memory the way dC0 is.
 *
 *   appliesTo() is bulk, not bulkAndEnvelope, because with the surface gate on
 *   this processor reads face neighbours of the mask, and a neighbour of an
 *   envelope cell is outside the block. The output box never includes the
 *   envelope, so nothing is lost.
 * ============================================================================= */
template<typename T, template<typename U> class Descriptor>
class sample_reaction_rates : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    sample_reaction_rates(plint nx_, plint ny_, plint nz_,
                          plint subsNum_, plint bioNum_, plint nch_,
                          plint solid_, plint bb_,
                          bool doBio_, bool doAbio_, int surfaceOnly_,
                          plint nrb_, plint nra_,
                          plint offBioC_, plint offBioM_, plint offBioR_,
                          plint offAbioC_, plint offAbioR_)
      : nx(nx_), ny(ny_), nz(nz_), subsNum(subsNum_), bioNum(bioNum_), nch(nch_),
        solid(solid_), bb(bb_), doBio(doBio_), doAbio(doAbio_), surfaceOnly(surfaceOnly_),
        nrb(nrb_), nra(nra_),
        offBioC(offBioC_), offBioM(offBioM_), offBioR(offBioR_),
        offAbioC(offAbioC_), offAbioR(offAbioR_),
        rateLoc(subsNum_ + bioNum_), maskLloc(subsNum_ + bioNum_ + nch_)
    {}

    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices)
    {
        Dot3D absoluteOffset = lattices[0]->getLocation();
        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);

        // face neighbours, used only by the surface gate
        const plint dloc[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

        // displacements of every block relative to the first, computed once
        std::vector<Dot3D> off;
        for (plint iT = 0; iT < maskLloc; ++iT) {
            off.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));
        }

        std::vector<T> conc(subsNum, T()), bmass(bioNum, T());
        std::vector<T> subs_rate(subsNum, T()), bio_rate(bioNum, T());
        std::vector<T> abio_rate(subsNum, T());
        std::vector<double> rxn_rate((size_t) (nrb > 0 ? nrb : 1), 0.0);
        std::vector<double> arxn_rate((size_t) (nra > 0 ? nra : 1), 0.0);

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint absY = iY + absoluteOffset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint absZ = iZ + absoluteOffset.z;

                    // Clear every channel at this voxel first, so a voxel that
                    // does not react reads as zero rather than as whatever the
                    // previous rate snapshot left there.
                    for (plint k = 0; k < nch; ++k) setChannel(lattices, off, iX, iY, iZ, k, T());

                    if (!(absX > 0 && absX < nx-1)) continue;

                    plint iXm = iX + maskOffset.x, iYm = iY + maskOffset.y, iZm = iZ + maskOffset.z;
                    plint mask = util::roundToInt(lattices[maskLloc]->get(iXm,iYm,iZm).computeDensity());
                    if (mask == solid || mask == bb) continue;

                    // ---- read the state, exactly as run_kinetics reads it ----
                    for (plint iS = 0; iS < subsNum; ++iS) {
                        T c0 = lattices[iS]->get(iX+off[iS].x, iY+off[iS].y, iZ+off[iS].z).computeDensity();
                        if (c0 < thrd) c0 = 0;
                        conc[(size_t) iS] = c0;
                    }
                    for (plint iM = 0; iM < bioNum; ++iM) {
                        T b0 = lattices[subsNum+iM]->get(iX+off[subsNum+iM].x, iY+off[subsNum+iM].y,
                                                         iZ+off[subsNum+iM].z).computeDensity();
                        if (b0 < thrd) b0 = 0;
                        bmass[(size_t) iM] = b0;
                    }

                    // ---- biotic ----
                    if (doBio) {
                        for (plint i = 0; i < subsNum; ++i) subs_rate[(size_t) i] = T();
                        for (plint i = 0; i < bioNum;  ++i) bio_rate[(size_t) i]  = T();
                        for (size_t i = 0; i < rxn_rate.size(); ++i) rxn_rate[i] = 0.0;
#ifdef COMPLAB_HAS_RXN_RATES
                        // rxnR is always non null here, even when the header
                        // declares no reactions, so the header can use it as the
                        // signal that this call is a report and not a step.
                        defineRxnKinetics(bmass, conc, subs_rate, bio_rate, mask, &rxn_rate);
#else
                        defineRxnKinetics(bmass, conc, subs_rate, bio_rate, mask);
#endif
                        for (plint i = 0; i < subsNum; ++i)
                            setChannel(lattices, off, iX, iY, iZ, offBioC + i, subs_rate[(size_t) i]);
                        for (plint i = 0; i < bioNum; ++i)
                            setChannel(lattices, off, iX, iY, iZ, offBioM + i, bio_rate[(size_t) i]);
                        for (plint i = 0; i < nrb; ++i)
                            setChannel(lattices, off, iX, iY, iZ, offBioR + i, (T) rxn_rate[(size_t) i]);
                    }

                    // ---- abiotic, with the same surface gate the solver uses ----
                    if (doAbio) {
                        bool isInterface = (surfaceOnly == 0);
                        if (surfaceOnly != 0) {
                            for (int d = 0; d < 6 && !isInterface; ++d) {
                                plint ax = absX+dloc[d][0], ay = absY+dloc[d][1], az = absZ+dloc[d][2];
                                if (ax<0||ax>nx-1||ay<0||ay>ny-1||az<0||az>nz-1) continue;
                                plint nb = util::roundToInt(lattices[maskLloc]->get(
                                    iXm+dloc[d][0], iYm+dloc[d][1], iZm+dloc[d][2]).computeDensity());
                                if (nb == solid || nb == bb) isInterface = true;
                            }
                        }
                        if (isInterface) {
                            for (plint i = 0; i < subsNum; ++i) abio_rate[(size_t) i] = T();
                            for (size_t i = 0; i < arxn_rate.size(); ++i) arxn_rate[i] = 0.0;
#ifdef COMPLAB_HAS_ABIO_RXN_RATES
                            defineAbioticRxnKinetics(conc, abio_rate, mask, &arxn_rate);
#else
                            defineAbioticRxnKinetics(conc, abio_rate, mask);
#endif
                            for (plint i = 0; i < subsNum; ++i)
                                setChannel(lattices, off, iX, iY, iZ, offAbioC + i, abio_rate[(size_t) i]);
                            for (plint i = 0; i < nra; ++i)
                                setChannel(lattices, off, iX, iY, iZ, offAbioR + i, (T) arxn_rate[(size_t) i]);
                        }
                    }
                }
            }
        }
    }

    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }

    virtual sample_reaction_rates<T,Descriptor>* clone() const {
        return new sample_reaction_rates<T,Descriptor>(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iT = 0; iT < (plint) modified.size(); ++iT) modified[iT] = modif::nothing;
        for (plint k = 0; k < nch; ++k) modified[rateLoc + k] = modif::staticVariables;
    }

private:
    /* Store one value so that computeDensity() of that lattice returns it.
     * The 1/4 and 1/8 split is the D3Q7 weighting the rest of the code uses for
     * the same purpose, and the seven weights sum to one. */
    void setChannel(std::vector<BlockLattice3D<T,Descriptor>*>& lattices,
                    std::vector<Dot3D> const& off,
                    plint iX, plint iY, plint iZ, plint k, T value)
    {
        plint L = rateLoc + k;
        plint jX = iX + off[L].x, jY = iY + off[L].y, jZ = iZ + off[L].z;
        Array<T,7> g;
        g[0] = value/(T)4;
        g[1] = value/(T)8; g[2] = value/(T)8; g[3] = value/(T)8;
        g[4] = value/(T)8; g[5] = value/(T)8; g[6] = value/(T)8;
        lattices[L]->get(jX,jY,jZ).setPopulations(g);
    }

    plint nx, ny, nz, subsNum, bioNum, nch, solid, bb;
    bool  doBio, doAbio;
    int   surfaceOnly;
    plint nrb, nra;
    plint offBioC, offBioM, offBioR, offAbioC, offAbioR;
    plint rateLoc, maskLloc;
};

#endif // COMPLAB3D_RATES_HH
