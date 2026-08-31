/* ============================================================================
 * surrogateModel.hh  —  USER-EDITABLE SURROGATE (MACHINE-LEARNED) GROWTH MODEL
 *                        shipped example: Geobacter metallireducens GS-15 (iAF987)
 * ============================================================================
 *
 * HOW TO READ THIS FILE (for readers who do not program):
 *
 *   - This is a "recipe" file, the sibling of defineKinetics.hh. It is not a
 *     program you run on its own. The CompLaB3D solver opens it and calls the
 *     one routine inside (defineSurrogateModel) ONE time for every little cube
 *     of the domain (each cube is a "voxel"), on every time step, for every
 *     microbe that is present there. The solver hands the routine an estimate
 *     of how fast that microbe could take up each dissolved substrate in that
 *     voxel, and the routine hands back how fast the microbe grows, and
 *     (optionally) a corrected substrate uptake.
 *
 *   - Anything after "//", or between the block-comment markers, is a note for
 *     humans. The computer ignores it. Read it as the explanation of the line.
 *
 *   - EDITING THIS FILE REQUIRES RECOMPILING. Unlike CompLaB.xml, which the
 *     solver reads at start-up, this file is baked into the executable when it
 *     is built. After any change here you must rebuild (./comp.sh, or your
 *     cmake/make step) before the change has any effect. If you edit this file
 *     and the results do not change, you almost certainly forgot to rebuild.
 *
 * ----------------------------------------------------------------------------
 * WHAT IS A "SURROGATE MODEL", AND WHEN WOULD YOU USE ONE?
 * ----------------------------------------------------------------------------
 *   CompLaB3D can compute microbial growth in three ways:
 *
 *     1. KINETICS  (defineKinetics.hh) — you write the rate laws by hand
 *        (Monod terms, yields, decay). Cheapest, fully transparent, but you
 *        must know the rate laws.
 *
 *     2. FBA — flux balance analysis. At every voxel and every time step the
 *        solver builds and solves a linear program over the organism's whole
 *        genome-scale metabolic network (hundreds to thousands of reactions;
 *        see extractMM.py, which prepares that network). This is mechanistic
 *        and needs no fitted rate law, but it is EXPENSIVE: one linear program
 *        per microbe per active voxel per step. In a 3D domain of a few million
 *        voxels this quickly dominates the entire run time.
 *
 *     3. SURROGATE (this file) — you run the FBA offline, once, over the range
 *        of substrate uptake rates the simulation will visit, and you fit a
 *        cheap function (here a small neural network) to the answer:
 *
 *              substrate uptake rates  ->  predicted growth rate
 *
 *        At run time the solver then just evaluates that fitted function, which
 *        costs a few hundred multiplications instead of a linear program. Speed-
 *        ups of 100-1000x versus in-line FBA are typical.
 *
 *   USE A SURROGATE WHEN: you want FBA-quality growth predictions but the FBA
 *   cost is prohibitive; the number of substrates that actually limit growth is
 *   small (1-3 inputs is where these fits behave well); and the simulation stays
 *   inside the range of uptake rates the network was TRAINED on.
 *
 *   DO NOT USE A SURROGATE WHEN: you need the individual internal fluxes of the
 *   metabolic network (a surrogate only reproduces the objective, i.e. growth);
 *   or when your simulation wanders outside the training range — a neural
 *   network EXTRAPOLATES BADLY and will confidently return nonsense. Check the
 *   training range of the shipped net below before trusting it on new problems.
 *
 * ----------------------------------------------------------------------------
 * THE SHIPPED EXAMPLE
 * ----------------------------------------------------------------------------
 *   Organism : Geobacter metallireducens GS-15, genome-scale model iAF987.
 *   Trained  : offline FBA sweep, exported from MATLAB's Deep Learning Toolbox
 *              (the "genFunction"/"gensim" output of a trained feed-forward net).
 *   Shape    : 2 inputs
 *              -> 4 hidden layers of 10 neurons each, tansig ("S-shaped") units
 *              -> 1 linear output neuron
 *              with MATLAB "mapminmax" rescaling on the inputs and on the output.
 *   Inputs   : Fin[microbeId][0] and Fin[microbeId][1] — the Michaelis-Menten
 *              uptake-rate estimates for SUBSTRATE 0 and SUBSTRATE 1, in the order
 *              the substrates are listed in CompLaB.xml. For the shipped net these
 *              are the electron donor (acetate) and the electron acceptor.
 *   Output   : the specific growth rate, 1/h (equivalently gDW/gDW/h).
 *
 *   THE WEIGHTS AND BIASES BELOW ARE THE TRAINED MODEL. They are reproduced to
 *   full double precision exactly as MATLAB exported them; do not round them.
 *   To install YOUR OWN network, replace the weight blocks, the two mapminmax
 *   scaling blocks and the layer calls — the arithmetic below is the standard
 *   MATLAB feed-forward evaluation and needs no other change.
 *
 * ----------------------------------------------------------------------------
 * UNITS (get these wrong and everything else is wrong)
 * ----------------------------------------------------------------------------
 *   Fin, Fout : substrate uptake fluxes, mmol / gDW / h
 *               (millimoles of substrate per gram dry weight of cells per hour;
 *                this is the usual unit of an FBA exchange flux). POSITIVE means
 *                the microbe CONSUMES the substrate.
 *   bioR      : specific growth rate, 1/h. POSITIVE = growth. The solver turns
 *               it into a biomass change with  dB = bioR * biomass * dt/3600,
 *               so it converts from per-hour to the solver's per-second step.
 *   Concentrations elsewhere in the solver are mmol/L for substrates and
 *   kgDW/m3 (= gDW/L) for biomass.
 *   mask      : an integer label of the voxel material (pore / biofilm / solid /
 *               a user-defined region), same numbering as everywhere else in
 *               CompLaB.xml. The solver never calls this routine in solid or
 *               bounce-back voxels.
 *
 * ----------------------------------------------------------------------------
 * WHAT CHANGED RELATIVE TO THE 2D CompLaB VERSION  (read this if you are porting)
 * ----------------------------------------------------------------------------
 *   [FIX-3D] THE SIGNATURE. The 2D routine was
 *
 *       void defineSurrogateModel(std::vector<std::vector<double>> Fin,
 *                                 std::vector<std::vector<double>> &Fout,
 *                                 std::vector<double> &bioR);
 *
 *   and the shipped body read Fin[0][0], Fin[0][1] and wrote bioR[0]. In 2D the
 *   solver COMPACTED those arrays: it built them only for the microbes that had
 *   biomass above threshold IN THAT VOXEL, so "0" meant "the first microbe that
 *   happens to be present here" — a different organism in different voxels. With
 *   more than one surrogate microbe the network was therefore applied to an
 *   arbitrary one of them, and every other organism got bioR left at 0, which
 *   the solver reads as "no growth" and replaces with pure first-order decay.
 *   Whole populations silently died for no biological reason.
 *
 *   The 3D routine instead receives:
 *       * microbeId — the GLOBAL microbe index, in the order the microbes are
 *         listed in CompLaB.xml, so you always know WHICH organism you are
 *         being asked about and can branch on it;
 *       * Fin/Fout/bioR indexed by that same GLOBAL index, i.e. full length
 *         [nMicrobes], not compacted;
 *       * mask, so rates can differ by region.
 *     The routine is called once per microbe present in the voxel, and each call
 *     must write bioR[microbeId] (see the stoichiometry note below for Fout).
 *
 *   [FIX-3D] The helper routines below are now "static inline" (so several
 *   translation units can include this header without duplicate-symbol errors,
 *   and the compiler can inline the small loops), and a dimension mismatch now
 *   PRINTS A WARNING AND RETURNS ZEROS instead of calling exit(). Killing an MPI
 *   job from inside a data processor leaves the other ranks hanging and destroys
 *   the run; returning zeros makes the affected voxel grow at zero rate, which is
 *   visible in the output and recoverable.
 * ============================================================================ */
#ifndef SURROGATE_MODEL_HH   // include guard: the three #ifndef/#define/#endif lines
#define SURROGATE_MODEL_HH   //   make sure this file is read only once at build time

#include <vector>       // the "list of numbers" type used to pass fluxes in and rates out
#include <cmath>        // std::exp, used by the tansig activation function
#include <cstddef>      // std::size_t, the type used to count list entries
#include <iostream>     // lets us print a warning to the screen if something goes wrong

/* ---------------------------------------------------------------------------
 * WARNING PRINTER
 *   Prints a short message about a dimension mismatch. Capped at a few messages
 *   so that a mistake cannot flood the log with one line per voxel per step.
 *   Nothing here changes the simulation; it only tells you something is wrong.
 * --------------------------------------------------------------------------- */
static inline void surrogateDimWarning(const char* where, std::size_t got, std::size_t want)
{
    static int printed = 0;              // "static" = remembered between calls
    const int maxPrints = 10;
    if (printed < maxPrints) {
        ++printed;
        std::cout << "[surrogateModel] dimension mismatch in " << where
                  << ": got " << got << ", expected " << want
                  << ". Returning zeros (this voxel will not grow)."
                  << (printed == maxPrints ? " [further messages suppressed]" : "")
                  << std::endl;
    }
}

/* ---------------------------------------------------------------------------
 * SMALL MATRIX HELPERS
 *   These four little routines are all the arithmetic a feed-forward neural
 *   network needs. They were carried over from the 2D version unchanged in
 *   meaning; only the error handling and the "static inline" changed.
 * --------------------------------------------------------------------------- */

// matrix x vector:  result[i] = sum_j mat1[i][j] * mat2[j]
static inline std::vector<double> matMultp(const std::vector<std::vector<double>>& mat1,
                                           const std::vector<double>& mat2)
{
    const std::size_t r3 = mat1.size();
    if (r3 == 0) return std::vector<double>();
    const std::size_t c3 = mat1[0].size();
    const std::size_t m3 = mat2.size();
    if (m3 != c3) {                                  // [FIX-3D] warn, do not exit()
        surrogateDimWarning("matMultp", m3, c3);
        return std::vector<double>(r3, 0.0);
    }
    std::vector<double> mat3(r3);
    for (std::size_t i = 0; i < r3; ++i) {
        mat3[i] = 0.0;
        for (std::size_t j = 0; j < m3; ++j) {
            mat3[i] += mat1[i][j] * mat2[j];
        }
    }
    return mat3;
}

// dot product of two lists:  result = sum_i mat1[i] * mat2[i]
static inline double vecMultp(const std::vector<double>& mat1, const std::vector<double>& mat2)
{
    const std::size_t r3 = mat1.size();
    if (r3 != mat2.size()) {                         // [FIX-3D] warn, do not exit()
        surrogateDimWarning("vecMultp", mat2.size(), r3);
        return 0.0;
    }
    double d3 = 0.0;
    for (std::size_t i = 0; i < r3; ++i) {
        d3 += mat1[i] * mat2[i];
    }
    return d3;
}

// element-by-element sum of two lists:  result[i] = mat1[i] + mat2[i]
static inline std::vector<double> matSum(const std::vector<double>& mat1,
                                         const std::vector<double>& mat2)
{
    const std::size_t r3 = mat1.size();
    if (r3 != mat2.size()) {                         // [FIX-3D] warn, do not exit()
        surrogateDimWarning("matSum", mat2.size(), r3);
        return std::vector<double>(r3, 0.0);
    }
    std::vector<double> mat3(r3);
    for (std::size_t i = 0; i < r3; ++i) {
        mat3[i] = mat1[i] + mat2[i];
    }
    return mat3;
}

// tansig: MATLAB's "hyperbolic tangent sigmoid" neuron. It squashes any number
// into the range -1 .. +1 following an S-shaped curve. 2/(1+exp(-2n))-1 is
// exactly tanh(n); this is the form MATLAB exports, kept for comparability.
static inline std::vector<double> tansig(const std::vector<double>& n)
{
    const std::size_t row = n.size();
    std::vector<double> mat3(row);
    for (std::size_t i = 0; i < row; ++i) {
        mat3[i] = 2.0 / (1.0 + std::exp(-2.0 * n[i])) - 1.0;
    }
    return mat3;
}

/* ===========================================================================
 * THE MAIN SURROGATE ROUTINE
 *   Runs once per microbe per voxel per step.
 *
 * microbeId : global microbe index (as ordered in CompLaB.xml), so a user can
 *             branch on WHICH organism this call is for.
 * Fin       : [nMicrobes][nSubstrates] Michaelis-Menten uptake flux estimates, mmol/gDW/h
 * Fout      : same shape, PRE-FILLED with Fin; the user overwrites entries to
 *             change substrate consumption. Must be written or consumption stays at Fin.
 * bioR      : [nMicrobes] growth rate, 1/h. Write bioR[microbeId].
 * mask      : local material number, so rates can differ by region.
 * =========================================================================== */
void defineSurrogateModel(plb::plint microbeId,
                          const std::vector<std::vector<double>>& Fin,
                          std::vector<std::vector<double>>& Fout,
                          std::vector<double>& bioR,
                          plb::plint mask)
{
    /* -----------------------------------------------------------------------
     * SAFETY CHECKS. If anything is the wrong size we simply return without
     * writing: bioR[microbeId] then keeps the 0 the solver put there, which the
     * solver interprets as "no growth" and replaces with first-order decay.
     * ----------------------------------------------------------------------- */
    const std::size_t iM = static_cast<std::size_t>(microbeId);
    if (microbeId < 0) return;
    if (iM >= bioR.size() || iM >= Fin.size() || iM >= Fout.size()) {
        surrogateDimWarning("defineSurrogateModel (microbeId out of range)", iM, bioR.size());
        return;
    }

    // "mask" is available so that you can make growth depend on the region, for
    // example switching the network off inside a low-permeability inclusion:
    //     if (mask == 5) { bioR[microbeId] = 0.0; return; }
    // The shipped example uses the same network everywhere, so we simply mark
    // the argument as deliberately unused to keep the compiler quiet.
    (void) mask;

    /* =======================================================================
     * ORGANISM 0 : Geobacter metallireducens GS-15 (iAF987)
     *
     * [FIX-3D] Note the explicit "microbeId == 0" test and the write to
     * bioR[microbeId] rather than bioR[0]. In the 2D version index 0 meant
     * "whichever microbe happened to be first in this voxel", so with more than
     * one surrogate organism the network was applied to an arbitrary one and
     * all the others were left at zero growth (= pure decay). Here every call
     * knows exactly which organism it is for.
     *
     * TO ADD A SECOND ORGANISM, copy this whole block and change the test and
     * the weights, for example:
     *
     *     else if (microbeId == 1) {          // e.g. Rhodoferax ferrireducens
     *         static const std::vector<double> b1_Rf = { ... its own weights ... };
     *         ...
     *         bioR[microbeId] = g_Rf;
     *     }
     *
     * Microbes that are NOT handled by any branch fall through with bioR left at
     * 0, i.e. they decay. If that is not what you want, give them a branch (or
     * an "else" default) of their own.
     * ======================================================================= */
    if (microbeId == 0) {

        /* -------------------------------------------------------------------
         * THE TRAINED NETWORK: weights and biases, exactly as MATLAB exported
         * them. "IW1_1" = Input Weights into layer 1; "LWk_j" = Layer Weights
         * from layer j into layer k; "bk" = the biases of layer k.
         *
         * These are declared "static const": they are built once, the first
         * time this routine runs, and then reused for every voxel and every
         * step (the routine is called millions of times, so rebuilding ~500
         * numbers on every call would be pure waste). "const" means nothing can
         * modify them afterwards, which also makes the sharing safe.
         * You may edit the numbers freely — remember to rebuild afterwards.
         * ------------------------------------------------------------------- */
    // IW1_1 : input weights, layer 1. 10 rows (one per neuron) x 2 columns (one per input).
    static const std::vector<std::vector<double>> IW1_1 = {
        {-2.4353481594548980205,1.2453128078005795132},
        {3.0447436282826410014,-0.2160224443494934421},
        {3.3652253034249217656,-2.8426260160350889095},
        {-3.1293074164349627964,-0.15144313699638009552},
        {-4.2466279892081377767,0.048067073437338070363},
        {-2.0644520776813846119,0.20579224152417988081},
        {-2.7118141917732496715,-1.1832463926401328713},
        {0.79348921033598163177,0.61294054741632619798},
        {-0.97016886502106802759,-0.44174306611501107378},
        {-2.6909608875897919056,0.41133381228525300877}
    };
    // b1 : biases of layer 1 (10 neurons).
    static const std::vector<double> b1 = {0.67792810869733710621,-0.9886167579067058897,-0.37003642622919635796,1.187986042338397219,1.4928311621218286476,0.74486892316637776101,1.2366922395968444892,0.58038069834857120011,0.69148356138717936847,0.76749740434725033378};
    // LW2_1 : weights from layer 1 to layer 2 (10 x 10).
    static const std::vector<std::vector<double>> LW2_1 = {
        {-0.13021009621502233067,0.46442254820360506784,-0.16254619072135081947,0.70479133981267150233,-0.80601268077157350866,0.35014543180417939672,0.56774008699689559876,-0.3906759623889550781,-0.898143664504920336,-0.24754261405271288377},
        {2.518750594262290754,-3.3890342415710121848,-0.71120619706686016848,3.8668651434643219744,4.4018552127422019282,2.2613443499182683816,3.9169814443830763828,-0.60161775818946705563,2.7068183968055965494,2.9795796512994212613},
        {0.42557023743950256334,-0.42910835299689670252,0.006623678487498242326,0.32585250139515931078,0.17252820512663838426,-0.33120755333000762022,0.15016415539439734173,0.13336646976018576294,0.23895690885079068355,-0.22793121065459601149},
        {0.20078429524850449628,-0.38157473464206032032,-0.099047076698002681217,0.070102567431937559683,0.35507260532355533478,-0.54461934713816939624,-0.1968711534564925314,0.60250507440798961589,-0.44426905569738106561,0.24315149459187768155},
        {0.31629747518418027674,-0.3118319044532789075,-0.9123417313631246861,-0.2573673873336732032,0.31267570421009849291,-0.13358228571689653719,-0.55898333022648172275,0.10666045246443203731,-0.66973344485640051715,0.24014261400320235929},
        {0.6924282332619320357,-0.078382125152097775755,-3.4840802350832711376,-0.18172819551184554721,0.1515799541514252502,-0.015265525306761567811,-0.32289229015509640641,0.16159579535671905748,-0.66863753416093285598,0.33647512057550205133},
        {-0.25886547528813286245,0.13633458336918899412,0.2947218735301452841,-0.0023282054626124577476,0.054456311793002872002,0.095008845918454917778,0.28959484682779945697,-0.029137712697683795793,0.49594911858347640043,0.039736674169332318607},
        {0.51073247336451210732,-0.14384984072356604701,-0.15639684424700089904,-0.11088610710943468118,-0.022872385465158150131,0.083135507534787433936,-0.21613475120211508851,0.16091041916938766954,-0.18378098674317042138,0.20554819668917792552},
        {0.25237798811669914789,-0.36105550555846621652,0.011381940728288644782,0.56407634030788489365,0.85647370945278011867,-0.14480711335240797899,0.3046788419351742494,0.0013516875590215597994,0.2757582991239274639,0.38411940065103516995},
        {0.33884964013140711492,-0.23090126343401515263,-0.96122609731222552476,-0.33304000614680201453,-0.14417899498289740712,0.33737575440946965255,-0.4533261036858540205,0.053629664722837928903,-0.19685127443443936612,0.093581538729661503662}
    };
    // b2 : biases of layer 2.
    static const std::vector<double> b2 = {-1.1913964702126025319,0.5315756456217934911,0.51812888918392441262,0.78473922822423625156,-0.04780782970920898628,-0.074872203343969892519,0.14691264129258957416,0.28025089037945333237,-0.25030251589138291513,-0.35615209218765614407};
    // LW3_2 : weights from layer 2 to layer 3 (10 x 10).
    static const std::vector<std::vector<double>> LW3_2 = {
        {-0.82091872069622850994,0.22438251557365621047,0.5265595653322431291,-0.69973974122533599829,0.37042549249387185517,-0.19308018099405999113,-0.10580395159475251832,-0.35336236161202799755,0.1916518290153929327,0.35086123362869919839},
        {0.31091264379149929908,-0.41007743706175431297,-0.48046580580431064167,0.58664915847258081172,0.22131366753307862849,0.15710960104902832457,0.12075542322889883107,-0.028059878601090128963,0.18240355990986434342,0.12440912996266705048},
        {0.96026291928512541585,-0.15062914356935228066,-0.30695289709264289568,0.26381162551885606327,0.38617277157996904302,-0.38170038724703242439,-0.0058172276994241274573,0.012571174361457283439,-0.38188721451308255128,0.28983069504579134223},
        {-0.43721611603198734519,-0.38226614450786755572,-0.064579274030037991938,0.68636814976872972949,0.29061976100115061161,-0.19957639893847745061,-0.11680693987323419181,0.048298821524873372657,-0.27169183799227214493,0.11861785588166218197},
        {-0.35502458178955448309,5.5157577586010413384,0.61302041852483513118,0.42279811612204359905,-0.10291230816136444359,-0.22076224661081123024,0.31129365685323112656,-0.061542142486668921508,0.63168159727335326803,-0.71462674446375573645},
        {-0.72917896299636197899,6.1331129692495380823,0.94029284802146151367,-0.62536401033617916578,-0.22084173697346587417,-0.19512019308629499625,0.11926434833272134273,0.031419645512173362267,0.57682446333386916404,-0.98336946536787828155},
        {0.098646442334866371593,-0.37893382044506546125,-0.34244456237860165793,0.032650961554379520635,0.41080882671301532927,-0.37502061540207848322,-0.42317490376427258081,0.072047219209698989961,0.18093273528051156962,0.043773091105185450711},
        {-0.4841899677432765503,2.1348341673979480682,0.43061476346751970112,-0.90133856813014867626,-0.64908551956644688907,-0.63958696082321175869,0.7022518391949992278,-0.070022874034757501271,0.028640308837575328277,-0.45938347454213213084},
        {-0.64169883917626657777,1.4295258938535240212,0.58995300762194913258,-0.34490238270280276778,-0.057939343936153311909,-0.41764281575526895907,0.3855552452338121272,0.028627464540422968564,0.030846031053558486956,-0.8270502539059545466},
        {-0.073115059806730969827,0.52703172973593737094,-0.43719128364029113953,0.33410566336205554938,0.050351153505058941773,-0.00099089668102476154837,0.24476931026614226483,0.4556475612623314686,-0.24209632021610011376,-0.25351670878630366834}
    };
    // b3 : biases of layer 3.
    static const std::vector<double> b3 = {0.83095050244120038929,-0.3825593621459533189,-0.30574795305885948959,-0.19276192855578366814,-0.44085952624946872502,0.033006063256286348462,0.041652187406828918015,0.7305377886218068495,0.020819257857229409719,0.4305093214655514311};
    // LW4_3 : weights from layer 3 to layer 4 (10 x 10).
    static const std::vector<std::vector<double>> LW4_3 = {
        {0.059376900871753603151,0.065211221136606573046,0.1271323327647514434,0.052177527427534793614,-0.091507903248997726764,-0.083387337494353508394,-0.12177727144466253539,0.09500991491455457183,-0.19593155730913847101,-0.088671200217663559418},
        {-0.089769771890704092021,0.76643193439072965223,0.60761302769193525908,0.12116809552328869359,-0.75750278302614015846,-0.6418405634060270204,0.098320409066735614534,-0.33103336379812819956,0.026766804490196069444,0.11232913636935391855},
        {0.037855696090291823808,0.014062954034362624631,-0.17573398872942558313,0.13200451401746191027,-0.078003996544937156954,-0.036978318665778164842,-0.060842383502702331033,0.027339729725764701923,0.039318031935295143231,-0.15346890150974940026},
        {0.57367122054253938401,-0.37859251209054933796,0.14403494026234808789,-0.61829591949860285283,1.4075355971549599055,1.4353399340805503837,-0.51715996925715501664,1.1563131883992683324,-0.12644624692200742699,-0.44830438379016496198},
        {-0.092141549145137904842,0.92593116399977504205,0.16547819777257360974,0.45047699556308307134,-0.25849350700304479789,-2.6775146159887666109,0.013845636158351385531,-0.45942665476506833189,-0.50885720864770955796,0.17826599281715602152},
        {0.98560263600769482117,-0.6068022719561978473,-0.14405888941267674941,-0.37402020425467558118,4.8060884258287162041,4.6669221352427765481,-0.31396656433810349318,2.6511657592973332243,1.7851072292153427057,0.26026441566828961705},
        {0.037965905954689092849,0.078946346238600167977,-0.051274048622614552817,-0.16525330221492048888,0.2445666393614064904,0.91982039311092222977,0.028837994974132648285,0.24866123979819970691,0.22687031011158714788,0.48090175189056399985},
        {0.25266185203144747584,-0.22712153817671429379,-0.15841465368015558712,-0.066760356633889891831,0.089591994516210693433,-0.10791148810566909833,-0.055028367611793319036,-0.15914571815020828183,0.081741948688057744499,-0.10092094147119363978},
        {1.3040893149970238518,-0.55801226086847832697,-0.98044396835311597993,0.23425980957510209035,0.18021286832592975369,0.2874211884096304348,0.40338674326675910686,0.33142284410054684285,-0.089475538864739537215,-0.32480745373934855058},
        {0.0096018183614487231936,0.062714743404617773193,0.027392093261605989646,0.049642895566454105227,0.037226016709552951778,-0.15834473600114601366,-0.18006539194424373007,0.041681844469807542708,-0.065903487355015166749,-0.096443697730371300003}
    };
    // b4 : biases of layer 4.
    static const std::vector<double> b4 = {-0.026055378788094950976,-0.2945281106715050834,0.085090100040697905226,0.29779965595809959611,0.19527693799985568202,0.12917058855744362189,0.48370574551530176599,0.098492259186305872176,0.47705363420121221774,-0.19315483138355993287};
    // LW5_4 : weights from layer 4 to the single output neuron (1 x 10).
    static const std::vector<double> LW5_4 = {0.082646451299537612711,0.37850905959073699591,0.059437553346362914652,-1.0173724340607672723,-1.1624070560978583266,-0.84840045855708201561,0.50599701115930084683,0.044083723718335611486,-1.3130720140869056589,-0.10527344255602824608};
    // b5 : bias of the output neuron (a single number, not a list).
    static const double b5 = 0.85641685816650836571;
        /* -------------------------------------------------------------------
         * INPUT SCALING ("mapminmax" in MATLAB)
         *   MATLAB trains on rescaled inputs: each input is mapped linearly onto
         *   the range -1 .. +1 using the minimum and maximum seen in training.
         *       xp = (x - offset) * gain + ymin
         *   The gains below therefore also tell you the TRAINING RANGE of the
         *   shipped network, which is the range in which you may trust it:
         *       input 0 : offset .. offset + 2/gain  =  8.9016e-04 .. ~10.0
         *       input 1 : offset .. offset + 2/gain  =  3.5130e-05 .. ~0.4998
         *   (both in mmol/gDW/h). Outside that range the network EXTRAPOLATES
         *   and its answer is not meaningful — see the header note.
         * ------------------------------------------------------------------- */
        static const double x_offset0 = 0.000890159834356918;
        static const double x_offset1 = 3.51303254819135e-05;
        static const double x_gain0   = 0.200059307889652;
        static const double x_gain1   = 4.00231590228802;
        static const double x_ymin    = -1;

        /* OUTPUT SCALING: the reverse map, turning the network's internal
         * -1 .. +1 answer back into a growth rate in 1/h. */
        static const double y_ymin   = -1;
        static const double y_gain   = 37.9181676275123;
        static const double y_offset = 0;

        // the two inputs for THIS microbe, rescaled ready for the network
        if (Fin[iM].size() < 2) {                    // needs at least two substrates
            surrogateDimWarning("defineSurrogateModel (needs 2 substrates)", Fin[iM].size(), 2);
            return;
        }
        std::vector<double> Xp0(2);
        Xp0[0] = (Fin[iM][0] - x_offset0) * x_gain0 + x_ymin;   // substrate 0 (electron donor)
        Xp0[1] = (Fin[iM][1] - x_offset1) * x_gain1 + x_ymin;   // substrate 1 (electron acceptor)

        /* -------------------------------------------------------------------
         * THE LAYERS. Each line is one layer: multiply by that layer's weight
         * matrix, add its biases, then squash through tansig. The last layer is
         * LINEAR (no tansig), which is what makes the net able to output an
         * unbounded number.
         * ------------------------------------------------------------------- */
        const std::vector<double> a1 = tansig( matSum(b1, matMultp(IW1_1, Xp0)) );
        const std::vector<double> a2 = tansig( matSum(b2, matMultp(LW2_1, a1 )) );
        const std::vector<double> a3 = tansig( matSum(b3, matMultp(LW3_2, a2 )) );
        const std::vector<double> a4 = tansig( matSum(b4, matMultp(LW4_3, a3 )) );
        const double a5 = b5 + vecMultp(LW5_4, a4);             // linear output neuron

        /* -------------------------------------------------------------------
         * BACK TO PHYSICAL UNITS, and the final answer.
         * ------------------------------------------------------------------- */
        double g = (a5 - y_ymin) / y_gain + y_offset;           // specific growth rate, 1/h
        if (g < 1e-8) { g = 0.0; }                              // numerical dust -> exactly zero
        if (std::isnan(g) || std::isinf(g)) { g = 0.0; }        // safety net (see defineKinetics.hh)

        bioR[microbeId] = g;                                    // [FIX-3D] was bioR[0]

        /* ###################################################################
         * ### IMPORTANT: STOICHIOMETRIC CONSISTENCY IS *NOT* ENFORCED HERE ###
         * ###################################################################
         *
         * WHAT THE PROBLEM IS
         *   The solver gives us Fout PRE-FILLED with Fin, and whatever is in
         *   Fout when we return is what actually gets removed from the pore
         *   water. The shipped example (like the 2D original) never touches
         *   Fout. So the amount of substrate consumed stays at the raw
         *   Michaelis-Menten estimate Fin, while the biomass produced comes from
         *   the neural network — TWO INDEPENDENT NUMBERS. Nothing makes them
         *   agree.
         *
         *   Consequences you should be aware of:
         *     * carbon / electron balance is not closed: the model can consume
         *       substrate and produce a growth rate that no yield coefficient
         *       could reconcile, in either direction;
         *     * if the net predicts zero growth (for instance because the inputs
         *       fell outside its training range) the substrate is STILL consumed
         *       at the full Michaelis-Menten rate — food disappears for nothing;
         *     * conversely, capping the uptake elsewhere does not cap growth.
         *
         *   For the shipped single-organism demonstration this is tolerable,
         *   because the network was trained on FBA solutions driven by exactly
         *   these Michaelis-Menten uptake rates, so the two are consistent BY
         *   CONSTRUCTION as long as you stay inside the training range. As soon
         *   as you deviate from that setup, you should close the balance
         *   yourself.
         *
         * HOW TO CLOSE IT — worked example, tie consumption to predicted growth
         * through a yield coefficient Y (gDW of cells per mmol of substrate; it
         * is the inverse of the usual FBA slope of growth versus uptake, and you
         * read it straight off the training data):
         *
         *     // Yield of THIS organism on each substrate, gDW per mmol.
         *     // Y[0] = donor (e.g. acetate), Y[1] = acceptor.
         *     static const std::vector<double> Y = { 0.0435, 0.0091 };
         *
         *     for (std::size_t iS = 0; iS < Fin[iM].size() && iS < Y.size(); ++iS) {
         *         if (Y[iS] <= 0.0) continue;              // substrate not consumed for growth
         *         double needed = g / Y[iS];               // mmol/gDW/h REQUIRED by the growth
         *                                                  //   the network just predicted
         *         if (needed > Fin[iM][iS]) {              // cannot take more than is available:
         *             needed = Fin[iM][iS];                //   the transport estimate is the cap
         *         }
         *         Fout[iM][iS] = needed;                   // consume exactly what the growth needs
         *     }
         *
         * A stricter variant also feeds the cap back into the growth rate, so
         * that a substrate shortage really does slow the organism down:
         *
         *     double gLimited = g;
         *     for (std::size_t iS = 0; iS < Fin[iM].size() && iS < Y.size(); ++iS) {
         *         if (Y[iS] <= 0.0) continue;
         *         double gMax = Fin[iM][iS] * Y[iS];       // most growth this substrate allows
         *         if (gMax < gLimited) { gLimited = gMax; } // Liebig's law of the minimum
         *     }
         *     bioR[microbeId] = gLimited;
         *     for (std::size_t iS = 0; iS < Fin[iM].size() && iS < Y.size(); ++iS) {
         *         if (Y[iS] > 0.0) { Fout[iM][iS] = gLimited / Y[iS]; }
         *     }
         *
         * Products (substrates the organism RELEASES rather than consumes) use
         * the same idea with a negative entry: Fout[iM][iS] = -p * gLimited,
         * with p in mmol released per gDW of new cells. Remember the sign
         * convention: POSITIVE Fout = consumed, NEGATIVE Fout = produced.
         * ################################################################### */
    }
    // else if (microbeId == 1) { ... a second organism's network goes here ... }
}

#endif // SURROGATE_MODEL_HH   // closes the include guard opened at the top
