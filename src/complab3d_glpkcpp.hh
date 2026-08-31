/* This file is a part of the CompLaB program.
 *
 * The CompLaB3D software (3-D pore-scale extension) is developed since 2024
 * by the University of Georgia (United States, Meile Lab, Department of
 * Marine Sciences). The original 2-D CompLaB v1.0 was a collaboration of
 * the University of Georgia and Chungnam National University (South Korea).
 *
 * Contact:
 * Shahram Asgari, Christof Meile
 * Department of Marine Sciences (Meile Lab)
 * University of Georgia, Athens, GA 30602, USA
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

/* ==========================================================================
 *  PROVENANCE / LICENCE NOTICE  --  READ BEFORE REDISTRIBUTING
 * ==========================================================================
 *
 *  This file is a derivative work of `glpkcc.cpp`, the C++ MEX/oct-file
 *  interface to the GNU Linear Programming Kit known as GLPKMEX / the
 *  Octave-Forge `glpk` interface:
 *
 *      Original author : Nicolo' Giorgetti  (GLPKMEX, glpkcc.cpp)
 *      Later authors   : Niels Klitgord (COBRA Toolbox maintenance),
 *                        Jonathan Currie (Windows/OPTI maintenance)
 *      The same file is redistributed with the COBRA Toolbox.
 *
 *  The structure kept here from that upstream is essentially verbatim in
 *  spirit: the `glpIntParam[NIntP]` / `glpRealParam[NRealP]` parameter block,
 *  the `IParam[]` / `RParam[]` name tables, the `ctype`/`vartype` character
 *  encodings, the `freeLB`/`freeUB` handling, the `save_pb` export switch,
 *  the `+100` offset on the simplex/MIP return codes, and the ordering of the
 *  control-parameter remapping.  The 2-D CompLaB file `complab_glpkcpp.hh`
 *  is the intermediate adaptation; this file is the 3-D port of that.
 *
 *  *** LICENCE COMPATIBILITY ***
 *  GLPK itself and glpkcc.cpp/GLPKMEX are distributed under the GNU General
 *  Public License (GPL).  CompLaB is distributed under the GNU Affero General
 *  Public License (AGPL).  Combining GPL-only code with AGPL code, and
 *  redistributing the result, is NOT automatically permitted: AGPLv3 and
 *  GPLv3 are individually compatible with each other only through the
 *  explicit provisions of GPLv3 section 13 / AGPLv3 section 13, and the
 *  upstream file does not carry an "or (at your option) any later version"
 *  grant in every distributed copy.  This is a real, unresolved question and
 *  it is left to the CompLaB maintainer to resolve before any binary or
 *  source release that includes this file.  Options are: (a) obtain
 *  relicensing permission from the upstream authors, (b) rewrite the
 *  parameter-marshalling from the GLPK reference manual without reference to
 *  glpkcc.cpp, or (c) ship this module as a separately-distributed optional
 *  plug-in.  Until then this file is compiled only when the user explicitly
 *  defines COMPLAB_ENABLE_GLPK.
 *
 *  Build usage:
 *      add_definitions(-DCOMPLAB_ENABLE_GLPK)
 *      target_link_libraries(${EXECUTABLE_NAME} glpk)
 *  and, from a translation unit that has already included palabos3D.h /
 *  palabos3D.hh, done `using namespace plb;` and `typedef double T;`:
 *      #include "complab3d_glpkcpp.hh"
 *  With COMPLAB_ENABLE_GLPK undefined the whole file is an empty no-op, so it
 *  can sit in the tree unconditionally and CompLB3D still builds without GLPK
 *  installed.
 * ========================================================================== */

#ifndef COMPLAB3D_GLPKCPP_HH
#define COMPLAB3D_GLPKCPP_HH

#ifdef COMPLAB_ENABLE_GLPK

#include <glpk.h>

#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <cfloat>   // DBL_MAX
#include <climits>  // INT_MAX
#include <ctime>    // clock(), CLOCKS_PER_SEC
#include <cmath>    // std::isinf

// NOTE: palabos3D.h / palabos3D.hh, `using namespace plb;` and
// `typedef double T;` are assumed to be in scope already (they come from
// complab_functions.hh in the CompLB3D tree).  `pcout` below is plb::pcout.

// -------------------------------------------------------------------------
// Size of the GLPKMEX parameter block (kept identical to 2-D / upstream).
// -------------------------------------------------------------------------
#define NIntP  21
#define NRealP 11

// Control-parameter identifiers that are not declared in glpk.h (upstream
// GLPKMEX invents these two).  Only consumed by the legacy lpx_* path below.
#ifndef LPX_K_PREPROCESS
#define LPX_K_PREPROCESS 401   /* preprocessing */
#endif
#ifndef LPX_K_RATIO_TEST
#define LPX_K_RATIO_TEST 402   /* ratio test */
#endif

// -------------------------------------------------------------------------
// [FIX-3D] Legacy-API detection.
// The whole lpx_* family (lpx_set_int_parm, lpx_set_real_parm,
// lpx_print_prob, the LPX_K_* / LPX_E_* identifiers) and the 128-bit
// `glp_long` type used by the 4-argument glp_mem_usage() were removed from
// GLPK after ~4.45 and do not exist at all in GLPK 5.x.  The 2-D file uses
// them unconditionally, so it simply does not compile against any GLPK a user
// is likely to have installed today.  Everything that depends on them is
// gated on this macro.
// -------------------------------------------------------------------------
#if defined(GLP_MAJOR_VERSION) && (GLP_MAJOR_VERSION == 4) && defined(GLP_MINOR_VERSION) && (GLP_MINOR_VERSION < 45)
#define COMPLAB3D_GLPK_LEGACY_API 1
#else
#define COMPLAB3D_GLPK_LEGACY_API 0
#endif

// GLPK release the CompLaB parameter defaults were originally tuned against.
#ifndef COMPLAB3D_GLPK_REFERENCE_VERSION
#define COMPLAB3D_GLPK_REFERENCE_VERSION "4.39"
#endif

// -------------------------------------------------------------------------
// WARM-START / PRESOLVE NOTE
// -------------------------------------------------------------------------
// CompLaB keeps one persistent glp_prob* per metabolic model and, per lattice
// voxel, only rewrites the bounds of the exchange columns before re-solving
// (see run_glpk).  The intent is a warm start: the simplex basis left over
// from the previous voxel should make the next solve a handful of dual
// simplex iterations instead of a full solve.
//
// That warm start does NOT actually happen with the CompLaB default.  The
// default parameter block sets glpIntParam[16] = 1, i.e. sParam.presolve =
// GLP_ON, and with the presolver enabled glp_simplex() builds its own
// presolved problem, solves that from a fresh basis, and post-solves; the
// basis carried in the problem object is discarded on every call.  (GLPK also
// only writes a basis back when the presolved problem solves to optimality.)
//
// Setting this knob to 0 (sParam.presolve = GLP_OFF) is the recommended
// tuning option for the voxel loop: the first solve is set up by
// glp_adv_basis() -- which set_glpk already calls when presolve is off -- and
// every later solve starts from the previous voxel's optimal basis.  Expect a
// large speedup on genome-scale models.  The trade-off is that with presolve
// off GLPK can return GLP_ENOPFS/GLP_ENODFS on a badly-scaled model instead
// of silently recovering, so it is left OFF-by-default-unchanged here and is
// exposed as an opt-in rather than flipped outright.
//
// It is routed through the existing glpIntParam[16] slot, so nothing about
// the GLPKMEX parameter-block idiom changes.
#ifndef COMPLAB3D_GLPK_PRESOLVE
#define COMPLAB3D_GLPK_PRESOLVE 1   // 1 = GLP_ON (2-D behaviour), 0 = GLP_OFF (warm start)
#endif

// -------------------------------------------------------------------------
// [FIX-3D] Error codes returned in place of exit(EXIT_FAILURE).
//
// The 2-D file calls exit(EXIT_FAILURE) in 15 places inside set_glpk /
// run_glpk / initialize_glpk.  run_glpk is called from inside a Palabos data
// processor, i.e. from inside the per-voxel loop of every MPI rank.  A data
// processor must never call exit(): the ranks do not agree on when a voxel is
// bad, so one rank leaving the job while the others are inside the next
// collective (a communication step, or the MPI_Finalize barrier) deadlocks
// the whole run -- the job then burns its full wall-clock allocation doing
// nothing instead of failing fast.  Every one of those exits is replaced by a
// pcout diagnostic plus a non-zero return, and the caller decides (CompLaB's
// existing convention is to zero the fluxes for that voxel and carry on).
//
// The codes are negative on purpose.  complab3d_processors decodes solver
// failures as `glpkerr - 100 == GLP_E*`, and the GLP_E* codes are the small
// positive integers 1..12; negative codes therefore fall through that decode
// chain harmlessly instead of being mis-reported as e.g. "singular matrix".
// -------------------------------------------------------------------------
static const int CLB3D_GLPK_ERR_WRITE     = -101; // problem export failed
static const int CLB3D_GLPK_ERR_SCALING   = -102; // unrecognised scaling option
static const int CLB3D_GLPK_ERR_PRIMDUAL  = -103; // unrecognised primal/dual method
static const int CLB3D_GLPK_ERR_MSGLEV    = -104; // bad message level
static const int CLB3D_GLPK_ERR_BRANCH    = -105; // bad branching rule
static const int CLB3D_GLPK_ERR_BTRACK    = -106; // bad backtracking rule
static const int CLB3D_GLPK_ERR_TOLINT    = -107; // tolint out of (0,1)
static const int CLB3D_GLPK_ERR_CUTSET    = -108; // bad cut-set selector
static const int CLB3D_GLPK_ERR_PPROCESS  = -109; // bad MIP pre-processing selector
static const int CLB3D_GLPK_ERR_LPSOLVER  = -110; // lpsolver not in {1,2,3}
static const int CLB3D_GLPK_ERR_CTYPE     = -111; // unrecognised row sense character
static const int CLB3D_GLPK_ERR_SAVEPB    = -112; // unrecognised/unsupported save_pb
static const int CLB3D_GLPK_ERR_DIM       = -113; // empty or ragged stoichiometric matrix
static const int CLB3D_GLPK_ERR_METHOD    = -114; // run_glpk called with an unknown method
static const int CLB3D_GLPK_ERR_SETUP     = -115; // set_glpk reported a failure


// =========================================================================
// set_glpk
// -------------------------------------------------------------------------
// Populates an *empty* glp_prob with a complete LP/MILP: adds `ncol`
// structural columns (bounds from freeLB/lb/freeUB/ub, objective row from c,
// kind from vartype when isMIP), adds `nrow` auxiliary rows (senses from
// ctype, right-hand sides from b), loads the `nz` constraint-matrix triplets
// (rn/cn/a, 1-based, index 0 unused), optionally scales the problem, and
// optionally builds an advanced initial basis.
//
// MUTATES: *lp (columns, rows, matrix, scaling, basis), sParam and iParam
// (the simplex and integer-optimiser control blocks are filled in from the
// glpIntParam/glpRealParam block), `method` (set to 'S', 'T', 'E' or 'I') and
// `isMIP` (forced consistent with `method`).  Reads everything else.
//
// Called once per metabolic model at start-up.  Returns 0 on success and a
// negative CLB3D_GLPK_ERR_* code on a bad parameter -- it never exits.
// =========================================================================
inline int set_glpk (glp_prob *lp, int sense, int nrow, int ncol, int nz, int &method,
            double *c, int *rn, int *cn, double *a, double *b, char *ctype,
            bool *freeLB, double *lb, bool *freeUB, double *ub, int *vartype, int &isMIP, int lpsolver, glp_smcp &sParam, glp_iocp &iParam,
            int save_pb, char *save_filename, int *glpIntParam, int *IParam, double *glpRealParam, int *RParam)
{
  // Redirect standard output.  NOTE (preserved from 2-D/upstream): passing a
  // NULL hook *restores* GLPK's default terminal output, it does not silence
  // it.  The actual silencing comes from glpIntParam[0] == 0, which maps to
  // msg_lev = GLP_MSG_OFF below.  Real silencing would be
  // glp_term_out(GLP_OFF); left as-is so behaviour matches 2-D exactly.
  glp_term_hook (NULL, NULL); // default

  //-- Set the sense of optimization (the optimization direction flag)
  if (sense == 1) { // sense is the first argument of the glpk function
    glp_set_obj_dir (lp, GLP_MIN); // minimization
  }
  else {
    glp_set_obj_dir (lp, GLP_MAX); // maximization
  }

  //-- Define the number of unknowns and their domains.

  // column (structural) variables
  glp_add_cols (lp, ncol);
  for (int iT = 0; iT < ncol; ++iT) {
    //-- Define type of the structural variables
    if (!freeLB[iT] && !freeUB[iT]) {
      if ( lb[iT] == ub[iT] ) {
        glp_set_col_bnds (lp, iT+1, GLP_FX, lb[iT], ub[iT]); // GLP_FX: Fixed variable
      }
      else {
        glp_set_col_bnds (lp, iT+1, GLP_DB, lb[iT], ub[iT]); // GLP_DB: Double-bounded variable
      }
    }
    else {
      if (!freeLB[iT] && freeUB[iT]) {
        glp_set_col_bnds (lp, iT+1, GLP_LO, lb[iT], ub[iT]); // Variable with lower bound
      }
      else {
        if (freeLB[iT] && !freeUB[iT]) {
          glp_set_col_bnds (lp, iT+1, GLP_UP, lb[iT], ub[iT]); // Variable with upper bound
        }
        else {
          glp_set_col_bnds (lp, iT+1, GLP_FR, lb[iT], ub[iT]); // Free (unbounded) variable
        }
      }
    }

    // -- Set the objective coefficient of the corresponding
    // -- structural variable. No constant term is assumed.
    glp_set_obj_coef(lp, iT+1, c[iT]);

    if (isMIP) {
      glp_set_col_kind (lp, iT+1, vartype[iT]);
    }
  }

  // row (auxiliary) variables
  int typx = 0;
  glp_add_rows (lp, nrow);
  for (int iT = 0; iT < nrow; ++iT)
  {
  /*  If the i-th row has no lower bound (types F,U),
      the corrispondent parameter will be ignored.
      If the i-th row has no upper bound (types F,L),
      the corrispondent parameter will be ignored.
      If the i-th row is of S type,
      the i-th LB is used, but the i-th UB is ignored. */
    switch (ctype[iT])
    {
      // Free bound
      case 'F': typx = GLP_FR; break;
      // upper bound
      case 'U': typx = GLP_UP; break;
      // lower bound
      case 'L': typx = GLP_LO; break;
      // fixed constraint
      case 'S': typx = GLP_FX; break;
      // double-bounded variable
      case 'D': typx = GLP_DB; break;
      // [FIX-3D] 2-D left typx at its previous value (0 on the first row) for
      // an unrecognised sense character and handed that to GLPK, which reacts
      // by calling its own error handler -> abort().  Fail cleanly instead.
      default:
        pcout << "complab3d_glpkcpp: unrecognized row sense '" << ctype[iT]
              << "' at row " << iT << " (expected one of F,U,L,S,D)." << std::endl;
        return CLB3D_GLPK_ERR_CTYPE;
    }

    // NOTE (preserved from 2-D): b[iT] is passed as BOTH the lower and the
    // upper bound.  For CompLaB this is always right because every row is of
    // type 'S' (S*v = 0), where GLPK ignores the upper argument anyway.  It
    // is wrong for the 'D' case, where GLPKMEX semantics are [-b, b].  Not
    // changed, because changing it would silently alter results for anyone
    // who does use 'D'.
    glp_set_row_bnds (lp, iT+1, typx, b[iT], b[iT]);
  }

  // Load the vectorized constraint matrix a (1-based triplets, [0] unused)
  glp_load_matrix (lp, nz, rn, cn, a);

  // Save problem -- OFF by default (save_pb == 0 at every CompLaB call site).
  if ( save_pb == 1 ){
    if (glp_write_lp (lp, NULL, save_filename) != 0) {
      pcout << "complab3d_glpkcpp: unable to write the problem." << std::endl;
      return CLB3D_GLPK_ERR_WRITE;   // [FIX-3D] was exit(EXIT_FAILURE)
    }
  }
  else if ( save_pb == 2) {
    if (glp_write_mps (lp, GLP_MPS_DECK, NULL, save_filename) != 0) {
      pcout << "complab3d_glpkcpp: unable to write the problem." << std::endl;
      return CLB3D_GLPK_ERR_WRITE;   // [FIX-3D] was exit(EXIT_FAILURE)
    }
  }
  else if ( save_pb == 3) {
    if (glp_write_mps (lp, GLP_MPS_FILE, NULL, save_filename) != 0) {
      pcout << "complab3d_glpkcpp: unable to write the problem." << std::endl;
      return CLB3D_GLPK_ERR_WRITE;   // [FIX-3D] was exit(EXIT_FAILURE)
    }
  }
  else if ( save_pb == 4 ) {
#if COMPLAB3D_GLPK_LEGACY_API
    if (lpx_print_prob (lp, save_filename) != 0) {
      pcout << "complab3d_glpkcpp: unable to write the problem." << std::endl;
      return CLB3D_GLPK_ERR_WRITE;   // [FIX-3D] was exit(EXIT_FAILURE)
    }
#else
    // [FIX-3D] 2-D line ~154: lpx_print_prob() was removed from GLPK after
    // ~4.45; there is no plain-text problem dump in modern GLPK.  Use
    // save_pb = 1 (CPLEX LP) or 2/3 (MPS) instead.
    pcout << "complab3d_glpkcpp: save_pb = 4 (plain text) needs lpx_print_prob(), "
          << "which was removed from GLPK after 4.45. Use save_pb = 1/2/3." << std::endl;
    return CLB3D_GLPK_ERR_SAVEPB;
#endif
  }

  //-- scale the problem data (if required) (see manual page 33)
  // In GLPK the scaling means a linear transformation applied to the
  // constraint matrix to improve its numerical properties.
  if (glpIntParam[1] && (! glpIntParam[16] || lpsolver != 1)) {
    switch ( glpIntParam[1] ) {
      case ( 1 ): glp_scale_prob( lp, GLP_SF_GM ); break; // geometric mean scaling
      case ( 2 ): glp_scale_prob( lp, GLP_SF_EQ ); break; // equilibration scaling
      case ( 3 ): glp_scale_prob( lp, GLP_SF_GM | GLP_SF_EQ ); break;
      case ( 4 ): glp_scale_prob( lp, GLP_SF_2N ); break; // round scale factors to nearest power of two
      default :
        pcout << "complab3d_glpkcpp: unrecognized scaling option" << std::endl;
        return CLB3D_GLPK_ERR_SCALING;   // [FIX-3D] was exit(EXIT_FAILURE)
    }
  }
  else {
    /* do nothing? or unscale?
        glp_unscale_prob( lp );
    */
  }

  //-- build advanced initial LP basis (if required).  Only meaningful when
  //   the presolver is off -- see the WARM-START note at the top of the file.
  if (lpsolver == 1 && ! glpIntParam[16])
    glp_adv_basis (lp, 0);

  //-- set control parameters for simplex (lpsolver == 1) / exact (lpsolver == 3) method (see manual page 41)
  if (lpsolver == 1 || lpsolver == 3) {
    //remap of control parameters for simplex method
    sParam.msg_lev=glpIntParam[0];  // message level (0 == GLP_MSG_OFF)

    // simplex method: primal/dual
    switch ( glpIntParam[2] ) {
      case 0: sParam.meth=GLP_PRIMAL; break;
      case 1: sParam.meth=GLP_DUAL;   break;
      case 2: sParam.meth=GLP_DUALP;  break;
      default:
        pcout << "complab3d_glpkcpp: unrecognized primal/dual method." << std::endl;
        return CLB3D_GLPK_ERR_PRIMDUAL;  // [FIX-3D] was exit(EXIT_FAILURE)
    }

    // pricing technique
    if (glpIntParam[3]==0) sParam.pricing=GLP_PT_STD;
    else sParam.pricing=GLP_PT_PSE;

    // ratio test
    if (glpIntParam[20]==0) sParam.r_test = GLP_RT_STD;
    else sParam.r_test=GLP_RT_HAR;

    //tollerances
    sParam.tol_bnd=glpRealParam[1]; // primal feasible tollerance
    sParam.tol_dj=glpRealParam[2];  // dual feasible tollerance
    sParam.tol_piv=glpRealParam[3]; // pivot tollerance
    sParam.obj_ll=glpRealParam[4];  // lower limit
    sParam.obj_ul=glpRealParam[5];  // upper limit

    // iteration limit
    if (glpIntParam[5]==-1) sParam.it_lim=INT_MAX;
    else sParam.it_lim=glpIntParam[5];

    // time limit
    if (glpRealParam[6]==-1) sParam.tm_lim=INT_MAX;
    else sParam.tm_lim=(int) glpRealParam[6];
    sParam.out_frq=glpIntParam[7];  // output frequency
    sParam.out_dly=(int) glpRealParam[7]; // output delay
    // presolver -- see the WARM-START note at the top of this file
    if (glpIntParam[16]) sParam.presolve=GLP_ON;
    else sParam.presolve=GLP_OFF;
  }
  else {
#if COMPLAB3D_GLPK_LEGACY_API
    // [FIX-3D] 2-D lines ~227-238: the lpx_* control-parameter setters are
    // gone in GLPK >= ~4.45 and absent entirely in 5.x, so this loop broke
    // the build on every modern GLPK.  Kept for old installs only.
    for (int iT = 0; iT < NIntP; ++iT) {
      // skip assinging ratio test
      if ( iT == 18 || iT == 20) {
        continue;
      }
      lpx_set_int_parm (lp, IParam[iT], glpIntParam[iT]);
    }

    for (int iT = 0; iT < NRealP; ++iT) {
      lpx_set_real_parm (lp, RParam[iT], glpRealParam[iT]);
    }
#else
    // Modern GLPK: this branch is only reached for lpsolver == 2 (interior
    // point); run_glpk builds its own glp_iptcp so it can silence the solver log,
    // defaults.  There is nothing to set, so the IParam/RParam tables are
    // simply unused here.  They are still carried through the signature so
    // the GLPKMEX parameter-block idiom -- and every existing call site --
    // stays unchanged.
    (void) IParam; (void) RParam; (void) glpRealParam;
#endif
  }

  // set MIP params if MIP....
  if (isMIP) {
    method = 'I';

    switch (glpIntParam[0]) { //message level
      case 0:  iParam.msg_lev = GLP_MSG_OFF;   break;
      case 1:  iParam.msg_lev = GLP_MSG_ERR;   break;
      case 2:  iParam.msg_lev = GLP_MSG_ON;    break;
      case 3:  iParam.msg_lev = GLP_MSG_ALL;   break;
      default: pcout << "complab3d_glpkcpp: msg_lev bad param" << std::endl;
               return CLB3D_GLPK_ERR_MSGLEV;   // [FIX-3D] was exit(EXIT_FAILURE)
    }
    switch (glpIntParam[14]) { //branching param
      case 0:  iParam.br_tech = GLP_BR_FFV;    break;
      case 1:  iParam.br_tech = GLP_BR_LFV;    break;
      case 2:  iParam.br_tech = GLP_BR_MFV;    break;
      case 3:  iParam.br_tech = GLP_BR_DTH;    break;
      default: pcout << "complab3d_glpkcpp: branch bad param" << std::endl;
               return CLB3D_GLPK_ERR_BRANCH;   // [FIX-3D] was exit(EXIT_FAILURE)
    }
    switch (glpIntParam[15]) { //backtracking heuristic
      case 0:  iParam.bt_tech = GLP_BT_DFS;    break;
      case 1:  iParam.bt_tech = GLP_BT_BFS;    break;
      case 2:  iParam.bt_tech = GLP_BT_BLB;    break;
      case 3:  iParam.bt_tech = GLP_BT_BPH;    break;
      default: pcout << "complab3d_glpkcpp: backtrack bad param" << std::endl;
               return CLB3D_GLPK_ERR_BTRACK;   // [FIX-3D] was exit(EXIT_FAILURE)
    }

    if ( glpRealParam[8] > 0.0 && glpRealParam[8] < 1.0 ) {
      iParam.tol_int = glpRealParam[8];  // absolute tolorence
    }
    else {
      pcout << "complab3d_glpkcpp: tolint must be between 0 and 1" << std::endl;
      return CLB3D_GLPK_ERR_TOLINT;      // [FIX-3D] was exit(EXIT_FAILURE)
    }

    iParam.tol_obj = glpRealParam[9];  // relative tolarence
    iParam.mip_gap = glpRealParam[10]; // realative gap tolerance

    // set time limit for mip
    if ( glpRealParam[6] < 0.0 || glpRealParam[6] > 1e6 ) {
      iParam.tm_lim = INT_MAX;
    }
    else {
      iParam.tm_lim = (int)(1000.0 * glpRealParam[6] );
    }

    // Choose Cutsets for mip
    // shut all cuts off, then start over....
    iParam.gmi_cuts = GLP_OFF;
    iParam.mir_cuts = GLP_OFF;
    iParam.cov_cuts = GLP_OFF;
    iParam.clq_cuts = GLP_OFF;

    switch( glpIntParam[17] ) {
      case 0: break;
      case 1: iParam.gmi_cuts = GLP_ON; break;
      case 2: iParam.mir_cuts = GLP_ON; break;
      case 3: iParam.cov_cuts = GLP_ON; break;
      case 4: iParam.clq_cuts = GLP_ON; break;
      case 5: iParam.clq_cuts = GLP_ON;
              iParam.gmi_cuts = GLP_ON;
              iParam.mir_cuts = GLP_ON;
              iParam.cov_cuts = GLP_ON;
              iParam.clq_cuts = GLP_ON; break;
      default: pcout << "complab3d_glpkcpp: cutset bad param" << std::endl;
               return CLB3D_GLPK_ERR_CUTSET;   // [FIX-3D] was exit(EXIT_FAILURE)
    }

    switch( glpIntParam[18] ) { // pre-processing for mip
        case 0: iParam.pp_tech = GLP_PP_NONE; break;
        case 1: iParam.pp_tech = GLP_PP_ROOT; break;
        case 2: iParam.pp_tech = GLP_PP_ALL;  break;
        default: pcout << "complab3d_glpkcpp: pprocess bad param" << std::endl;
                 return CLB3D_GLPK_ERR_PPROCESS;  // [FIX-3D] was exit(EXIT_FAILURE)
    }

    if (glpIntParam[16])  iParam.presolve=GLP_ON;
    else                  iParam.presolve=GLP_OFF;

    if (glpIntParam[19])  iParam.binarize = GLP_ON;
    else                  iParam.binarize = GLP_OFF;

  }
  else {
     /* Choose simplex method ('S')
     or interior point method ('T')
     or Exact method          ('E')
     to solve the problem  */
    switch (lpsolver) {
      case 1: method = 'S'; break;
      case 2: method = 'T'; break;
      case 3: method = 'E'; break;
      default:
            pcout << "complab3d_glpkcpp: lpsolver must be 1 (simplex), 2 (interior) or 3 (exact); got "
                  << lpsolver << std::endl;
            return CLB3D_GLPK_ERR_LPSOLVER;   // [FIX-3D] was exit(EXIT_FAILURE)
    }
  }

  return 0;
}


// =========================================================================
// run_glpk
// -------------------------------------------------------------------------
// Solves ONE flux-balance problem for ONE lattice voxel.
//
// This does NOT rebuild the problem.  `lp` is the persistent glp_prob* that
// initialize_glpk built once for this metabolic model at start-up and that is
// shared by every voxel; all this function does is reset the bounds of the
// exchange columns listed in `location` (column index `location[iL]`, 0-based
// on the CompLaB side, gets GLP_DB with [lb[iL], ub[iL]]) and re-solve.  The
// stoichiometric matrix, the objective, the internal columns and every row
// are untouched from call to call.  Entries with location[iL] < 0 mean "this
// substrate has no exchange reaction in this model" and are skipped.
//
// MUTATES: *lp (exchange-column bounds, and whatever the solver leaves in the
// basis / interior-point solution), and the outputs xmin (primal flux vector,
// ncol entries), fmin (objective = growth rate), status (glp_*_status), lambda
// (nrow row duals), redcosts (ncol reduced costs), time (seconds), mem (kB).
// sParam/iParam are taken by value exactly as in 2-D, so a solve cannot
// perturb the shared control blocks.
//
// Returns 0 on success, the (offset) GLPK return code on a solver failure, or
// a negative CLB3D_GLPK_ERR_* code.  It never exits -- see the note on the
// error codes above; this function runs inside a Palabos data processor.
// =========================================================================
inline int run_glpk (glp_prob *lp, int method, int isMIP, int nrow, int ncol, int lpsolver, glp_smcp sParam, glp_iocp iParam,
             std::vector<double> &xmin, double &fmin, int &status, std::vector<int> location, std::vector<double> lb,
             std::vector<double> ub, std::vector<double> &lambda, std::vector<double> &redcosts, double &time, double &mem)

{
  clock_t t_start = clock();

  // [FIX-3D] CRITICAL.  This function writes xmin[0..ncol-1], lambda[0..nrow-1]
  //   and redcosts[0..ncol-1] with operator[], but never sized them.  The 2-D
  //   caller happened to size them at its call site (complab_processors.hh:144);
  //   the 3-D caller cleared them instead, so the very first solve wrote through
  //   a null data pointer and segfaulted.  Sizing belongs HERE, with the code
  //   that decides how much it writes -- a caller cannot be expected to know.
  //   Also guarantees a defined value on every early-return path below.
  xmin.assign((size_t) (ncol > 0 ? ncol : 0), 0.0);
  lambda.assign((size_t) (nrow > 0 ? nrow : 0), 0.0);
  redcosts.assign((size_t) (ncol > 0 ? ncol : 0), 0.0);
  fmin = 0.0; status = 0; time = 0.0; mem = 0.0;

  // Only the exchange columns are re-bounded; the problem itself persists.
  // [FIX-3D] GLP_DB with lb == ub makes glp_simplex return GLP_EBOUND.  That is
  //   reachable through <equate_bounds>, which sets hi = lo + 1e-12: for any
  //   |lo| above about 4.5e3 that sum is exactly lo in double precision.  Use
  //   GLP_FX when the two coincide, which is what GLPK wants for a fixed column.
  for (size_t iL = 0; iL < location.size(); ++iL) {
    if (location[iL] >= 0 ) {
      if (lb[iL] == ub[iL]) glp_set_col_bnds (lp, location[iL]+1, GLP_FX, lb[iL], ub[iL]);
      else                  glp_set_col_bnds (lp, location[iL]+1, GLP_DB, lb[iL], ub[iL]);
    }
  }

  // now run the problem...
  int errnum;
  switch (method) {
    case 'I': errnum = glp_intopt( lp, &iParam ); // solve MIP problem with the branch-and-cut method
              errnum += 100; //this is to avoid ambiguity in the return codes.
              break;

    case 'S': errnum = glp_simplex(lp, &sParam); // solve LP problem with the primal or dual simplex method
              errnum += 100; //this is to avoid ambiguity in the return codes.
              break;

    case 'T': {
                // [FIX-3D] glp_interior(lp, NULL) uses the glp_iptcp defaults, and the
                //   default msg_lev is GLP_MSG_ALL.  glpIntParam[0] (message level)
                //   only reaches the simplex and exact paths, so the interior-point
                //   solver printed a full Cholesky/iteration log for EVERY voxel of
                //   EVERY step on EVERY rank.  Silence it explicitly.
                glp_iptcp ipParam;
                glp_init_iptcp(&ipParam);
                ipParam.msg_lev = GLP_MSG_OFF;
                errnum = glp_interior(lp, &ipParam);
              } break; // solve LP problem with the interior-point method

    case 'E': errnum = glp_exact(lp, &sParam); break; // solve LP problem in exact arithmetic

    // [FIX-3D] 2-D line ~375: exit(EXIT_FAILURE) from inside a data processor.
    default:
        pcout << "complab3d_glpkcpp: run_glpk called with unknown method code " << method
              << " (expected 'S', 'T', 'E' or 'I')." << std::endl;
        status = CLB3D_GLPK_ERR_METHOD;
        return CLB3D_GLPK_ERR_METHOD;
  }

  /* Which return codes count as "solved".
     'S' (glp_simplex) and 'I' (glp_intopt) have been offset by +100, so
        100 == GLP_OK, 108 == GLP_EITLIM + 100, 109 == GLP_ETMLIM + 100
        (the last two are accepted because GLPK still leaves a usable
        incumbent solution in the problem object).
     'T' (glp_interior) and 'E' (glp_exact) are NOT offset, so success is 0.

     [FIX-3D] 2-D line ~382 tested
         errnum == LPX_E_OK || errnum == 100 || errnum == 109 || errnum == 108
     and LPX_E_OK is 200, not 0 -- so plain 0, the value glp_interior() and
     glp_exact() actually return on success, was never in the accepted set.
     Every successful lpsolver == 2 and lpsolver == 3 solve therefore fell
     into the failure branch, returned errnum == 0 (which the caller reads as
     "no error" but with fmin/xmin left completely unwritten), and the voxel
     silently got whatever was in the uninitialised output vectors.  In other
     words lpsolver 2 and 3 have never worked.  Accept 0 for 'T'/'E'.
     LPX_E_OK is dropped: no glp_* entry point can return 200, and the symbol
     does not exist in modern GLPK headers. */
  bool solved;
  if (method == 'S' || method == 'I') {
    solved = (errnum == 100 || errnum == 108 || errnum == 109);
  }
  else { // 'T' or 'E'
    solved = (errnum == 0);
  }

  if (solved) {
    // Get status and object value
    if (isMIP) {
      status = glp_mip_status (lp);
      fmin = glp_mip_obj_val (lp);
    }
    else {
      if (lpsolver == 1 || lpsolver == 3) {
        status = glp_get_status (lp);
        fmin = glp_get_obj_val (lp);
      }
      else {
        status = glp_ipt_status (lp);
        fmin = glp_ipt_obj_val (lp);
      }
    }

    // Get optimal solution (if exists)
    if (isMIP) {
      for (int iT = 0; iT < ncol; ++iT) {
        xmin[iT] = glp_mip_col_val (lp, iT+1);
      }
    }
    else {
      /* Primal values */
      for (int iT = 0; iT < ncol; ++iT) {
        if (lpsolver == 1 || lpsolver == 3) {
          xmin[iT] = glp_get_col_prim (lp, iT+1);
        }
        else {
          xmin[iT] = glp_ipt_col_prim (lp, iT+1);
        }
      }
      /* Dual values */
      for (int iT = 0; iT < nrow; ++iT) {
        if (lpsolver == 1 || lpsolver == 3) {
          lambda[iT] = glp_get_row_dual (lp, iT+1);
        }
        else {
          lambda[iT] = glp_ipt_row_dual (lp, iT+1);
        }
      }
      /* Reduced costs */
      for (int iT = 0; iT < ncol; ++iT) {
        if (lpsolver == 1 || lpsolver == 3) {
          redcosts[iT] = glp_get_col_dual (lp, iT+1);
        }
        else {
          redcosts[iT] = glp_ipt_col_dual (lp, iT+1);
        }
      }
    }

    // [FIX-3D] 2-D line ~435 divided clock_t by CLOCKS_PER_SEC in integer
    // arithmetic, so `time` was 0 for anything shorter than a second, i.e.
    // always.  Cast to double first.
    time = (double)(clock () - t_start) / (double) CLOCKS_PER_SEC;

#if COMPLAB3D_GLPK_LEGACY_API
    // [FIX-3D] 2-D lines ~437-439: `glp_long` and the glp_long form of
    // glp_mem_usage() were removed from GLPK after ~4.45 (and the type does
    // not exist at all in 5.x), so this block alone stopped the 2-D file
    // compiling against any current GLPK.  Kept for legacy installs.
    glp_long tpeak;
    glp_mem_usage(NULL, NULL, NULL, &tpeak);
    mem = (double)(4294967296.0 * tpeak.hi + tpeak.lo) / (1024);
#else
    // Modern GLPK still has glp_mem_usage(), but with a size_t 4th argument.
    // Peak allocator usage is a diagnostic CompLaB never consumes, and
    // probing it per voxel is pure overhead, so report 0 here.  A maintainer
    // who wants it back can use:
    //     size_t tpeak; glp_mem_usage(NULL, NULL, NULL, &tpeak);
    //     mem = (double) tpeak / 1024.0;
    mem = 0.0;
#endif

    return 0;
  }
  else {
    // Solver failure.
    // [FIX-3D] The 'T'/'E' paths used to return the raw GLP_E* code and copy it
    //   into `status`.  GLP_EFAIL is 5 and GLP_OPT is also 5, so a failed
    //   interior-point or exact solve was indistinguishable from an optimal one:
    //   the caller's `status == GLP_OPT` test passed and it read a growth rate
    //   and flux vector that the solver never wrote.  Offset these the same way
    //   'S'/'I' are offset, so no failure code can ever land in the 1..6 status
    //   range, and set `status` to a value that cannot be mistaken for success.
    if (method == 'T' || method == 'E') errnum += 100;
    status = -1;                 // never a valid glp_*_status() value
    return errnum;
  }
  status = errnum;

  return errnum;
}


// =========================================================================
// initialize_glpk
// -------------------------------------------------------------------------
// Builds the persistent LP for one metabolic model, once, at start-up.
// Converts the CompLaB std::vector representation of the model (dense
// stoichiometric matrix S, rhs vec_b, objective vec_c, bounds vec_lb/vec_ub,
// row senses ctype, column kinds vtype) into the flat GLPKMEX-style arrays,
// fills the default glpIntParam/glpRealParam block, and hands everything to
// set_glpk.
//
// MUTATES: *lp (which must be an empty glp_prob* from glp_create_prob(); it
// comes back fully populated and is then owned by the caller for the whole
// run), `method`, `isMIP`, `sParam` and `iParam` (initialised with
// glp_init_smcp/glp_init_iocp and then filled in by set_glpk).  Every input
// vector is read only.
//
// [FIX-3D] Returns int instead of 2-D's void: every failure path here used to
// be exit(EXIT_FAILURE).  Existing call sites compile unchanged (the return
// value is simply discarded), but callers should now check it -- 0 means the
// model is ready, negative means it is not and *lp must not be used.
// =========================================================================
inline int initialize_glpk(int iter, glp_prob *lp, std::vector< std::vector<double> > S, std::vector<double> vec_b, std::vector<double> vec_c,
                    std::vector<double> vec_lb, std::vector<double> vec_ub, char *ctype, char *vtype,
                    int sense, int lpsolver, int save_pb, int &method, int &isMIP, glp_smcp &sParam, glp_iocp &iParam) {

  std::string siter = std::to_string(iter+1);
  if (iter == 0) siter += "st ";
  else if (iter == 1) siter += "nd ";
  else if (iter == 2) siter += "rd ";
  else siter += "th ";

  // [FIX-3D] 2-D lines ~462-470: the version check was broken twice over.
  //   (a) `sizeof(glp_version())` is sizeof(const char*) == 8, not the length
  //       of the string, so the loop copied a fixed 8 bytes out of a 4- or
  //       5-character literal and read past its end.
  //   (b) std::string::compare() returns 0 when the strings are EQUAL, so
  //       `if (!glpk_ver.compare("4.39"))` printed "compatibility is not
  //       guaranteed" precisely when the version DID match, and stayed quiet
  //       for every version that did not.
  // Compare the whole string and warn only on a genuine mismatch, once per
  // run rather than once per model.
  {
    static bool versionWarned = false;
    const char *rawVer = glp_version();
    std::string glpk_ver = (rawVer != NULL) ? std::string(rawVer) : std::string("unknown");
    if (!versionWarned && glpk_ver != std::string(COMPLAB3D_GLPK_REFERENCE_VERSION)) {
      pcout << "complab3d_glpkcpp: running against GLPK " << glpk_ver
            << "; the CompLaB solver defaults were tuned for GLPK "
            << COMPLAB3D_GLPK_REFERENCE_VERSION
            << ". Compatibility is not guaranteed." << std::endl;
      versionWarned = true;
    }
  }

  // Integer/Real Param Defaults.  Slot 16 is the presolver -- see the
  // WARM-START note at the top of this file; override at compile time with
  // -DCOMPLAB3D_GLPK_PRESOLVE=0 to get a real warm start across voxels.
  int glpIntParam[NIntP] = { 0,1,0,1,0,INT_MAX,INT_MAX,200,1,2,0,1,0,0,3,2,
                             COMPLAB3D_GLPK_PRESOLVE, 0,2,0,1 };
  double glpRealParam[NRealP] = {0.07,1e-7,1e-7,1e-10,-DBL_MAX,DBL_MAX,INT_MAX,0.0,1e-5,1e-7,0.0};

  // Integer/Real Param Names (GLPKMEX parameter-block idiom, kept intact).
  // Only consumed by the legacy lpx_* path in set_glpk.
#if COMPLAB3D_GLPK_LEGACY_API
  int IParam[NIntP] = {LPX_K_MSGLEV,LPX_K_SCALE,LPX_K_DUAL,LPX_K_PRICE,LPX_K_ROUND,LPX_K_ITLIM,LPX_K_ITCNT,
    LPX_K_OUTFRQ,LPX_K_MPSINFO,LPX_K_MPSOBJ,LPX_K_MPSORIG,LPX_K_MPSWIDE,LPX_K_MPSFREE,LPX_K_MPSSKIP,
    LPX_K_BRANCH,LPX_K_BTRACK,LPX_K_PRESOL,LPX_K_USECUTS,LPX_K_PREPROCESS,LPX_K_BINARIZE,LPX_K_RATIO_TEST};
  int RParam[NRealP] = {LPX_K_RELAX,LPX_K_TOLBND,LPX_K_TOLDJ,LPX_K_TOLPIV,LPX_K_OBJLL,LPX_K_OBJUL,LPX_K_TMLIM,
    LPX_K_OUTDLY,LPX_K_TOLINT,LPX_K_TOLOBJ,LPX_K_MIPGAP};
#else
  // [FIX-3D] The LPX_K_* identifiers vanished with the lpx_* API.  The tables
  // are still declared (and still passed to set_glpk unchanged) so that the
  // signature and the parameter-block idiom are identical on both GLPK
  // generations; they are never read on this path.
  int IParam[NIntP];
  int RParam[NRealP];
  for (int iT = 0; iT < NIntP;  ++iT) IParam[iT] = -1;
  for (int iT = 0; iT < NRealP; ++iT) RParam[iT] = -1;
#endif

  // [FIX-3D] guard against an empty/ragged model instead of dereferencing
  // S[0] on an empty S.
  if (S.empty() || S[0].empty()) {
    pcout << "complab3d_glpkcpp: " << siter << "metabolic model has an empty stoichiometric matrix." << std::endl;
    return CLB3D_GLPK_ERR_DIM;
  }
  const int nrow = (int) S.size();    // number of substrates (metabolites)
  const int ncol = (int) S[0].size(); // number of reactions
  if ((int) vec_b.size() < nrow || (int) vec_c.size() < ncol ||
      (int) vec_lb.size() < ncol || (int) vec_ub.size() < ncol) {
    pcout << "complab3d_glpkcpp: " << siter << "metabolic model has inconsistent vector sizes "
          << "(nrow = " << nrow << ", ncol = " << ncol << ")." << std::endl;
    return CLB3D_GLPK_ERR_DIM;
  }

  // [FIX-3D] 2-D lines ~487-490 declared `double c[ncol], lb[ncol], ub[ncol];`
  // -- variable-length arrays, which are not standard C++ and put a
  // genome-scale model's worth of doubles on the stack -- and then zeroed
  // them with `memset(c, 0, ncol*sizeof(int))`.  sizeof(int) is 4 and
  // sizeof(double) is 8, so exactly the first half of each buffer was
  // cleared and the second half kept whatever was on the stack.  It only
  // escaped notice because the loop below overwrites all ncol entries anyway.
  // std::vector is heap-allocated, value-initialised, correctly sized and
  // self-freeing.
  std::vector<double> c(ncol, 0.0), lb(ncol, 0.0), ub(ncol, 0.0);

  for (int iT = 0; iT < ncol; ++iT) {
     c[iT] =  vec_c[iT];
    lb[iT] = vec_lb[iT];
    ub[iT] = vec_ub[iT];
  }

  // [FIX-3D] 2-D lines ~498-503 calloc'd rn/cn/a at the DENSE worst case,
  // nrow*ncol+1 entries each, and never freed them.  A genome-scale model
  // (say 1600 metabolites x 2600 reactions, ~5000 nonzeros) burns
  // 1600*2600*(4+4+8) = ~66 MB to hold ~80 kB of data, per model, for the
  // whole run.  Count the nonzeros in a first pass and allocate exactly
  // nz+1 (index 0 is unused: glp_load_matrix takes 1-based triplets).  The
  // vectors free themselves once glp_load_matrix has copied the data in.
  int nz = 0;
  for (int iT0 = 0; iT0 < nrow; ++iT0) {
    if ((int) S[iT0].size() != ncol) {
      pcout << "complab3d_glpkcpp: " << siter << "metabolic model has a ragged stoichiometric matrix "
            << "(row " << iT0 << ")." << std::endl;
      return CLB3D_GLPK_ERR_DIM;
    }
    for (int iT1 = 0; iT1 < ncol; ++iT1) {
      if (S[iT0][iT1] != 0) { ++nz; }
    }
  }

  std::vector<int>    rn(nz+1, 0), cn(nz+1, 0);
  std::vector<double> a(nz+1, 0.0);

  int iz = 0;
  for (int iT0 = 0; iT0 < nrow; ++iT0) {
    for (int iT1 = 0; iT1 < ncol; ++iT1) {
      if (S[iT0][iT1] != 0) {
        ++iz;
        rn[iz] = iT0 + 1;
        cn[iz] = iT1 + 1;
         a[iz] = S[iT0][iT1]; // a matrix containing the constraints coefficients.
      }
    }
  }

  std::vector<double> b(nrow, 0.0);
  for (int iT0 = 0; iT0 < nrow; ++iT0) {
    b[iT0] = vec_b[iT0];
  }

  //-- freeLB/UB arguments, default: Free.
  // std::vector<bool> is bit-packed and has no contiguous bool* storage, so a
  // plain owning array is used here; unique_ptr frees it on every return path
  // (the 2-D calloc'd versions leaked).
  std::unique_ptr<bool[]> freeLB(new bool[ncol]());
  std::unique_ptr<bool[]> freeUB(new bool[ncol]());
  for (int iT = 0; iT < ncol; iT++) {
    freeLB[iT] = (std::isinf(lb[iT]) || lb[iT] <= -1e30) ? true : false;   // [FIX-3D] 1e30 sentinel
    freeUB[iT] = (std::isinf(ub[iT]) || ub[iT] >=  1e30) ? true : false;   // [FIX-3D] 1e30 sentinel
  }

  //-- vartype arguments, default: Continuous.
  // NOTE (preserved from 2-D): the loop runs to ncol+1, i.e. it reads one
  // element past the reactions.  Every CompLaB call site allocates
  // vtype[nrxns+1] and memsets the whole thing to 'C', so this is in bounds,
  // and set_glpk only ever reads vartype[0..ncol-1].  Left as-is so that a
  // caller passing a genuinely (ncol+1)-long vtype keeps the same isMIP
  // detection.
  isMIP = 0;
  std::vector<int> vartype(ncol+1, GLP_CV);
  for (int iT = 0; iT < (ncol+1) ; ++iT) {
    switch ( vtype[iT] ) {
      case 'I': vartype[iT] = GLP_IV; isMIP = 1; break;
      case 'B': vartype[iT] = GLP_BV; isMIP = 1; break;
      default : vartype[iT] = GLP_CV;
    }
  }

  //-- Save option (off by default: every CompLaB call site passes save_pb = 0)
  char save_filename[32];
  strcpy(save_filename,"glpk_output");
  if (save_pb > 0) {
    char save_filetype[8]; // .mps; .txt; .lp
    if (save_pb == 1) {
      strcpy(save_filetype,".lp");
      if (isMIP == 1) {
        pcout << "WARNING: save file type incompatible with the input vtype" << std::endl;
      }
    }
    else if (save_pb == 2 || save_pb == 3) {
      strcpy(save_filetype,".mps");
      if (isMIP == 0) {
        pcout << "WARNING: save file type incompatible with the input vtype" << std::endl;
      }
    }
    else if (save_pb == 4) {
      strcpy(save_filetype,".mps");
    }
    else {
      pcout << "complab3d_glpkcpp: unspecified save_pb (" << save_pb << ")" << std::endl;
      return CLB3D_GLPK_ERR_SAVEPB;   // [FIX-3D] was exit(EXIT_FAILURE)
    }
    strcat(save_filename,save_filetype);
  }

  glp_init_iocp(&iParam);
  glp_init_smcp(&sParam);

  int errchk = set_glpk (lp, sense, nrow, ncol, nz, method, &c[0], &rn[0], &cn[0], &a[0], &b[0],
          ctype, freeLB.get(), &lb[0], freeUB.get(), &ub[0], &vartype[0], isMIP, lpsolver, sParam, iParam,
          save_pb, save_filename, glpIntParam, IParam, glpRealParam, RParam);

  if (!errchk) {
    pcout << siter << "metabolic model has been initialized successfully. "
          << "(" << nrow << " metabolites, " << ncol << " reactions, " << nz << " nonzeros)" << std::endl;
    return 0;
  }
  else {
    // [FIX-3D] was exit(EXIT_FAILURE).  Report and let the caller abort
    // cleanly and collectively.
    pcout << "complab3d_glpkcpp: ERROR in set_glpk for the " << siter
          << "metabolic model (code " << errchk << "). This model is unusable." << std::endl;
    return (errchk != 0) ? errchk : CLB3D_GLPK_ERR_SETUP;
  }
}

#endif  // COMPLAB_ENABLE_GLPK

#endif  // COMPLAB3D_GLPKCPP_HH
