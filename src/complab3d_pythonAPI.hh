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

#ifndef COMPLAB3D_PYTHONAPI_HH
#define COMPLAB3D_PYTHONAPI_HH

/* ===============================================================================================================
   ============================ OPTIONAL EMBEDDED-PYTHON / COBRApy BRIDGE (CompLaB3D) ============================
   ===============================================================================================================

   This whole file is a no-op unless COMPLAB_ENABLE_COBRAPY is defined at compile time. Builds without a
   Python development environment simply include this header and get nothing, so the COBRApy path is a
   strictly optional feature of CompLaB3D:

       g++ -DCOMPLAB_ENABLE_COBRAPY $(python3-config --includes) ... $(python3-config --ldflags) -lpython3.x

   The header assumes that Palabos has already been included, that `using namespace plb;` is in effect and
   that `typedef double T;` is in scope (both come from complab3d_processors.hh). It is meant to be included
   in exactly ONE translation unit (as all other CompLaB *.hh files are); the module/function cache below is
   file-static and relies on that.

   Ported from the 2-D reference implementation src/complab_pythonAPI.hh. Every behavioural change with
   respect to the 2-D code is marked with a `// [FIX-3D]` comment that names the offending 2-D line.
*/

#ifdef COMPLAB_ENABLE_COBRAPY

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <vector>
#include <string>
#include <cmath>
#include <cstdio>

/* ===============================================================================================================
   ============================== REMOVED FROM THE 2-D VERSION (DELIBERATELY NOT PORTED) =========================
   ===============================================================================================================

   1) `check_model(char*, std::vector<PyObject*>, int, int)`   -- 2D complab_pythonAPI.hh:199-221
      Removed because it is dead AND harmful. It looks up a Python function called "check_model", which does
      not exist in complab_cobrapy.py (nor in complab3d_cobrapy.py), so the PyCallable_Check branch is never
      taken. Worse, had it ever been taken, it hands the cached model objects to PyTuple_SetItem (which STEALS
      the reference) and then calls Py_DECREF(pArgs) -- destroying every cached cobra model on the spot. It is
      never called from complab.cpp. Deleted rather than repaired.

   2) `load_metabolic_model(char*, char*, std::vector<const char*>, ...)` -- 2D complab_pythonAPI.hh:307-378
      Removed because it is never called: CompLaB loads stoichiometry from the model XML through
      load_metabolic_models() in complab_functions.hh, and only the resulting dense arrays are shipped to
      COBRApy by prep_cobrapy(). Keeping a second, unexercised model loader alive is a maintenance trap.
      (Its only unique feature -- getDoubleArray3D -- is kept below, since the helper itself is still useful.)

   ===============================================================================================================
   ================================== GIL / MPI / THREADING -- READ BEFORE REUSING ================================
   ===============================================================================================================

   There is NO PyGILState_Ensure()/PyGILState_Release() anywhere in this file, and none in the 2-D original
   either. That is not an oversight that happens to work; it is a hard constraint on how this code may be used:

     * Palabos executes data processors serially inside a single rank (one processor at a time, one thread per
       MPI rank). Only that serialization makes it legal to touch the interpreter from a data processor without
       ever acquiring the GIL: the calling thread is the thread that ran Py_Initialize(), and it is the only
       thread that ever enters the interpreter.

     * Consequently this bridge is UNSUITABLE for any thread-parallel build (OpenMP-inside-a-rank, a threaded
       Palabos backend, or any user code that calls optimize_cobrapy() from more than one thread). Two threads
       entering CPython without the GIL corrupt reference counts and crash or, worse, silently produce garbage
       fluxes. If CompLaB3D ever grows intra-rank threading, every Py* call below must be wrapped in
       PyGILState_Ensure()/PyGILState_Release() and the model cache must become thread-local.

     * Under MPI each rank runs its OWN interpreter (Py_Initialize() is called per process) and therefore holds
       its OWN full copy of every metabolic model. Memory scales as
             (size of all cobra models) x (MPI ranks per node),
       which for genome-scale models is the dominant memory term on fat nodes. There is no sharing and no
       copy-on-write benefit, because the models are built (not forked) inside each interpreter. Size jobs
       accordingly, or use the GLPK path, which keeps one compact glp_prob per model per rank.

     * Py_Initialize() must be called by the host program before prep_cobrapy(), and Py_FinalizeEx() after
       finalize_cobrapy().

   ===============================================================================================================
   ================================ FLUX-VECTOR CONTRACT (C++ <-> complab3d_cobrapy.py) ==========================
   ===============================================================================================================

   THIS BLOCK IS DUPLICATED VERBATIM AT THE TOP OF complab3d_cobrapy.py. If you change one, change both.

   Argument vector handed to run_FBA() by optimize_cobrapy() (all values are Python floats, models appended
   last as opaque handles):

       args[0]                                    = M            number of microbes present in this voxel
       args[1 + i]                                = n_i          number of USED exchange reactions of microbe i
       args[1 + M + off_i + k]                    = loc[i][k]    reaction INDEX inside microbe i's own model
       args[1 + M + N + off_i + k]                = lb[i][k]     lower bound for reaction loc[i][k]
       args[1 + M + 2N + off_i + k]               = ub[i][k]     upper bound for reaction loc[i][k]
       args[len(args) - M + i]                    = model handle of microbe i (returned by feed_cobraModel)

       with N = sum_i n_i,  off_i = sum_{j<i} n_j,  k = 0 .. n_i-1.

   Return value of run_FBA(): the Python list [grate, flux, status] with

       grate[i]     float   objective value of microbe i, or 0.0 when the LP did not solve
       flux[i]      list of EXACTLY n_i floats.  flux[i][k] is the flux of reaction loc[i][k].
                    ==> flux[i] is indexed by the PACKED index k, in the SAME ORDER as loc[i] was supplied.
                        It is NOT indexed by the global substrate id iS of the calling processor.
       status[i]    str     solver status of microbe i ('optimal' when the LP solved, anything else = failure)

   [FIX-3D] 2D complab_pythonAPI.hh:412 / complab_processors.hh:395 -- the flux-vector index mismatch.
       The 2-D processor reads `vec2_flux[iB][iS]` with iS running over ALL substrates (0..subsNum-1), while
       run_FBA returns a vector packed over that microbe's USED exchange reactions only (length n_i <= subsNum).
       As soon as one microbe does not use one substrate, every subsequent flux is read from the wrong slot and
       silently attributed to the wrong solute -- and for iS >= n_i the read is out of bounds. The 2-D GLPK path
       does the right thing (it pads with 0.0 for unused substrates, complab_processors.hh:151-159), so the two
       solver back-ends disagree.
       Fixed here by (a) pinning the contract above: length(flux[i]) == n_i == length(loc[i]), same order,
       (b) VERIFYING that length at run time in optimize_cobrapy() against n_i taken from args[1+i], and
       (c) providing expandFluxToSubstrates() below, which converts the packed vector into the full
       subsNum-long, substrate-indexed vector the processor wants -- padding unused substrates with 0.0 exactly
       like the GLPK path. Processors must call expandFluxToSubstrates(); they must never index the packed
       vector with a substrate id.
*/

/* ---------------------------------------------------------------------------------------------------------------
   Error codes returned by the helpers and by the public entry points.
   [FIX-3D] 2D complab_pythonAPI.hh:50,56,62,68,88,94,99,114,120,135,141,161,167,173,189,194 -- every
   unmarshalling helper called exit(EXIT_FAILURE) on malformed input. A single rank calling exit() under MPI
   tears down that process while its peers sit in a collective forever: the job hangs until the wall clock
   kills it, and no output is flushed. All of them now report through pcout and return an error code.
   --------------------------------------------------------------------------------------------------------------- */
enum CobrapyError {
    COBRAPY_OK            = 0,
    COBRAPY_NOT_A_LIST    = 1,  // expected a Python list at this nesting level, got something else
    COBRAPY_BAD_NESTING   = 2,  // found a nested list where a scalar was expected
    COBRAPY_PY_EXCEPTION  = 3,  // a Python exception was raised (traceback already printed)
    COBRAPY_NULL_ARG      = 4,  // NULL PyObject handed to a helper
    COBRAPY_NOT_PREPPED   = 5,  // optimize_cobrapy() called before prep_cobrapy() succeeded
    COBRAPY_BAD_LENGTH    = 6,  // returned vector violates the flux-vector contract above
    COBRAPY_CALL_FAILED   = 7   // the Python call itself failed
};

/* Report and clear a pending Python exception. Returns true if there was one.
   [FIX-3D] the 2-D code checked PyErr_Occurred() in only two places, so a Python-side exception typically
   surfaced as a NULL dereference or as a silent 0.0 from PyFloat_AsDouble instead of a traceback. */
inline bool cobrapyCheckError(const char *where)
{
    if (PyErr_Occurred()) {
        pcout << "COBRApy: Python exception in " << where << " -- traceback follows." << std::endl;
        PyErr_Print();   // prints to stderr AND clears the error indicator
        return true;
    }
    return false;
}

/* ===============================================================================================================
   ==================================== PyObject -> std::vector UNMARSHALLING =====================================
   Names and nesting semantics are unchanged from the 2-D version; the signatures now return an error code and
   deliver the data through an output reference instead of returning it by value and exit()ing on failure.
   =============================================================================================================== */

// PyObject -> Vector (DOUBLE)
inline int getDoubleArray3D(PyObject *inputArg, std::vector< std::vector< std::vector<double> > > &data)
{
    if (inputArg == NULL) { pcout << "getDoubleArray3D ERROR: NULL argument." << std::endl; return COBRAPY_NULL_ARG; }
    if (!PyList_Check(inputArg)) {
        pcout << "getDoubleArray3D ERROR: 1st dimension is not a list." << std::endl;
        return COBRAPY_NOT_A_LIST;
    }
    data.assign(PyList_Size(inputArg), std::vector< std::vector<double> >());
    for (Py_ssize_t iT2 = 0; iT2 < PyList_Size(inputArg); ++iT2) {
        PyObject *value2 = PyList_GetItem(inputArg, iT2);   // borrowed reference
        if (!PyList_Check(value2)) {
            pcout << "getDoubleArray3D ERROR: 2nd dimension of the stoichiometric matrix is not a list." << std::endl;
            return COBRAPY_NOT_A_LIST;
        }
        data[iT2].assign(PyList_Size(value2), std::vector<double>());
        for (Py_ssize_t iT0 = 0; iT0 < PyList_Size(value2); ++iT0) {
            PyObject *value0 = PyList_GetItem(value2, iT0);
            if (!PyList_Check(value0)) {
                pcout << "getDoubleArray3D ERROR: 3rd dimension of the stoichiometric matrix is not a list." << std::endl;
                return COBRAPY_NOT_A_LIST;
            }
            data[iT2][iT0].assign(PyList_Size(value0), 0.);
            for (Py_ssize_t iT1 = 0; iT1 < PyList_Size(value0); ++iT1) {
                PyObject *value1 = PyList_GetItem(value0, iT1);
                if (PyList_Check(value1)) {
                    pcout << "getDoubleArray3D ERROR: 4th dimension of the stoichiometric matrix is still a list." << std::endl;
                    return COBRAPY_BAD_NESTING;
                }
                data[iT2][iT0][iT1] = PyFloat_AsDouble(value1);
                if (cobrapyCheckError("getDoubleArray3D")) { return COBRAPY_PY_EXCEPTION; }
            }
        }
    }
    return COBRAPY_OK;
}

inline int getDoubleArray2D(PyObject *inputArg, std::vector< std::vector<double> > &data)
{
    if (inputArg == NULL) { pcout << "getDoubleArray2D ERROR: NULL argument." << std::endl; return COBRAPY_NULL_ARG; }
    if (!PyList_Check(inputArg)) {
        pcout << "getDoubleArray2D ERROR: 1st dimension is not a list." << std::endl;
        return COBRAPY_NOT_A_LIST;
    }
    data.assign(PyList_Size(inputArg), std::vector<double>());
    for (Py_ssize_t iT0 = 0; iT0 < PyList_Size(inputArg); ++iT0) {
        PyObject *value0 = PyList_GetItem(inputArg, iT0);
        if (!PyList_Check(value0)) {
            pcout << "getDoubleArray2D ERROR: 2nd dimension is not a list." << std::endl;
            return COBRAPY_NOT_A_LIST;
        }
        data[iT0].assign(PyList_Size(value0), 0.);
        for (Py_ssize_t iT1 = 0; iT1 < PyList_Size(value0); ++iT1) {
            PyObject *value1 = PyList_GetItem(value0, iT1);
            if (PyList_Check(value1)) {
                pcout << "getDoubleArray2D ERROR: 3rd dimension is still a list." << std::endl;
                return COBRAPY_BAD_NESTING;
            }
            data[iT0][iT1] = PyFloat_AsDouble(value1);
            if (cobrapyCheckError("getDoubleArray2D")) { return COBRAPY_PY_EXCEPTION; }
        }
    }
    return COBRAPY_OK;
}

inline int getDoubleArray1D(PyObject *inputArg, std::vector<double> &data)
{
    if (inputArg == NULL) { pcout << "getDoubleArray1D ERROR: NULL argument." << std::endl; return COBRAPY_NULL_ARG; }
    if (!PyList_Check(inputArg)) {
        pcout << "getDoubleArray1D ERROR: 1st dimension of the array is not a list." << std::endl;
        return COBRAPY_NOT_A_LIST;
    }
    data.assign(PyList_Size(inputArg), 0.);
    for (Py_ssize_t iT0 = 0; iT0 < PyList_Size(inputArg); ++iT0) {
        PyObject *value0 = PyList_GetItem(inputArg, iT0);
        if (PyList_Check(value0)) {
            pcout << "getDoubleArray1D ERROR: 2nd dimension of the array shouldn't be a list." << std::endl;
            return COBRAPY_BAD_NESTING;
        }
        data[iT0] = PyFloat_AsDouble(value0);
        if (cobrapyCheckError("getDoubleArray1D")) { return COBRAPY_PY_EXCEPTION; }
    }
    return COBRAPY_OK;
}

inline int getPyObject1D(PyObject *inputArg, std::vector<PyObject *> &data)
{
    if (inputArg == NULL) { pcout << "getPyObject1D ERROR: NULL argument." << std::endl; return COBRAPY_NULL_ARG; }
    if (!PyList_Check(inputArg)) {
        pcout << "getPyObject1D ERROR: 1st dimension of the array is not a list." << std::endl;
        return COBRAPY_NOT_A_LIST;
    }
    data.assign(PyList_Size(inputArg), (PyObject *) NULL);
    for (Py_ssize_t iT0 = 0; iT0 < PyList_Size(inputArg); ++iT0) {
        PyObject *value0 = PyList_GetItem(inputArg, iT0);   // borrowed reference, as in the 2-D version
        if (PyList_Check(value0)) {
            pcout << "getPyObject1D ERROR: 2nd dimension of the array shouldn't be a list." << std::endl;
            return COBRAPY_BAD_NESTING;
        }
        data[iT0] = value0;
    }
    return COBRAPY_OK;
}

// PyObject -> Vector (INT)
inline int getIntArray2D(PyObject *inputArg, std::vector< std::vector<int> > &data)
{
    if (inputArg == NULL) { pcout << "getIntArray2D ERROR: NULL argument." << std::endl; return COBRAPY_NULL_ARG; }
    if (!PyList_Check(inputArg)) {
        pcout << "getIntArray2D ERROR: 1st dimension is not a list." << std::endl;
        return COBRAPY_NOT_A_LIST;
    }
    data.assign(PyList_Size(inputArg), std::vector<int>());
    for (Py_ssize_t iT0 = 0; iT0 < PyList_Size(inputArg); ++iT0) {
        PyObject *value0 = PyList_GetItem(inputArg, iT0);
        if (!PyList_Check(value0)) {
            pcout << "getIntArray2D ERROR: 2nd dimension is not a list." << std::endl;
            return COBRAPY_NOT_A_LIST;
        }
        data[iT0].assign(PyList_Size(value0), 0);
        for (Py_ssize_t iT1 = 0; iT1 < PyList_Size(value0); ++iT1) {
            PyObject *value1 = PyList_GetItem(value0, iT1);
            if (PyList_Check(value1)) {
                pcout << "getIntArray2D ERROR: 3rd dimension is still a list." << std::endl;
                return COBRAPY_BAD_NESTING;
            }
            data[iT0][iT1] = (int) PyLong_AsLong(value1);
            if (cobrapyCheckError("getIntArray2D")) { return COBRAPY_PY_EXCEPTION; }
        }
    }
    return COBRAPY_OK;
}

inline int getIntArray1D(PyObject *inputArg, std::vector<int> &data)
{
    if (inputArg == NULL) { pcout << "getIntArray1D ERROR: NULL argument." << std::endl; return COBRAPY_NULL_ARG; }
    if (!PyList_Check(inputArg)) {
        pcout << "getIntArray1D ERROR: 1st dimension of the array is not a list." << std::endl;
        return COBRAPY_NOT_A_LIST;
    }
    data.assign(PyList_Size(inputArg), 0);
    for (Py_ssize_t iT0 = 0; iT0 < PyList_Size(inputArg); ++iT0) {
        PyObject *value0 = PyList_GetItem(inputArg, iT0);
        if (PyList_Check(value0)) {
            pcout << "getIntArray1D ERROR: 2nd dimension of the array shouldn't be a list." << std::endl;
            return COBRAPY_BAD_NESTING;
        }
        data[iT0] = (int) PyLong_AsLong(value0);
        if (cobrapyCheckError("getIntArray1D")) { return COBRAPY_PY_EXCEPTION; }
    }
    return COBRAPY_OK;
}

/* Solver-status list (list of Python strings) -> per-microbe integer flag.
   1 = 'optimal' (usable solution), 0 = anything else (infeasible / unbounded / solver error).
   New in 3-D: the 2-D bridge had no notion of solver status at all. */
inline int getStatusArray1D(PyObject *inputArg, std::vector<int> &data)
{
    if (inputArg == NULL) { pcout << "getStatusArray1D ERROR: NULL argument." << std::endl; return COBRAPY_NULL_ARG; }
    if (!PyList_Check(inputArg)) {
        pcout << "getStatusArray1D ERROR: the status entry is not a list." << std::endl;
        return COBRAPY_NOT_A_LIST;
    }
    data.assign(PyList_Size(inputArg), 0);
    for (Py_ssize_t iT0 = 0; iT0 < PyList_Size(inputArg); ++iT0) {
        PyObject *value0 = PyList_GetItem(inputArg, iT0);
        if (value0 == NULL || !PyUnicode_Check(value0)) {
            pcout << "getStatusArray1D ERROR: status entry " << (plint) iT0 << " is not a string." << std::endl;
            return COBRAPY_NOT_A_LIST;
        }
        const char *s = PyUnicode_AsUTF8(value0);           // borrowed, valid while value0 lives
        if (cobrapyCheckError("getStatusArray1D")) { return COBRAPY_PY_EXCEPTION; }
        data[iT0] = (s != NULL && std::string(s) == "optimal") ? 1 : 0;
    }
    return COBRAPY_OK;
}

/* ===============================================================================================================
   ================================= CACHED PYTHON MODULE / FUNCTION OBJECTS =====================================
   [FIX-3D] 2D complab_pythonAPI.hh:383-387 -- optimize_cobrapy() called PyUnicode_FromString + PyImport_Import
   + PyObject_GetAttrString on EVERY invocation, i.e. once per voxel per time step. PyImport_Import hits
   sys.modules so it does not re-execute the module, but it still builds a unicode object, takes the import
   lock, walks the module dict and creates/destroys the bound function object every single call. Measured
   against the LP itself this is the single largest remaining per-voxel overhead. The module and the run_FBA
   callable are now looked up ONCE in prep_cobrapy() and held as strong references until finalize_cobrapy().

   These are file-static on purpose: this header is included in exactly one translation unit. They are also
   per-process, which is exactly right under MPI (one interpreter per rank) and exactly wrong under threads
   (see the GIL note at the top).
   =============================================================================================================== */
static PyObject *cobrapy_pModule  = NULL;   // strong reference to the complab3d_cobrapy module
static PyObject *cobrapy_pRunFBA  = NULL;   // strong reference to complab3d_cobrapy.run_FBA

// The Python module shipped with CompLaB3D. Renamed from the 2-D "complab_cobrapy".
static const char *COMPLAB3D_PY_MODULE = "complab3d_cobrapy";

/* Release the cached module/function (and, optionally, the cached cobra models).
   Call before Py_FinalizeEx(). Safe to call more than once and safe to call if prep_cobrapy() never ran. */
inline void finalize_cobrapy(std::vector<PyObject *> *vec_model = NULL)
{
    if (vec_model != NULL) {
        for (size_t iM = 0; iM < vec_model->size(); ++iM) {
            Py_XDECREF((*vec_model)[iM]);       // strong references handed out by prep_cobrapy()
            (*vec_model)[iM] = NULL;
        }
    }
    Py_XDECREF(cobrapy_pRunFBA); cobrapy_pRunFBA = NULL;
    Py_XDECREF(cobrapy_pModule); cobrapy_pModule = NULL;
}

/* ===============================================================================================================
   ============================================== prep_cobrapy ===================================================
   Builds one cobra model per microbe and stores a strong reference to it in vec_model, then caches the module
   and the run_FBA callable for optimize_cobrapy(). Signature and parameter order are unchanged from 2-D
   (arrays are now taken by const reference -- they are read-only and can be several hundred MB).

   Wire format for feed_cobraModel (unchanged from 2-D, dense):
       args[0] = nmets, args[1] = nrxns,
       args[2 .. 2+nmets*nrxns-1]            = S, metabolite-major: S[i*nrxns+j] is metabolite i in reaction j
       args[2+ nmets   *nrxns .. +nrxns-1]   = lower bounds
       args[2+(nmets+1)*nrxns .. +nrxns-1]   = upper bounds
       args[2+(nmets+2)*nrxns .. +nrxns-1]   = objective coefficients
   NOTE this transfers the FULL dense stoichiometric matrix as boxed Python floats. It is a one-off start-up
   cost (not per voxel), so it is kept as-is for wire compatibility, but it is the reason model loading is slow
   for genome-scale networks; the Python side compresses it to a sparse representation immediately.
   Returns 0 on success, non-zero on failure.
   =============================================================================================================== */
inline int prep_cobrapy(char *pyFileName, char *src_path, std::vector<PyObject *> &vec_model,
                        const std::vector< std::vector<T> > &S, const std::vector< std::vector<T> > &c,
                        const std::vector< std::vector<T> > &lb, const std::vector< std::vector<T> > &ub,
                        const std::vector<int> &vec_nrxns, const std::vector<int> &vec_nmets,
                        const std::vector<plint> &vec_objLoc,
                        const std::vector<int> &vec_sense)
{
    (void) vec_objLoc;   // kept for signature compatibility with the 2-D API; the objective is carried in c

    /* [FIX-3D] vec_sense is new.  <objective_direction> reached GLPK (as `sense`,
     * -1 maximize / +1 minimize) but never reached COBRApy, so a microbe asking
     * to MINIMISE was silently maximised on that back end -- and since the two
     * back ends cannot be mixed in one run, nothing would ever have caught it by
     * cross-checking.  It is appended as one extra element on the wire so the
     * rest of the format is unchanged. */

    pcout << "Feeding metabolic models to COBRApy ..." << std::endl;
    const char funName[] = "feed_cobraModel";     // python function name inside COMPLAB3D_PY_MODULE
    plint num_of_microbes = (plint) S.size();

    // Make the CompLaB3D src directory importable.
    PyObject *sysPath = PySys_GetObject((char *) "path");
    if (sysPath == NULL) {
        pcout << "COBRApy ERROR: cannot access sys.path -- was Py_Initialize() called?" << std::endl;
        return COBRAPY_PY_EXCEPTION;
    }
    PyObject *pPath = PyUnicode_FromString(src_path);
    PyList_Append(sysPath, pPath);      // PyList_Append does NOT steal
    Py_XDECREF(pPath);                  // [FIX-3D] 2D complab_pythonAPI.hh:230/324 leaked this unicode object

    // ------------------------------------------------------------------------------------------------------
    // Import the module ONCE and cache it together with run_FBA.
    // ------------------------------------------------------------------------------------------------------
    PyObject *pName = PyUnicode_FromString(COMPLAB3D_PY_MODULE);
    PyObject *pModule = PyImport_Import(pName);
    Py_XDECREF(pName);
    if (pModule == NULL) {
        cobrapyCheckError("prep_cobrapy (import)");
        pcout << "COBRApy ERROR: failed to import module \"" << COMPLAB3D_PY_MODULE
              << "\" (host passed \"" << (pyFileName ? pyFileName : "(null)") << "\"). "
              << "Check that " << COMPLAB3D_PY_MODULE << ".py is in " << src_path << "." << std::endl;
        return COBRAPY_PY_EXCEPTION;
    }

    PyObject *pFeed = PyObject_GetAttrString(pModule, funName);
    if (pFeed == NULL || !PyCallable_Check(pFeed)) {
        cobrapyCheckError("prep_cobrapy (getattr feed_cobraModel)");
        pcout << "COBRApy ERROR: cannot find callable \"" << funName << "\" in " << COMPLAB3D_PY_MODULE << "." << std::endl;
        Py_XDECREF(pFeed);
        Py_DECREF(pModule);
        return COBRAPY_PY_EXCEPTION;
    }

    int erck = COBRAPY_OK;
    for (plint j = 0; j < num_of_microbes; ++j) {
        int nmets = vec_nmets[j];
        int nrxns = vec_nrxns[j];
        int argsize = (nmets + 3) * nrxns + 3;   // +1 for the trailing objective sense

        std::vector<double> args(argsize);
        args[0] = double(nmets);
        args[1] = double(nrxns);
        /* trailing element: -1 maximize, +1 minimize (default maximize) */
        args[argsize - 1] = (j < (plint) vec_sense.size()) ? double(vec_sense[j]) : -1.0;
        for (int i = 0; i < nmets * nrxns; ++i) { args[2 + i] = S[j][i]; }
        for (int i = 0; i < nrxns; ++i) {
            args[2 +  nmets      * nrxns + i] = lb[j][i];
            args[2 + (nmets + 1) * nrxns + i] = ub[j][i];
            args[2 + (nmets + 2) * nrxns + i] = c[j][i];
        }

        PyObject *pArgs = PyTuple_New(argsize);
        if (pArgs == NULL) {
            cobrapyCheckError("prep_cobrapy (tuple alloc)");
            erck = COBRAPY_PY_EXCEPTION;
            break;
        }
        bool convOK = true;
        for (int i = 0; i < argsize; ++i) {
            PyObject *pValue = PyFloat_FromDouble(args[i]);
            if (pValue == NULL) {
                pcout << "COBRApy ERROR: cannot convert argument " << i << " of microbe " << j << "." << std::endl;
                cobrapyCheckError("prep_cobrapy (PyFloat_FromDouble)");
                convOK = false;
                break;
            }
            PyTuple_SetItem(pArgs, i, pValue);   // steals pValue
        }
        if (!convOK) { Py_DECREF(pArgs); erck = COBRAPY_PY_EXCEPTION; break; }

        vec_model[j] = PyObject_CallObject(pFeed, pArgs);   // strong reference, kept for the whole run
        Py_DECREF(pArgs);                                   // safe: pArgs owns only the boxed floats it stole
        if (vec_model[j] == NULL) {
            cobrapyCheckError("prep_cobrapy (feed_cobraModel)");
            pcout << "COBRApy ERROR: feed_cobraModel failed for microbe " << j << "." << std::endl;
            erck = COBRAPY_CALL_FAILED;
            break;
        }
    }

    Py_DECREF(pFeed);
    /* [FIX-3D] 2D complab_pythonAPI.hh:296 -- Py_DECREF(pModule) sat INSIDE the per-microbe loop although the
       import happened once outside it. With two or more microbes the module's reference count was decremented
       once per microbe, so from the second microbe on the code was dropping references it did not own; with
       enough microbes the module object is freed while pModule still points at it (use-after-free), and the
       loop then calls PyObject_GetAttrString on dead memory. It only ever "worked" because sys.modules keeps
       an extra reference alive for the common 1-2 microbe cases. The DECREF is now outside the loop -- and in
       fact the reference is not dropped at all here, it is transferred to the module cache below and released
       by finalize_cobrapy(). */
    if (erck != COBRAPY_OK) {
        Py_DECREF(pModule);
        return erck;
    }

    // Cache module + run_FBA for optimize_cobrapy().
    finalize_cobrapy();                 // drop anything a previous prep_cobrapy() left behind
    cobrapy_pModule = pModule;          // ownership transferred (no DECREF here)
    cobrapy_pRunFBA = PyObject_GetAttrString(cobrapy_pModule, "run_FBA");
    if (cobrapy_pRunFBA == NULL || !PyCallable_Check(cobrapy_pRunFBA)) {
        cobrapyCheckError("prep_cobrapy (getattr run_FBA)");
        pcout << "COBRApy ERROR: cannot find callable \"run_FBA\" in " << COMPLAB3D_PY_MODULE << "." << std::endl;
        finalize_cobrapy();
        return COBRAPY_PY_EXCEPTION;
    }

    pcout << "COBRApy: " << (plint) num_of_microbes << " metabolic model(s) loaded; module \""
          << COMPLAB3D_PY_MODULE << "\" and run_FBA cached." << std::endl;
    return COBRAPY_OK;
}

/* ===============================================================================================================
   ============================================ optimize_cobrapy =================================================
   One FBA solve per microbe present in the current voxel. Called from the runFBA_cobrapy data processor.

   Parameter order follows the 2-D API; `status` is appended (per-microbe flag: 1 = optimal, 0 = otherwise) and
   the function now returns an error code instead of void. `fileName` is retained for source compatibility with
   the 2-D call site and is used for diagnostics only -- the module is resolved from the cache built by
   prep_cobrapy().

   `flux` obeys the FLUX-VECTOR CONTRACT documented at the top of this file: flux[iB] has EXACTLY n_iB entries,
   ordered like the loc list of microbe iB. Use expandFluxToSubstrates() to obtain a substrate-indexed vector.

   On a non-optimal status the growth rate and all fluxes of that microbe are forced to zero, mirroring the
   GLPK path (complab_processors.hh:156-159 in 2-D). Returns 0 on success.
   =============================================================================================================== */
inline int optimize_cobrapy(char *fileName, const std::vector<double> &args, const std::vector<PyObject *> &model,
                            std::vector<double> &grate, std::vector< std::vector<double> > &flux,
                            std::vector<int> &status)
{
    grate.clear(); flux.clear(); status.clear();

    if (cobrapy_pRunFBA == NULL) {
        pcout << "COBRApy ERROR: optimize_cobrapy() called before a successful prep_cobrapy() (module \""
              << (fileName ? fileName : COMPLAB3D_PY_MODULE) << "\")." << std::endl;
        return COBRAPY_NOT_PREPPED;
    }
    if (args.empty()) {
        pcout << "COBRApy ERROR: empty argument vector handed to optimize_cobrapy()." << std::endl;
        return COBRAPY_BAD_LENGTH;
    }

    const int len             = (int) args.size();
    const int num_of_microbes = (int) std::floor(args[0] + 0.5);
    if (num_of_microbes <= 0 || (int) model.size() < num_of_microbes || len < 1 + num_of_microbes) {
        pcout << "COBRApy ERROR: inconsistent argument vector (M = " << num_of_microbes
              << ", models = " << (plint) model.size() << ", args = " << len << ")." << std::endl;
        return COBRAPY_BAD_LENGTH;
    }

    // Expected packed flux length per microbe, straight from the contract: n_i = args[1+i].
    std::vector<int> nsubs(num_of_microbes);
    for (int i = 0; i < num_of_microbes; ++i) { nsubs[i] = (int) std::floor(args[1 + i] + 0.5); }

    PyObject *pArgs = PyTuple_New(len + num_of_microbes);
    if (pArgs == NULL) {
        cobrapyCheckError("optimize_cobrapy (tuple alloc)");
        return COBRAPY_PY_EXCEPTION;
    }
    for (int i = 0; i < len; ++i) {
        PyObject *pValue = PyFloat_FromDouble(args[i]);
        if (pValue == NULL) {
            pcout << "COBRApy ERROR: cannot convert argument " << i << " (optimize_cobrapy)." << std::endl;
            cobrapyCheckError("optimize_cobrapy (PyFloat_FromDouble)");
            Py_DECREF(pArgs);
            return COBRAPY_PY_EXCEPTION;
        }
        PyTuple_SetItem(pArgs, i, pValue);          // steals pValue
    }
    for (int i = 0; i < num_of_microbes; ++i) {
        /* [FIX-3D] 2D complab_pythonAPI.hh:405+409 -- PyTuple_SetItem STEALS the reference it is given, so
           handing it the cached model without owning an extra reference made the tuple co-owner of the model.
           The 2-D author noticed that Py_DECREF(pArgs) then destroyed the cached models and "fixed" it by
           commenting the DECREF out (line 409), which leaks the tuple AND every boxed float inside it on every
           solve, at every voxel, at every time step -- tens of megabytes per hundred thousand solves, growing
           without bound over a run.
           The correct fix is to give the tuple its own reference: INCREF before the steal, so the model
           survives Py_DECREF(pArgs) with its original reference count intact and the tuple is freed normally. */
        Py_INCREF(model[i]);
        PyTuple_SetItem(pArgs, len + i, model[i]);  // steals the reference we just added
    }

    PyObject *pValue = PyObject_CallObject(cobrapy_pRunFBA, pArgs);
    Py_DECREF(pArgs);                               // [FIX-3D] now correct and mandatory: no leak, no loss
    if (pValue == NULL) {
        cobrapyCheckError("optimize_cobrapy (run_FBA)");
        pcout << "COBRApy ERROR: run_FBA call failed." << std::endl;
        return COBRAPY_CALL_FAILED;
    }
    if (cobrapyCheckError("optimize_cobrapy (post-call)")) { Py_DECREF(pValue); return COBRAPY_PY_EXCEPTION; }

    if (!PyList_Check(pValue) || PyList_Size(pValue) < 3) {
        pcout << "COBRApy ERROR: run_FBA must return the list [grate, flux, status]." << std::endl;
        Py_DECREF(pValue);
        return COBRAPY_BAD_LENGTH;
    }

    int erck = getDoubleArray1D(PyList_GetItem(pValue, 0), grate);
    if (erck == COBRAPY_OK) { erck = getDoubleArray2D(PyList_GetItem(pValue, 1), flux); }
    if (erck == COBRAPY_OK) { erck = getStatusArray1D(PyList_GetItem(pValue, 2), status); }
    Py_DECREF(pValue);
    if (erck != COBRAPY_OK) {
        pcout << "COBRApy ERROR: could not unmarshal the run_FBA result (code " << erck << ")." << std::endl;
        return erck;
    }

    if ((int) grate.size() != num_of_microbes || (int) flux.size() != num_of_microbes ||
        (int) status.size() != num_of_microbes) {
        pcout << "COBRApy ERROR: run_FBA returned " << (plint) grate.size() << " growth rates, "
              << (plint) flux.size() << " flux vectors and " << (plint) status.size()
              << " status flags for " << num_of_microbes << " microbes." << std::endl;
        return COBRAPY_BAD_LENGTH;
    }

    for (int iB = 0; iB < num_of_microbes; ++iB) {
        // Contract check: the packed flux vector must match the loc list this microbe was given.
        if ((int) flux[iB].size() != nsubs[iB]) {
            pcout << "COBRApy ERROR: flux-vector contract violated for microbe " << iB << ": got "
                  << (plint) flux[iB].size() << " fluxes, expected " << nsubs[iB]
                  << " (one per used exchange reaction, in loc order)." << std::endl;
            return COBRAPY_BAD_LENGTH;
        }
        /* Infeasibility handling. cobra's model.optimize() yields objective_value = None for an infeasible LP;
           in the 2-D code that None travelled straight into PyFloat_AsDouble, which returns -1.0 and sets a
           Python exception that nobody ever checked -- so an infeasible voxel silently produced a growth rate
           of -1.0 h^-1 and garbage fluxes. The Python side now maps that case to 0.0 and reports a status; here
           we additionally hard-zero the whole microbe, exactly like the GLPK branch does for a bad glp status. */
        if (status[iB] != 1) {
            grate[iB] = 0.;
            for (size_t k = 0; k < flux[iB].size(); ++k) { flux[iB][k] = 0.; }
        }
    }
    return COBRAPY_OK;
}

/* ===============================================================================================================
   ========================================= expandFluxToSubstrates ==============================================
   Converts the PACKED flux vector returned for one microbe (length n_i, ordered like that microbe's loc list)
   into the FULL substrate-indexed vector (length subsNum) that the data processors consume, padding substrates
   the microbe does not exchange with 0.0.

   subsLocRow is the vec_EX_loc / vec2_subsLoc row of that microbe: subsLocRow[iS] is the reaction index of
   substrate iS in the microbe's model, or a negative value if the microbe does not use substrate iS. The packed
   order is defined by the caller that built `args`: it walks iS = 0..subsNum-1 and appends every entry with
   subsLocRow[iS] >= 0, so walking it the same way here reproduces the mapping exactly.

   This is the C++ half of the [FIX-3D] flux-vector index-mismatch repair. Processors must use it instead of
   indexing the packed vector with a substrate id.  Returns 0 on success.
   =============================================================================================================== */
inline int expandFluxToSubstrates(const std::vector<double> &packedFlux, const std::vector<int> &subsLocRow,
                                  plint subsNum, std::vector<T> &fullFlux)
{
    fullFlux.assign(subsNum, T());
    if ((plint) subsLocRow.size() < subsNum) {
        pcout << "expandFluxToSubstrates ERROR: substrate-location row has " << (plint) subsLocRow.size()
              << " entries, expected " << (plint) subsNum << "." << std::endl;
        return COBRAPY_BAD_LENGTH;
    }
    size_t k = 0;
    for (plint iS = 0; iS < subsNum; ++iS) {
        if (subsLocRow[iS] >= 0) {
            if (k >= packedFlux.size()) {
                pcout << "expandFluxToSubstrates ERROR: packed flux vector has only " << (plint) packedFlux.size()
                      << " entries but the microbe uses more substrates -- flux-vector contract violated."
                      << std::endl;
                return COBRAPY_BAD_LENGTH;
            }
            fullFlux[iS] = (T) packedFlux[k];
            ++k;
        }
        // else: substrate not exchanged by this microbe -> stays 0, as in the GLPK path.
    }
    if (k != packedFlux.size()) {
        pcout << "expandFluxToSubstrates ERROR: consumed " << (plint) k << " of " << (plint) packedFlux.size()
              << " packed fluxes -- flux-vector contract violated." << std::endl;
        return COBRAPY_BAD_LENGTH;
    }
    return COBRAPY_OK;
}

#endif  // COMPLAB_ENABLE_COBRAPY

#endif  // COMPLAB3D_PYTHONAPI_HH
