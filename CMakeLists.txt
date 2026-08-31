cmake_minimum_required(VERSION 3.5)

#set the project name
project(complab)
enable_language(CXX)

set(EXECUTABLE_NAME "complab")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

#if(NOT CMAKE_BUILD_TYPE)
    #set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
#endif()




if(NOT CMAKE_BUILD_TYPE)
    #set(CMAKE_BUILD_TYPE "Debug")
    set(CMAKE_BUILD_TYPE "Release")
endif()



message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")
message(${CMAKE_BUILD_TYPE})

# Compiler flags
# Append flags: set(CMAKE_XXX_FLAGS "${CMAKE_XXX_FLAGS} ...")
if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    message("GCC.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    message("Clang.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0 -faddress-sanitizer -O1 -fno-omit-frame-pointer")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("MSVC.")
    set(CMAKE_CXX_FLAGS_RELEASE "/Ox /Ot /GS- /GL /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/DPLB_DEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else()
    message( FATAL_ERROR "CXX compiler not recognized. CMake will quit." )
endif()

option(ENABLE_MPI "Enable MPI" ON)
if(ENABLE_MPI)
    message("Enabling MPI")
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
        #set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
        include_directories(${MPI_CXX_INCLUDE_PATH})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        add_definitions(-DPLB_MPI_PARALLEL)
    endif()
endif()

if(WIN32)
    option(ENABLE_POSIX "Enable POSIX" OFF)
else()
    option(ENABLE_POSIX "Enable POSIX" ON)
endif()

if(ENABLE_POSIX)
    message("Enabling POSIX")
    add_definitions(-DPLB_USE_POSIX)
endif()

if(APPLE)
    add_definitions(-DPLB_MAC_OS_X)
endif()

if(WIN32 OR CYGWIN)
    add_definitions(-DPLB_WINDOWS)
endif()
###############################################################################
# Palabos Library
###############################################################################


include_directories("/usr/local/include")
include_directories("/scratch/sa01687/Tahoma/versionControl/palabos-v2.3.0/src")
include_directories("/scratch/sa01687/Tahoma/versionControl/palabos-v2.3.0/externalLibraries")
include_directories("/scratch/sa01687/Tahoma/versionControl/palabos-v2.3.0/externalLibraries/Eigen3")
#include_directories("/apps/eb/Anaconda3/5.0.1/include/python3.6m") # location for Python.h
link_directories("/usr/local/lib")
#link_directories("/apps/eb/Anaconda3/5.0.1/lib")

file(GLOB_RECURSE PALABOS_SRC "/scratch/sa01687/Tahoma/versionControl/palabos-v2.3.0/src/*.cpp")
file(GLOB_RECURSE EXT_SRC "/scratch/sa01687/Tahoma/versionControl/palabos-v2.3.0/externalLibraries/tinyxml/*.cpp")

add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})


###############################################################################

add_executable(${EXECUTABLE_NAME} "./src/${EXECUTABLE_NAME}.cpp")

# Link with the following libraries
target_link_libraries(${EXECUTABLE_NAME} palabos)

###############################################################################
# OPTIONAL METABOLIC LAYER
#
#   CompLB3D runs perfectly well with hand-written kinetics alone.  These two
#   options add the flux balance analysis solvers back.  Both default to OFF so
#   that the code builds on a machine that has neither GLPK nor Python -- turn
#   on only what you actually have installed and intend to use.
#
#     cmake -DENABLE_GLPK=ON ..            flux balance analysis through GLPK
#     cmake -DENABLE_COBRAPY=ON ..         flux balance analysis through COBRApy
#
#   The surrogate-model solver needs neither and is always available.
#
#   Turning an option ON only makes the solver AVAILABLE.  You still switch it
#   on for a given run in CompLaB.xml, under <simulation_mode>.  If the XML asks
#   for a solver that was not built in, the program stops at start-up with a
#   message naming the option to enable, rather than ignoring you silently.
###############################################################################

option(ENABLE_GLPK "Build with GLPK flux balance analysis support" OFF)
if(ENABLE_GLPK)
    message("Enabling GLPK flux balance analysis")
    # Point these at your GLPK installation if it is not on the default path.
    # NOTE: complab3d_glpkcpp.hh supports GLPK 4.39 through 5.x; the legacy
    #       lpx_* / glp_long code paths are compiled only for 4.x below 4.45.
    find_path(GLPK_INCLUDE_DIR glpk.h)
    find_library(GLPK_LIBRARY NAMES glpk)
    if(GLPK_INCLUDE_DIR AND GLPK_LIBRARY)
        message("  found glpk.h in ${GLPK_INCLUDE_DIR}")
        message("  found libglpk at ${GLPK_LIBRARY}")
        include_directories(${GLPK_INCLUDE_DIR})
        target_link_libraries(${EXECUTABLE_NAME} ${GLPK_LIBRARY})
        target_compile_definitions(${EXECUTABLE_NAME} PRIVATE COMPLAB_ENABLE_GLPK)
    else()
        message(FATAL_ERROR
            "ENABLE_GLPK is ON but GLPK was not found.\n"
            "  Install it (Debian/Ubuntu: libglpk-dev; macOS: brew install glpk)\n"
            "  or point CMake at it:  -DGLPK_INCLUDE_DIR=/path/include -DGLPK_LIBRARY=/path/lib/libglpk.so")
    endif()
endif()

option(ENABLE_COBRAPY "Build with COBRApy flux balance analysis support" OFF)
if(ENABLE_COBRAPY)
    message("Enabling COBRApy flux balance analysis")
    # Needs the Python development headers (python3-dev / python-devel) AND the
    # cobra package importable by that same interpreter at run time.
    # NOTE: FindPython3 and Python3_Development_FOUND arrived in CMake 3.12.  Many
    #       HPC module trees still ship 3.5-3.11, where the call below fails with
    #       "could not find a package configuration file provided by Python3" --
    #       which reads as a missing dependency rather than an old CMake.  Fall back
    #       to the legacy module so the message is at least truthful.
    if(CMAKE_VERSION VERSION_LESS 3.12)
        message(WARNING "CMake ${CMAKE_VERSION} is older than 3.12; using the legacy FindPythonLibs.")
        find_package(PythonLibs)
        set(Python3_Development_FOUND ${PYTHONLIBS_FOUND})
        set(Python3_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
        set(Python3_LIBRARIES ${PYTHON_LIBRARIES})
        set(Python3_VERSION ${PYTHONLIBS_VERSION_STRING})
    else()
        find_package(Python3 COMPONENTS Development)
    endif()
    if(Python3_Development_FOUND)
        message("  found Python ${Python3_VERSION} headers in ${Python3_INCLUDE_DIRS}")
        include_directories(${Python3_INCLUDE_DIRS})
        target_link_libraries(${EXECUTABLE_NAME} ${Python3_LIBRARIES})
        target_compile_definitions(${EXECUTABLE_NAME} PRIVATE COMPLAB_ENABLE_COBRAPY)
    else()
        message(FATAL_ERROR
            "ENABLE_COBRAPY is ON but the Python development headers were not found.\n"
            "  Install them (Debian/Ubuntu: python3-dev; conda: the python package ships them)\n"
            "  and make sure `python3 -c \"import cobra\"` works for that interpreter.")
    endif()
endif()

# Speed option for the FBA processors: solve only on the bulk and let Palabos
# communicate the increments, instead of also solving a linear program for every
# envelope cell.  Typically 20-40% faster.  Verify against a default run on your
# own case before trusting it.
option(FBA_BULK_ONLY "Restrict the FBA processors to the bulk domain" OFF)
if(FBA_BULK_ONLY)
    message("FBA processors restricted to the bulk domain")
    target_compile_definitions(${EXECUTABLE_NAME} PRIVATE COMPLAB3D_FBA_BULK_ONLY)
endif()

if(ENABLE_MPI)
    target_link_libraries(${EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
endif()
