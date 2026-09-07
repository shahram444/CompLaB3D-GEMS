/* ================================================================================================
 *  complab3d_configpath.hh  --  WHICH CONFIGURATION FILE THIS RUN READS
 *
 *  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
 *  Meile Lab, University of Georgia.  shahram.asgari@uga.edu
 *
 *  Six places in the program open the input file: initialize_complab(), initialize_metabolic() and
 *  four readers in complab.cpp.  All six used to open the literal string "CompLaB.xml", while
 *  main() accepted argv[1] and did nothing with it.  So
 *
 *      ./complab variant.xml
 *
 *  ran the case in CompLaB.xml and said nothing about it.  That is worse than refusing the
 *  argument, because the obvious way to test what a setting does -- run the case twice, changing
 *  one line -- returns two byte-identical results, and the natural reading of that is that the
 *  setting does nothing.
 *
 *  One string, set once from argv at the top of main(), read by all six.  A header of its own
 *  rather than a line in complab_functions.hh so that complab3d_metabolic.hh, which is compiled on
 *  its own by the tests, does not have to pull in Palabos to find it.
 * ================================================================================================ */

#ifndef COMPLAB3D_CONFIGPATH_HH
#define COMPLAB3D_CONFIGPATH_HH

#include <string>

namespace complab_input {

/* The path every reader opens.  Assign to it, once, before anything reads it. */
inline std::string &configPath()
{
    static std::string p("CompLaB.xml");
    return p;
}

}  /* namespace complab_input */

#endif  /* COMPLAB3D_CONFIGPATH_HH */
