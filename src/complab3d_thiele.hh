/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =================================================================================================
 * complab3d_thiele.hh  --  the classical effectiveness factor, and the convention it belongs to
 *
 * One function, split out of complab3d_upscale.hh so that it can be tested on its own: the rest of
 * that header needs Palabos, this does not, and tests/test_upscale.cpp checks this file with
 * nothing but <cmath>.
 *
 * WHY THE CONVENTION IS WRITTEN ON THE TIN
 *
 *   The sphere with first-order kinetics has one solution and two common spellings, and they look
 *   nothing like each other because they are written for two different characteristic lengths:
 *
 *       on the RADIUS            phi = R sqrt(k/D)         eta = (3/phi^2)(phi coth phi - 1)
 *       on volume-over-surface   phi = (R/3) sqrt(k/D)     eta = (1/phi)(coth 3phi - 1/(3phi))
 *
 *   The second is often called the generalized or Aris modulus. Substituting one modulus into the
 *   other's formula evaluates the curve at three times the right argument, and nothing about the
 *   result looks wrong: it is still between 0 and 1, still falls monotonically, still tends to the
 *   right shape. It is simply too low, by about 0.15 over the range a pore-scale aggregate lives
 *   in, and only a comparison against an independent calculation catches it.
 *
 *   [v1.3] That is exactly what had happened. complab.cpp computes phi on the RADIUS, every README
 *   prints it that way, and this function was the generalized spelling. Every eta_classical the
 *   solver reported was eta(3 phi), which inverted the headline conclusion of example 20: the
 *   measured factor appeared to sit above the classical curve when it sits below it.
 *
 *   So: this file is the RADIUS convention, it says so three times, and the test pins it by the
 *   small-phi series (1 - phi^2/15 here against 1 - 0.6 phi^2 there, a factor of nine in the first
 *   correction) and by the large-phi asymptote (3/phi here against 1/phi there).
 * ============================================================================================== */

#ifndef COMPLAB3D_THIELE_HH
#define COMPLAB3D_THIELE_HH

#include <cmath>

namespace complab_upscale {

/* eta = (3/phi^2)(phi coth phi - 1), for phi = R sqrt(k/D).
 *
 * Below phi = 1e-3 the closed form is the difference of two nearly equal large numbers and loses
 * every digit it has, so the series is taken explicitly there. The two branches agree to better
 * than 1e-8 across the switch, which the test checks. */
inline double classicalEta(double phi)
{
    if (!(phi > 0)) return 1.0;
    if (phi < 1e-3) return 1.0 - phi * phi / 15.0;     /* the series, good to 1e-12 there */
    return (3.0 / (phi * phi)) * (phi / std::tanh(phi) - 1.0);
}

}  /* namespace complab_upscale */

#endif
