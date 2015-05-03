/* This file is part of SMART+.
 
 SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 SMART+ is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with SMART+.  If not, see <http://www.gnu.org/licenses/>.
 
 */

///@file stats.hpp
///@brief Usefull statistical functions
///@version 1.0

#pragma once

#include <armadillo>

using namespace std;
using namespace arma;

namespace smart{

///Approximation of a normal distribution
double normal_distrib(const double &, const double &, const double &);

//tri_sum of a and b
int tri_sum(const int &, const int &);

///1. Classic ODF: a1 * cos(Theta)^(2*p1) + a2 * cos(Theta)^(2*p2 + 1) * sin(Theta)^(2*p2) 
double ODF_sd(const double &, const double &, const double &, const double &, const double &, const double &);

///2. Classic ODF - hardening-like
double ODF_hard(const double &, const double &, const double &, const double &);

///3. Gaussian
double Gaussian(const double &, const double &, const double &, const double & = 1.);

///30. Several Gaussian
double Mult_Gaussian(const double &, const int &, const vec &, const vec &, const vec &);

///4. Lorentzian
double Lorentzian(const double &, const double &, const double &, const double & = 1.);

///40. Several Lorentzian
double Mult_Lorentzian(const double &, const int &, const vec &, const vec &, const vec &);

///5. Pseudo-Voigt
double PseudoVoigt(const double &, const double &, const double &, const double &, const double &, const double & = 1.);

///50. Several Pseudo-Voigt
double Mult_PseudoVoigt(const double &, const int &, const vec &, const vec &, const vec &, const vec &, const vec &);

///6. Pearson VII
double Pearson7(const double &, const double &, const double &, const double &, const double & = 1.);

///60. Several Pearson VII
double Mult_Pearson7(const double &, const int &, const vec &, const vec &, const vec &, const vec &);

} //namespace smart
