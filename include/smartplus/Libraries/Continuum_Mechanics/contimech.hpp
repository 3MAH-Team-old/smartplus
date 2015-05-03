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

///@file contimech.hpp
///@brief Functions that computes Mises stress/strains, directions, etc
///@version 1.0

#pragma once
#include <armadillo>

using namespace std;
using namespace arma;

namespace smart{

//This function returns the trace of the tensor v
double tr(const vec &);

//This function returns the deviatoric part of v
vec dev(const vec &);

//This function determines the Mises equivalent of a stress tensor, according to the Voigt convention for stress 
double Mises_stress(const vec &);

//This function determines the strain flow (direction) from a stress tensor, according to the Voigt convention for strains
vec eta_stress(const vec &);

//This function determines the Mises equivalent of a strain tensor, according to the Voigt convention for strains 
double Mises_strain(const vec &);

//This function determines the strain flow (direction) from a strain tensor, according to the Voigt convention for strains
vec eta_strain(const vec &);

//This function transforms the strain Voigt vector into a 3*3 strain matrix
mat v2t_strain(const vec &v);

//This function transforms a 3*3 strain matrix into a strain Voigt vector
vec t2v_strain (const mat &);

//This function transforms the stress Voigt vector into a 3*3 stress matrix
mat v2t_stress(const vec &);

//This function transforms a 3*3 stress matrix into a stress Voigt vector
vec t2v_stress (const mat &);

//Returns the second invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J2_stress(const vec &);

//Returns the second invariant of the deviatoric part of a second order strain tensor written as a Voigt vector
double J2_strain(const vec &);

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_stress(const vec &);

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_strain(const vec &);

//Returns the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3. The direction of the normalized vector is set by angles u
vec normal_ellipsoid(const double &, const double &, const double &, const double &, const double &);

//Returns the normal and tangent components of the stress vector in the normal direction n to an ellipsoid with axes a1, a2, a3. The direction of the normalized vector is set by angles u
vec sigma_int(const vec &, const double &, const double &, const double &, const double &, const double &);

///This computes the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)
mat p_ikjl(const vec &);

} //namespace smart
