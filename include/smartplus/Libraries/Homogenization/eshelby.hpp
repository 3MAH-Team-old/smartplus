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

///@file eshelby.hpp
///@brief Definition of the Eshelby tensor for ellipsoidal inclusions with
// Parts of this methods are copyrighted by Gavazzi & Lagoudas 1992 - Fair use only
///@version 1.0

#pragma once

#include <math.h>
#include <armadillo>
#include "../../parameter.hpp"

using namespace std;
using namespace arma;

namespace smart{

//Eshelby tensor for a sphere
mat Eshelby_sphere(const double &);

//	Eshelby tensor determination. The cylinder is oriented in such a way that the axis direction is the 1 direction. a2=a3 here
mat Eshelby_cylinder(const double &);

//	Eshelby tensor determination. The prolate shape is oriented in such a way that the axis direction is the 1 direction. a1>a2=a3 here
mat Eshelby_prolate(const double &, const double &);

//	Eshelby tensor determination. The oblate shape is oriented in such a way that the axis direction is the 1 direction. a1<a2=a3 here
mat Eshelby_oblate(const double &, const double &);

//This methods is using the Voigt notations for the tensors.
void calG(const double &, const double &, const double &, const double &, const double &, const Mat<int> &, const mat &, mat &);

//Weighted Gauss integration over a sphere to represent the integration over the ellipsoid
void Gauss(Mat<int> &, const mat &, mat &, const double &, const double &, const double &, const vec &, const vec &, const vec &, const vec &, const int &, const int &);

//Numerical Eshelby tensor determination
mat Eshelby(const mat &, const double &, const double &, const double &, const vec &, const vec &, const vec &, const vec &, const int &mp, const int &np);

//mat T_II_sphere(const double &, const double &); {

//Numerical Hill Interaction tensor determination
mat T_II(const mat &Lt, const double &, const double &, const double &, const vec &, const vec &, const vec &, const vec &, const int &, const int &);

//This function computes the integration points and weights
void points(vec &x, vec &wx, vec &y, vec &wy, const int &mp, const int &np);

} //namespace smart
