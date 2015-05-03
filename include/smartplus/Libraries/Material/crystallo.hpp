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

///@file crystallo.hpp
///@brief Some definitions coming from crystallography
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

namespace smart{

//This function returns the Shmidt tensor (3x3 matrix)
mat Schmid(const vec &n, const vec &m);

//This function returns the Shmidt tensor (6 vector), with the convention of strain
vec Schmid_v(const vec &n, const vec &m);

//This function returns a matrix utilized for the Hill interfacial operator
mat F_nm(const vec &N);

//This function returns the Hill interfacial operator for an isotropic material
mat Q_nm(const vec &N, const double &mu, const double &lambda);

} //namespace smart
