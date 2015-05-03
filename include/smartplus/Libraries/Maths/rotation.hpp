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

///@file rotation.hpp
///@brief rotation of a Voigt tensor
///@version 1.0

#pragma once

#include <armadillo>

using namespace std;
using namespace arma;

namespace smart{

void Rot_strain(vec &, const mat &);

void Rot_stress(vec &, const mat &);

mat fillQS(const double &, const int &);

mat fillQE(const double &, const int &);

//To rotate a stiffness matrix (6,6)
mat rotateL(const mat &, const double &, const int &);

//To rotate a compliance matrix (6,6)
mat rotateM(const mat &, const double &, const int &);

//To rotate a interaction matrix (6,6)
mat rotateA(const mat &, const double &, const int &);

//To rotate a stress vector (6)
vec rotate_stress(const vec &, const double &, const int &);

//To rotate a strain vector (6)
vec rotate_strain(const vec &, const double &, const int &);

//To rotate from local to global a stiffness matrix (6,6)
mat rotate_l2g_L(const mat &, const double &, const double &, const double &);

//To rotate from global to local a stiffness matrix (6,6)
mat rotate_g2l_L(const mat &, const double &, const double &, const double &);

//To rotate from local to global a localisation matrix (6,6)
mat rotate_l2g_A(const mat &, const double &, const double &, const double &);

//To rotate from global to local a localisation matrix (6,6)
mat rotate_g2l_A(const mat &, const double &, const double &, const double &);

} //namespace smart
