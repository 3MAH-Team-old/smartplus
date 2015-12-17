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

namespace smart{

void Rot_strain(arma::vec &, const arma::mat &);

void Rot_stress(arma::vec &, const arma::mat &);

arma::mat fillQS(const double &, const int &);

arma::mat fillQE(const double &, const int &);

//To rotate a stiffness matrix (6,6)
arma::mat rotateL(const arma::mat &, const double &, const int &);

//To rotate a compliance matrix (6,6)
arma::mat rotateM(const arma::mat &, const double &, const int &);

//To rotate a interaction matrix (6,6)
arma::mat rotateA(const arma::mat &, const double &, const int &);

//To rotate a stress vector (6)
arma::vec rotate_stress(const arma::vec &, const double &, const int &);

//To rotate a strain vector (6)
arma::vec rotate_strain(const arma::vec &, const double &, const int &);

//To rotate from local to global a stiffness matrix (6,6)
arma::mat rotate_l2g_L(const arma::mat &, const double &, const double &, const double &);

//To rotate from global to local a stiffness matrix (6,6)
arma::mat rotate_g2l_L(const arma::mat &, const double &, const double &, const double &);

//To rotate from local to global a localisation matrix (6,6)
arma::mat rotate_l2g_A(const arma::mat &, const double &, const double &, const double &);

//To rotate from global to local a localisation matrix (6,6)
arma::mat rotate_g2l_A(const arma::mat &, const double &, const double &, const double &);

} //namespace smart
