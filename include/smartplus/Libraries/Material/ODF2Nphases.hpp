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

///@file ODF2Nphases.hpp
///@brief ODF2Nphases discretization of ODFs
///@version 1.0

#pragma once

#include <iostream>
#include <string.h>
#include <armadillo>

namespace smart{

//This function computes the ODF of the selected angle, according to different methods (Lorentzian, Pearson...)
double ODF(const double&, const int&, const arma::vec&, const bool&, const double& = 0.);

//Writes the Nphases.dat file for multiphase modeling, according to specific ODFs
void ODF2Nphases(const arma::Col<int> &, const arma::Col<int> &, const arma::Col<int> &, const std::vector<std::string> &, const arma::mat &, const bool& = false, const double& = 0.);

} //namespace smart
