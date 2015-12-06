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

///@file Self-Consistent.hpp
///@brief User subroutine for non-linear N-phases heterogeneous materials using
///@brief Self-Consistent scheme
///@version 1.0

#pragma once

#include <armadillo>

using namespace std;
using namespace arma;

namespace smart{

///@brief props[0] : Number of phases
///@brief props[1] : Number of the file NPhase[i].dat utilized
///@brief props[2] : Number of integration points in the 1 direction
///@brief props[3] : Number of integration points in the 2 direction

    void umat_SC_N(const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, const int &, const int &, const bool &, double &);

} //namespace smart
