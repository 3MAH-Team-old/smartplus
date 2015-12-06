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

///@file elastic_transversely_isotropic.hpp
///@brief User subroutine for transversely isotropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <armadillo>

namespace smart{
    
    ///@brief The elastic UMAT requires 2 constants:
    ///@brief props[0] : Young modulus
    ///@brief props[1] : Poisson ratio
    ///@brief props[2] : CTE
    
    ///@brief No statev is required for thermoelastic constitutive law
    
    void umat_elasticity_trans_iso(const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, const int &, const int &, const bool &, double &);
    
} //namespace smart
