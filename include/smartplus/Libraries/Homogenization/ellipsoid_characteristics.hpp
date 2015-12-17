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

///@file ellipsoid_characteristics.hpp
///@brief characteristics of an ellipsoid
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include "phase_characteristics.hpp"

namespace smart{

//======================================
class ellipsoid_characteristics : public phase_characteristics
//======================================
{
private:
    
protected:
    
	public :
    
    int coatingof;
    int coatedby;
    
    double a1;
    double a2;
    double a3;
    
    double psi_geom;
    double theta_geom;
    double phi_geom;
    
    arma::mat S;
    arma::mat T;
    
    ellipsoid_characteristics(); 	//default constructor
    ellipsoid_characteristics(int, int, bool=true, double=0.);	//constructor - allocates memory for statev
    ellipsoid_characteristics(int, int, int, std::string, double, double, double, double, double, double, double, double, double, double, int, const arma::vec&, int, const arma::vec&, const state_variables&, const state_variables&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&); //Constructor with parameters
    ellipsoid_characteristics(const ellipsoid_characteristics&);	//Copy constructor
    ~ellipsoid_characteristics();
    
    void fillS(const arma::mat&, const arma::vec&, const arma::vec&, const arma::vec&, const arma::vec&, const int&, const int&); //need the L_global of the matrix
    void fillT(const arma::mat&, const arma::vec&, const arma::vec&, const arma::vec&, const arma::vec&, const int&, const int&); //need the L_global of the matrix
    
	virtual ellipsoid_characteristics& operator = (const ellipsoid_characteristics&);
    
    friend std::ostream& operator << (std::ostream&, const ellipsoid_characteristics&);
    
};

} //namespace smart
