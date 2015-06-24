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

///@file state_variables.hpp
///@brief State variables of a phase, in a defined coordinate system:
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

namespace smart{

//======================================
class state_variables
//======================================
{
	private:

	protected:

	public :
		
		arma::vec Etot;
		arma::vec DEtot;
		arma::vec sigma;
		arma::mat L;
		arma::mat Lt;
		
		state_variables(); 	//default constructor
		state_variables(arma::vec, arma::vec, arma::vec, arma::mat, arma::mat, arma::mat); //Constructor with parameters
		state_variables(const state_variables &);	//Copy constructor
		~state_variables();
		
		virtual state_variables& operator = (const state_variables&);
		
        virtual state_variables& rotate_l2g(const state_variables&, const double&, const double&, const double&);
        virtual state_variables& rotate_g2l(const state_variables&, const double&, const double&, const double&);
    
        friend std::ostream& operator << (std::ostream&, const state_variables&);
};

} //namespace smart
