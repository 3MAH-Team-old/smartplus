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

using namespace std;
using namespace arma;

namespace smart{

//======================================
class state_variables
//======================================
{
	private:

	protected:

	public :
		
		vec Etot;
		vec DEtot;
		vec sigma;
		mat L;
		mat Lt;
		
		state_variables(); 	//default constructor
		state_variables(vec, vec, vec, mat, mat, mat); //Constructor with parameters
		state_variables(const state_variables &);	//Copy constructor
		~state_variables();
		
		virtual state_variables& operator = (const state_variables&);
		
        virtual state_variables& rotate_l2g(const state_variables&, const double&, const double&, const double&);
        virtual state_variables& rotate_g2l(const state_variables&, const double&, const double&, const double&);
    
		friend ostream& operator << (ostream&, const state_variables&);
};

} //namespace smart
