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

///@file layer_characteristics.hpp
///@brief Characteristics of a layer, similar to a phase in this version
///@version 1.0
//
#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include "phase_characteristics.hpp"

using namespace std;
using namespace arma;

namespace smart{

//======================================
class layer_characteristics : public phase_characteristics
//======================================
{
	private:

	protected:
    
	public :

        layer_characteristics(); 	//default constructor
        layer_characteristics(int, int, bool=true, double=0.);	//constructor - allocates memory for statev
        layer_characteristics(int, string, double, double, double, double, int, const vec&, int, const vec&, const state_variables&, const state_variables&, const mat&, const mat&); //Constructor with parameters
        layer_characteristics(const layer_characteristics&);	//Copy constructor
        ~layer_characteristics();
    
    	virtual layer_characteristics& operator = (const layer_characteristics&);
    
        friend ostream& operator << (ostream&, const layer_characteristics&);
    
};

} //namespace smart
