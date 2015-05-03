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

///@file general_characteristics.hpp
///@brief general characteristics of a phase, including the fact that it is a coating or coated by another phase
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include "phase_characteristics.hpp"

using namespace std;
using namespace arma;

namespace smart{

//======================================
class general_characteristics : public phase_characteristics
//======================================
{
private:
    
protected:
    
	public :
    
    int coatingof;
    int coatedby;

    general_characteristics(); 	//default constructor
    general_characteristics(int, int, bool=true, double=0.);	//constructor - allocates memory for statev
    general_characteristics(int, int, int, string, double, double, double, double, int, const vec&, int, const vec&, const state_variables&, const state_variables&, const mat&, const mat&); //Constructor with parameters
    general_characteristics(const general_characteristics&);	//Copy constructor
    ~general_characteristics();
    
	virtual general_characteristics& operator = (const general_characteristics&);
    
    friend ostream& operator << (ostream&, const general_characteristics&);
    
};

} //namespace smart
