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

///@file output.hpp
///@brief object that defines the output
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

namespace smart{

//======================================
class solver_output
//======================================
{
private:
    
protected:
    
public :

    //output values
    int o_nb_meca;
    Col<int> o_meca;
    int o_nb_T;
    
	int o_nw_statev;
	Col<int> o_wanted_statev;
	Col<int> o_range_statev;
    
    Col<int> o_type;
    Col<int> o_nfreq;
    vec o_tfreq;
    
    solver_output(); 	//default constructor
    solver_output(const int&);	//Constructor with parameters
    solver_output(const solver_output &);	//Copy constructor
    ~solver_output();
    
    virtual solver_output& operator = (const solver_output&);
    
    friend  ostream& operator << (ostream&, const solver_output&);
};

} //namespace smart
