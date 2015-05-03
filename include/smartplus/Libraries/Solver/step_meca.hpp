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

///@file step_meca.hpp
///@brief object that defines a mechanical step
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "step.hpp"

using namespace std;
using namespace arma;

namespace smart{

//======================================
class step_meca : public step
//======================================
{
private:
    
protected:
    
	public :

    Col<int> cBC_meca; //True is for stress (flux), false if for strain (state)
    vec BC_meca;
    mat mecas;
    double BC_T;
    int cBC_T;
    vec Ts;
    
    vec Etot;
    vec DEtot;
    vec sigma;
    double T;
 
    step_meca(); 	//default constructor
    step_meca(int, int, int, const Col<int>&, const vec&, const mat&, const double&, const int&, const vec&, const vec&, const vec&, const vec&, const double&); //Constructor with parameters
    step_meca(const step_meca&);	//Copy constructor
    ~step_meca();
    
    virtual void generate(const double&, const vec&, const vec&, const double&);
    
    virtual step_meca& operator = (const step_meca&);
    
    virtual void output(ostream&, const solver_output&, const int&, const int&, const int&, const vec&);
    
    friend  ostream& operator << (ostream&, const step_meca&);
};

} //namespace smart
