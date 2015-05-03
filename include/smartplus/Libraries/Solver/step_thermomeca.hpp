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

///@file step_thermomeca.hpp
///@brief object that defines a thermomechanical step
///@version 1.0

#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>

#include <armadillo>
#include "step.hpp"

using namespace std;
using namespace arma;

namespace smart{

//======================================
class step_thermomeca : public step
//======================================
{
private:
    
protected:
    
	public :

    Col<int> cBC_meca; //True is for stress (flux), false if for strain (state)
    vec BC_meca;
    mat mecas;
    double BC_T;
    int cBC_T;         //True (1) is for a heat flux entering in a material point, 0 is for fixed temperature
    vec Ts;
    
    vec Etot;
    vec DEtot;
    vec sigma;
    double T;
    double Q;
    
    step_thermomeca(); 	//default constructor
    step_thermomeca(int, int, int, const Col<int>&, const vec&, const mat&, const double&, const int&, const vec&, const vec&, const vec&, const vec&, const double&, const double&); //Constructor with parameters
    step_thermomeca(const step_thermomeca&);	//Copy constructor
    ~step_thermomeca();
    
    virtual void generate(const double&, const vec&, const vec&, const double&);
    
    virtual step_thermomeca& operator = (const step_thermomeca&);
    
    virtual void output(ostream&, const solver_output&, const int&, const int&, const int&, const vec&);
    
    friend  ostream& operator << (ostream&, const step_thermomeca&);
};

} //namespace smart
