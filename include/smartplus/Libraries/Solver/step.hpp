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

///@file step.hpp
///@brief object that defines an step
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "output.hpp"

namespace smart{

//======================================
class step
//======================================
{
private:
    
protected:
    
	public :
    int number;
    int ninc;
    int mode;
    
    double Time;
    arma::vec times;
    double BC_Time;
    
    std::string file; //  It is used for input/output values of the loading path
    
    step(); 	//default constructor
    step(int, int, int);	//Constructor with parameters
    step(const step &);	//Copy constructor
    ~step();
   
    virtual void generate(const double &, const arma::vec &, const arma::vec &, const double &);
    
    virtual step& operator = (const step&);

    //This function serves to output
    virtual void output(std::ostream&, const solver_output &, const int &, const int &, const int&, const arma::vec&);
    
    friend  std::ostream& operator << (std::ostream&, const step&);
};

} //namespace smart
