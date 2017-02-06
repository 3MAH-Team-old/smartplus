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

///@file material_characteristics.hpp
///@brief Characteristics of a material, the parent class of:
// - phase_characteristics
// - ellipsoid_characteristics
// - layer_characteristics
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include "../Solver/step.hpp"
#include "../Solver/step_meca.hpp"

namespace smart{
    
//======================================
class aba_step_meca : public step_meca
//======================================
{
private:
    
protected:
    
    public :
    
    std::string name;     //name of the step
    bool nlgeom;
    int type; //Automatic or fixed
        
    aba_step_meca(); 	//default constructor
    aba_step_meca(const std::string &, const bool &, const int &, const int &, const double &, const double &, const double &, const int &, const arma::Col<int>&, const arma::vec&, const arma::mat&, const double&, const int&, const arma::vec&); //Constructor with parameters
    
    aba_step_meca(const aba_step_meca&);	//Copy constructor
    virtual ~aba_step_meca();

    using step::generate;
    virtual void update(const step_meca&, const std::string &, const bool &, const int &);
    
    virtual aba_step_meca& operator = (const aba_step_meca&);
    
    virtual void write(const std::string &, const std::string &inputfile, const bool & = true);
    
    friend  std::ostream& operator << (std::ostream&, const aba_step_meca&);
};
    
} //namespace smart