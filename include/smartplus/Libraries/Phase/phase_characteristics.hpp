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

///@file phase_characteristics.hpp
///@brief Characteristics of a phase, which hereditates from:
///-material characteristics
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include "../Geometry/geometry.hpp"
#include "../Homogenization/phase_multi.hpp"
#include "material_characteristics.hpp"
#include "state_variables.hpp"
#include "../Solver/output.hpp"

namespace smart{

//======================================
class phase_characteristics
//======================================
{
	private:

	protected:

	public :

        int shape_type;
        int sv_type;
        std::shared_ptr<geometry> sptr_shape; //The shape of the phase
        std::shared_ptr<phase_multi> sptr_multi; //The multiscale information of the phase
        std::shared_ptr<material_characteristics> sptr_matprops;
        std::shared_ptr<state_variables> sptr_sv_global;
        std::shared_ptr<state_variables> sptr_sv_local;
        std::shared_ptr<std::ofstream> sptr_out;
    
        std::vector<phase_characteristics> sub_phases;
    
		phase_characteristics(); 	//default constructor
    
        phase_characteristics(const int &, const int &, const std::shared_ptr<geometry> &, const std::shared_ptr<phase_multi> &, const std::shared_ptr<material_characteristics> &, const std::shared_ptr<state_variables> &, const std::shared_ptr<state_variables> &, const std::shared_ptr<std::ofstream> &);

		phase_characteristics(const phase_characteristics&);	//Copy constructor
        virtual ~phase_characteristics();
    
        virtual void construct(const int &, const int &);
        virtual void sub_phases_construct(const int &, const int &, const int &);
        virtual void to_start();
        virtual void set_start();
        virtual void local2global();
        virtual void global2local();
    
		virtual phase_characteristics& operator = (const phase_characteristics&);
    
        virtual void define_output(const std::string &);
        virtual void output(const solver_output &, const int &, const int &, const int &, const int &, const double &, const std::string & = "global");
    
    
        friend std::ostream& operator << (std::ostream&, const phase_characteristics&);
};

} //namespace smart
