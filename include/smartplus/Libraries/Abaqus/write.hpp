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

///@file write
///@brief To write the necessary Abaqus input files
///@version 1.0

#pragma once
#include <armadillo>
#include <string>
#include "materials.hpp"
#include "../Phase/phase_characteristics.hpp"
#include "section_characteristics.hpp"
#include "steps.hpp"
#include "../Solver/block.hpp"

namespace smart{
    
/// Function that reads the output parameters
void update_sections(section_characteristics &, const int &, const int &, const int &, const std::string & = "data");
void write_sections(section_characteristics &, const std::string & = "data", const std::string & = "Nmat_0.inp");
void update_materials(std::vector<aba_material> &, const int &, const int &, const int &, const std::string & = "data");
void write_materials(std::vector<aba_material> &, const std::string & = "data", const std::string & = "Nmat_0.inp");
void update_steps(std::vector<aba_step_meca> &, const std::vector<block> &, const bool & = false);
void write_steps(std::vector<aba_step_meca> &, const double &, const std::string & = "data", const std::string & = "Nstep_0.inp");

} //namespace smart
