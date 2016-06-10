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

///@file individual.hpp
///@brief individual for genetic algorithm (among others)
///@author Chemisky & Despringre
///@version 1.0

#pragma once
#include <armadillo>
#include "generation.hpp"

namespace smart{
    
//Genetic method
void genetic(generation &, generation &, int &, const double &, const double &, const std::vector<parameters> &);

///Genrun creation
void to_run(generation &, generation &, generation &, const double &, const std::vector<parameters> &);

} //namespace smart