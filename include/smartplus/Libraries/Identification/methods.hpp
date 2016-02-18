/* This file is part of SMART+ private.
 
 Only part of SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
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