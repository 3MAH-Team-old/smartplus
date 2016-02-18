/* This file is part of SMART+ private.
 
 Only part of SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file script.hpp
///@brief Scripts that allows to run identification algorithms based on Smart+ Control functions
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "parameters.hpp"
#include "opti_data.hpp"
#include "generation.hpp"

namespace smart{

//This function will copy the parameters files
void copy_parameters(const std::vector<parameters> &, const std::string &, const std::string &);

//This function will copy the parameters files
void copy_constants(const std::vector<constants> &, const std::string &, const std::string &);
    
//This function will replace the keys by the parameters
void apply_parameters(const std::vector<parameters> &, const std::string &);

//This function will replace the keys by the parameters
void apply_constants(const std::vector<constants> &, const std::string &);
    
    
//Read the control parameters of the optimization algorithm
void launch_solver(const generation &, const int &, std::vector<parameters> &, std::vector<constants> &, const std::string &, const std::string &, const std::string &);

} //namespace smart
