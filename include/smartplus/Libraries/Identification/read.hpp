/* This file is part of SMART+ private.
 
 Only part of SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file factories
///@brief Generation of the complex objects for identification library
///@author Chemisky
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "parameters.hpp"
#include "constants.hpp"
#include "opti_data.hpp"

namespace smart{

//Generation of the parameters from the parameters file
void read_parameters(const int &, std::vector<parameters> &);
    
void read_constants(const int &, std::vector<constants> &, const int &);
    
void read_data_exp(const int &, std::vector<opti_data> &);

void read_data_weights(const int &, arma::Col<int> &, arma::vec &, std::vector<arma::vec> &, std::vector<opti_data> &, const std::vector<opti_data> &);
    
void read_data_num(const int &, const std::vector<opti_data> &, std::vector<opti_data> &);

//Read the control parameters of the optimization algorithm
void ident_control(int &, int &, int &, int &, int &, int &, int &, int &, int &, double &, double &, double &, double &, double &, double &);

void read_gen(int &, arma::mat &, const int &);
    
} //namespace smart
