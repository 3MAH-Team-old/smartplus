/* This file is part of SMART+ private.
 
 Only part of SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file damage_LLD_m.hpp
///@brief Modified Ladeveze-Le Dantec model that accounts for uniaxial damage and anisotopic plasticity
///@brief Damage evolution is considered as a function of time
///@author Y. Chemisky

#pragma once

#include <iostream>
#include <armadillo>

namespace smart {
    
    void umat_damage_LLD_0(const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, arma::mat &, arma::vec &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, const int &, const int &, const bool &, const int &, double &);
                          
} //namespace smart
