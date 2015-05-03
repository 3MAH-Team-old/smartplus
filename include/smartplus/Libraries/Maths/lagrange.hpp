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

///@file lagrange.cpp
///@brief Function that are utilized in constrained problems
///@version 1.0

#pragma once

#include <math.h>
#include <armadillo>

using namespace arma;

namespace smart{

//This function is used to determine an exponential Lagrange Multiplier (like contact in Abaqus)
double lagrange_exp(const double &, const double &, const double &);

//This function is used to determine the first derivative of an exponential Lagrange Multiplier
double dlagrange_exp(const double &, const double &, const double &);

//This function is used to determine a power-law Lagrange Multiplier for problem such x >= 0
double lagrange_pow_0(const double &, const double &, const double &, const double &, const double &);

//This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x >= 0
double dlagrange_pow_0(const double &, const double &, const double &, const double &, const double &);

//This function is used to determine a power-law Lagrange Multiplier for problem such x <= 1
double lagrange_pow_1(const double &, const double &, const double &, const double &, const double &);

//This function is used to determine the first derivative of a power-law Lagrange Multiplier for problem such x <= 1
double dlagrange_pow_1(const double &, const double &, const double &, const double &, const double &);
    
} //namespace smart

