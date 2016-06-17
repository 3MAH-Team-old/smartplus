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

///@file read.hpp
///@brief To read from material.dat and path.dat
///@version 1.0

#pragma once
#include <armadillo>
#include <string>
#include "block.hpp"
#include "output.hpp"

namespace smart{

arma::Col<int> subdiag2vec();

/// Function that fills the matrix Tdsde for mix strain/stress conditions
void Lt_2_K(const arma::mat &, arma::mat &, const arma::Col<int> &, const double &);

/// Function that fills the matrix Tdsde for mix strain/stress conditions
void Lth_2_K(const arma::mat &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::Col<int> &, const int &, const double &);

/// Function that reads the material properties
void read_matprops(std::string &, int &, arma::vec &, int &, double &, double &, double &, double &, double &, const string &);
    
/// Function that reads the output parameters
void read_output(solver_output &, const int &, const int &);

/// Function that checks the coherency between the path and the step increments provided
void check_path_output(const std::vector<block> &, const solver_output &);
    
/// Function that reads the loading path
void read_path(std::vector<block> &, double &, const std::string & = "path.txt");

} //namespace smart
