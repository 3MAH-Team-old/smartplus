///@file solve.hpp
///@brief random number generators
///@author Chemisky & Anagnostou
///@version 1.0
///@date 10/23/2014
#pragma once

#include <armadillo>

namespace smart{

arma::vec quadratic(const double &, const double &, const double &);

arma::cx_vec cx_quadratic(const double &, const double &, const double &);

arma::cx_vec cx_quadratic(const cx_double &, const cx_double &, const cx_double &);

} //namespace smart    
