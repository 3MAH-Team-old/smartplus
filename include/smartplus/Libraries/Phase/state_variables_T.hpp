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

///@file state_variables.hpp
///@brief State variables of a phase, in a defined coordinate system:
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "state_variables.hpp"

namespace smart{

//======================================
    class state_variables_T : public state_variables
//======================================
{
	private:

	protected:

	public :
		
        arma::mat dSdE;
        arma::mat dSdEt;
        arma::mat dSdT;
        double Q;
        double rpl;
        arma::mat drpldE;
        arma::mat drpldT;

		state_variables_T(); 	//default constructor
		state_variables_T(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const double &, const double &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const double &, const double &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &); //Constructor with parameters
		state_variables_T(const state_variables_T &);	//Copy constructor
		virtual ~state_variables_T();
		
		virtual state_variables_T& operator = (const state_variables_T&);
		
        using state_variables::update;
        virtual void update(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const double &, const double &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const double &, const double &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &);
    
        using state_variables::rotate_l2g;
        virtual state_variables_T& rotate_l2g(const state_variables_T&, const double&, const double&, const double&);
        using state_variables::rotate_g2l
    ;
        virtual state_variables_T& rotate_g2l(const state_variables_T&, const double&, const double&, const double&);
    
        friend std::ostream& operator << (std::ostream&, const state_variables_T&);
};

} //namespace smart
