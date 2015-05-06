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

///@file phase_characteristics.hpp
///@brief Characteristics of a phase, the parent class of:
// - general_characteristics
// - ellipsoid_characteristics
// - layer_characteristics
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include "state_variables.hpp"

namespace smart{

//======================================
class phase_characteristics
//======================================
{
	private:

	protected:

	public :

		int number;
        std::string umat_name;
		double concentration;
        double psi_mat;
        double theta_mat;
        double phi_mat;
    
		int nprops;
		arma::vec props;
		int nstatev;
		arma::vec statev;
		
		state_variables local;
		state_variables global;
		
		arma::mat A;	//Concentration tensor (strain)
		arma::mat B;	//Concentration tensor (stress)
		
		phase_characteristics(); 	//default constructor
		phase_characteristics(int, int, bool=true, double=0.);	//constructor - allocates memory for statev
        phase_characteristics(int, std::string, double, double, double, double, int, const arma::vec&, int, const arma::vec&, const state_variables&, const state_variables&, const arma::mat&, const arma::mat&); //Constructor with parameters
		phase_characteristics(const phase_characteristics&);	//Copy constructor
        ~phase_characteristics();
		
		void resize(int, int, bool=true, double=0.);
		int dimprops () const {return nprops;}       // returns the number of props, nprops
		int dimstatev () const {return nstatev;}       // returns the number of statev, nstatev

        void local2global();
        void global2local();
    
		virtual phase_characteristics& operator = (const phase_characteristics&);
		
        friend std::ostream& operator << (std::ostream&, const phase_characteristics&);
};

} //namespace smart
