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

using namespace std;
using namespace arma;

namespace smart{

//======================================
class phase_characteristics
//======================================
{
	private:

	protected:

	public :

		int number;
		string umat_name;
		double concentration;
        double psi_mat;
        double theta_mat;
        double phi_mat;
    
		int nprops;
		vec props;
		int nstatev;
		vec statev;
		
		state_variables local;
		state_variables global;
		
		mat A;	//Concentration tensor (strain)
		mat B;	//Concentration tensor (stress)	
		
		phase_characteristics(); 	//default constructor
		phase_characteristics(int, int, bool=true, double=0.);	//constructor - allocates memory for statev
		phase_characteristics(int, string, double, double, double, double, int, const vec&, int, const vec&, const state_variables&, const state_variables&, const mat&, const mat&); //Constructor with parameters
		phase_characteristics(const phase_characteristics&);	//Copy constructor
        ~phase_characteristics();
		
		void resize(int, int, bool=true, double=0.);
		int dimprops () const {return nprops;}       // returns the number of props, nprops
		int dimstatev () const {return nstatev;}       // returns the number of statev, nstatev

        void local2global();
        void global2local();
    
		virtual phase_characteristics& operator = (const phase_characteristics&);
		
		friend ostream& operator << (ostream&, const phase_characteristics&);
};

} //namespace smart
