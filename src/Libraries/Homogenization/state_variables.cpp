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

///@file state_variables.cpp
///@brief State variables of a phase, in a defined coordinate system:
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <smartplus/Libraries/Homogenization/state_variables.hpp>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Maths/rotation.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for state_variables===================================

//=====Public methods for state_variables============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
state_variables::state_variables() : Etot(6), DEtot(6), sigma(6), L(6,6), Lt(6,6)
//-------------------------------------------------------------
{

	Etot = zeros(6);
	DEtot = zeros(6);
	sigma = zeros(6);
	L = zeros(6,6);
	Lt = zeros(6,6);
}

/*!
  \brief Constructor with parameters
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
state_variables::state_variables(vec mEtot, vec mDEtot, vec msigma, mat mL, mat mLt) : Etot(6), DEtot(6), sigma(6), L(6,6), Lt(6,6)
//-------------------------------------------------------------
{	
	assert (mEtot.size() == 6);
	assert (mDEtot.size() == 6);
	assert (msigma.size() == 6);
	assert (mL.n_rows == 6);
	assert (mL.n_cols == 6);
	assert (mLt.n_rows == 6);
	assert (mLt.n_cols == 6);	
	
	Etot = mEtot;
	DEtot = mDEtot;
	sigma = msigma;
	L = mL;	
	Lt = mLt;
}

/*!
  \brief Copy constructor
  \param s state_variables object to duplicate
*/

//------------------------------------------------------
state_variables::state_variables(const state_variables& sv) : Etot(6), DEtot(6), sigma(6), L(6,6), Lt(6,6)
//------------------------------------------------------
{
	Etot = sv.Etot;
	DEtot = sv.DEtot;
	sigma = sv.sigma;
	L = sv.L;
	Lt = sv.Lt;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
state_variables::~state_variables()
//-------------------------------------
{

}

/*!
  \brief Standard operator = for state_variables
*/

//----------------------------------------------------------------------
state_variables& state_variables::operator = (const state_variables& sv)
//----------------------------------------------------------------------
{
	Etot = sv.Etot;
	DEtot = sv.DEtot;
	sigma = sv.sigma;
	L = sv.L;
	Lt = sv.Lt;

	return *this;

}

//----------------------------------------------------------------------
state_variables& state_variables::rotate_l2g(const state_variables& sv, const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{
    
	Etot = sv.Etot;
	DEtot = sv.DEtot;
	sigma = sv.sigma;
	L = sv.L;
	Lt = sv.Lt;
    
  	if(fabs(phi) > iota) {
		Etot = rotate_strain(Etot, -phi, axis_phi);
		DEtot = rotate_strain(DEtot, -phi, axis_phi);
		sigma = rotate_stress(sigma, -phi, axis_phi);
		Lt = rotateL(Lt, -phi, axis_phi);
		L = rotateL(L, -phi, axis_phi);
	}
  	if(fabs(theta) > iota) {
		Etot = rotate_strain(Etot, -theta, axis_theta);
		DEtot = rotate_strain(DEtot, -theta, axis_theta);
		sigma = rotate_stress(sigma, -theta, axis_theta);				
		L = rotateL(L, -theta, axis_theta);		
		Lt = rotateL(Lt, -theta, axis_theta);
	}
	if(fabs(psi) > iota) {
		Etot = rotate_strain(Etot, -psi, axis_psi);
		DEtot = rotate_strain(DEtot, -psi, axis_psi);
		sigma = rotate_stress(sigma, -psi, axis_psi);
		L = rotateL(L, -psi, axis_psi);
		Lt = rotateL(Lt, -psi, axis_psi);
	}
    
	return *this;
}

//----------------------------------------------------------------------
state_variables& state_variables::rotate_g2l(const state_variables& sv, const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{

	Etot = sv.Etot;
	DEtot = sv.DEtot;
	sigma = sv.sigma;
	L = sv.L;
	Lt = sv.Lt;
    
  	if(fabs(psi) > iota) {
		Etot = rotate_strain(Etot, psi, axis_psi);
		DEtot = rotate_strain(DEtot, psi, axis_psi);
		sigma = rotate_stress(sigma, psi, axis_psi);
		L = rotateL(L, psi, axis_psi);
		Lt = rotateL(Lt, psi, axis_psi);
	}			
	if(fabs(theta) > iota) {
		Etot = rotate_strain(Etot, theta, axis_theta);
		DEtot = rotate_strain(DEtot, theta, axis_theta);
		sigma = rotate_stress(sigma, theta, axis_theta);
		L = rotateL(L, theta, axis_theta);	
		Lt = rotateL(Lt, theta, axis_theta);
	}
	if(fabs(phi) > iota) {
		Etot = rotate_strain(Etot, phi, axis_phi);
		DEtot = rotate_strain(DEtot, phi, axis_phi);
		sigma = rotate_stress(sigma, phi, axis_phi);
		L = rotateL(L, phi, axis_phi);
		Lt = rotateL(Lt, phi, axis_phi);
    }
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const state_variables& sv)
//--------------------------------------------------------------------------
{
	s << "Etot: \n" << sv.Etot << "\n";
	s << "DEtot: \n" << sv.DEtot << "\n";
	s << "sigma: \n" << sv.sigma << "\n";
	s << "L: \n" << sv.L << "\n";
	s << "Lt: \n" << sv.Lt << "\n";
	s << "\n";

	return s;
}

} //namespace smart
