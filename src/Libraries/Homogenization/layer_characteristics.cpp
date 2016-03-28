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

///@file layer_characteristics.cpp
///@brief Characteristics of a layer, similar to a phase in this version
///@version 1.0
//

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <smartplus/Libraries/Homogenization/phase_characteristics.hpp>
#include <smartplus/Libraries/Homogenization/layer_characteristics.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for general_characteristics===================================

//=====Public methods for general_characteristics============================================

//-------------------------------------------------------------
layer_characteristics::layer_characteristics() : phase_characteristics()
//-------------------------------------------------------------
{ }

/*!
 \brief Constructor
 \param nprops : size of the table statev
 \param nstatev : size of the table statev
 \param init boolean that indicates if the constructor has to initialize the state variables vectors and statev table (default value is true) \n
 \n\n
 \f$ \textbf{Examples :} \f$ \n
 */

//-------------------------------------------------------------
layer_characteristics::layer_characteristics(int n, int m, bool init, double value) : phase_characteristics(n,m,init,value)
//-------------------------------------------------------------
{ }

/*!
 \brief Constructor with parameters
 \param nstatev : size of the table statev
 \param init boolean that indicates if the constructor has to initialize the state variables vectors and statev table (default value is true) \n
 \n\n
 \f$ \textbf{Examples :} \f$ \n
 */

//-------------------------------------------------------------
layer_characteristics::layer_characteristics(int mnumber, string mumat_name, double mconcentration, double mpsi_mat, double mtheta_mat, double mphi_mat, int mnprops, const vec &mprops, int mnstatev, const vec &mstatev, const state_variables &mlocal, const state_variables &mglobal, const mat &mA, const mat &mB) : phase_characteristics(mnumber, mumat_name, mconcentration, mpsi_mat, mtheta_mat, mphi_mat, mnprops, mprops, mnstatev, mstatev, mlocal, mglobal, mA, mB)
//-------------------------------------------------------------
{ }

/*!
 \brief Copy constructor
 \param s phases_characteristics object to duplicate
 */

//------------------------------------------------------
layer_characteristics::layer_characteristics(const layer_characteristics& lv) : phase_characteristics(lv)
//------------------------------------------------------
{ }

/*!
 \brief Destructor
 
 Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
 */

//-------------------------------------
layer_characteristics::~layer_characteristics()
//-------------------------------------
{ }

//----------------------------------------------------------------------
layer_characteristics& layer_characteristics::operator = (const layer_characteristics& lv)
//----------------------------------------------------------------------
{
	assert(lv.nprops);
    
	number = lv.number;
	umat_name = lv.umat_name;
	concentration = lv.concentration;
    
    psi_mat = lv.psi_mat;
	theta_mat = lv.theta_mat;
	phi_mat = lv.phi_mat;
    
	local = lv.local;
	global = lv.global;
	
	A = lv.A;
	B = lv.B;
	
	nprops = lv.nprops;
	props = lv.props;
	nstatev = lv.nstatev;
	statev = lv.statev;
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const layer_characteristics& lv)
//--------------------------------------------------------------------------
{
	assert(lv.nprops);
    
	s << "Display state variables\n";
	s << "Number of the phase: " << lv.number << "\n";
	s << "Name of the material = " << lv.umat_name << "\n";
	s << "concentration: " << lv.concentration << "\n";
	s << "local material orientation: psi = " << lv.psi_mat << "\t theta = " << lv.theta_mat << "\t phi = " << lv.phi_mat << "\n";
	s << "nprops: \n" << lv.nprops << "\n";
	s << "props: \n";
    s << lv.props.t();
	s << "\n";
    
	s << "local state variables: \n" << lv.local << "\n";
	s << "global state variables: \n" << lv.global << "\n";
	
	s << "A: \n" << lv.A << "\n";
	s << "B: \n" << lv.B << "\n";
    
	s << "nstatev: \n" << lv.nstatev << "\n";
	if (lv.nstatev) {
		s << "statev: \n";
		s << lv.statev.t();
		s << "\n";
	}
    
	s << "\n\n";
    
	return s;
}

} //namespace smart
