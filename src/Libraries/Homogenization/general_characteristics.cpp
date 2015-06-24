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

///@file general_characteristics.cpp
///@brief general characteristics of a phase, including the fact that it is a coating or coated by another phase
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <smartplus/Libraries/Homogenization/phase_characteristics.hpp>
#include <smartplus/Libraries/Homogenization/general_characteristics.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for general_characteristics===================================

//=====Public methods for general_characteristics============================================


/*!
  \brief default constructor
*/

//-------------------------------------------------------------
general_characteristics::general_characteristics() : phase_characteristics()
//-------------------------------------------------------------
{
    coatingof = 0;
    coatedby = 0;
}

/*!
  \brief Constructor
  \param nprops : size of the table statev 
  \param nstatev : size of the table statev 
  \param init boolean that indicates if the constructor has to initialize the state variables vectors and statev table (default value is true) \n
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
general_characteristics::general_characteristics(int n, int m, bool init, double value) : phase_characteristics(n,m,init,value)
//-------------------------------------------------------------
{
    coatingof = 0;
    coatedby = 0;
}

/*!
  \brief Constructor with parameters
  \param nstatev : size of the table statev 
  \param init boolean that indicates if the constructor has to initialize the state variables vectors and statev table (default value is true) \n
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
general_characteristics::general_characteristics(int mnumber, int mcoatingof, int mcoatedby, string mumat_name, double mconcentration, double mpsi_mat, double mtheta_mat, double mphi_mat, int mnprops, const vec &mprops, int mnstatev, const vec &mstatev, const state_variables &mlocal, const state_variables &mglobal, const mat &mA, const mat &mB) : phase_characteristics(mnumber, mumat_name, mconcentration, mpsi_mat, mtheta_mat, mphi_mat, mnprops, mprops, mnstatev, mstatev, mlocal, mglobal, mA, mB)
//-------------------------------------------------------------
{
	coatingof = mcoatingof;
    coatedby = mcoatedby;
}

/*!
  \brief Copy constructor
  \param s phases_characteristics object to duplicate
*/

//------------------------------------------------------
general_characteristics::general_characteristics(const general_characteristics& gv) : phase_characteristics(gv)
//------------------------------------------------------
{
	coatingof = gv.coatingof;
    coatedby = gv.coatedby;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
general_characteristics::~general_characteristics()
//-------------------------------------
{ }

/*!
  \brief Standard operator = for phases_characteristics
*/

//----------------------------------------------------------------------
general_characteristics& general_characteristics::operator = (const general_characteristics& gv)
//----------------------------------------------------------------------
{
	assert(gv.nprops);

	number = gv.number;
	coatingof = gv.coatingof;
    coatedby = gv.coatedby;
	umat_name = gv.umat_name;

	concentration = gv.concentration;
	psi_mat = gv.psi_mat;
	theta_mat = gv.theta_mat;
	phi_mat = gv.phi_mat;
    
	local = gv.local;
	global = gv.global;
	
	A = gv.A;
	B = gv.B;
	
	nprops = gv.nprops;
	props = gv.props;
	nstatev = gv.nstatev;
	statev = gv.statev;

	return *this;

}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const general_characteristics& gv)
//--------------------------------------------------------------------------
{
	assert(gv.nprops);

	s << "Display state variables\n";
	s << "Number of the phase: " << gv.number << "\n";
	s << "is a coating of the phase: " << gv.coatingof << "\n";
	s << "is coated by the phase: " << gv.coatedby << "\n";
	s << "Name of the material = " << gv.umat_name << "\n";
	s << "concentration: " << gv.concentration << "\n";
	s << "local material orientation: psi = " << gv.psi_mat << "\t theta = " << gv.theta_mat << "\t phi = " << gv.phi_mat << "\n";
	s << "nprops: \n" << gv.nprops << "\n";
	s << "props: \n";
    s << gv.props.t();
    s << "\n";

	s << "local state variables: \n" << gv.local << "\n";
	s << "global state variables: \n" << gv.global << "\n";
	
	s << "A: \n" << gv.A << "\n";
	s << "B: \n" << gv.B << "\n";	
	s << "nstatev: \n" << gv.nstatev << "\n";
	if (gv.nstatev) {
		s << "statev: \n";
		s << gv.statev.t();
		s << "\n";
	}

	s << "\n\n";
    
	return s;
}

} //namespace smart
