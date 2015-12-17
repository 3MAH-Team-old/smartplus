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

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <smartplus/Libraries/Homogenization/phase_characteristics.hpp>
#include <smartplus/Libraries/Homogenization/state_variables.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for phase_characteristics===================================

//=====Public methods for phase_characteristics============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
phase_characteristics::phase_characteristics() : A(6,6), B(6,6)
//-------------------------------------------------------------
{
	number=0;
	concentration=0.;
    
	psi_mat=0.;
	theta_mat=0.;
	phi_mat=0.;
	
	nprops=0;
	nstatev=0;
	
	A = zeros(6,6);
	B = zeros(6,6);
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
phase_characteristics::phase_characteristics(int n, int m, bool init, double value) : A(6,6), B(6,6)
//-------------------------------------------------------------
{

	assert(n>0);
	assert(m>=0);  
  
	number = 0;
	concentration = 0.;
    
	psi_mat=0.;
	theta_mat=0.;
	phi_mat=0.;
	
    nprops = n;
    if (init) {
        props = value*ones(n);
        
    }
    else{
        props = zeros(n);
        statev = zeros(m);
    }
    
    nstatev = m;
	if (m>0) {
		if (init)
            statev = value*ones(m);
        else
            statev = zeros(m);
	}
}

/*!
  \brief Constructor with parameters
  \param nstatev : size of the table statev 
  \param init boolean that indicates if the constructor has to initialize the state variables vectors and statev table (default value is true) \n
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
phase_characteristics::phase_characteristics(int mnumber, string mumat_name, double mconcentration, double mpsi_mat, double mtheta_mat, double mphi_mat, int mnprops, const vec &mprops, int mnstatev, const vec &mstatev, const state_variables &mlocal, const state_variables &mglobal, const mat &mA, const mat &mB) : A(6,6), B(6,6)
//-------------------------------------------------------------
{	
	assert(mnprops);

	assert (mA.n_rows == 6);
	assert (mA.n_cols == 6);	
	assert (mB.n_rows == 6);
	assert (mB.n_cols == 6);
	
	number = mnumber;
	umat_name = mumat_name;
	concentration = mconcentration;
    
    psi_mat = mpsi_mat;
	theta_mat = mtheta_mat;
	phi_mat = mphi_mat;
	
    local = mlocal;
    global = mglobal;
    
	A = mA;
	B = mB;		
	
	nprops = mnprops;
	props = mprops;
	nstatev = mnstatev;
	statev = mstatev;
}

/*!
  \brief Copy constructor
  \param s phase_characteristics object to duplicate
*/

//------------------------------------------------------
phase_characteristics::phase_characteristics(const phase_characteristics& sv) : A(6,6), B(6,6)
//------------------------------------------------------
{
	number = sv.number;
	umat_name = sv.umat_name;
	concentration = sv.concentration;

    psi_mat = sv.psi_mat;
	theta_mat = sv.theta_mat;
	phi_mat = sv.phi_mat;
    
	local = sv.local;
	global=sv.global;
	
	A = sv.A;
	B = sv.B;
	
	nprops = sv.nprops;
	props = sv.props;
	nstatev = sv.nstatev;
	statev = sv.statev;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
phase_characteristics::~phase_characteristics() {}
//-------------------------------------

//-------------------------------------------------------------
void phase_characteristics::resize(int n, int m, bool init, double value)
//-------------------------------------------------------------
{
	
    nprops = n;
    if (init) {
        props = value*ones(n);
        
    }
    else{
        props = zeros(n);
        statev = zeros(m);
    }
    
    nstatev = m;
	if (m>0) {
		if (init)
            statev = value*ones(m);
        else
            statev = zeros(m);
	}

}

/*!
  \brief Standard operator = for phase_characteristics
*/

//----------------------------------------------------------------------
phase_characteristics& phase_characteristics::operator = (const phase_characteristics& sv)
//----------------------------------------------------------------------
{
	assert(sv.nprops);

	number = sv.number;
	umat_name = sv.umat_name;
	concentration = sv.concentration;
    
    psi_mat = sv.psi_mat;
	theta_mat = sv.theta_mat;
	phi_mat = sv.phi_mat;

	local = sv.local;
	global = sv.global;
	
	A = sv.A;
	B = sv.B;		
	
	nprops = sv.nprops;
	props = sv.props;
	nstatev = sv.nstatev;
	statev = sv.statev;
    
	return *this;
}

//-------------------------------------------------------------
void phase_characteristics::local2global()
//-------------------------------------------------------------
{
	global.rotate_l2g(local, psi_mat, theta_mat, phi_mat);
}

//-------------------------------------------------------------
void phase_characteristics::global2local()
//-------------------------------------------------------------
{
	local.rotate_g2l(global, psi_mat, theta_mat, phi_mat);
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const phase_characteristics& sv)
//--------------------------------------------------------------------------
{
	assert(sv.nprops);

	s << "Display state variables\n";
	s << "Number of the phase: " << sv.number << "\n";
	s << "Name of the material = " << sv.umat_name << "\n";
	s << "concentration: " << sv.concentration << "\n";
	s << "local material orientation: psi = " << sv.psi_mat << "\t theta = " << sv.theta_mat << "\t phi = " << sv.phi_mat << "\n";
	s << "nprops: \n" << sv.nprops << "\n";
	s << "props: \n";
    s << sv.props.t();
	s << "\n";	

	s << "local state variables: \n" << sv.local << "\n";
	s << "global state variables: \n" << sv.global << "\n";
	
	s << "A: \n" << sv.A << "\n";
	s << "B: \n" << sv.B << "\n";
    
	s << "nstatev: \n" << sv.nstatev << "\n";
	if (sv.nstatev) {
		s << "statev: \n";
		s << sv.statev.t();
		s << "\n";
	}

	s << "\n\n";

	return s;
}

} //namespace smart
