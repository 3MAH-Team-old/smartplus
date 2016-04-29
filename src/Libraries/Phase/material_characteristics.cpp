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

///@file material_characteristics.cpp
///@brief Characteristics of a material,
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <smartplus/Libraries/Phase/material_characteristics.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for material_characteristics===================================

//=====Public methods for material_characteristics============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
material_characteristics::material_characteristics()
//-------------------------------------------------------------
{
	number=-1;
    save = 0;
    
	psi_mat=0.;
	theta_mat=0.;
	phi_mat=0.;
	
	nprops=0;
    
    rho = 0.;
    c_p = 0.;
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
material_characteristics::material_characteristics(const int &n, const bool &init, const double &value)
//-------------------------------------------------------------
{

	assert(n>0);
  
	number = -1;
    save = 0;
    
	psi_mat=0.;
	theta_mat=0.;
	phi_mat=0.;
	
    nprops = n;
    if (init) {
        props = value*ones(n);
    }
    else{
        props = zeros(n);
    }

    rho = 0.;
    c_p = 0.;
}

/*!
  \brief Constructor with parameters
  \param nstatev : size of the table statev 
  \param init boolean that indicates if the constructor has to initialize the state variables vectors and statev table (default value is true) \n
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
material_characteristics::material_characteristics(const int &mnumber, const string &mumat_name, const int &msave, const double &mpsi_mat, const double &mtheta_mat, const double &mphi_mat, const int &mnprops, const vec &mprops, const double &mrho, const double &mc_p)
//-------------------------------------------------------------
{	
	assert(mnprops);
	
	number = mnumber;
	umat_name = mumat_name;
    save = msave;
    
    psi_mat = mpsi_mat;
	theta_mat = mtheta_mat;
	phi_mat = mphi_mat;
    
	nprops = mnprops;
	props = mprops;
    
    rho = mrho;
    c_p = mc_p;
}

/*!
  \brief Copy constructor
  \param s phase_characteristics object to duplicate
*/

//------------------------------------------------------
material_characteristics::material_characteristics(const material_characteristics& sv)
//------------------------------------------------------
{
	number = sv.number;
	umat_name = sv.umat_name;
    save = sv.save;

    psi_mat = sv.psi_mat;
	theta_mat = sv.theta_mat;
	phi_mat = sv.phi_mat;
    
	nprops = sv.nprops;
	props = sv.props;
    
    rho = sv.rho;
    c_p = sv.c_p;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
material_characteristics::~material_characteristics() {}
//-------------------------------------

//-------------------------------------------------------------
void material_characteristics::resize()
//-------------------------------------------------------------
{
    assert(nprops > 0);
    props = zeros(nprops);
}
    
//-------------------------------------------------------------
void material_characteristics::resize(const int &n, const bool &init, const double &value)
//-------------------------------------------------------------
{
    assert(n > 0);
    nprops = n;
    if (init) {
        props = value*ones(n);
        
    }
    else{
        props = zeros(n);
    }
}

/*!
  \brief Standard operator = for phase_characteristics
*/

//-------------------------------------------------------------
void material_characteristics::update(const int &mnumber, const string &mumat_name, const int &msave, const double &mpsi_mat, const double &mtheta_mat, const double &mphi_mat, const int &mnprops, const vec &mprops, const double &mrho, const double &mc_p)
//-------------------------------------------------------------
{
    assert(mnprops);
    
    number = mnumber;
    umat_name = mumat_name;
    save = msave;
    
    psi_mat = mpsi_mat;
    theta_mat = mtheta_mat;
    phi_mat = mphi_mat;
    
    nprops = mnprops;
    props = mprops;
    
    rho = mrho;
    c_p = mc_p;
}
    
//----------------------------------------------------------------------
material_characteristics& material_characteristics::operator = (const material_characteristics& sv)
//----------------------------------------------------------------------
{
	assert(sv.nprops);

	number = sv.number;
	umat_name = sv.umat_name;
    save = sv.save;
    
    psi_mat = sv.psi_mat;
	theta_mat = sv.theta_mat;
	phi_mat = sv.phi_mat;
		
	nprops = sv.nprops;
	props = sv.props;
    
    rho = sv.rho;
    c_p = sv.c_p;
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const material_characteristics& sv)
//--------------------------------------------------------------------------
{
	assert(sv.nprops);

	s << "Display state variables\n";
	s << "Number of the phase: " << sv.number << "\n";
	s << "Name of the material = " << sv.umat_name << "\n";
	s << "local material orientation: psi = " << sv.psi_mat << "\t theta = " << sv.theta_mat << "\t phi = " << sv.phi_mat << "\n";
	s << "nprops: \n" << sv.nprops << "\n";
	s << "props: \n";
    s << sv.props.t();
	s << "\n";	
        
    s << "rho: \n" << sv.rho << "\n";
    s << "c_p: \n" << sv.c_p << "\n";

	s << "\n\n";

	return s;
}

} //namespace smart
