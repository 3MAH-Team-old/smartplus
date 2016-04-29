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

///@file state_variables_T.cpp
///@brief State variables of a phase, in a defined coordinate system:
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <smartplus/Libraries/Phase/state_variables.hpp>
#include <smartplus/Libraries/Phase/state_variables_T.hpp>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Maths/rotation.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for state_variables_T===================================

//=====Public methods for state_variables_T============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
state_variables_T::state_variables_T() : state_variables(), dSdE(6,6), dSdEt(6,6), dSdT(1,6), drpldE(1,6), drpldT(1,1)
//-------------------------------------------------------------
{
    rpl = 0.;
    dSdE = zeros(6,6);
    dSdEt = zeros(6,6);
    dSdT = zeros(1,6);
    drpldE = zeros(1,6);
    drpldT = zeros(1,1);
}

/*!
  \brief Constructor with parameters
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
state_variables_T::state_variables_T(const vec &mEtot, const vec &mDEtot, const vec &msigma, const vec &msigma_start, const double &mT, const double &mDT, const double &msse, const double &mspd, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const double &mQ, const double &mrpl, const mat &mdSdE, const mat &mdSdEt, const mat &mdSdT, const mat &mdrpldE, const mat &mdrpldT) : state_variables(mEtot, mDEtot, msigma, msigma_start, mT, mDT, msse, mspd, mnstatev, mstatev, mstatev_start), dSdE(6,6), dSdEt(6,6), dSdT(1,6), drpldE(1,6), drpldT(1,1)
//-------------------------------------------------------------
{	
    
	assert (mdSdE.n_rows == 6);
	assert (mdSdE.n_cols == 6);

    assert (mdSdEt.n_rows == 6);
	assert (mdSdEt.n_cols == 6);
    
    assert (mdSdT.n_rows == 1);
    assert (mdSdT.n_cols == 6);
    
    assert (mdrpldE.n_rows == 1);
    assert (mdrpldE.n_cols == 6);
    
    assert (mdrpldT.n_rows == 1);
    assert (mdrpldT.n_cols == 1);

    Q = mQ;
    rpl = mrpl;
    dSdE = mdSdE;
	dSdEt = mdSdEt;
    dSdT = mdSdT;
    drpldE = mdrpldE;
    drpldT = mdrpldT;
}

/*!
  \brief Copy constructor
  \param s state_variables_T object to duplicate
*/

//------------------------------------------------------
state_variables_T::state_variables_T(const state_variables_T& sv) : state_variables(sv), dSdE(6,6), dSdEt(6,6), dSdT(1,6), drpldE(1,6), drpldT(1,1)
//------------------------------------------------------
{
    Q = sv.Q;
    rpl = sv.rpl;
    dSdE = sv.dSdE;
	dSdEt = sv.dSdEt;
    dSdT = sv.dSdT;
    drpldE = sv.drpldE;
    drpldT = sv.drpldT;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
state_variables_T::~state_variables_T()
//-------------------------------------
{

}

/*!
  \brief Standard operator = for state_variables
*/

//----------------------------------------------------------------------
state_variables_T& state_variables_T::operator = (const state_variables_T& sv)
//----------------------------------------------------------------------
{
    Etot = sv.Etot;
    DEtot = sv.DEtot;
    sigma = sv.sigma;
    sigma_start = sv.sigma_start;
    dSdE = sv.dSdE;
    dSdEt = sv.dSdEt;
    dSdT = sv.dSdT;
    drpldE = sv.drpldE;
    drpldT = sv.drpldT;
    sse = sv.sse;
    spd = sv.spd;
    Q = sv.Q;
    rpl = sv.rpl;
    T = sv.T;
    DT = sv.DT;

	return *this;
}

//-------------------------------------------------------------
void state_variables_T::update(const vec &mEtot, const vec &mDEtot, const vec &msigma, const vec &msigma_start, const double &mT, const double &mDT, const double &msse, const double &mspd, const int &mnstatev, const vec &mstatev, const vec &mstatev_start, const double &mQ, const double &mrpl, const mat &mdSdE, const mat &mdSdEt, const mat &mdSdT, const mat &mdrpldE, const mat &mdrpldT)
//-------------------------------------------------------------
{
    state_variables::update(mEtot, mDEtot, msigma, msigma_start, mT, mDT, msse, mspd, mnstatev, mstatev, mstatev_start);
    
    assert (mdSdE.n_rows == 6);
    assert (mdSdE.n_cols == 6);
    
    assert (mdSdEt.n_rows == 6);
    assert (mdSdEt.n_cols == 6);
    
    assert (mdSdT.n_rows == 1);
    assert (mdSdT.n_cols == 6);
    
    assert (mdrpldE.n_rows == 1);
    assert (mdrpldE.n_cols == 6);
    
    assert (mdrpldT.n_rows == 1);
    assert (mdrpldT.n_cols == 1);
    
    Q = mQ;
    rpl = mrpl;
    dSdE = mdSdE;
    dSdEt = mdSdEt;
    dSdT = mdSdT;
    drpldE = mdrpldE;
    drpldT = mdrpldT;
}
    
//----------------------------------------------------------------------
state_variables_T& state_variables_T::rotate_l2g(const state_variables_T& sv, const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{

    state_variables::rotate_l2g(sv, psi, theta, phi);

    dSdE = sv.dSdE;
    dSdEt = sv.dSdEt;
    dSdT = sv.dSdT;
    drpldE = sv.drpldE;
    drpldT = sv.drpldT;
    Q = sv.Q;
    rpl = sv.rpl;
    
  	if(fabs(phi) > iota) {
		dSdE = rotateL(dSdE, -phi, axis_phi);
		dSdEt = rotateL(dSdEt, -phi, axis_phi);
        dSdT = rotate_stress(dSdT, -phi, axis_phi);
		drpldE = rotate_strain(drpldE, -phi, axis_phi);
	}
  	if(fabs(theta) > iota) {
		dSdE = rotateL(dSdE, -theta, axis_theta);
		dSdEt = rotateL(dSdEt, -theta, axis_theta);
        dSdT = rotate_stress(dSdT, -theta, axis_theta);
		drpldE = rotate_strain(drpldE, -theta, axis_theta);
	}
	if(fabs(psi) > iota) {
		dSdE = rotateL(dSdE, -psi, axis_psi);
		dSdEt = rotateL(dSdEt, -psi, axis_psi);
        dSdT = rotate_stress(dSdT, -psi, axis_psi);
		drpldE = rotate_strain(drpldE, -psi, axis_psi);
	}
    
	return *this;
}

//----------------------------------------------------------------------
state_variables_T& state_variables_T::rotate_g2l(const state_variables_T& sv, const double &psi, const double &theta, const double &phi)
//----------------------------------------------------------------------
{

    state_variables::rotate_g2l(sv, psi, theta, phi);
    
    dSdE = sv.dSdE;
    dSdEt = sv.dSdEt;
    dSdT = sv.dSdT;
    drpldE = sv.drpldE;
    drpldT = sv.drpldT;
    Q = sv.Q;
    rpl = sv.rpl;
    
  	if(fabs(psi) > iota) {
		dSdE = rotateL(dSdE, psi, axis_psi);
		dSdEt = rotateL(dSdEt, psi, axis_psi);
        dSdT = rotate_stress(dSdT, psi, axis_psi);
        drpldE = rotate_strain(drpldE, psi, axis_psi);
        
	}			
	if(fabs(theta) > iota) {
		dSdE = rotateL(dSdE, theta, axis_theta);
		dSdEt = rotateL(dSdEt, theta, axis_theta);
        dSdT = rotate_stress(dSdT, theta, axis_theta);
        drpldE = rotate_strain(drpldE, theta, axis_theta);
	}
	if(fabs(phi) > iota) {
		dSdE = rotateL(dSdE, phi, axis_phi);
		dSdEt = rotateL(dSdEt, phi, axis_phi);
        dSdT = rotate_stress(dSdT, phi, axis_phi);
        drpldE = rotate_strain(drpldE, phi, axis_phi);
    }
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const state_variables_T& sv)
//--------------------------------------------------------------------------
{
	s << "Etot: \n" << sv.Etot << "\n";
	s << "DEtot: \n" << sv.DEtot << "\n";
	s << "sigma: \n" << sv.sigma << "\n";
	s << "sigma_start: \n" << sv.sigma_start << "\n";
    s << "T: \n" << sv.T << "\n";
    s << "DT: \n" << sv.DT << "\n";
    s << "Q: \n" << sv.Q << "\n";
    s << "rpl: \n" << sv.rpl << "\n";
	s << "dSdE: \n" << sv.dSdE << "\n";
	s << "dSdEt: \n" << sv.dSdEt << "\n";
	s << "dSdT: \n" << sv.dSdT << "\n";
	s << "drpldE: \n" << sv.drpldE << "\n";
	s << "drpldT: \n" << sv.drpldT << "\n";
    
    s << "nstatev: \n" << sv.nstatev << "\n";
    if (sv.nstatev) {
        s << "statev: \n";
        s << sv.statev.t();
        s << "\n";
        s << sv.statev_start.t();
        s << "\n";
    }    
    
	s << "\n";

	return s;
}

} //namespace smart
