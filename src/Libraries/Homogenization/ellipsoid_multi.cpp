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

///@file ellipsoid_multi.cpp
///@brief Micromechanical characteristics of a phase
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <smartplus/Libraries/Maths/rotation.hpp>
#include <smartplus/Libraries/Geometry/ellipsoid.hpp>
#include <smartplus/Libraries/Homogenization/ellipsoid_multi.hpp>
#include <smartplus/Libraries/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;

namespace smart{

//Definition of the static variables
int ellipsoid_multi::mp;
int ellipsoid_multi::np;
vec ellipsoid_multi::x;
vec ellipsoid_multi::wx;
vec ellipsoid_multi::y;
vec ellipsoid_multi::wy;
    
    
//=====Private methods for ellipsoid_multi===================================

//=====Public methods for ellipsoid_multi====================================

/*!
  \brief default constructor
*/
    
//-------------------------------------------------------------
ellipsoid_multi::ellipsoid_multi() : phase_multi(), S_loc(6,6), P_loc(6,6), T_loc(6,6), T(6,6)
//-------------------------------------------------------------
{
    //This calls only the constructor of the two matrix A & B
    
}

/*!
  \brief Constructor with parameters
*/

//-------------------------------------------------------------
ellipsoid_multi::ellipsoid_multi(const mat &mA, const mat &mA_start, const mat &mB, const mat &mB_start, const mat &mS_loc, const mat &mP_loc, const mat &mT_loc, const mat &mT) : phase_multi(mA, mA_start, mB, mB_start), S_loc(6,6), P_loc(6,6), T_loc(6,6), T(6,6)
//-------------------------------------------------------------
{
    S_loc = mS_loc;
    P_loc = mP_loc;
    T_loc = mT_loc;
    T = mT;
}

/*!
  \brief Copy constructor
  \param s ellipsoid_multi object to duplicate
*/
    
//------------------------------------------------------
ellipsoid_multi::ellipsoid_multi(const ellipsoid_multi& pc) : phase_multi(pc), S_loc(6,6), P_loc(6,6), T_loc(6,6), T(6,6)
//------------------------------------------------------
{
    S_loc = pc.S_loc;
    P_loc = pc.P_loc;
    T_loc = pc.T_loc;
    T = pc.T;
}

/*!
  \brief Destructor

  Deletes phase_multi (the arma::mat).
*/

//-------------------------------------
ellipsoid_multi::~ellipsoid_multi() {}
//-------------------------------------

/*!
  \brief Standard operator = for phase_multi
*/

//-------------------------------------
void ellipsoid_multi::fillS_loc(const mat& Lt_m, const ellipsoid &ell)
//-------------------------------------
{
    mat Ltm_local_geom = rotate_g2l_L(Lt_m, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    S_loc = Eshelby(Ltm_local_geom, ell.a1, ell.a2, ell.a3, x, wx, y, wy, mp, np);
}
    
//-------------------------------------
void ellipsoid_multi::fillP_loc(const mat& Lt_m, const ellipsoid &ell)
//-------------------------------------
{
    mat Ltm_local_geom = rotate_g2l_L(Lt_m, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    P_loc = T_II(Ltm_local_geom, ell.a1, ell.a2, ell.a3, x, wx, y, wy, mp, np);
}
    

//-------------------------------------
void ellipsoid_multi::fillT(const mat& Lt_m, const mat& Lt, const ellipsoid &ell)
//This method correspond to the classical Eshelby method
//-------------------------------------
{
    mat Lt_m_local_geom = rotate_g2l_L(Lt_m, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    S_loc = Eshelby(Lt_m_local_geom, ell.a1, ell.a2, ell.a3, x, wx, y, wy, mp, np);
    mat Lt_local_geom = rotate_g2l_L(Lt, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    
    T_loc = inv(eye(6,6) + S_loc*inv(Lt_m_local_geom)*(Lt_local_geom - Lt_m_local_geom));
    
    T = rotate_l2g_A(T_loc, ell.psi_geom, ell.theta_geom, ell.phi_geom);
}

//-------------------------------------
void ellipsoid_multi::fillT_m(const mat& Lt_m, const mat& Lt, const ellipsoid &ell)
//This method corresponf to the Ponte-astenada and Willis method
//-------------------------------------
{
    mat Lt_m_local_geom = rotate_g2l_L(Lt_m, ell.psi_geom, ell.theta_geom, ell.phi_geom);    
    P_loc = T_II(Lt_m_local_geom, ell.a1, ell.a2, ell.a3, x, wx, y, wy, mp, np);
    mat P = rotate_l2g_M(P_loc, ell.psi_geom, ell.theta_geom, ell.phi_geom);
    
    T = inv(inv(Lt - Lt_m) + P);
}
        
//----------------------------------------------------------------------
ellipsoid_multi& ellipsoid_multi::operator = (const ellipsoid_multi& pc)
//----------------------------------------------------------------------
{
    A = pc.A;
    B = pc.B;

    A_start = pc.A_start;
    B_start = pc.B_start;
    
    S_loc = pc.S_loc;
    P_loc = pc.P_loc;
    T_loc = pc.T_loc;
    T = pc.T;
    
	return *this;
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const ellipsoid_multi& pc)
//--------------------------------------------------------------------------
{
	s << "Display phase multi:\n";
	s << "Display strain concentration tensor:\n";
    s << pc.A;
    s << "Display stress concentration tensor:\n";
    s << pc.B;

    s << "Display Eshelby tensor (local coordinates):\n";
    s << pc.S_loc;
    s << "Display Polarization tensor (local coordinates):\n";
    s << pc.P_loc;
    s << "Display Interaction concentration tensor (global coordinates):\n";
    s << pc.T_loc;
    s << "Display Interaction concentration tensor (global coordinates):\n";
    s << pc.T;
    
    s << "\n\n";

	return s;
}

} //namespace smart
