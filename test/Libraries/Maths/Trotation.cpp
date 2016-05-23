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

///@file Tconstitutive.cpp
///@brief Test for Constitutive tensors in Voigt notation
///@version 1.0

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "rotation"
#include <boost/test/unit_test.hpp>

#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Maths/rotation.hpp>
#include <smartplus/Libraries/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;
using namespace smart;

BOOST_AUTO_TEST_CASE( rotation )
{

    vec test = zeros(6);
    test(0) = 4.;
    test(1) = 2.;
    test(2) = 6.;
    test(3) = 8.;
    test(4) = 3.;
    test(5) = 7.;

    double psi = 12.5;
    double theta = 32.;
    double phi = -4.5;
    
    psi *= (pi/180.);
    theta *= (theta/180.);
    phi *= (psi/180.);
    
    mat R1 = { {cos(psi),sin(psi),0}, {-sin(psi),cos(psi), 0}, {0,0,1}};
    mat R2 = { {1,0,0}, {0,cos(theta),sin(theta)}, {0,-sin(theta),cos(theta)} };
    mat R3 = { {cos(phi),sin(phi),0}, {-sin(phi),cos(phi), 0}, {0,0,1}};
    
    mat R = R3*R2*R1;

    //test of rotate A
    mat S_c = Eshelby_cylinder(0.12);
    
    mat S_c2 = rotateA(S_c, psi, 3);
    S_c2 = rotateA(S_c2, theta, 1);
    S_c2 = rotateA(S_c2, phi, 3);
        
    mat S_c3 = rotateA(S_c,R);
    
    BOOST_CHECK( norm(S_c3-S_c2,2) < 1.E-9 );
    
}
