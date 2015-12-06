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
#define BOOST_TEST_MODULE "constitutive"
#include <boost/test/unit_test.hpp>

#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Continuum_Mechanics/constitutive.hpp>

using namespace std;
using namespace arma;
using namespace smart;

BOOST_AUTO_TEST_CASE( L_iso_M_iso )
{
    double E = 70000.;
    double nu = 0.3;
    
    double mu = E/(2.*(1+nu));
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    
    mat Lt = zeros(6,6);
    Lt(0,0) = 2.*mu+lambda;
    Lt(0,1) = lambda;
    Lt(0,2) = lambda;
    Lt(1,0) = lambda;
    Lt(1,1) = 2.*mu+lambda;
    Lt(1,2) = lambda;
    Lt(2,0) = lambda;
    Lt(2,1) = lambda;
    Lt(2,2) = 2.*mu+lambda;
    Lt(3,3) = mu;
    Lt(4,4) = mu;
    Lt(5,5) = mu;
    
    mat Mt = zeros(6,6);
    Mt(0,0) = 1./E;
    Mt(0,1) = -nu/E;
    Mt(0,2) = -nu/E;
    Mt(1,0) = -nu/E;
    Mt(1,1) = 1./E;
    Mt(1,2) = -nu/E;
    Mt(2,0) = -nu/E;
    Mt(2,1) = -nu/E;
    Mt(2,2) = 1./E;
    Mt(3,3) = 2.*(1.+nu)/E;
    Mt(4,4) = 2.*(1.+nu)/E;
    Mt(5,5) = 2.*(1.+nu)/E;
    
    //Test of L_iso function
    mat Ltest = L_iso(E, nu, "Enu");
    BOOST_CHECK( norm(Ltest - Lt,2) < 1.E-9 );
	Ltest = L_iso(lambda, mu, "lambdamu");
	BOOST_CHECK( norm(Ltest - Lt,2) < 1.E-9 );
    
	//Test of M_iso function
    mat Mtest = M_iso(E, nu, "Enu");
	BOOST_CHECK( norm(Mtest - Mt,2) < 1.E-9 );
	Mtest = M_iso(mu, lambda, "mulambda");
	BOOST_CHECK( norm(Mtest - Mt,2) < 1.E-9 );
}

BOOST_AUTO_TEST_CASE( L_cubic_M_cubic )
{
    double C11 = 1000;
    double C12 = 400;
    double C44 = 500;
    
    mat Lcub = zeros(6,6);
    Lcub(0,0) = C11;
    Lcub(0,1) = C12;
    Lcub(0,2) = C12;
    Lcub(1,0) = C12;
    Lcub(1,1) = C11;
    Lcub(1,2) = C12;
    Lcub(2,0) = C12;
    Lcub(2,1) = C12;
    Lcub(2,2) = C11;
    Lcub(3,3) = C44;
    Lcub(4,4) = C44;
    Lcub(5,5) = C44;
    
    //Test of L_cubic function
    mat Ltest = L_cubic(C11,C12,C44);
    BOOST_CHECK( norm(Ltest - Lcub,2) < 1.E-9 );
    //Test of M_cubic function
    mat Mtest = M_cubic(C11,C12,C44);
    BOOST_CHECK( norm(Mtest - inv(Lcub),2) < 1.E-9 );
}

BOOST_AUTO_TEST_CASE( L_isotrans_M_isotrans )
{
    double EL = 10000.;
    double ET = 20000.;
    double nuTL = 0.3;
    double nuTT = 0.3;
    double GLT = 8000.;
    
    mat Mtrans = zeros(6,6);
    Mtrans(0,0) = 1./EL;
    Mtrans(0,1) = -nuTL/ET;
    Mtrans(0,2) = -nuTL/ET;
    Mtrans(1,0) = -nuTL/ET;
    Mtrans(1,1) = 1./ET;
    Mtrans(1,2) = -nuTT/ET;
    Mtrans(2,0) = -nuTL/ET;
    Mtrans(2,1) = -nuTT/ET;
    Mtrans(2,2) = 1./ET;
    Mtrans(3,3) = 1/GLT;
    Mtrans(4,4) = 1/GLT;
    Mtrans(5,5) = (2.*(1.+nuTT))/ET;
    
	//Test of L_isotrans function axis 1
    mat Ltest = L_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
	BOOST_CHECK( norm(Ltest - inv(Mtrans),2) < 1.E-9 );
    //Test of M_isotrans function axis 1
    mat Mtest = M_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
    BOOST_CHECK( norm(Mtest - Mtrans,2) < 1.E-9 );
    
}