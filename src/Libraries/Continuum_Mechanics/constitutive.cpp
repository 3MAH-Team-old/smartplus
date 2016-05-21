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

///@file constitutive.cpp
///@brief Constitutive tensors in Voigt notation
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <armadillo>
#include <smartplus/Libraries/Continuum_Mechanics/constitutive.hpp>

using namespace std;
using namespace arma;

namespace smart{

//Returns the fourth order identity tensor written in Voigt notation Ireal
mat Ireal(){
	
	mat Ireal = eye(6,6);
	for (int i=3; i<6; i++)
		Ireal(i,i) = 0.5;
		
	return Ireal;
}

//Returns the volumic of the identity tensor Ireal written in Voigt notation
mat Ivol(){
	
	mat Ivol = zeros(6,6);
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			Ivol(i,j) = (1./3.);
		}
	}
	return Ivol;
}

//Returns the deviatoric of the identity tensor Ireal written in Voigt notation
mat Idev(){
	
	return	Ireal() - Ivol();
}

//Returns the fourth order identity tensor Iˆ written in Voigt notation
mat Ireal2(){
	
	mat Ireal2 = eye(6,6);
	for (int i=3; i<6; i++)
		Ireal2(i,i) = 2.;
		
	return Ireal2;
}

//Returns the deviatoric of the identity tensor Iˆ written in Voigt notation
mat Idev2(){
	
	return	Ireal2() - Ivol();
}

//Returns the expansion vector
vec Ith(){
	
	vec Ith = zeros(6);	
	for (int i=0; i<3; i++)
		Ith(i) = 1.;
		
	return Ith;
}

//Returns the stress 2 strain operator
vec Ir2(){
	
	vec Ir2 = ones(6);	
	for (int i=3; i<6; i++)
		Ir2(i) = 2.;
	
	return Ir2;
}

//Returns the strain 2 stress operator
vec Ir05(){
	
	vec Ir05 = ones(6);	
	for (int i=3; i<6; i++)
		Ir05(i) = 0.5;
	
	return Ir05;
}

//Provides the elastic stiffness tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
// ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
mat L_iso(const double &C1, const double &C2, const std::string &conv) {
	
	double K = 0.;
	double mu = 0.;	
	
	if(conv == "Enu") {
		assert(C1 > 0.);
		assert(C2 < 0.5);
		K = C1/(3.*(1.-2*C2));
		mu = C1/(2.*(1.+C2));			
	}
	else if(conv == "nuE") {	
		assert(C2 > 0.);
		assert(C1 < 0.5);
		K = C2/(3.*(1.-2*C1));
		mu = C2/(2.*(1.+C1));
	}
	else if((conv == "Kmu")||(conv == "KG")) {
		assert(C1 > 0.);
		assert(C2 > 0.);
		K = C1;
		mu = C2;
	}
	else if((conv == "muK")||(conv == "GK")) {
		assert(C1 > 0.);
		assert(C2 > 0.);
		K = C2;
		mu = C1;
	}
	else if((conv == "lambdamu")||(conv == "lambdaG")) {	
		assert(C1 > 0.);
		assert(C2 > 0.);		
		K = C1+(2./3.)*C2;
		mu = C2;
	}
	else if((conv == "mulambda")||(conv == "Glambda")) {	
		assert(C1 > 0.);
		assert(C2 > 0.);		
		K = C2+(2./3.)*C1;
		mu = C1;
	}	
	else {
		cout << "ERROR : Please use a valid couple of elastic constants";
	}
	
	return 3.*K*Ivol() + 2.*mu*Idev();
}

//Provides the elastic compliance tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
//‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
mat M_iso(const double &C1, const double &C2, const string &conv) {
	
	double K = 0.;
	double mu = 0.;	
	
	if(conv == "Enu") {
		assert(C1 > 0.);
		assert(C2 < 0.5);
		K = C1/(3.*(1.-2*C2));
		mu = C1/(2.*(1.+C2));		
	}
	else if(conv == "nuE") {	
		assert(C2 > 0.);
		assert(C1 < 0.5);
		K = C2/(3.*(1.-2*C1));
		mu = C2/(2.*(1.+C1));
	}
	else if((conv == "Kmu")||(conv == "KG")) {
		assert(C1 > 0.);
		assert(C2 > 0.);
		K = C1;
		mu = C2;
	}
	else if((conv == "muK")||(conv == "GK")) {
		assert(C1 > 0.);
		assert(C2 > 0.);
		K = C2;
		mu = C1;
	}
	else if((conv == "lambdamu")||(conv == "lambdaG")) {	
		assert(C1 > 0.);
		assert(C2 > 0.);		
		K = C1+(2./3.)*C2;
		mu = C2;
	}
	else if((conv == "mulambda")||(conv == "Glambda")) {	
		assert(C1 > 0.);
		assert(C2 > 0.);		
		K = C2+(2./3.)*C1;
		mu = C1;
	}	
	else {
		cout << "ERROR : Please use a valid couple of elastic constants";
        exit(0);
	}
	
	return 1/(3.*K)*Ivol() + 1/(2.*mu)*Idev2();
}

//Returns the elastic stiffness tensor for a cubic material.
//Arguments are the stiffness coefficients C11, C12 and C44, or E, nu and G
mat L_cubic(const double &C1, const double &C2, const double &C4, const string &conv){
	
    mat L = zeros(6,6);
	if (conv == "Cii") {
        L(0,0) = C1;
        L(0,1) = C2;
        L(0,2) = C2;
        L(1,0) = C2;
        L(1,1) = C1;
        L(1,2) = C2;
        L(2,0) = C2;
        L(2,1) = C2;
        L(2,2) = C1;
        L(3,3) = C4;
        L(4,4) = C4;
        L(5,5) = C4;
    }
    else if(conv == "EnuG") {
        double E = C1;
        double nu = C2;
        double G = C4;
    
        double C11 = E*( 1-pow(nu,2))/(1 - 3*pow(nu,2) - 2*pow(nu,3) );
        double C12 = C11/(1./nu - 1.);
        double C44 = G;

        L(0,0) = C11;
        L(0,1) = C12;
        L(0,2) = C12;
        L(1,0) = C12;
        L(1,1) = C11;
        L(1,2) = C12;
        L(2,0) = C12;
        L(2,1) = C12;
        L(2,2) = C11;
        L(3,3) = C44;
        L(4,4) = C44;
        L(5,5) = C44;
    }
    else {
        cout << "ERROR : Please use a valid couple of elastic constants";
        exit(0);
    }
    
	return L;
}

//Returns the elastic compliance tensor for an isotropic material.
//Arguments are the stiffness coefficients C11, C12 and C44, or E, nu and G.
mat M_cubic(const double &C1, const double &C2, const double &C4, const string &conv){
	
    mat M = zeros(6,6);
	if (conv == "Cii") {
        double delta = (C1-C2)*(C1+2.*C2);
        M(0,0) = (C1+C2)/delta;
        M(0,1) = -C2/delta;
        M(0,2) = -C2/delta;
        M(1,0) = -C2/delta;
        M(1,1) = (C1+C2)/delta;
        M(1,2) = -C2/delta;
        M(2,0) = -C2/delta;
        M(2,1) = -C2/delta;
        M(2,2) = (C1+C2)/delta;
        M(3,3) = 1./C4;
        M(4,4) = 1./C4;
        M(5,5) = 1./C4;
    }
    else if (conv == "EnuG") {
        double E = C1;
        double nu = C2;
        double G = C4;
        
        M(0,0) = 1./E;
        M(0,1) = -nu/E;
        M(0,2) = -nu/E;
        M(1,0) = -nu/E;
        M(1,1) = 1./E;
        M(1,2) = -nu/E;
        M(2,0) = -nu/E;
        M(2,1) = -nu/E;
        M(2,2) = 1./E;
        M(3,3) = 1./G;
        M(4,4) = 1./G;
        M(5,5) = 1./G;
    }
    else {
        cout << "ERROR : Please use a valid couple of elastic constants";
        exit(0);
    }
	
	return M;
}

//Returns the elastic stiffness tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
mat L_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const string &conv){
	mat L = zeros(6,6);
	
	if (conv == "Cii") {
	    L(0,0) = C11;
	    L(0,1) = C12;
	    L(0,2) = C13;
	    L(1,0) = C12;
	    L(1,1) = C22;
	    L(1,2) = C23;
	    L(2,0) = C13;
	    L(2,1) = C23;
	    L(2,2) = C33;
	    L(3,3) = C44;
	    L(4,4) = C55;
	    L(5,5) = C66;
		
	}
	else if (conv == "EnuG") {
        double Ex = C11;
        double Ey = C12;
        double Ez = C13;
        double nuxy = C22;
        double nuxz = C23;
        double nuyz = C33;
        double Gxy = C44;
        double Gxz = C55;
        double Gyz = C66;
        
        double nuyx=nuxy*(Ey/Ex);
        double nuzy=nuyz*(Ez/Ey);
        double nuzx=nuxz*(Ez/Ex);

        double delta =1./(1-nuxy*nuyx-nuyz*nuzy-nuzx*nuxz-2.*nuyx*nuzy*nuxz);
         
        L(0,0) = Ex*(1-nuyz*nuzy)*delta;
        L(0,1) = Ex*(nuyx+nuzx*nuyz)*delta;
        L(0,2) = Ex*(nuzx+nuyx*nuzy)*delta;
        L(1,0) = L(0,1);
        L(1,1) = Ey*(1-nuxz*nuzx)*delta;
        L(1,2) = Ey*(nuzy+nuxy*nuzx)*delta;
        L(2,0) = L(0,2);
        L(2,1) = L(1,2);
        L(2,2) = Ez*(1-nuxy*nuyx)*delta;
        L(3,3) = Gxy;
        L(4,4) = Gxz;
        L(5,5) = Gyz;
    }
    else {
        cout << "ERROR : Please use a valid couple of elastic constants";
        exit(0);
    }
	
	return L;
}

//Returns the elastic compliance tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
mat M_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const string &conv){
	
    mat M = zeros(6,6);
    mat L = zeros(6,6);
    
    if (conv == "Cii") {
	    L(0,0) = C11;
	    L(0,1) = C12;
	    L(0,2) = C13;
	    L(1,0) = C12;
	    L(1,1) = C22;
	    L(1,2) = C23;
	    L(2,0) = C13;
	    L(2,1) = C23;
	    L(2,2) = C33;
	    L(3,3) = C44;
	    L(4,4) = C55;
	    L(5,5) = C66;
		
        M = inv(L);
	}
    
    else if (conv == "EnuG") {
		double Ex = C11;
		double Ey = C12;
		double Ez = C13;
		double nuxy = C22;
		double nuxz = C23;
		double nuyz = C33;
		double Gxy = C44;
		double Gxz = C55;
		double Gyz = C66;
        
		M(0,0) = 1/Ex;
		M(0,1) = -nuxy/Ex;
		M(0,2) = -nuxz/Ex;
		M(1,0) = -nuxy/Ex;
		M(1,1) = 1/Ey;
		M(1,2) = -nuyz/Ey;
		M(2,0) = -nuxz/Ex;
		M(2,1) = -nuyz/Ey;
		M(2,2) = 1/Ez;
		M(3,3) = 1/Gxy;
		M(4,4) = 1/Gxz;
		M(5,5) = 1/Gyz;
	}
    else {
        cout << "ERROR : Please use a valid couple of elastic constants";
        exit(0);
    }
	
	return M;
}

//Returns the elastic stiffness tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
mat L_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis){
	
    mat L = zeros(6,6);
    double delta = (1.+nuTT)*(2*EL*nuTL*nuTL+ET*(nuTT-1.))/ET;

	switch(axis) {
		case 1: {	
		    L(0,0) = EL*(nuTT*nuTT-1.)/delta;
		    L(0,1) = -EL*nuTL*(1+nuTT)/delta;
		    L(0,2) = -EL*nuTL*(1+nuTT)/delta;
		    L(1,0) = -EL*nuTL*(1+nuTT)/delta;
		    L(1,1) = (EL*nuTL*nuTL-ET)/delta;
		    L(1,2) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
		    L(2,0) = -EL*nuTL*(1+nuTT)/delta;
		    L(2,1) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
		    L(2,2) = (EL*nuTL*nuTL-ET)/delta;
		    L(3,3) = GLT;
		    L(4,4) = GLT;
		    L(5,5) = ET/(2.*(1.+nuTT));
		    break;
		}
		case 2: {
		    L(0,0) = (EL*nuTL*nuTL-ET)/delta;
		    L(0,1) = -EL*nuTL*(1+nuTT)/delta;
		    L(0,2) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
		    L(1,0) = -EL*nuTL*(1+nuTT)/delta;
		    L(1,1) = EL*(nuTT*nuTT-1.)/delta;
		    L(1,2) = -EL*nuTL*(1+nuTT)/delta;
		    L(2,0) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
		    L(2,1) = -EL*nuTL*(1+nuTT)/delta;
		    L(2,2) = (EL*nuTL*nuTL-ET)/delta;
		    L(3,3) = GLT;
		    L(4,4) = ET/(2.*(1.+nuTT));
		    L(5,5) = GLT;
			break;
		}
		case 3: {
			L(0,0) = (EL*nuTL*nuTL-ET)/delta;
			L(0,1) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
			L(0,2) = -EL*nuTL*(1+nuTT)/delta;
			L(1,0) = -(EL*nuTL*nuTL+ET*nuTT)/delta;
			L(1,1) = (EL*nuTL*nuTL-ET)/delta;
			L(1,2) = -EL*nuTL*(1+nuTT)/delta;
			L(2,0) = -EL*nuTL*(1+nuTT)/delta;
			L(2,1) = -EL*nuTL*(1+nuTT)/delta;
			L(2,2) = EL*(nuTT*nuTT-1.)/delta;
			L(3,3) = ET/(2.*(1.+nuTT));
			L(4,4) = GLT;
			L(5,5) = GLT;
			break;	
		}
        default : {
            cout << "ERROR : Please use a valid couple of elastic constants";
            exit(0);
        }
            
	}

    
	return L;
}

//Returns the elastic compliance tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
mat M_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis){
	
	mat M = zeros(6,6);
	switch(axis) {
		case 1: {	
			M(0,0) = 1./EL;
			M(0,1) = -nuTL/ET;
			M(0,2) = -nuTL/ET;
			M(1,0) = -nuTL/ET;
			M(1,1) = 1./ET;
			M(1,2) = -nuTT/ET;
			M(2,0) = -nuTL/ET;
			M(2,1) = -nuTT/ET;
			M(2,2) = 1./ET;
			M(3,3) = 1/GLT;
			M(4,4) = 1/GLT;
			M(5,5) = (2.*(1.+nuTT))/ET;
			break;
		}
		case 2: {	
			M(0,0) = 1./ET;
			M(0,1) = -nuTL/ET;
			M(0,2) = -nuTT/ET;
			M(1,0) = -nuTL/ET;
			M(1,1) = 1./EL;
			M(1,2) = -nuTL/ET;
			M(2,0) = -nuTT/ET;
			M(2,1) = -nuTL/ET;
			M(2,2) = 1./ET;
			M(3,3) = 1/GLT;
			M(4,4) = (2.*(1.+nuTT))/ET;
			M(5,5) = 1/GLT;
			break;
		}
		case 3: {	
			M(0,0) = 1./ET;
			M(0,1) = -nuTT/ET;
			M(0,2) = -nuTL/ET;
			M(1,0) = -nuTT/ET;
			M(1,1) = 1./ET;
			M(1,2) = -nuTL/ET;
			M(2,0) = -nuTL/ET;
			M(2,1) = -nuTL/ET;
			M(2,2) = 1./EL;
			M(3,3) = (2.*(1.+nuTT))/ET;
			M(4,4) = 1/GLT;
			M(5,5) = 1/GLT;
			break;
		}
        default : {
            cout << "ERROR : Please use a valid couple of elastic constants";
            exit(0);
        }
	}			
	return M;
}

//Provides the viscous tensor H an isotropic material.
//The two first arguments are a couple of viscous coefficients (the first is bulk, the second is shear).
mat H_iso(const double &etaB, const double &etaS) {
		
	return 3.*etaB*Ivol() + 2.*etaS*Idev();
}
    
} //namespace smart
