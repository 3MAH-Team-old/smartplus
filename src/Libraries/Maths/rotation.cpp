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

///@file rotation.cpp
///@brief rotation of a Voigt tensor
///@version 1.0

#include <math.h>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Continuum_Mechanics/contimech.hpp>

using namespace std;
using namespace arma;

namespace smart{

void Rot_strain(vec &E, const mat &DR) {
	
	double a = DR(0,0);
	double b = DR(0,1);
	double c = DR(0,2);
	double d = DR(1,0);
	double e = DR(1,1);
	double f = DR(1,2);
	double g = DR(2,0);
	double h = DR(2,1);
	double i = DR(2,2);
	
	mat QE= zeros(6,6);
	QE(0,0) = a*a;
	QE(0,1) = b*b;
	QE(0,2) = c*c;
	QE(0,3) = a*b;
	QE(0,4) = a*c;
	QE(0,5) = b*c;
	QE(1,0) = d*d;
	QE(1,1) = e*e;
	QE(1,2) = f*f;
	QE(1,3) = d*e;
	QE(1,4) = d*f;
	QE(1,5) = e*f;
	QE(2,0) = g*g; 
	QE(2,1) = h*h;
	QE(2,2) = i*i;
	QE(2,3) = g*h;
	QE(2,4) = g*i;
	QE(2,5) = h*i;
	QE(3,0) = 2.*a*d;
	QE(3,1) = 2.*b*e;
	QE(3,2) = 2.*c*f;
	QE(3,3) = b*d+a*e;
	QE(3,4) = c*d+a*f;
	QE(3,5) = c*e+b*f;
	QE(4,0) = 2.*a*g;
	QE(4,1) = 2.*b*h;
	QE(4,2) = 2.*c*i;
	QE(4,3) = b*g+a*h;
	QE(4,4) = c*g+a*i;
	QE(4,5) = c*h+b*i;
	QE(5,0) = 2.*d*g;
	QE(5,1) = 2.*e*h;
	QE(5,2) = 2.*f*i;
	QE(5,3) = e*g+d*h;
	QE(5,4) = f*g+d*i;
	QE(5,5) = f*h+e*i;

	E = QE * E;
}

void Rot_stress(vec &S, const mat &DR) {

	double a = DR(0,0);
	double b = DR(0,1);
	double c = DR(0,2);
	double d = DR(1,0);
	double e = DR(1,1);
	double f = DR(1,2);
	double g = DR(2,0);
	double h = DR(2,1);
	double i = DR(2,2);
	
	mat QS= zeros(6,6);
	QS(0,0) = a*a;
	QS(0,1) = b*b;
	QS(0,2) = c*c;
	QS(0,3) = 2.*a*b;
	QS(0,4) = 2.*a*c;
	QS(0,5) = 2.*b*c;
	QS(1,0) = d*d;
	QS(1,1) = e*e;
	QS(1,2) = f*f;
	QS(1,3) = 2.*d*e;
	QS(1,4) = 2.*d*f;
	QS(1,5) = 2.*e*f;
	QS(2,0) = g*g; 
	QS(2,1) = h*h;
	QS(2,2) = i*i;
	QS(2,3) = 2.*g*h;
	QS(2,4) = 2.*g*i;
	QS(2,5) = 2.*h*i;
	QS(3,0) = a*d;
	QS(3,1) = b*e;
	QS(3,2) = c*f;
	QS(3,3) = b*d+a*e;
	QS(3,4) = c*d+a*f;
	QS(3,5) = c*e+b*f;
	QS(4,0) = a*g;
	QS(4,1) = b*h;
	QS(4,2) = c*i;
	QS(4,3) = b*g+a*h;
	QS(4,4) = c*g+a*i;
	QS(4,5) = c*h+b*i;
	QS(5,0) = d*g;
	QS(5,1) = e*h;
	QS(5,2) = f*i;
	QS(5,3) = e*g+d*h;
	QS(5,4) = f*g+d*i;
	QS(5,5) = f*h+e*i;
	
	S = QS * S;
}

mat fillQS(const double &alpha, const int &axis) {

	mat QS = zeros(6,6);
	double c = cos(alpha);
	double s = sin(alpha);
	
	switch(axis) {
		case 1: {
			QS(0,0) = 1.;
			QS(0,1) = 0.;
			QS(0,2) = 0.;
			QS(0,3) = 0.;
			QS(0,4) = 0.;
			QS(0,5) = 0.;
			QS(1,0) = 0.;
			QS(1,1) = c*c;
			QS(1,2) = s*s;
			QS(1,3) = 0.;
			QS(1,4) = 0.;
			QS(1,5) = 2*s*c;
			QS(2,0) = 0.;
			QS(2,1) = s*s;
			QS(2,2) = c*c;
			QS(2,3) = 0.;
			QS(2,4) = 0.;
			QS(2,5) = -2*s*c;
			QS(3,0) = 0.;
			QS(3,1) = 0.;
			QS(3,2) = 0.;
			QS(3,3) = c;
			QS(3,4) = s;
			QS(3,5) = 0.;
			QS(4,0) = 0.;
			QS(4,1) = 0.;
			QS(4,2) = 0.;
			QS(4,3) = -s;
			QS(4,4) = c;
			QS(4,5) = 0.;
			QS(5,0) = 0.;
			QS(5,1) = -s*c;
			QS(5,2) = s*c;
			QS(5,3) = 0.;
			QS(5,4) = 0.;
			QS(5,5) = c*c-s*s;
			
			break;
		}
		case 2: {
			QS(0,0) = c*c;
			QS(0,1) = 0.;
			QS(0,2) = s*s;
			QS(0,3) = 0.;
			QS(0,4) = -2*c*s;
			QS(0,5) = 0.;
			QS(1,0) = 0.;
			QS(1,1) = 1.;
			QS(1,2) = 0.;
			QS(1,3) = 0.;
			QS(1,4) = 0.;
			QS(1,5) = 0.;
			QS(2,0) = s*s;
			QS(2,1) = 0.;
			QS(2,2) = c*c;
			QS(2,3) = 0.;
			QS(2,4) = 2*c*s;
			QS(2,5) = 0.;
			QS(3,0) = 0.;
			QS(3,1) = 0.;
			QS(3,2) = 0.;
			QS(3,3) = c;
			QS(3,4) = 0.;
			QS(3,5) = -s;
			QS(4,0) = c*s;
			QS(4,1) = 0.;
			QS(4,2) = -c*s;
			QS(4,3) = 0.;
			QS(4,4) = c*c-s*s;
			QS(4,5) = 0.;
			QS(5,0) = 0.;
			QS(5,1) = 0.;
			QS(5,2) = 0.;
			QS(5,3) = s;
			QS(5,4) = 0.;
			QS(5,5) = c;
			
			break;
			
		}
		case 3: {
			QS(0,0) = c*c;
			QS(0,1) = s*s;
			QS(0,2) = 0.;
			QS(0,3) = 2*s*c;
			QS(0,4) = 0.;
			QS(0,5) = 0.;
			QS(1,0) = s*s;
			QS(1,1) = c*c;
			QS(1,2) = 0.;
			QS(1,3) = -2*s*c;
			QS(1,4) = 0.;
			QS(1,5) = 0.;
			QS(2,0) = 0.;
			QS(2,1) = 0.;
			QS(2,2) = 1.;
			QS(2,3) = 0.;
			QS(2,4) = 0.;
			QS(2,5) = 0.;
			QS(3,0) = -s*c;
			QS(3,1) = s*c;
			QS(3,2) = 0.;
			QS(3,3) = c*c-s*s;
			QS(3,4) = 0.;
			QS(3,5) = 0.;
			QS(4,0) = 0.;
			QS(4,1) = 0.;
			QS(4,2) = 0.;
			QS(4,3) = 0.;
			QS(4,4) = c;
			QS(4,5) = s;
			QS(5,0) = 0.;
			QS(5,1) = 0.;
			QS(5,2) = 0.;
			QS(5,3) = 0.;
			QS(5,4) = -s;
			QS(5,5) = c;
						
			break;
		}
		default: {
			cout << "Please choose the axis 1,2 or 3\n";
		}
	}
	return QS;		
}

mat fillQE(const double &alpha, const int &axis) {

	mat QE = zeros(6,6);
	double c = cos(alpha);
	double s = sin(alpha);
	
	switch(axis) {
		case 1: {
			QE(0,0) = 1.;
			QE(0,1) = 0.;
			QE(0,2) = 0.;
			QE(0,3) = 0.;
			QE(0,4) = 0.;
			QE(0,5) = 0.;
			QE(1,0) = 0.;
			QE(1,1) = c*c;
			QE(1,2) = s*s;
			QE(1,3) = 0.;
			QE(1,4) = 0.;
			QE(1,5) = s*c;
			QE(2,0) = 0.;
			QE(2,1) = s*s;
			QE(2,2) = c*c;
			QE(2,3) = 0.;
			QE(2,4) = 0.;
			QE(2,5) = -s*c;
			QE(3,0) = 0.;
			QE(3,1) = 0.;
			QE(3,2) = 0.;
			QE(3,3) = c;
			QE(3,4) = s;
			QE(3,5) = 0.;
			QE(4,0) = 0.;
			QE(4,1) = 0.;
			QE(4,2) = 0.;
			QE(4,3) = -s;
			QE(4,4) = c;
			QE(4,5) = 0.;
			QE(5,0) = 0.;
			QE(5,1) = -2.*s*c;
			QE(5,2) = 2.*s*c;
			QE(5,3) = 0.;
			QE(5,4) = 0.;
			QE(5,5) = c*c-s*s;
			break;
		}
		case 2: {
			QE(0,0) = c*c;
			QE(0,1) = 0.;
			QE(0,2) = s*s;
			QE(0,3) = 0.;
			QE(0,4) = -c*s;
			QE(0,5) = 0.;
			QE(1,0) = 0.;
			QE(1,1) = 1.;
			QE(1,2) = 0.;
			QE(1,3) = 0.;
			QE(1,4) = 0.;
			QE(1,5) = 0.;
			QE(2,0) = s*s;
			QE(2,1) = 0.;
			QE(2,2) = c*c;
			QE(2,3) = 0.;
			QE(2,4) = c*s;
			QE(2,5) = 0.;
			QE(3,0) = 0.;
			QE(3,1) = 0.;
			QE(3,2) = 0.;
			QE(3,3) = c;
			QE(3,4) = 0.;
			QE(3,5) = -s;
			QE(4,0) = 2.*c*s;
			QE(4,1) = 0.;
			QE(4,2) = -2.*c*s;
			QE(4,3) = 0.;
			QE(4,4) = c*c-s*s;
			QE(4,5) = 0.;
			QE(5,0) = 0.;
			QE(5,1) = 0.;
			QE(5,2) = 0.;
			QE(5,3) = s;
			QE(5,4) = 0.;
			QE(5,5) = c;
			break;
		}
		case 3: {
			QE(0,0) = c*c;
			QE(0,1) = s*s;
			QE(0,2) = 0.;
			QE(0,3) = s*c;
			QE(0,4) = 0.;
			QE(0,5) = 0.;
			QE(1,0) = s*s;
			QE(1,1) = c*c;
			QE(1,2) = 0.;
			QE(1,3) = -s*c;
			QE(1,4) = 0.;
			QE(1,5) = 0.;
			QE(2,0) = 0.;
			QE(2,1) = 0.;
			QE(2,2) = 1.;
			QE(2,3) = 0.;
			QE(2,4) = 0.;
			QE(2,5) = 0.;
			QE(3,0) = -2.*s*c;
			QE(3,1) = 2.*s*c;
			QE(3,2) = 0.;
			QE(3,3) = c*c-s*s;
			QE(3,4) = 0.;
			QE(3,5) = 0.;
			QE(4,0) = 0.;
			QE(4,1) = 0.;
			QE(4,2) = 0.;
			QE(4,3) = 0.;
			QE(4,4) = c;
			QE(4,5) = s;
			QE(5,0) = 0.;
			QE(5,1) = 0.;
			QE(5,2) = 0.;
			QE(5,3) = 0.;
			QE(5,4) = -s;
			QE(5,5) = c;
			break;
		}
		default: {
			cout << "Please choose the axis 1,2 or 3\n";
		}
	}
	return QE;		
}


//To rotate a stiffness matrix (6,6)
mat rotateL(const mat &L, const double &alpha, const int &axis) {

	mat QS = fillQS(alpha, axis);
	return QS*(L*trans(QS));
}

//To rotate a compliance matrix (6,6)
mat rotateM(const mat &M, const double &alpha, const int &axis) {

	mat QE = fillQE(alpha, axis);
	return QE*(M*trans(QE));
}

//To rotate a interaction matrix (6,6)
mat rotateA(const mat &A, const double &alpha, const int &axis) {

	mat QS = fillQS(alpha, axis);
	mat Ainit = zeros(6,6);	
	mat Arot = zeros(6,6);	

	Ainit = A;
	
	for(int i=3; i<6; i++) {
		for(int j=0; j<6; j++) {
			Ainit(i,j)*=0.5;
		}
	}		

	Arot = QS*(Ainit*trans(QS));

	for(int i=3; i<6; i++) {
		for(int j=0; j<6; j++) {
			Arot(i,j)*=2.;
		}
	}

	return Arot;
}

//To rotate a stress vector (6)
vec rotate_stress(const vec &V, const double &alpha, const int &axis) {

	mat QS = fillQS(alpha, axis);
	return QS*V;
}

//To rotate a strain vector (6)
vec rotate_strain(const vec &V, const double &alpha, const int &axis) {

	mat QE = fillQE(alpha, axis);
	return QE*V;
}

//To rotate from local to global a stiffness matrix (6,6)
mat rotate_l2g_L(const mat &Lt, const double &psi, const double &theta, const double &phi) {
    
    mat Lt_temp = Lt;
  	if(fabs(phi) > iota) {
		Lt_temp = rotateL(Lt_temp, -phi, axis_phi);
	}
  	if(fabs(theta) > iota) {
		Lt_temp = rotateL(Lt_temp, -theta, axis_theta);
	}
	if(fabs(psi) > iota) {
		Lt_temp = rotateL(Lt_temp, -psi, axis_psi);
	}
    
	return Lt_temp;
}

//To rotate from global to local a stiffness matrix (6,6)
mat rotate_g2l_L(const mat &Lt, const double &psi, const double &theta, const double &phi) {
    
    mat Lt_temp = Lt;
  	if(fabs(psi) > iota) {
		Lt_temp = rotateL(Lt_temp, psi, axis_psi);
	}
	if(fabs(theta) > iota) {
		Lt_temp = rotateL(Lt_temp, theta, axis_theta);
	}
	if(fabs(phi) > iota) {
		Lt_temp = rotateL(Lt_temp, phi, axis_phi);
    }
    
	return Lt_temp;
}


//To rotate from local to global a localisation matrix (6,6)
mat rotate_l2g_A(const mat &A, const double &psi, const double &theta, const double &phi) {
    
    mat A_temp = A;
  	if(fabs(phi) > iota) {
		A_temp = rotateA(A_temp, -phi, axis_phi);
	}
  	if(fabs(theta) > iota) {
		A_temp = rotateA(A_temp, -theta, axis_theta);
	}
	if(fabs(psi) > iota) {
		A_temp = rotateA(A_temp, -psi, axis_psi);
	}
    
	return A_temp;
}

//To rotate from global to local a localisation matrix (6,6)
mat rotate_g2l_A(const mat &A, const double &psi, const double &theta, const double &phi) {
    
    mat A_temp = A;
  	if(fabs(psi) > iota) {
		A_temp = rotateA(A_temp, psi, axis_psi);
	}
	if(fabs(theta) > iota) {
		A_temp = rotateA(A_temp, theta, axis_theta);
	}
	if(fabs(phi) > iota) {
		A_temp = rotateA(A_temp, phi, axis_phi);
    }
    
	return A_temp;
}

} //namespace smart
