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

vec rotate_vec(const vec &v, const mat &DR) {
    return DR*v;
}

vec rotate_mat(const mat &m, const mat &DR) {
    return trans(DR)*m*DR;
}

mat fillR(const double &alpha, const int &axis) {
    
    mat R = zeros(3,3);
    double c = cos(alpha);
    double s = sin(alpha);
    
    switch(axis) {
        case 1: {
            R = { {1,0,0}, {0,c,-s}, {0,s,c} };
            break;
        }
        case 2: {
            R = mat {{c, 0, s}, {0,1,0}, {-s,0,c}};
            break;
        }
        case 3: {
            R = {{c,-s,0}, {s,c, 0}, {0,0,1}};
            break;
        }
        default: {
            cout << "Please choose the axis 1,2 or 3\n";
        }
    }
    return R;
}
    
mat fillR(const double &psi, const double &theta, const double &phi) {
    
    mat R = zeros(3,3);
    
    mat R1 = fillR(psi, axis_psi);
    mat R2 = fillR(theta, axis_theta);
    mat R3 = fillR(phi, axis_phi);
    
    R = R3*R2*R1;
    return R;
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

mat fillQS(const mat &DR) {
    
    double a = DR(0,0);
    double d = DR(0,1);
    double g = DR(0,2);
    double b = DR(1,0);
    double e = DR(1,1);
    double h = DR(1,2);
    double c = DR(2,0);
    double f = DR(2,1);
    double i = DR(2,2);
    
    mat QS= zeros(6,6);
    QS(0,0) = a*a;
    QS(0,1) = d*d;
    QS(0,2) = g*g;
    QS(0,3) = 2.*a*d;
    QS(0,4) = 2.*a*g;
    QS(0,5) = 2.*d*g;
    QS(1,0) = b*b;
    QS(1,1) = e*e;
    QS(1,2) = h*h;
    QS(1,3) = 2.*b*e;
    QS(1,4) = 2.*b*h;
    QS(1,5) = 2.*e*h;
    QS(2,0) = c*c;
    QS(2,1) = f*f;
    QS(2,2) = i*i;
    QS(2,3) = 2.*c*f;
    QS(2,4) = 2.*c*i;
    QS(2,5) = 2.*f*i;
    QS(3,0) = a*b;
    QS(3,1) = d*e;
    QS(3,2) = g*h;
    QS(3,3) = b*d+a*e;
    QS(3,4) = b*g+a*h;
    QS(3,5) = e*g+d*h;
    QS(4,0) = a*c;
    QS(4,1) = d*f;
    QS(4,2) = g*i;
    QS(4,3) = c*d+a*f;
    QS(4,4) = c*g+a*i;
    QS(4,5) = f*g+d*i;
    QS(5,0) = b*c;
    QS(5,1) = e*f;
    QS(5,2) = h*i;
    QS(5,3) = c*e+b*f;
    QS(5,4) = c*h+b*i;
    QS(5,5) = f*h+e*i;
    
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
    
mat fillQE(const mat &DR) {
    
    double a = DR(0,0);
    double d = DR(0,1);
    double g = DR(0,2);
    double b = DR(1,0);
    double e = DR(1,1);
    double h = DR(1,2);
    double c = DR(2,0);
    double f = DR(2,1);
    double i = DR(2,2);
    
    mat QE= zeros(6,6);
    QE(0,0) = a*a;
    QE(0,1) = d*d;
    QE(0,2) = g*g;
    QE(0,3) = a*d;
    QE(0,4) = a*g;
    QE(0,5) = d*g;
    QE(1,0) = b*b;
    QE(1,1) = e*e;
    QE(1,2) = h*h;
    QE(1,3) = b*e;
    QE(1,4) = b*h;
    QE(1,5) = e*h;
    QE(2,0) = c*c;
    QE(2,1) = f*f;
    QE(2,2) = i*i;
    QE(2,3) = c*f;
    QE(2,4) = c*i;
    QE(2,5) = f*i;
    QE(3,0) = 2.*a*b;
    QE(3,1) = 2.*d*e;
    QE(3,2) = 2.*g*h;
    QE(3,3) = b*d+a*e;
    QE(3,4) = b*g+a*h;
    QE(3,5) = e*g+d*h;
    QE(4,0) = 2.*a*c;
    QE(4,1) = 2.*d*f;
    QE(4,2) = 2.*g*i;
    QE(4,3) = c*d+a*f;
    QE(4,4) = c*g+a*i;
    QE(4,5) = f*g+d*i;
    QE(5,0) = 2.*b*c;
    QE(5,1) = 2.*e*f;
    QE(5,2) = 2.*h*i;
    QE(5,3) = c*e+b*f;
    QE(5,4) = c*h+b*i;
    QE(5,5) = f*h+e*i;
    
    return QE;
}


//To rotate a stiffness matrix (6,6)
mat rotateL(const mat &L, const double &alpha, const int &axis) {

	mat QS = fillQS(alpha, axis);
	return QS*(L*trans(QS));
}

mat rotateL(const mat &L, const mat &DR) {
    
    mat QS = fillQS(DR);
    return QS*(L*trans(QS));
}
    
    
//To rotate a compliance matrix (6,6)
mat rotateM(const mat &M, const double &alpha, const int &axis) {

	mat QE = fillQE(alpha, axis);
	return QE*(M*trans(QE));
}

mat rotateM(const mat &M, const mat &DR) {
    
    mat QE = fillQE(DR);
    return QE*(M*trans(QE));
}
    
    
//To rotate an interaction matrix of type A (strain) (6,6)
mat rotateA(const mat &A, const double &alpha, const int &axis) {
        
    mat QE = fillQE(alpha, axis);
    mat QS = fillQS(alpha, axis);
    
    return QE*(A*trans(QS));
}

mat rotateA(const mat &A, const mat &DR) {
    
    mat QE = fillQE(DR);
    mat QS = fillQS(DR);
    
    return QE*(A*trans(QS));
}
    
//To rotate an interaction matrix of type B (stress) (6,6)
mat rotateB(const mat &B, const double &alpha, const int &axis) {
    
    mat QE = fillQE(alpha, axis);
    mat QS = fillQS(alpha, axis);
    
    return QS*(B*trans(QE));
}

mat rotateB(const mat &B, const mat &DR) {
    
    mat QE = fillQE(DR);
    mat QS = fillQS(DR);
    
    return QS*(B*trans(QE));
}
    
//To rotate a stress vector (6)
vec rotate_stress(const vec &V, const double &alpha, const int &axis) {

	mat QS = fillQS(alpha, axis);
	return QS*V;
}

//To rotate a stress vector (6)
vec rotate_stress(const vec &V, const mat &DR) {
    
    mat QS = fillQS(DR);
    return QS*V;
}
    
//To rotate a strain vector (6)
vec rotate_strain(const vec &V, const double &alpha, const int &axis) {

	mat QE = fillQE(alpha, axis);
	return QE*V;
}

vec rotate_strain(const vec &V, const mat &DR) {
    
    mat QE = fillQE(DR);
    return QE*V;
}
    

//To rotate from local to global a strain tensor (6) using Euler angles
mat rotate_l2g_strain(const vec &E, const double &psi, const double &theta, const double &phi) {
    
    mat E_temp = E;
    if(fabs(phi) > iota) {
        E_temp = rotate_strain(E_temp, -phi, axis_phi);
    }
    if(fabs(theta) > iota) {
        E_temp = rotate_strain(E_temp, -theta, axis_theta);
    }
    if(fabs(psi) > iota) {
        E_temp = rotate_strain(E_temp, -psi, axis_psi);
    }
    
    return E_temp;
}
    
//To rotate from global to local a strain tensor (6) using Euler angles
mat rotate_g2l_strain(const vec &E, const double &psi, const double &theta, const double &phi) {
    
    mat E_temp = E;
    if(fabs(psi) > iota) {
        E_temp = rotate_strain(E_temp, psi, axis_psi);
    }
    if(fabs(theta) > iota) {
        E_temp = rotate_strain(E_temp, theta, axis_theta);
    }
    if(fabs(phi) > iota) {
        E_temp = rotate_strain(E_temp, phi, axis_phi);
    }
    
    return E_temp;
}

//To rotate from local to global a stress tensor (6)
mat rotate_l2g_stress(const vec &S, const double &psi, const double &theta, const double &phi) {
    
    mat S_temp = S;
    if(fabs(phi) > iota) {
        S_temp = rotate_stress(S_temp, -phi, axis_phi);
    }
    if(fabs(theta) > iota) {
        S_temp = rotate_stress(S_temp, -theta, axis_theta);
    }
    if(fabs(psi) > iota) {
        S_temp = rotate_stress(S_temp, -psi, axis_psi);
    }
    
    return S_temp;
}

//To rotate from global to local a stress tensor (6)
mat rotate_g2l_stress(const vec &S, const double &psi, const double &theta, const double &phi) {
    
    mat S_temp = S;
    if(fabs(psi) > iota) {
        S_temp = rotate_stress(S_temp, psi, axis_psi);
    }
    if(fabs(theta) > iota) {
        S_temp = rotate_stress(S_temp, theta, axis_theta);
    }
    if(fabs(phi) > iota) {
        S_temp = rotate_stress(S_temp, phi, axis_phi);
    }
    
    return S_temp;
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

//To rotate from local to global a  localisation matrix (6,6)
mat rotate_l2g_M(const mat &M, const double &psi, const double &theta, const double &phi) {
    
    mat M_temp = M;
    if(fabs(phi) > iota) {
        M_temp = rotateM(M_temp, -phi, axis_phi);
    }
    if(fabs(theta) > iota) {
        M_temp = rotateM(M_temp, -theta, axis_theta);
    }
    if(fabs(psi) > iota) {
        M_temp = rotateM(M_temp, -psi, axis_psi);
    }
    
    return M_temp;
}

//To rotate from global to local a localisation matrix (6,6)
mat rotate_g2l_M(const mat &M, const double &psi, const double &theta, const double &phi) {
    
    mat M_temp = M;
    if(fabs(psi) > iota) {
        M_temp = rotateM(M_temp, psi, axis_psi);
    }
    if(fabs(theta) > iota) {
        M_temp = rotateM(M_temp, theta, axis_theta);
    }
    if(fabs(phi) > iota) {
        M_temp = rotateM(M_temp, phi, axis_phi);
    }
    
    return M_temp;
}

//To rotate from local to global a strain localisation matrix (6,6)
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

//To rotate from global to local a strain localisation matrix (6,6)
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

//To rotate from local to global a stress localisation matrix (6,6)
mat rotate_l2g_B(const mat &B, const double &psi, const double &theta, const double &phi) {
    
    mat B_temp = B;
  	if(fabs(phi) > iota) {
		B_temp = rotateB(B_temp, -phi, axis_phi);
	}
  	if(fabs(theta) > iota) {
		B_temp = rotateB(B_temp, -theta, axis_theta);
	}
	if(fabs(psi) > iota) {
		B_temp = rotateB(B_temp, -psi, axis_psi);
	}
    
	return B_temp;
}

//To rotate from global to local a stress localisation matrix (6,6)
mat rotate_g2l_B(const mat &B, const double &psi, const double &theta, const double &phi) {
    
    mat B_temp = B;
  	if(fabs(psi) > iota) {
		B_temp = rotateB(B_temp, psi, axis_psi);
	}
	if(fabs(theta) > iota) {
		B_temp = rotateB(B_temp, theta, axis_theta);
	}
	if(fabs(phi) > iota) {
		B_temp = rotateB(B_temp, phi, axis_phi);
    }
    
	return B_temp;
}
    
} //namespace smart
