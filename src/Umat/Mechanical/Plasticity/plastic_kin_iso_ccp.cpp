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

///@file plastic_kin_iso_ccp.cpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief This subroutines uses a convex cutting plane algorithm
///@brief Linear Kinematical hardening coupled with a power-law hardenig is considered
///@version 1.0

#include <iostream>
#include <fstream>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Continuum_Mechanics/contimech.hpp>
#include <smartplus/Libraries/Continuum_Mechanics/constitutive.hpp>
#include <smartplus/Libraries/Maths/rotation.hpp>

using namespace std;
using namespace arma;

namespace smart {

///@brief The elastic-plastic UMAT with isotropic hardening requires 6 constants:

///@brief props[1] : Young modulus
///@brief props[2] : Poisson ratio
///@brief props[3] : CTE
///@brief props[4] : J2 equivalent yield stress limit : sigmaY
///@brief props[5] : hardening parameter k
///@brief props[6] : exponent m
///@brief props[7] : linear kinematical hardening h

///@brief The elastic-plastic UMAT with isotropic hardening requires 8 statev:
///@brief statev[0] : T_init : Initial temperature
///@brief statev[1] : Accumulative plastic parameter: p
///@brief statev[2] : Plastic strain 11: EP(0,0)
///@brief statev[3] : Plastic strain 22: EP(1,1)
///@brief statev[4] : Plastic strain 33: EP(2,2)
///@brief statev[5] : Plastic strain 12: EP(0,1) (*2)
///@brief statev[6] : Plastic strain 13: EP(0,2) (*2)
///@brief statev[7] : Plastic strain 23: EP(1,2) (*2)
///@brief statev[8] : Backstress 11: X(0,0)
///@brief statev[8] : Backstress 11: X(1,1)
///@brief statev[8] : Backstress 11: X(2,2)
///@brief statev[8] : Backstress 11: X(0,1)
///@brief statev[8] : Backstress 11: X(0,2)
///@brief statev[8] : Backstress 11: X(1,2)
    

void umat_plasticity_kin_iso_CCP(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &sse, double &spd, const int &ndi, const int &nshr, const bool &start)
{
	vec sigma_start = sigma;

	///@brief Temperature initialization
	double Tinit = statev(0);

	//From the statev to the internal variables
	double p = statev(1);
	double p_temp = p;
	double p_start = p;
        
	vec EP = zeros(6);
	EP(0) = statev(2);
	EP(1) = statev(3);
	EP(2) = statev(4);
	EP(3) = statev(5);
	EP(4) = statev(6);
	EP(5) = statev(7);

    ///@brief X is the backstress associated with reorientation
    vec X = zeros(6);
    X(0) = statev(8);
    X(1) = statev(9);
    X(2) = statev(10);
    X(3) = statev(11);
    X(4) = statev(12);
    X(5) = statev(13);
    
	double h = 0.;

	//Rotation of internal variables (tensors)
	Rot_strain(EP, DR);
	Rot_stress(X, DR);

	//From the props to the material properties
	double E = props(0);
	double nu= props(1);
	double alpha = props(2);
	double sigmaY = props(3);
	double k = props(4);
	double m = props(5);
    double kX = props(6);
	
	// ######################  Elastic compliance and stiffness #################################			
	//defines L
	mat L = L_iso(E, nu, "Enu");
	mat M = M_iso(E, nu, "Enu");
    
	///@brief Initialization
	if(start)
	{
		Tinit = T;
		vec vide = zeros(6);
		sigma = vide;
		sigma_start = vide;
        X = vide;
		EP = vide;
		p = 0.;
		p_temp = 0.;
		p_start = 0.;
		h = 0.;
	}
	
	///Elastic prediction - Accounting for the thermal prediction	
	
	vec Eelstart = Etot - (alpha*Ith()*T -  alpha*Ith()*Tinit) - EP;
	vec DEtot2 = DEtot - alpha*Ith()*DT;

	if (ndi == 1) {
        sigma(0) = sigma_start(0) + E*(DEtot2(0));
    }
	else if (ndi == 2) {
        sigma(0) = sigma_start(0) + E/(1. - (nu*nu))*(DEtot2(0)) + nu*(DEtot2(1));
        sigma(1) = sigma_start(1) + E/(1. - (nu*nu))*(DEtot2(1)) + nu*(DEtot2(0));
        sigma(3) = sigma_start(3) + E/(1.+nu)*0.5*DEtot2(3);
    }
    else
        sigma = sigma_start + (L*DEtot2);
        
    //Compute the explicit flow direction
	vec Lambdap = eta_stress(sigma - X);
	vec lambdaX = -1.*eta_stress(sigma-X);
	
	//Define the plastic function and the stress
	double Phi=0.;
	double dPhi = 0.;
	
	double dPhidp=0.;
	vec dPhidsigma = zeros(6);
	vec dPhidEP = zeros(6);
	vec dPhidX = zeros(6);
    
	double Hp=0.;
	double dHpdp=0.;

	vec dEP(6);
	
	vec dsigmadp = zeros(6);
	
	if (p > 0.)
        Hp = k*pow(p, m);
	else
        Hp = 0.;
    
    //Initialization of the function that defines the elastic domain
	Phi = Mises_stress(sigma - X) - Hp - sigmaY;
	vec diff = zeros(6);
	int compteur = 0;
    
    if(fabs(Phi)/sigmaY > precision_umat) {
		for(compteur = 0 ; (compteur < maxiter_umat) && (fabs(Phi)/sigmaY > precision_umat) ; compteur++) {
	        
	        Lambdap = eta_stress(sigma - X);            
            lambdaX = kX*(eta_stress(sigma - X)%Ir05());
	        dPhidsigma = eta_stress(sigma - X);
            
			Hp = k*pow(p, m);
			Phi = Mises_stress(sigma - X) - Hp - sigmaY;
			
			if (p > 1.E-12)	{		
				dHpdp = m*k*pow(p, m-1);
				Hp = k*pow(p, m);
			}
			else {
				dHpdp = 0.;
				Hp = 0.;
			}
			
			dPhidp = -dHpdp;
            dPhidX = -1.*eta_stress(sigma - X);
            dPhi = sum((-1.*L*Lambdap)%dPhidsigma) + dPhidp + sum(dPhidX%lambdaX);
			
			p_temp = p;
			
			if(fabs(dPhi) > 0.)
			    p = p - 1.*Phi/dPhi;
			else {
			    //pnewdt = 0.1;
			    Phi = 0.;
			}
			
			if (p < p_start) {
				//Phi = 0.;
				p = p_start+1E-6;
			}		
			
			EP = EP + (p - p_temp)*Lambdap;
            X = X + (p - p_temp)*lambdaX;
            
			//the stress is now computed using the relationship sigma = L(E-Ep)
			diff = Etot + DEtot - alpha*Ith()*(T + DT - Tinit) - EP;
			
			if(ndi == 1) {						// 1D
				double Eeff = 1./M(0,0);
	            sigma(0) = Eeff*(diff(0));
	        }
	        else if(ndi == 2){					// 2D Plane Stress
				double Eeff = 1./M(0,0);
	            double nueff = -M(0,1)/M(0,0);
	            sigma(0) = Eeff/(1. - nueff*nueff)*(diff(0)) + nueff*(diff(1));
	            sigma(1) = Eeff/(1. - nueff*nueff)*(diff(1)) + nueff*(diff(0));
	            sigma(3) = Eeff/(1. - nueff*nueff)*(1. - nueff)*0.5*diff(3);
	        }
	        else
				sigma = (L*diff); // 2D Generalized Plane Strain (Plane Strain, Axisymetric) && 3D			
			
		}
			
		// Computation of the elastoplastic tangent moduli
        if (p > 1E-12)
            dHpdp = m*k*pow(p, m-1);
        else
            dHpdp =  0.;

        dPhidp = -dHpdp;
        dPhidX = -1.*eta_stress(sigma - X);
        dPhi = dPhidp + sum(dPhidX%lambdaX);
        dPhidsigma = eta_stress(sigma - X);
        
        //Computation of the tangent modulus !
        vec B1(6);
        vec B2(6);
        double At2;
        
        B1 = L*Lambdap;
        B2 = L*dPhidsigma;
        At2 = (dPhi - sum(dPhidsigma%B1));
        
        Lt = L + (B1*trans(B2))/At2;
	}
	else {
		diff = Etot + DEtot - alpha*Ith()*(T + DT - Tinit) - EP;
		
		if(ndi == 1){						// 1D
			double Eeff = 1./M(0,0);
			sigma(0) = Eeff*(diff(0));
		}
		else if(ndi == 2){					// 2D Plane Stress
			double Eeff = 1./M(0,0);
			double nueff = -M(0,1)/M(0,0);
			sigma(0) = Eeff/(1. - nueff*nueff)*(diff(0)) + nueff*(diff(1));
			sigma(1) = Eeff/(1. - nueff*nueff)*(diff(1)) + nueff*(diff(0));
			sigma(3) = Eeff/(1. - nueff*nueff)*(1. - nueff)*0.5*diff(3);
		}
		else
			sigma = (L*diff); // 2D Generalized Plane Strain (Plane Strain, Axisymetric) && 3D
		
		Lt = L;
	}
	
	///@brief statev evolving variables
	//statev
	statev(0) = Tinit;
	statev(1) = p;
    
	statev(2) = EP(0);
	statev(3) = EP(1);
	statev(4) = EP(2);
	statev(5) = EP(3);
	statev(6) = EP(4);
	statev(7) = EP(5);
    
    statev(8) = X(0);
    statev(9) = X(1);
    statev(10) = X(2);
    statev(11) = X(3);
    statev(12) = X(4);
    statev(13) = X(5);
    
	//Returning the energy
	vec Eel = Etot + DEtot - (alpha*Ith()*(T + DT) -  alpha*Ith()*Tinit) - EP;
	vec DEel = Eel - Eelstart;
	vec Dsigma = sigma - sigma_start;
    
	double Dtde = 0.5*sum((sigma_start+sigma)%DEtot);
	double Dsse = sum(sigma_start%DEel) + 0.5*sum(Dsigma%DEel);
    
	sse += Dsse;
	spd += Dtde - Dsse;	
}
    
} //namespace smart
