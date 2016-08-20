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

///@file plastic_isotropic_ccp.cpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief This subroutines uses a convex cutting plane algorithm
///@brief Isotropic hardening with a power-law hardenig is considered
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

///@brief The elastic-plastic UMAT with isotropic hardening requires 8 statev:
///@brief statev[0] : T_init : Initial temperature
///@brief statev[1] : Accumulative plastic parameter: p
///@brief statev[2] : Plastic strain 11: EP(0,0)
///@brief statev[3] : Plastic strain 22: EP(1,1)
///@brief statev[4] : Plastic strain 33: EP(2,2)
///@brief statev[5] : Plastic strain 12: EP(0,1) (*2)
///@brief statev[6] : Plastic strain 13: EP(0,2) (*2)
///@brief statev[7] : Plastic strain 23: EP(1,2) (*2)

void umat_plasticity_iso_CCP(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);    
    
	//From the props to the material properties
	double E = props(0);
	double nu= props(1);
	double alpha = props(2);
	double sigmaY = props(3);
	double k=props(4);
	double m=props(5);
	
	// ######################  Elastic compliance and stiffness #################################			
	//defines L
	mat L = L_iso(E, nu, "Enu");
	mat M = M_iso(E, nu, "Enu");
    
    ///@brief Temperature initialization
    double Tinit = statev(0);
    //From the statev to the internal variables
    double p = statev(1);
    vec EP(6);
    EP(0) = statev(2);
    EP(1) = statev(3);
    EP(2) = statev(4);
    EP(3) = statev(5);
    EP(4) = statev(6);
    EP(5) = statev(7);
    
    //Rotation of internal variables (tensors)
    EP = rotate_strain(EP, DR);
    
	///@brief Initialization
	if(start)
	{
		Tinit = T;
		vec vide = zeros(6);
		sigma = vide;
		EP = vide;
		p = 0.;
	}
    vec sigma_start = sigma;
    vec EP_start = EP;
    
    double p_temp = p;
    double p_start = p;
    double h = 0.;
	
	///Elastic prediction - Accounting for the thermal prediction	
	
	vec Eel_start = Etot - (alpha*Ith()*T -  alpha*Ith()*Tinit) - EP;
	vec DEel_start = DEtot - alpha*Ith()*DT;

	if (ndi == 1) {
        sigma(0) = sigma_start(0) + E*(DEel_start(0));
    }
	else if (ndi == 2) {
        sigma(0) = sigma_start(0) + E/(1. - (nu*nu))*(DEel_start(0)) + nu*(DEel_start(1));
        sigma(1) = sigma_start(1) + E/(1. - (nu*nu))*(DEel_start(1)) + nu*(DEel_start(0));
        sigma(3) = sigma_start(3) + E/(1.+nu)*0.5*DEel_start(3);
    }
    else
        sigma = sigma_start + (L*DEel_start);
        
    //Compute the explicit flow direction
	vec Lambdap = eta_stress(sigma);
	
	//Define the plastic function and the stress
	double Phi=0.;
	double dPhi = 0.;
	
	double dPhidp=0.;
	vec dPhidsigma = zeros(6);
	vec dPhidEP = zeros(6);

	double Hp=0.;
	double dHpdp=0.;

	vec dEP(6);
	vec dsigmadp = zeros(6);
	
	if (p > 0.)
        Hp = k*pow(p, m);
	else
        Hp = 0.;
    double Hp_start = Hp;
    
    //Initialization of the function that defines the elastic domain
	Phi = Mises_stress(sigma) - Hp - sigmaY;
	vec Eel = zeros(6);
	int compteur = 0;
    
    if(fabs(Phi)/sigmaY > precision_umat) {
		for(compteur = 0 ; (compteur < maxiter_umat) && (fabs(Phi)/sigmaY > precision_umat) ; compteur++) {
	        
	        Lambdap = eta_stress(sigma);
	        dPhidsigma = eta_stress(sigma);
				
			Hp = k*pow(p, m);
			Phi = Mises_stress(sigma) - Hp - sigmaY;
			
			if (p > 1.E-12)	{		
				dHpdp = m*k*pow(p, m-1);
				Hp = k*pow(p, m);
			}
			else {
				dHpdp = 0.;
				Hp = 0.;
			}
			
			dPhidp = -dHpdp;
			dPhi = sum((-1.*L*Lambdap)%dPhidsigma) + dPhidp;
			
			p_temp = p;
			
			if(fabs(dPhi) > 0.)
			    p = p - 1.*Phi/dPhi;
			else {
			    //pnewdt = 0.1;
			    Phi = 0.;
			}
			
			if (p < p_start) {
				//Phi = 0.;
                p = p_start+precision_umat;
			}		
			
			EP = EP + (p - p_temp)*Lambdap;
			
			//the stress is now computed using the relationship sigma = L(E-Ep)
			Eel = Etot + DEtot - alpha*Ith()*(T + DT - Tinit) - EP;
			
			if(ndi == 1) {						// 1D
				double Eeff = 1./M(0,0);
	            sigma(0) = Eeff*(Eel(0));
	        }
	        else if(ndi == 2){					// 2D Plane Stress
				double Eeff = 1./M(0,0);
	            double nueff = -M(0,1)/M(0,0);
	            sigma(0) = Eeff/(1. - nueff*nueff)*(Eel(0)) + nueff*(Eel(1));
	            sigma(1) = Eeff/(1. - nueff*nueff)*(Eel(1)) + nueff*(Eel(0));
	            sigma(3) = Eeff/(1. - nueff*nueff)*(1. - nueff)*0.5*Eel(3);
	        }
	        else
				sigma = (L*Eel); // 2D Generalized Plane Strain (Plane Strain, Axisymetric) && 3D
			
		}
			
		// Computation of the elastoplastic tangent moduli
        if (p > precision_umat)
            dHpdp = m*k*pow(p, m-1);
        else
            dHpdp =  0.;
		
		double mu = E/(2.*(1+nu));
		h = 3.*mu+dHpdp;

        if (p-p_start > 1.1*precision_umat) {
            Lambdap = Lambdap%Ir05();
            Lt = (4.*pow(mu, 2.)/h)*(Lambdap*trans(Lambdap));
            Lt = L - Lt;
        }
        else
            Lt = L;
	}
	else {
		Eel = Etot + DEtot - alpha*Ith()*(T + DT - Tinit) - EP;
		
		if(ndi == 1){						// 1D
			double Eeff = 1./M(0,0);
			sigma(0) = Eeff*(Eel(0));
		}
		else if(ndi == 2){					// 2D Plane Stress
			double Eeff = 1./M(0,0);
			double nueff = -M(0,1)/M(0,0);
			sigma(0) = Eeff/(1. - nueff*nueff)*(Eel(0)) + nueff*(Eel(1));
			sigma(1) = Eeff/(1. - nueff*nueff)*(Eel(1)) + nueff*(Eel(0));
			sigma(3) = Eeff/(1. - nueff*nueff)*(1. - nueff)*0.5*Eel(3);
		}
		else
			sigma = (L*Eel); // 2D Generalized Plane Strain (Plane Strain, Axisymetric) && 3D
		
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
    
    //Computation of the increments of variables
    vec Dsigma = sigma - sigma_start;
    vec DEP = EP - EP_start;
    double Dp = p-p_start;
    
    double gamma_loc = 0.5*sum((sigma_start+sigma)%DEP) - 0.5*(Hp_start + Hp)*Dp;
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DEP));
    Wm_ir += 0.5*(Hp_start + Hp)*Dp;
    Wm_d += gamma_loc;
}
    
} //namespace smart
