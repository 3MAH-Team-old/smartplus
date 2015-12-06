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

///@file elastic_isotropic.cpp
///@brief User thermomechanical subroutine for Isotropic elastic materials in 3D case
///@version 1.0

#include <iostream>
#include <fstream>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Continuum_Mechanics/constitutive.hpp>

using namespace std;
using namespace arma;

namespace smart{

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_elasticity_iso_T(const vec &Etot, const vec &DEtot, vec &sigma, double &rpl, mat &dSdE, mat &dSdT, mat &drpldE, mat &drpldT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &sse, double &spd, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{  	

    UNUSED(Etot);
    UNUSED(DR);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(statev);
    UNUSED(T);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);    
    
	//From the props to the material properties
	double E = props(0);
	double nu = props(1);
	double alpha = props(2);
	
	// ######################  Elastic compliance and stiffness #################################			
	//defines L
	dSdE = L_iso(E, nu, "Enu");
	
	if(start) { //Initialization
		sigma = zeros(6);
	}	
	
	vec sigma_start = sigma;
	
	//Compute the elastic strain and the related stress	
	vec DEel = DEtot - alpha*Ith()*DT;

	if (ndi == 1) {
        sigma(0) = sigma_start(0) + E*(DEel(0));
    }
	else if (ndi == 2) {
        sigma(0) = sigma_start(0) + E/(1. - (nu*nu))*(DEel(0)) + nu*(DEel(1));
        sigma(1) = sigma_start(1) + E/(1. - (nu*nu))*(DEel(1)) + nu*(DEel(0));
        sigma(3) = sigma_start(3) + E/(1.+nu)*0.5*DEel(3);
    }
    else
        sigma = sigma_start + (dSdE*DEel);
    
    //Returning the energy
    vec Dsigma = sigma - sigma_start;
    
    double Dtde = 0.5*sum((sigma_start+sigma)%DEtot);
    double Dsse = sum(sigma_start%DEel) + 0.5*sum(Dsigma%DEel);
    
    sse += Dsse;
    spd += Dtde - Dsse;
    
    
    //Computation of the Thermal coupling variables
    vec Lth = -1.*(dSdE*Ith())*alpha;
    
    if(DTime < 1.E-12) {
        rpl = 0.;
        dSdT = Lth;
        drpldE = zeros(6);
        drpldT = 0.;
    }
    else {
        rpl = -alpha/DTime*(T + DT)*sum(Ith()%Dsigma);
        dSdT = Lth;
        drpldE = 1./DTime*(T + DT)*Lth;
        drpldT = -1./DTime*alpha*(sum(Ith()%Dsigma)+(T + DT)*sum(Lth%Ith()));
    }
                
}
    
} //namespace smart
