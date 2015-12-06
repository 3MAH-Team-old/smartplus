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

///@file stats.cpp
///@brief Usefull statistical functions
///@version 1.0

#include <math.h>
#include <assert.h>
#include <armadillo>
#include <smartplus/parameter.hpp>

using namespace std;
using namespace arma;

namespace smart{
    
///Approximation of a normal distribution
double normal_distrib(const double &x, const double &mean, const double &dev){
	
	double x_norm = (x-mean)/dev;
	
	if (x_norm < 0.){
		return (1. - normal_distrib(-1.*x_norm, 0., 1.));
	}
	else {
		double k = 1.0/(1.0 + 0.2316419*x_norm);
		double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
		return (1.0 - (1.0/(sqrt(2*pi)))*exp(-0.5*x_norm*x_norm) * k_sum);
	}	
}

//tri_sum of a and b
int tri_sum(const int &a, const int &b) {
    assert(a>0);
    assert(b>0);
    
    return b*(2*a+(b-1))/2;
}

///1. Classic ODF: a1 * cos(Theta)^(2*p1) + a2 * cos(Theta)^(2*p2 + 1) * sin(Theta)^(2*p2) 
double ODF_sd(const double& Theta, const double& m, const double& alpha1, const double& alpha2, const double& pow1, const double& pow2){
	  
	if (fabs(Theta - m) < 1.E-6)
		return alpha1;  
	else if (fabs((Theta - m)-0.5*pi) < 1.E-6)
		return 0.;
	else
		return fabs((alpha1*pow(cos(Theta - m),2.*pow1) + alpha2*pow(cos(Theta - m),2.*pow2)*pow(sin(Theta - m),2.*pow2))*cos(Theta - m));
}

///2. Classic ODF - hardening-like
double ODF_hard(const double& Theta, const double& m, const double& k, const double& stdev) {
    assert(k>0);
    assert(stdev>0);
	  
	return k*exp(-0.5*pow(fabs(Theta - m)/stdev,2.));
}

///3. Gaussian
double Gaussian(const double& X, const double& mean, const double& sd, const double& ampl){
    assert(sd>0);
	
	return ampl/(sd * sqrt(2.*pi)) * exp(-1./2. * pow((X - mean)/sd, 2.));
}

///30. Several Gaussian
double Mult_Gaussian(const double& X, const int& Nb_peak, const vec& Mean, const vec& SD, const vec& Ampl){
	double result = 0.;
	
	for (int i = 0; i < Nb_peak; i++) {
		assert(SD(i)>0);
		result += Gaussian(X, Mean(i), SD(i), Ampl(i));
	}
	
	return result;
}

///4. Lorentzian
double Lorentzian(const double& X, const double& mean, const double& width, const double& ampl){
    assert(width>0);
	
	return ampl * width/(2.*pi*(pow(X - mean, 2.)+ pow(width / 2., 2.)));
}


///40. Several Lorentzian
double Mult_Lorentzian(const double& X, const int& Nb_peak, const vec& Mean, const vec& Width, const vec& Ampl){
	double result = 0.;
	
	for (int i = 0; i < Nb_peak; i++) {
		assert(Width(i)>0);
		result += Lorentzian(X, Mean(i), Width(i), Ampl(i));
	}
	return result;
}

///5. Pseudo-Voigt
double PseudoVoigt(const double& X, const double& eta, const double& mean, const double &width_Lor, const double& sd_Gau, const double& ampl){
    assert(width_Lor>0);
    assert(sd_Gau>0);
	  
	return eta * Lorentzian(X, mean, width_Lor, ampl) + (1.-eta) * Gaussian(X, mean, sd_Gau, ampl);
}

///50. Several Pseudo-Voigt
double Mult_PseudoVoigt(const double& X, const int& Nb_peak, const vec& Eta, const vec& Mean, const vec &Width_Lor, const vec& SD_Gau, const vec& Ampl){
	double result = 0.;
	  
	for (int i = 0; i < Nb_peak; i++) {
		assert(Width_Lor(i)>0);
		assert(SD_Gau(i)>0);
		result += PseudoVoigt(X, Eta(i), Mean(i), Width_Lor(i), SD_Gau(i), Ampl(i));
	}
	return result;
}

///6. Pearson VII
double Pearson7(const double& X, const double& mean, const double &inv_width, const double& shape, const double& max){
	assert(shape>0);
	      	
	return max * pow(1. + pow(inv_width*(X - mean), 2.)/shape, -1.*shape);
}

///60. Several Pearson VII
double Mult_Pearson7(const double& X, const int& Nb_peak, const vec& Mean, const vec &Inv_Width, const vec& Shape, const vec& Max){
	double result = 0.;
	  	
	for (int i = 0; i < Nb_peak; i++) {
		assert(Shape(i)>0);
		result += Pearson7(X, Mean(i), Inv_Width(i), Shape(i), Max(i));
	}
	return result;
}

} //namespace smart
