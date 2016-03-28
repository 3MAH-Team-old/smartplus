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

///@file criteria.cpp
///@brief Provide function for yield surfaces
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <smartplus/Libraries/Continuum_Mechanics/contimech.hpp>
#include <smartplus/Libraries/Continuum_Mechanics/constitutive.hpp>

using namespace std;
using namespace arma;

namespace smart{
    
//This function returns the Prager equivalent stress.
double Prager_stress(const vec &v, const double &b, const double &n)
{
     assert(v.size() == 6);
     assert(b >= 0.);
     assert(n > 0.);
     double temp;
     double m;
    
     if (Mises_stress(v) > 0.) {
         if (n < 10.)
         {
             m = 1. / n;
             temp = Mises_stress(v) * pow((1. + b * J3_stress(v) / pow(J2_stress(v), 1.5)), m);
         }
         else {
             m = 0.;
             temp = Mises_stress(v);
         }
     }
     else {
         m = 0.;
         temp = 0.;
     }
     return temp;
}

//This function returns the derivative of the Prager equivalent stress.
vec dPrager_stress(const vec &v, const double &b, const double &n)
{
     assert(v.size() == 6);
     assert(b >= 0.);
     assert(n > 0.);
     vec vdev = dev(v);
     mat devstress_t = v2t_stress(vdev);
     vec square_stressdev = t2v_stress(devstress_t * devstress_t);
     double m;
     
     vec temp;
    
     if (Mises_stress(v) > 0.)
     {
         if (n < 10.) {
             m = 1. / n;
             temp = sqrt(3.) * pow((1. + b * J3_stress(v) / pow(J2_stress(v), 1.5)), (m - 1.)) * (0.5 / sqrt(J2_stress(v)) * vdev + b * m / (6. * pow(J2_stress(v), 2.)) * (6. * J2_stress(v) * square_stressdev - 4. * pow(J2_stress(v), 2.) * Ith() + (3. / m - 9.) * J3_stress(v) * vdev));
             
             for (int i = 3; i < 6; i++)
             {
             temp(i) = 2. * temp(i);
             }
         }
         else {
             m = 0.;
             temp = eta_stress(v);
         }         
     }
     else {
         m = 0.;
         temp = zeros(6);
         
     }
     return temp;

}

//This function returns the Prager equivalent stress.
double Tresca_stress(const vec &v)
{
    mat sigma = v2t_stress(v);
    vec lambda = eig_sym(sigma);
    
    return lambda(0) - lambda(2);
}

//This function returns the Prager equivalent stress.
vec dTresca_stress(const vec &v)
{
    return eta_stress(v);
}
    
    
double Eq_stress(const vec &v, const string &eq_type, const vec &param)
{
    if(eq_type == "Mises") {
        return Mises_stress(v);
    }
    else if(eq_type == "Tresca") {
        return Tresca_stress(v);
    }
    else if(eq_type == "Prager") {
        return Prager_stress(v, param(0), param(1));
    }
    else {
        cout << "Error in Eq_stress : No valid arguement is given\n";
        exit(0);
    }
}

vec dEq_stress(const vec &v, const string &eq_type, const vec &param)
{
    if(eq_type == "Mises") {
        return eta_stress(v);
    }
    else if(eq_type == "Tresca") {
        return dTresca_stress(v);
    }
    else if(eq_type == "Prager") {
        return dPrager_stress(v, param(0), param(1));
    }
    else {
        cout << "Error in dEq_stress : No valid arguement is given\n";
        exit(0);
    }
}
    
} //namespace smart
