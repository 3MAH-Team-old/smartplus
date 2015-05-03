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

///@file constitutive.hpp
///@brief Constitutive tensors in Voigt notation
///@version 1.0

#pragma once
#include <string.h>
#include <armadillo>

using namespace std;
using namespace arma;

namespace smart{

//Returns the fourth order identity tensor written in Voigt notation Ireal
mat Ireal();

//Returns the volumic of the identity tensor Ireal written in Voigt notation
mat Ivol();

//Returns the deviatoric of the identity tensor Ireal written in Voigt notation
mat Idev();

//Returns the fourth order identity tensor Iˆ written in Voigt notation
mat Ireal2();

//Returns the deviatoric of the identity tensor Iˆ written in Voigt notation
mat Idev2();

//Returns the expansion vector
vec Ith();

//Returns the stress 2 strain operator
vec Ir2();

//Returns the strain 2 stress operator
vec Ir05();

//Provides the elastic stiffness tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
// ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
mat L_iso(const double &, const double &, string);

//Provides the elastic compliance tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
//‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
mat M_iso(const double &, const double &, string);

//Returns the elastic stiffness tensor for a cubic material.
//Arguments are the stiffness coefficients C11, C12 and C44.
mat L_cubic(const double &, const double &, const double &);

//Returns the elastic compliance tensor for an isotropic material.
//Arguments are the stiffness coefficients C11, C12 and C44.
mat M_cubic(const double &, const double &, const double &);

//Returns the elastic stiffness tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
mat L_ortho(const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, string);

//Returns the elastic compliance tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
mat M_ortho(const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, string);

//Returns the elastic stiffness tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
mat L_isotrans(const double &, const double &, const double &, const double &, const double &, const int &);

//Returns the elastic compliance tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
mat M_isotrans(const double &, const double &, const double &, const double &, const double &, const int &);

//Provides the viscous tensor H an isotropic material.
//The two first arguments are a couple of viscous coefficients (the first is bulk, the second is shear).
mat H_iso(const double &, const double &);

} //namespace smart
