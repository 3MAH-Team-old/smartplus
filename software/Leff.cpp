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
///@brief solver: solve the mechanical thermomechanical equilibrium			//
//	for a homogeneous loading path, allowing repeatable steps
///@version 1.9

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Umat/umat_smart.hpp>

#include <smartplus/Libraries/Solver/read.hpp>
#include <smartplus/Libraries/Solver/block.hpp>
#include <smartplus/Libraries/Solver/step.hpp>
#include <smartplus/Libraries/Solver/step_meca.hpp>
#include <smartplus/Libraries/Solver/step_thermomeca.hpp>

using namespace std;
using namespace arma;
using namespace smart;

ofstream output("L.txt");

int main() {
    
	///Material properties reading, use "material.dat" to specify parameters values
	string umat_name;
	int nprops = 0;
	int nstatev = 0;
	vec props;
	vec statev;
    
    double rho = 0.;
    double c_p = 0.;
    
	bool start = true;
	double Time = 0.;
	double DTime = 0.;
	double T = 0.;
	double DT = 0.;
    double tnew_dt = 1.;
	vec sigma = zeros(6);
	vec Etot = zeros(6);
	vec DEtot = zeros(6);
	mat Lt = zeros(6,6);
	mat DR = eye(3,3);
    
    double sse = 0.;
    double spd = 0.;
    
    int ndi = 3;
    int nshr = 3;
    
    //read the material properties
    read_matprops(umat_name, nprops, props, nstatev, statev, rho, c_p);
        
    select_umat(umat_name, Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
    
    output << Lt << "\n";
    
	return 0;
}
