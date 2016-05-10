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

///@file solver.cpp
///@brief solver: solve the mechanical thermomechanical equilibrium			//
//	for a homogeneous loading path, allowing repeatable steps
///@version 1.9

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Umat/umat_smart.hpp>
#include <smartplus/Libraries/Solver/read.hpp>
#include <smartplus/Libraries/Solver/block.hpp>
#include <smartplus/Libraries/Solver/step.hpp>
#include <smartplus/Libraries/Solver/step_meca.hpp>
#include <smartplus/Libraries/Solver/step_thermomeca.hpp>
#include <smartplus/Libraries/Solver/solver.hpp>

using namespace std;
using namespace arma;
using namespace smart;

int main() {

    string outputfile = "results_job.txt";
    string pathfile = "path.txt";
	string umat_name;
	int nprops = 0;
	int nstatev = 0;
	vec props;
    
    double rho = 0.;
    double c_p = 0.;
        
	double psi_rve = 0.;
	double theta_rve = 0.;
	double phi_rve = 0.;
    
    read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, rho, c_p);
    solver(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve, rho, c_p, pathfile, outputfile);

	return 0;
}
