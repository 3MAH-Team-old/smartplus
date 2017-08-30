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

    string path_data = "data";
    string path_results = "results";
    string outputfile = "results_job.txt";
    string pathfile = "path.txt";
	string materialfile = "material.dat";
    string sol_essentials = "solver_essentials.inp";
    string sol_control = "solver_control.inp";


    string umat_name;
	int nprops = 0;
	int nstatev = 0;
	vec props;
    
	double psi_rve = 0.;
	double theta_rve = 0.;
	double phi_rve = 0.;

    int solver_type = 0;
    double div_tnew_dt_solver = 0.;
    double mul_tnew_dt_solver = 0.;
    int miniter_solver = 0;
    int maxiter_solver = 0;
    int inforce_solver = 0;
    double precision_solver = 0.;
    double lambda_solver = 0.;
    
    solver_essentials(solver_type, path_data, sol_essentials);
    solver_control(div_tnew_dt_solver, mul_tnew_dt_solver, miniter_solver, maxiter_solver, inforce_solver, precision_solver, lambda_solver, path_data, sol_control);
    
    read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    solver(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve, solver_type, div_tnew_dt_solver, mul_tnew_dt_solver, miniter_solver, maxiter_solver, inforce_solver, precision_solver, lambda_solver, path_data, path_results, pathfile, outputfile);
    
	return 0;
}
