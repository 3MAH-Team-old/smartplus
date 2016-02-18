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

using namespace std;
using namespace arma;

namespace smart{

void solver(const string &umat_name, const vec &props, const double &nstatev, const double &rho, const double &c_p, const string &pathfile, const string &outputfile) {

    const char *car;
    car = outputfile.c_str();
    ofstream output(car);
    
	///Usefull UMAT variables
	int ndi = 3;
	int nshr = 3;
    
    mat dSdE = zeros(6,6);
    mat dSdT = zeros(1,6);
    mat dQdE = zeros(6,1);
    mat dQdT = zeros(1,1);
    
    vector<block> blocks;
    
	///Material properties reading, use "material.dat" to specify parameters values
	vec statev = zeros(nstatev);
    double nprops = props.n_elem;
    
	bool start = true;
	double Time = 0.;
	double DTime = 0.;
	double T = 0.;
	double DT = 0.;
    double tnew_dt = 1.;
	vec sigma = zeros(6);
    double Q = 0.;
    double rpl = 0.;            //The part of the heat linked with material's behavior (thermoelasticity, dissipation, phase change..)
	vec Etot = zeros(6);
	vec DEtot = zeros(6);
	mat Lt = zeros(6,6);
	mat DR = eye(3,3);
    
    double sse = 0.;
    double spd = 0.;
    
    //read the material properties
    //Read the loading path
    read_path(blocks, T, pathfile);
    
    //Output
    int o_ncount = 0;
    double o_tcount = 0.;
    solver_output so(blocks.size());
    read_output(so, blocks.size(), nstatev);
    
    //Chck output and step files
    check_path_output(blocks, so);
    
    vec statev_start = statev;
    
    double error = 0.;
    vec residual;
    vec Delta;
    int nK = 0; // The size of the problem to solve
    mat K;
    mat invK;
    int compteur = 0.;
    
    int inc = 0.;
    double tinc=0.;
    double Dtinc=0.;
    double Dtinc_cur=0.;
    double tau = 0.;        //Time constant for convexion thermal mechanical conditions
    double Tinit = T;       //Initial temperature
    
    /// Block loop
    for(unsigned int i = 0 ; i < blocks.size() ; i++){
        
//        cout << blocks[i];
        
        switch(blocks[i].type) {
            case 1: {
        
                /// resize the problem to solve
                residual = zeros(6);
                Delta = zeros(6);
                K = zeros(6, 6);
                invK = zeros(6, 6);
                
                vec sigma_start = sigma;
                Lt = zeros(6,6);
                
                DEtot = zeros(6);
                DR = eye(3,3);
                DTime = 0.;
                DT = 0.;
                
                //Run the umat for the first time in the block. So that we get the proper tangent properties
                run_umat(umat_name, Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
                
                statev_start = statev;
                start = false;
                mat Lt_start = Lt;
                
                /// Cycle loop
                for(int n = 0; n < blocks[i].ncycle; n++){
                    
                    /// Step loop
                    for(int j = 0; j < blocks[i].nstep; j++){
                        
                        shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[j]);
                        
                        sptr_meca->generate(Time, sigma, Etot, T);
                        
                        //Write the initial results
                        if ((n == 0)&&(j == 0)) {
                            //Write the results
                            sptr_meca->output(output, so, i, n, 0, statev);
                        }
                        
                        nK = sum(sptr_meca->cBC_meca);
                        
                        inc = 0;
                        while(inc < sptr_meca->ninc) {
                            
                            if(error > precision_solver) {
                                for(int k = 0 ; k < 6 ; k++)
                                {
                                    if(sptr_meca->cBC_meca(k)) {
                                        sptr_meca->mecas(inc,k) -= residual(k);
                                    }
                                }
                            }
                            
                            
                            while (tinc<1.) {
                                
                                sptr_meca->compute_inc(tnew_dt, inc, tinc, Dtinc, Dtinc_cur);
                                
                                if(nK == 0){
                                    
                                    DEtot = Dtinc*sptr_meca->mecas.row(inc).t();
                                    DT = Dtinc*sptr_meca->Ts(inc);
                                    DTime = Dtinc*sptr_meca->times(inc);
                                    
                                    run_umat(umat_name, Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
                                }
                                else{
                                    /// ********************** SOLVING THE MIXED PROBLEM NRSTRUCT ***********************************
                                    ///Saving stress and stress set point at the beginning of the loop
                                    
                                    error = 1.;
                                   
                                    DEtot = zeros(6);
                                    Lt = Lt_start;
                                    
                                    for(int k = 0 ; k < 6 ; k++)
                                    {
                                        if (sptr_meca->cBC_meca(k)) {
                                            residual(k) = sigma(k) - sigma_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                        }
                                        else {
                                            residual(k) = lambda_solver*(DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                        }
                                    }
                                    
                                    while((error > precision_solver)&&(compteur < maxiter_solver)) {
                                        
                                        ///Prediction of the strain increment using the tangent modulus given from the umat_ function
                                        //we use the ddsdde (Lt) from the previous increment
                                        Lt_2_K(Lt, K, sptr_meca->cBC_meca, lambda_solver);
                                        
                                        ///jacobian inversion
                                        invK = inv(K);
                                        
                                        /// Prediction of the component of the strain tensor
                                        Delta = -invK * residual;
                                        
                                        DEtot += Delta;
                                        DT = Dtinc*sptr_meca->Ts(inc);
                                        DTime = Dtinc*sptr_meca->times(inc);
                                        
                                        sigma = sigma_start;
                                        statev = statev_start;
                                        
                                        run_umat(umat_name, Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
                                        
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
                                                residual(k) = sigma(k) - sigma_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }

                                        compteur++;
                                        error = sqrt(norm(residual, 2.));

                                        if(tnew_dt < 1.) {
                                            if((fabs(Dtinc_cur - sptr_meca->Dn_mini) > iota)||(inforce_solver == 0)) {
                                                compteur = maxiter_solver;
                                            }
                                        }
                                        
                                    }
                                    
                                }
                                
                                if((fabs(Dtinc_cur - sptr_meca->Dn_mini) < iota)&&(inforce_solver == 1)&&(tnew_dt < 1.)) {
                                    cout << "The solver has been inforced to proceed (User subroutine issue) at step:" << sptr_meca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    tnew_dt = 1.;
                                }
                                
                                if(error > precision_solver) {
                                    if((Dtinc_cur == sptr_meca->Dn_mini)&&(inforce_solver == 1)) {
                                        cout << "The solver has been inforced to proceed (Solver issue) at step:" << sptr_meca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                        cout << "The next increment has integrated the error to avoid propagation\n";
                                        //The solver has been inforced!
                                        tnew_dt = 1.;

                                        if (inc+1<sptr_meca->ninc) {
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if(sptr_meca->cBC_meca(k)) {
                                                    sptr_meca->mecas(inc+1,k) -= residual(k);
                                                }
                                            }
                                        }
                                        
                                    }
                                    else {
                                        tnew_dt = div_tnew_dt_solver;
                                    }
                                }
                                
                                if((compteur < miniter_solver)&&(tnew_dt >= 1.)) {
                                    tnew_dt = mul_tnew_dt_solver;
                                }
                                compteur = 0;
                                
                                sptr_meca->assess_inc(tnew_dt, tinc, Dtinc, Etot, DEtot, T, DT, Time, DTime, sigma, sigma_start, statev, statev_start, Lt, Lt_start);
                            
                                
                            }
                            
                            //At the end of each increment, check if results should be written
                            if (so.o_type(i) == 1) {
                                o_ncount++;
                            }
                            if (so.o_type(i) == 2) {
                                o_tcount+=DTime;
                            }
                            
                            sptr_meca->Time = Time;
                            sptr_meca->Etot = Etot;
                            sptr_meca->sigma = sigma;
                            sptr_meca->T = T;
                            
                            //Write the results
                            if (((so.o_type(i) == 1)&&(o_ncount == so.o_nfreq(i)))||(((so.o_type(i) == 2)&&(fabs(o_tcount - so.o_tfreq(i)) < 1.E-12)))) {
                                
                                sptr_meca->output(output, so, i, n, inc, statev);
                                
                                if (so.o_type(i) == 1) {
                                    o_ncount = 0;
                                }
                                if (so.o_type(i) == 2) {
                                    o_tcount = 0.;
                                }
                            }
                            
                            tinc = 0.;
                            inc++;
                         }
                                                
                    }
                        
                }
                break;
            } //Mechanical
            case 2: {
                
                /// resize the problem to solve
                residual = zeros(7);
                Delta = zeros(7);
                K = zeros(7, 7);
                invK = zeros(7, 7);
                
                vec sigma_start = sigma;
                // double Q_start = Q;
                dSdE = zeros(6,6);
                dSdT = zeros(1,6);
                dQdE = zeros(6,1);
                dQdT = zeros(1,1);
                
                mat drpldE = dQdE.t();
                mat drpldT = dQdT;
                
                DEtot = zeros(6);
                DR = eye(3,3);
                DTime = 0.;
                DT = 0.;
                
                //Run the umat for the first time in the block. So that we get the proper tangent properties
                run_umat_T(umat_name, Etot, DEtot, sigma, rpl, dSdE, dSdT, drpldE, drpldT, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
                
                Q = -1.*rpl;    //Since DTime=0;
                
                dQdE = -drpldE.t();
                dQdT = lambda_solver;

                statev_start = statev;                
                start = false;
                mat dSdE_start = dSdE;
                mat dSdT_start = dSdT;
                mat dQdE_start = dQdE;
                mat dQdT_start = dQdT;
                
                
                /// Cycle loop
                for(int n = 0; n < blocks[i].ncycle; n++){
                    
                    /// Step loop
                    for(int j = 0; j < blocks[i].nstep; j++){
                        
                        shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[j]);
                        sptr_thermomeca->generate(Time, sigma, Etot, T);
                        
                        //Write the initial results
                        if ((n == 0)&&(j == 0)) {
                            //Write the results
                            sptr_thermomeca->output(output, so, i, n, 0, statev);
                        }
                        
                        nK = sum(sptr_thermomeca->cBC_meca);
                        
                        inc = 0;
                        if(sptr_thermomeca->cBC_T == 3)
                            tau = sptr_thermomeca->BC_T;
                            
                        
                        while(inc < sptr_thermomeca->ninc) {
                            
                            
                            if(error > precision_solver) {
                                for(int k = 0 ; k < 6 ; k++)
                                {
                                    if (sptr_thermomeca->cBC_meca(k)) {
                                        sptr_thermomeca->mecas(inc,k) -= residual(k);
                                    }
                                }
                                if (sptr_thermomeca->cBC_T) {
                                    sptr_thermomeca->Ts(inc) -= residual(6);
                                }
                            }
                            
                            while (tinc<1.) {

                                sptr_thermomeca->compute_inc(tnew_dt, inc, tinc, Dtinc, Dtinc_cur);
                                
                                if(nK + sptr_thermomeca->cBC_T == 0){
                                    
                                    DEtot = Dtinc*sptr_thermomeca->mecas.row(inc).t();
                                    DT = Dtinc*sptr_thermomeca->Ts(inc);
                                    DTime = Dtinc*sptr_thermomeca->times(inc);

                                    run_umat_T(umat_name, Etot, DEtot, sigma, rpl, dSdE, dSdT, drpldE, drpldT, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
                                    
                                    if (DTime < 1.E-12) {
                                        Q = -1.*rpl;    //Since DTime=0;
                                    }
                                    else{
                                        if(sptr_thermomeca->cBC_T == 3)
                                            Q = rho*c_p*DT*(1./DTime) - rpl;
                                        else
                                            Q = rho*c_p*DT/DTime - rpl;
                                    }
                                    
                                }
                                else{
                                    /// ********************** SOLVING THE MIXED PROBLEM NRSTRUCT ***********************************
                                    ///Saving stress and stress set point at the beginning of the loop
                                    //                                    Etot_start ?
                                    
                                    error = 1.;
                                    
                                    DEtot = zeros(6);
                                    dSdE = dSdE_start;
                                    dSdT = dSdT_start;
                                    dQdE = dQdE_start;
                                    dQdT = dQdT_start;
                                    
                                    DEtot = zeros(6);
                                    DT = 0.;
                                                                    
                                    //Construction of the initial residual
                                    for(int k = 0 ; k < 6 ; k++)
                                    {
                                        if (sptr_thermomeca->cBC_meca(k)) {
                                            residual(k) = sigma(k) - sigma_start(k) - Dtinc*sptr_thermomeca->mecas(inc,k);
                                        }
                                        else {
                                            residual(k) = lambda_solver*(DEtot(k) - Dtinc*sptr_thermomeca->mecas(inc,k));
                                        }
                                    }
                                    if (sptr_thermomeca->cBC_T == 1) {
                                        residual(6) = Q - sptr_thermomeca->Ts(inc);
                                    }
                                    else if(sptr_thermomeca->cBC_T == 0) {
                                        residual(6) = lambda_solver*(DT - Dtinc*sptr_thermomeca->Ts(inc));
                                    }
                                    else if(sptr_thermomeca->cBC_T == 3) { //Special case of 0D convexion that depends on temperature assumption
                                        residual(6) = Q + (rho*c_p/tau)*(T-Tinit);
                                    }
                                    else {
                                        cout << "error : The Thermal BC is not recognized\n";
                                        exit(0);
                                    }
                                    
                                    while((error > precision_solver)&&(compteur < maxiter_solver)) {
                                        
                                        ///Prediction of the strain increment using the tangent modulus given from the umat_ function
                                        //we use the ddsdde (Lt) from the previous increment
                                                                        
                                        Lth_2_K(dSdE, dSdT, dQdE, dQdT, K, sptr_thermomeca->cBC_meca, sptr_thermomeca->cBC_T, lambda_solver);
                                        
                                        ///jacobian inversion
                                        invK = inv(K);
                                        
                                        /// Prediction of the component of the strain tensor
                                        Delta = -invK * residual;
                                        
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            DEtot(k) += Delta(k);
                                        }
                                        DT += Delta(6);
                                        DTime = Dtinc*sptr_thermomeca->times(inc);
                                        
                                        sigma = sigma_start;
                                        statev = statev_start;
                                        
                                        run_umat_T(umat_name, Etot, DEtot, sigma, rpl, dSdE, dSdT, drpldE, drpldT, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);

                                        if (DTime < 1.E-12) {
                                            Q = -1.*rpl;    //Since DTime=0;
                                            
                                            dQdE = -drpldE.t();
                                            dQdT = lambda_solver;  //To avoid any singularity in the system
                                            
                                        }
                                        else{
                                            if(sptr_thermomeca->cBC_T == 3)
                                                Q = rho*c_p*DT*(1./DTime) + rho*c_p/tau  - rpl;
                                            else
                                                Q = rho*c_p*DT/DTime - rpl;
                        
                                            dQdE = -drpldE.t();
                                            
                                            if(sptr_thermomeca->cBC_T == 3)
                                                dQdT = rho*c_p*(1./DTime) - drpldT;
                                            else
                                                dQdT = rho*c_p/DTime - drpldT;
                                            
                                        }
                                        
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_thermomeca->cBC_meca(k)) {
                                                residual(k) = sigma(k) - sigma_start(k) - Dtinc*sptr_thermomeca->mecas(inc,k);
                                            }
                                            else{
                                                residual(k) = lambda_solver*(DEtot(k) - Dtinc*sptr_thermomeca->mecas(inc,k));
                                            }
                                        }
                                        if (sptr_thermomeca->cBC_T == 1) {
                                            residual(6) = Q - sptr_thermomeca->Ts(inc);
                                        }
                                        else if(sptr_thermomeca->cBC_T == 0) {
                                            residual(6) = lambda_solver*(DT - Dtinc*sptr_thermomeca->Ts(inc));
                                        }
                                        else if(sptr_thermomeca->cBC_T == 3) { //Special case of 0D convexion that depends on temperature assumption
                                            residual(6) = Q + (rho*c_p/tau)*(T-Tinit);
                                        }
                                        else {
                                            cout << "error : The Thermal BC is not recognized\n";
                                            exit(0);
                                        }

                                        
                                        error = sqrt(norm(residual, 2.));
                                        compteur++;
                                                                            
                                    }
                                    
                                }
                                
                                if((fabs(Dtinc_cur - sptr_thermomeca->Dn_mini) < iota)&&(inforce_solver == 1)&&(tnew_dt < 1.)) {
                                    cout << "The solver has been inforced to proceed (User subroutine issue) at step:" << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    tnew_dt = 1.;
                                }
                                
                                if(error > precision_solver) {
                                    if((Dtinc_cur == sptr_thermomeca->Dn_mini)&&(inforce_solver == 1)) {
                                        cout << "The solver has been inforced to proceed (Solver issue) at step:" << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                        cout << "The next increment has integrated the error to avoid propagation\n";
                                        //The solver has been inforced!
                                        tnew_dt = 1.;
                                        
                                    }
                                    else {
                                        tnew_dt = div_tnew_dt_solver;
                                    }
                                }
                                
                                if((compteur < miniter_solver)&&(tnew_dt >= 1.)) {
                                    tnew_dt = mul_tnew_dt_solver;
                                }
                                compteur = 0;
                                
                                sptr_thermomeca->assess_inc(tnew_dt, tinc, Dtinc, Etot, DEtot, T, DT, Time, DTime, sigma, sigma_start, statev, statev_start, dSdE, dSdE_start, dSdT, dSdT_start, dQdE, dQdE_start, dQdT, dQdT_start);
                                
                            }

                            //At the end of each increment, check if results should be written
                            if (so.o_type(i) == 1) {
                                o_ncount++;
                            }
                            if (so.o_type(i) == 2) {
                                o_tcount+=DTime;
                            }
                            
                            sptr_thermomeca->Time = Time;
                            sptr_thermomeca->Etot = Etot;
                            sptr_thermomeca->sigma = sigma;
                            sptr_thermomeca->T = T;
                            sptr_thermomeca->Q = Q;
                            
                            //Write the results
                            if (((so.o_type(i) == 1)&&(o_ncount == so.o_nfreq(i)))||(((so.o_type(i) == 2)&&(fabs(o_tcount - so.o_tfreq(i)) < 1.E-12)))) {
                                
                                sptr_thermomeca->output(output, so, i, n, inc, statev);
                                
                                if (so.o_type(i) == 1) {
                                    o_ncount = 0;
                                }
                                if (so.o_type(i) == 2) {
                                    o_tcount = 0.;
                                }
                            }
                            
                            tinc = 0.;
                            inc++;
                            
                        }
                        
                    }
                    
                }
                break;
            } //Thermomechanical
            default: {
                cout << "the block type is not defined!\n";
                break;
            }
        }
    //end of blocks loops
    }


        
}

} //namespace smart
