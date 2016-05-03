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
#include <memory>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Phase/material_characteristics.hpp>
#include <smartplus/Libraries/Phase/state_variables.hpp>
#include <smartplus/Libraries/Phase/state_variables_M.hpp>
#include <smartplus/Libraries/Phase/state_variables_T.hpp>
#include <smartplus/Umat/umat_smart.hpp>
#include <smartplus/Libraries/Solver/read.hpp>
#include <smartplus/Libraries/Solver/block.hpp>
#include <smartplus/Libraries/Solver/step.hpp>
#include <smartplus/Libraries/Solver/step_meca.hpp>
#include <smartplus/Libraries/Solver/step_thermomeca.hpp>

using namespace std;
using namespace arma;

namespace smart{

void solver(const string &umat_name, const vec &props, const double &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve,const double &rho, const double &c_p, const string &pathfile, const string &outputfile) {
    
//    ofstream output(outputfile);
    
	///Usefull UMAT variables
	int ndi = 3;
	int nshr = 3;    
    std::vector<block> blocks;  //loading blocks
    phase_characteristics rve;  // Representative volume element
    
	bool start = true;
	double Time = 0.;
	double DTime = 0.;
    double T_init = 0.;
//	double T = 0.;
//	double DT = 0.;
    double tnew_dt = 1.;
//    double Q = 0.;
//    double rpl = 0.;            //The part of the heat linked with material's behavior (thermoelasticity, dissipation, phase change..)
//    double sse = 0.;
//    double spd = 0.;
    
    //What follows constitutes the base variables to solve the mechanical/thermomechanical problem
//    vec sigma = zeros(6);
//    vec sigma_start = zeros(6);
//    vec Etot = zeros(6);
//    vec DEtot = zeros(6);
    mat L = zeros(6,6);
    mat Lt = zeros(6,6);
    mat DR = eye(3,3);
    
//    mat dSdE = zeros(6,6);
//    mat dSdT = zeros(1,6);
    mat dQdE = zeros(6,1);
    mat dQdT = zeros(1,1);
    
    //read the material properties
    //Read the loading path
    read_path(blocks, T_init, pathfile);
    
    for(auto b : blocks) {
        cout << "blocks = " << b << "\n";
    }
    
//    std::vector<shared_ptr<material_characteristics> > rve_save(1); //There is at least one material, the effective one
//    shared_ptr<material_characteristics> eff_material = rve_save[0]; //There is at least one material, the effective one
//    shared_ptr<material_characteristics> eff_material
    
    ///Material properties reading, use "material.dat" to specify parameters values
    rve.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props, rho, c_p);
    
    //Output
    int o_ncount = 0;
    double o_tcount = 0.;
    solver_output so(blocks.size());
//    std::vector<solver_output> rve_so;   //This is to output the response of the phases
    read_output(so, blocks.size(), nstatev);
    
    //Check output and step files
    check_path_output(blocks, so);
    
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
    
    /// Block loop
    for(unsigned int i = 0 ; i < blocks.size() ; i++){

        switch(blocks[i].type) {
            case 1: { //Mechanical
                
                /// resize the problem to solve
                residual = zeros(6);
                Delta = zeros(6);
                K = zeros(6,6);
                invK = zeros(6,6);
                
                shared_ptr<state_variables_M> sv_M;
                
                if(start) {
                    rve.construct(0,blocks[i].type);
                    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), T_init, 0., 0., 0., nstatev, zeros(nstatev), zeros(nstatev));
                    sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                }
                else {
                    //Dynamic cast from some other (possible state_variable_M/T)
                    /*sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                    rve.construct(0,blocks[i].type);
                    rve.sptr_sv_global->update(sv_M->Etot, sv_M->DEtot, sv_M->sigma, sv_M->sigma_start, sv_M->T, sv_M->DT, sv_M->sse, sv_M->spd, nstatev, sv_M->statev, sv_M->statev_start);*/
                    //sv_M is reassigned properly
                    sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                }
                
                sv_M->L = zeros(6,6);
                sv_M->Lt = zeros(6,6);
                
                DR = eye(3,3);
                DTime = 0.;
                sv_M->DEtot = zeros(6);
                sv_M->DT = 0.;
                
                //Run the umat for the first time in the block. So that we get the proper tangent properties
                run_umat_M(rve, DR, Time, DTime, ndi, nshr, start, tnew_dt);
                
                if(start) {
                    //Use the number of phases saved to define the files
                    rve.define_output(outputfile);
                    //Write the initial results
                    rve.output(so, -1, -1, -1, -1, Time);
                }
                //Set the start values of sigma_start=sigma and statev_start=statev for all phases
                rve.set_start(); //DEtot = 0 and DT = 0 so we can use it safely here
                start = false;
                
                /// Cycle loop
                for(int n = 0; n < blocks[i].ncycle; n++){
                    
                    /// Step loop
                    for(int j = 0; j < blocks[i].nstep; j++){
                        
                        shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[j]);
                        sptr_meca->generate(Time, sv_M->Etot, sv_M->sigma, sv_M->T);
                        
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
                                    
                                    sv_M->DEtot = Dtinc*sptr_meca->mecas.row(inc).t();
                                    sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                    DTime = Dtinc*sptr_meca->times(inc);
                                    
                                    run_umat_M(rve, DR, Time, DTime, ndi, nshr, start, tnew_dt);
                                }
                                else{
                                    /// ********************** SOLVING THE MIXED PROBLEM NRSTRUCT ***********************************
                                    ///Saving stress and stress set point at the beginning of the loop
                                    
                                    error = 1.;
                                   
                                    sv_M->DEtot = zeros(6);
                                    
                                    for(int k = 0 ; k < 6 ; k++)
                                    {
                                        if (sptr_meca->cBC_meca(k)) {
                                            residual(k) = sv_M->sigma(k) - sv_M->sigma_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                        }
                                        else {
                                            residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                        }
                                    }
                                    
                                    while((error > precision_solver)&&(compteur < maxiter_solver)) {
                                        
                                        ///Prediction of the strain increment using the tangent modulus given from the umat_ function
                                        //we use the ddsdde (Lt) from the previous increment
                                        Lt_2_K(sv_M->Lt, K, sptr_meca->cBC_meca, lambda_solver);
                                        
                                        ///jacobian inversion
                                        invK = inv(K);
                                        
                                        /// Prediction of the component of the strain tensor
                                        Delta = -invK * residual;
                                        
                                        sv_M->DEtot += Delta;
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        DTime = Dtinc*sptr_meca->times(inc);
                                        
                                        rve.to_start();
                                        run_umat_M(rve, DR, Time, DTime, ndi, nshr, start, tnew_dt);
                                        
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
                                                residual(k) = sv_M->sigma(k) - sv_M->sigma_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }
                                        
                                        compteur++;
                                        error = norm(residual, 2.);
                                        
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
                                        cout << "The solver has been inforced to proceed (Solver issue) at step:" << sptr_meca->number << " inc: " << inc << " and fraction:" << tinc << ", with the error: " << error << "\n";
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
                                
                                sptr_meca->assess_inc(tnew_dt, tinc, Dtinc, rve ,Time, DTime);
                                //start variables ready for the next increment
                                
                            }
                            
                            //At the end of each increment, check if results should be written
                            if (so.o_type(i) == 1) {
                                o_ncount++;
                            }
                            if (so.o_type(i) == 2) {
                                o_tcount+=DTime;
                            }
                            
                            //Write the results
                            if (((so.o_type(i) == 1)&&(o_ncount == so.o_nfreq(i)))||(((so.o_type(i) == 2)&&(fabs(o_tcount - so.o_tfreq(i)) < 1.E-12)))) {
                                
                                rve.output(so, i, n, j, inc, Time);

                                
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
            }
            case 2: { //Thermomechanical
                
                /// resize the problem to solve
                residual = zeros(7);
                Delta = zeros(7);
                K = zeros(7,7);
                invK = zeros(7,7);
                
                shared_ptr<state_variables_T> sv_T;
                
                if(start) {
                    rve.construct(0,blocks[i].type);
                    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), T_init, 0., 0., 0., nstatev, zeros(nstatev), zeros(nstatev));
                    sv_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
                }
                else {
                    //Dynamic cast from some other (possible state_variable_M/T)
                    /*sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                     rve.construct(0,blocks[i].type);
                     rve.sptr_sv_global->update(sv_M->Etot, sv_M->DEtot, sv_M->sigma, sv_M->sigma_start, sv_M->T, sv_M->DT, sv_M->sse, sv_M->spd, nstatev, sv_M->statev, sv_M->statev_start);*/
                    //sv_M is reassigned properly
                    sv_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
                }
                
                sv_T->dSdE = zeros(6,6);
                sv_T->dSdT = zeros(6,1);
                dQdE = zeros(1,6);
                dQdT = zeros(1,1);
                
                DR = eye(3,3);
                DTime = 0.;
                sv_T->DEtot = zeros(6);
                sv_T->DT = 0.;
                
                //Run the umat for the first time in the block. So that we get the proper tangent properties
                run_umat_T(rve, DR, Time, DTime, ndi, nshr, start, tnew_dt);
                
                sv_T->Q = -1.*sv_T->rpl;    //Since DTime=0;
                dQdT = lambda_solver;  //To avoid any singularity in the system                
                
                if(start) {
                    //Use the number of phases saved to define the files
                    rve.define_output(outputfile);
                    //Write the initial results
                    rve.output(so, -1, -1, -1, -1, Time);
                }
                //Set the start values of sigma_start=sigma and statev_start=statev for all phases
                rve.set_start(); //DEtot = 0 and DT = 0 so we can use it safely here
                start = false;
                
                /// Cycle loop
                for(int n = 0; n < blocks[i].ncycle; n++){
                    
                    /// Step loop
                    for(int j = 0; j < blocks[i].nstep; j++){
                        
                        
                        shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[j]);
                        sptr_thermomeca->generate(Time, sv_T->Etot, sv_T->sigma, sv_T->T);
                        
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
                                    
                                    sv_T->DEtot = Dtinc*sptr_thermomeca->mecas.row(inc).t();
                                    sv_T->DT = Dtinc*sptr_thermomeca->Ts(inc);
                                    DTime = Dtinc*sptr_thermomeca->times(inc);
                                    
                                    run_umat_T(rve, DR, Time, DTime, ndi, nshr, start, tnew_dt);
                                    
                                    if (DTime < 1.E-12) {
                                        sv_T->Q = -1.*sv_T->rpl;    //Since DTime=0;
                                    }
                                    else{
                                        sv_T->Q = rve.sptr_matprops->rho*rve.sptr_matprops->c_p*sv_T->DT*(1./DTime) - sv_T->rpl;
                                    }
                                    
                                }
                                else{
                                    /// ********************** SOLVING THE MIXED PROBLEM NRSTRUCT ***********************************
                                    ///Saving stress and stress set point at the beginning of the loop
                                    
                                    error = 1.;
                                    
                                    sv_T->DEtot = zeros(6);
                                    sv_T->DT = 0.;
                                    
                                    //Construction of the initial residual
                                    for(int k = 0 ; k < 6 ; k++)
                                    {
                                        if (sptr_thermomeca->cBC_meca(k)) {
                                            residual(k) = sv_T->sigma(k) - sv_T->sigma_start(k) - Dtinc*sptr_thermomeca->mecas(inc,k);
                                        }
                                        else {
                                            residual(k) = lambda_solver*(sv_T->DEtot(k) - Dtinc*sptr_thermomeca->mecas(inc,k));
                                        }
                                    }
                                    if (sptr_thermomeca->cBC_T == 1) {
                                        residual(6) = sv_T->Q - sptr_thermomeca->Ts(inc);
                                    }
                                    else if(sptr_thermomeca->cBC_T == 0) {
                                        residual(6) = lambda_solver*(sv_T->DT - Dtinc*sptr_thermomeca->Ts(inc));
                                    }
                                    else if(sptr_thermomeca->cBC_T == 3) { //Special case of 0D convexion that depends on temperature assumption
                                        residual(6) = sv_T->Q + (rve.sptr_matprops->rho*rve.sptr_matprops->c_p/tau)*(sv_T->T-T_init);
                                    }
                                    else {
                                        cout << "error : The Thermal BC is not recognized\n";
                                        exit(0);
                                    }
                                    
                                    while((error > precision_solver)&&(compteur < maxiter_solver)) {
                                        
                                        ///Prediction of the strain increment using the tangent modulus given from the umat_ function
                                        //we use the ddsdde (Lt) from the previous increment
                                        Lth_2_K(sv_T->dSdE, sv_T->dSdT, dQdE, dQdT, K, sptr_thermomeca->cBC_meca, sptr_thermomeca->cBC_T, lambda_solver);
                                        
                                        ///jacobian inversion
                                        invK = inv(K);
                                        
                                        /// Prediction of the component of the strain tensor
                                        Delta = -invK * residual;
                                        
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            sv_T->DEtot(k) += Delta(k);
                                        }
                                        sv_T->DT += Delta(6);
                                        DTime = Dtinc*sptr_thermomeca->times(inc);
                                        
                                        rve.to_start();
                                        run_umat_T(rve, DR, Time, DTime, ndi, nshr, start, tnew_dt);
                                        
                                        if (DTime < 1.E-12) {
                                            sv_T->Q = -1.*sv_T->rpl;    //Since DTime=0;
                                            
                                            dQdE = -sv_T->drpldE.t();
                                            dQdT = lambda_solver;  //To avoid any singularity in the system
                     
                                        }
                                        else{
                                            //Attention here, we solve Phi = Q - rho*c_p*(1./tau)*(sv_T->T-T_init) = rho*c_p*DT*(1./DTime) - rpl -(-rho*c_p*(1./tau)*(sv_T->T-T_init)) = 0
                                            //So the derivative / T has the extra rho*c_p*(1./tau)
                                            //It is actually icluded here in dQdT
                                            sv_T->Q = rve.sptr_matprops->rho*rve.sptr_matprops->c_p*sv_T->DT*(1./DTime) - sv_T->rpl;
                                            
                                            dQdE = -sv_T->drpldE.t();
                                            
                                            if (sptr_thermomeca->cBC_T < 3) {
                                                dQdT = rve.sptr_matprops->rho*rve.sptr_matprops->c_p*(1./DTime) - sv_T->drpldT;
                                            }
                                            else if(sptr_thermomeca->cBC_T == 3) {
                                                dQdT = rve.sptr_matprops->rho*rve.sptr_matprops->c_p*(1./DTime) - sv_T->drpldT + rve.sptr_matprops->rho*rve.sptr_matprops->c_p*(1./tau);

                                            }
                                            
                                        }
                                        
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_thermomeca->cBC_meca(k)) {
                                                residual(k) = sv_T->sigma(k) - sv_T->sigma_start(k) - Dtinc*sptr_thermomeca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_T->DEtot(k) - Dtinc*sptr_thermomeca->mecas(inc,k));
                                            }
                                        }
                                        if (sptr_thermomeca->cBC_T == 1) {
                                            residual(6) = sv_T->Q - sptr_thermomeca->Ts(inc);
                                        }
                                        else if(sptr_thermomeca->cBC_T == 0) {
                                            residual(6) = lambda_solver*(sv_T->DT - Dtinc*sptr_thermomeca->Ts(inc));
                                        }
                                        else if(sptr_thermomeca->cBC_T == 3) { //Special case of 0D convexion that depends on temperature assumption
                                            residual(6) = sv_T->Q + (rve.sptr_matprops->rho*rve.sptr_matprops->c_p/tau)*(sv_T->T-T_init);
                                        }
                                        else {
                                            cout << "error : The Thermal BC is not recognized\n";
                                            exit(0);
                                        }
                                        
                                        compteur++;
                                        error = norm(residual, 2.);
                                        
                                        if(tnew_dt < 1.) {
                                            if((fabs(Dtinc_cur - sptr_thermomeca->Dn_mini) > iota)||(inforce_solver == 0)) {
                                                compteur = maxiter_solver;
                                            }
                                        }
                                        
                                    }
                                    
                                }
                                
                                if((fabs(Dtinc_cur - sptr_thermomeca->Dn_mini) < iota)&&(inforce_solver == 1)&&(tnew_dt < 1.)) {
                                    cout << "The solver has been inforced to proceed (User subroutine issue) at step:" << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    tnew_dt = 1.;
                                }
                                
                                if(error > precision_solver) {
                                    if((Dtinc_cur == sptr_thermomeca->Dn_mini)&&(inforce_solver == 1)) {
                                        cout << "The solver has been inforced to proceed (Solver issue) at step:" << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << ", with the error: " << error << "\n";
                                        cout << "The next increment has integrated the error to avoid propagation\n";
                                        //The solver has been inforced!
                                        tnew_dt = 1.;
                                        
                                        if (inc+1<sptr_thermomeca->ninc) {
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if(sptr_thermomeca->cBC_meca(k)) {
                                                    sptr_thermomeca->mecas(inc+1,k) -= residual(k);
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
                                
                                sptr_thermomeca->assess_inc(tnew_dt, tinc, Dtinc, rve ,Time, DTime);
                                //start variables ready for the next increment
                                
                            }
                            
                            //At the end of each increment, check if results should be written
                            if (so.o_type(i) == 1) {
                                o_ncount++;
                            }
                            if (so.o_type(i) == 2) {
                                o_tcount+=DTime;
                            }
                            
                            //Write the results
                            if (((so.o_type(i) == 1)&&(o_ncount == so.o_nfreq(i)))||(((so.o_type(i) == 2)&&(fabs(o_tcount - so.o_tfreq(i)) < 1.E-12)))) {
                                
                                rve.output(so, i, n, j, inc, Time);
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
            }
            default: {
                cout << "the block type is not defined!\n";
                break;
            }
        }
    //end of blocks loops
    }


        
}

} //namespace smart
