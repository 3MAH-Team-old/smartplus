/* This file is part of SMART+ private.
 
 Only part of SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file damage_LLD_m.cpp
///@brief Modified Ladeveze-Le Dantec model that accounts for uniaxial damage and anisotopic plasticity
///@brief Damage evolution is considered as a function of time
///@author Y. Chemisky

#include <iostream>
#include <fstream>
#include <math.h>
#include <armadillo>

#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Maths/lagrange.hpp>
#include <smartplus/Libraries/Maths/rotation.hpp>
#include <smartplus/Libraries/Continuum_Mechanics/contimech.hpp>
#include <smartplus/Libraries/Continuum_Mechanics/constitutive.hpp>
#include <smartplus/Libraries/Maths/num_solve.hpp>

#include <smartplus/Umat/Mechanical/Damage/damage_LLD_0.hpp>

using namespace std;
using namespace arma;

namespace smart {
    
///@brief Material properties
///@param props(0) : E Young's modulus
///@param props(1) : nu Poisson's ration
///@param props(2) : alpha CTE
///@param props(3) : alphaD Damage evolution parameter alpha
///@param props(4) : betaD Damage evolution parameter beta
///@param props(5) : lambdaD Damage evolution parameter lambda
///@param props(6) : deltaD Damage evolution parameter delta
    
void umat_damage_LLD_0(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &sse, double &spd, const int &ndi, const int &nshr, const bool &start, double &tnew_dt) {
    
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    double Tinit = statev(0);
    double d_12 = statev(1);
    double d_22 = statev(2);
    double p_ts = statev(3);    //Accumulated plasticity
        
    vec EP(6);
    EP(0) = statev(4);
    EP(1) = statev(5);
    EP(2) = statev(6);
    EP(3) = statev(7);
    EP(4) = statev(8);
    EP(5) = statev(9);

    //Rotation of internal variables (tensors)
    EP = rotate_strain(EP, DR);
    
    ///@brief Initialization
    if (start) {
        Tinit = T;
        d_12 = 0.;
        d_22 = 0.;
        
        p_ts = 10.*iota;
        EP = zeros(6);
        sigma = zeros(6);
    }
    
    vec sigma_start = sigma;
    vec EP_start = EP;
    
    double EL = props(0);
    double ET = props(1);
    double nuTL = props(2);
    double nuTT = props(3);
    double GLT = props(4);
    double alphaL = props(5);
    double alphaT = props(6);
    
    //Shear damage
    double Y_12_0 = props(7);
    double Y_12_c = props(8);
    
    //Transverse damage
    double Y_22_0 = props(9);
    double Y_22_c = props(10);
    double Y_22_u = props(11);
    double b = props(12);
    
    //Coupled transverse-shear yield & plasticity
    double A_ts = props(13);
    double sigma_ts_0 = props(14);
    double alpha_ts = props(15);
    double beta_ts = props(16);
    
    //Set the Lagrange multipliers coefficients
    double c_lambda = 1.E-2;
    double p0_lambda = 1.E-1;
    double n_lambda = 1.0;
    double alpha_lambda = 0.;
    
    double E1 = EL;
    double E2_0 = ET;
    double E2 = ET*(1-d_22);
    double E3 = ET*(1-d_22);
    double nu12 = nuTL*(E1/E2_0);
    double nu13 = nuTL*(E1/E3);
    double nu32 = nuTT;
    double nu23 = nu32;
    double G12_0 = GLT;
    double G12 = G12_0*(1-d_12);
    double G13 = G12_0*(1-d_12);
    double G23 = E3/(2.*(1.+nu32));
    
    // ######################  Elastic stiffness #################################
    //defines L
    mat L = L_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
    
    //definition of the CTE tensor
    vec alpha = zeros(6);
    alpha = alphaT*Ith();
    alpha(0) += alphaL-alphaT;                              //WARNING CTE tensor is defined according to L is the direction 1
    
    //Compute the elastic strain and the related stress
    vec Eel_start = Etot - alpha*(T-Tinit) - EP;
    vec DEel_start = DEtot - alpha*DT;
    vec Eel = Eel_start + DEel_start;
    
    vec sigma_eff = zeros(6);
    vec sigma_eff_ts = zeros(6);
    
    ///Compute the Yd_12 and Yd_22 which are the criterium associated with damage
    double Yd_12 = 0.;
    double Yd_13 = 0.;
    double Yd_22 = 0.;
    double Yd_33 = 0.;
    
    //Compute Y_l, Y_t and Y_ts, which are used in the evolution equation of damage
    double Y_ts = 0.;              //transvers/shear coupling
    double Y_t = 0.;               //transverse

    double lambda_12 = lagrange_pow_1(d_12, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    double lambda_22 = lagrange_pow_1(d_22, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    double dlambda_12 = 0.;
    double dlambda_22 = 0.;
    
    //Compute the Plasticity functions
    //Compute the explicit flow direction
    vec Lambdap_ts = zeros(6);
    
    //Define the plastic function and the stress
    double Hp_ts = 0.;
    
    double Phi_p_ts = 0.;

    double dPhi_p_tsd_p = 0.;
    vec dPhi_p_tsd_sigma = zeros(6);
    
    //Note : The sup function are not required here, since we utilize Kuhn-Tucker conditions instead (dD >= 0)
    vec Phi_d = zeros(2);
    vec Phi_p = zeros(1);
    vec Dd = zeros(2);
    vec dd = zeros(2);
    vec Dp = zeros(1);
    vec dp = zeros(1);
    
    vec Y_dcrit = ones(2);
    vec Y_pcrit = ones(1);
    
    mat denom_d = zeros(2,2);
    mat denom_p = zeros(1,1);
    
    double dYd_12dd = 0.;
    double dYd_13dd = 0.;
    double dYd_22dd = 0.;
    double dYd_33dd = 0.;
    
    double dY_tsdd_12 = 0.;
    double dY_tsdd_22 = 0.;

    int compteur = 0;
    double error = 1.;
    
    mat Theta_ts = zeros(6,6);
    Theta_ts(1,1) = A_ts;
    Theta_ts(2,2) = A_ts;
    Theta_ts(3,3) = 1.;
    Theta_ts(4,4) = 1.;
    
    error = 1.;
    if (Y_t > Y_22_u) {
        d_12 = 0.99;
        d_22 = 0.99;
        error = 0.;
    }
    
    //First we find the plasticity
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        //Plasticity computations
                //Compute the hardening
        if (p_ts > iota)
            Hp_ts = beta_ts*pow(p_ts, alpha_ts);
        else
            Hp_ts = 0.;
        
        //effective stress
        sigma_eff = L*Eel;
        sigma_eff_ts = Theta_ts*sigma_eff;
        
        //Determine the Phi functions and their derivatives
        Phi_p_ts = Mises_stress(sigma_eff_ts) - Hp_ts - sigma_ts_0;

        //Compute the explicit flow direction
        Lambdap_ts = Theta_ts*eta_stress(sigma_eff_ts);
        
        if (p_ts > iota)
            dPhi_p_tsd_p = -1.*alpha_ts*beta_ts*pow(p_ts, alpha_ts-1.);
        else
            dPhi_p_tsd_p = 0.;
        
        dPhi_p_tsd_sigma = Theta_ts*eta_stress(sigma_eff_ts); //Here as well
        
        //compute Phi and the derivatives
        Phi_p(0) = Phi_p_ts;
        
        denom_p(0, 0) = -1.*sum(dPhi_p_tsd_sigma % (L*Lambdap_ts)) + dPhi_p_tsd_p;
        
        Y_pcrit(0) = sigma_ts_0;
        
        Fischer_Burmeister_m(Phi_p, Y_pcrit, denom_p, Dp, dp, error);
        
        p_ts += dp(0);
        
        EP = EP + dp(0)*Lambdap_ts;
        Eel = Etot + DEtot - alpha*(T+DT-Tinit) - EP;
    }
    
    error = 1.;
    if (Y_t > Y_22_u) {
        d_12 = 0.99;
        d_22 = 0.99;
        error = 0.;
    }
    
    //So it is forced to enter the damage loop once
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        mat L_tilde = L_ortho(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23, "EnuG");
        
        //damaged modulus
        //Compute the elastic strain and the related stress
        Eel = Etot + DEtot - alpha*(T+DT-Tinit) - EP;
        
        if (ndi == 1) {
            sigma(0) = sigma_start(0) + L_tilde(0,0)*(DEel_start(0));
        }
        else if (ndi == 2) {
            
            double Q11 = L_tilde(0,0)-L_tilde(0,2)*L_tilde(2,0)/L_tilde(2,2);
            double Q12 = L_tilde(0,1)-L_tilde(0,2)*L_tilde(2,1)/L_tilde(2,2);
            double Q14 = L_tilde(0,3)-L_tilde(0,2)*L_tilde(2,3)/L_tilde(2,2);
            double Q21 = L_tilde(1,0)-L_tilde(1,2)*L_tilde(2,0)/L_tilde(2,2);
            double Q22 = L_tilde(1,1)-L_tilde(1,2)*L_tilde(2,1)/L_tilde(2,2);
            double Q24 = L_tilde(1,3)-L_tilde(1,2)*L_tilde(2,3)/L_tilde(2,2);
            double Q41 = L_tilde(3,0)-L_tilde(3,2)*L_tilde(2,0)/L_tilde(2,2);
            double Q42 = L_tilde(3,1)-L_tilde(3,2)*L_tilde(2,1)/L_tilde(2,2);
            double Q44 = L_tilde(3,3)-L_tilde(3,2)*L_tilde(2,3)/L_tilde(2,2);
            
            sigma(0) = sigma_start(0) + Q11*DEel_start(0) + Q12*DEel_start(1) + Q14*DEel_start(3);
            sigma(1) = sigma_start(1) + Q21*DEel_start(0) + Q22*DEel_start(1) + Q24*DEel_start(3);
            sigma(3) = sigma_start(3) + Q41*DEel_start(0) + Q42*DEel_start(1) + Q44*DEel_start(3);
        }
        else
            sigma = L_tilde*Eel;
        
        if (fabs(sigma(3)) < iota )
            Yd_12 = 0.;
        else
            Yd_12 = 0.5*(pow(sigma(3),2.)/(G12_0*pow(1.-d_12,2.)));

        if (fabs(sigma(4)) < iota )
            Yd_13 = 0.;
        else
            Yd_13 = 0.5*(pow(sigma(4),2.)/(G12_0*pow(1.-d_12,2.)));
        
        if (fabs(sigma(1)) < iota )
            Yd_22 = 0.;
        else
            Yd_22 = 0.5*(pow(Macaulay_p(sigma(1)),2.)/(E2_0*pow(1.-d_22,2.)));

        if (fabs(sigma(2)) < iota )
            Yd_33 = 0.;
        else
            Yd_33 = 0.5*(pow(Macaulay_p(sigma(2)),2.)/(E2_0*pow(1.-d_22,2.)));
        
        
        //Compute Y_t and Y_ts, which are used in the evolution equation of damage
        Y_ts = sqrt(Yd_12 + Yd_13 + b*(Yd_22 + Yd_33));     //transvers/shear coupling
        Y_t = sqrt(Yd_22 + Yd_33);               //transverse
        
        lambda_12 = lagrange_pow_1(d_12, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        lambda_22 = lagrange_pow_1(d_22, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        
        //Preliminaries to compute damage the derivatives of Phi
        if (fabs(sigma(3)) < iota )
            dYd_12dd = 0.;
        else
            dYd_12dd = -0.25*(pow(sigma(3),2.)/(G12_0*pow(1-d_12,3.)));

        if (fabs(sigma(4)) < iota )
            dYd_13dd = 0.;
        else
            dYd_13dd = -0.25*(pow(sigma(4),2.)/(G12_0*pow(1-d_12,3.)));

                    
        if (fabs(sigma(1)) < iota )
            dYd_22dd = 0.;
        else
            dYd_22dd = -0.25*(pow(Macaulay_p(sigma(1)),2.)/(E2_0*pow(1-d_22,3.)));


        if (fabs(sigma(2)) < iota )
            dYd_33dd = 0.;
        else
            dYd_33dd = -0.25*(pow(Macaulay_p(sigma(2)),2.)/(E2_0*pow(1-d_22,3.)));
        
        if (fabs(Yd_12 + Yd_13 + b*(Yd_22 + Yd_33)) < iota ) {
            dY_tsdd_12 = 0.;
            dY_tsdd_22 = 0.;
        }
        else {
            dY_tsdd_12 = 0.5*(dYd_12dd + dYd_13dd)*pow(Yd_12 + Yd_13+ b*(Yd_22+Yd_33),-0.5);
            dY_tsdd_22 = 0.5*(b*(dYd_22dd+dYd_22dd))*pow(Yd_12 + Yd_13+ b*(Yd_22+Yd_33),-0.5);
        }
                    
        dlambda_12 = dlagrange_pow_1(d_12, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dlambda_22 = dlagrange_pow_1(d_22, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        
        //compute Phi and the derivatives
        Phi_d(0) = Macaulay_p(Y_ts - Y_12_0)/Y_12_c - lambda_12 - d_12;
        Phi_d(1) = Macaulay_p(Y_ts - Y_22_0)/Y_22_c - lambda_22 - d_22;
        
        denom_d(0, 0) = Macaulay_p(dY_tsdd_12)/Y_12_c - dlambda_12 - 1.;
        denom_d(0, 1) = Macaulay_p(dY_tsdd_22)/Y_12_c;
        
        denom_d(1, 0) = Macaulay_p(dY_tsdd_12)/Y_22_c;
        denom_d(1, 1) = Macaulay_p(dY_tsdd_22)/Y_22_c - dlambda_22 - 1.;
        
        Fischer_Burmeister_m(Phi_d, Y_dcrit, denom_d, Dd, dd, error);
        
        d_12 += dd(0);
        d_22 += dd(1);
        
        //Transverse Damage term
        if (Y_t > Y_22_u) {
            d_12 = 0.99;
            d_22 = 0.99;
            error = 0.;
        }
        
        //Update constitutive parameters
        E2 = ET*(1.-d_22);
        E3 = ET*(1.-d_22);
        G12 = G12_0*(1.-d_12);
        G13 = G12_0*(1.-d_12);
    }
    
    //Update constitutive parameters
    E2 = ET*(1.-d_22);
    E3 = ET*(1.-d_22);
    G12 = G12_0*(1.-d_12);
    G13 = G12_0*(1.-d_12);
    
    mat L_tilde = L_ortho(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23, "EnuG");
    
    //damaged modulus
    //Compute the elastic strain and the related stress
    Eel = Etot + DEtot - alpha*(T+DT-Tinit) - EP;
    
    if (ndi == 1) {
        sigma(0) = sigma_start(0) + L_tilde(0,0)*(DEel_start(0));
    }
    else if (ndi == 2) {
        
        double Q11 = L_tilde(0,0)-L_tilde(0,2)*L_tilde(2,0)/L_tilde(2,2);
        double Q12 = L_tilde(0,1)-L_tilde(0,2)*L_tilde(2,1)/L_tilde(2,2);
        double Q14 = L_tilde(0,3)-L_tilde(0,2)*L_tilde(2,3)/L_tilde(2,2);
        double Q21 = L_tilde(1,0)-L_tilde(1,2)*L_tilde(2,0)/L_tilde(2,2);
        double Q22 = L_tilde(1,1)-L_tilde(1,2)*L_tilde(2,1)/L_tilde(2,2);
        double Q24 = L_tilde(1,3)-L_tilde(1,2)*L_tilde(2,3)/L_tilde(2,2);
        double Q41 = L_tilde(3,0)-L_tilde(3,2)*L_tilde(2,0)/L_tilde(2,2);
        double Q42 = L_tilde(3,1)-L_tilde(3,2)*L_tilde(2,1)/L_tilde(2,2);
        double Q44 = L_tilde(3,3)-L_tilde(3,2)*L_tilde(2,3)/L_tilde(2,2);
        
        sigma(0) = sigma_start(0) + Q11*DEel_start(0) + Q12*DEel_start(1) + Q14*DEel_start(3);
        sigma(1) = sigma_start(1) + Q21*DEel_start(0) + Q22*DEel_start(1) + Q24*DEel_start(3);
        sigma(3) = sigma_start(3) + Q41*DEel_start(0) + Q42*DEel_start(1) + Q44*DEel_start(3);
    }
    else
        sigma = L_tilde*Eel;
    
    //Compute the explicit flow direction
    Lambdap_ts = Theta_ts*eta_stress(sigma_eff_ts);
    
    mat B = L*inv(L_tilde);         //stress "localization factor" in damage
    
    vec kappamat0 = zeros(6);
    kappamat0 = B*Lambdap_ts;
    
    dPhi_p_tsd_sigma = B*Theta_ts*eta_stress(sigma_eff_ts);
    
    double Bhat = -1.*sum((dPhi_p_tsd_sigma) % (L*kappamat0)) + dPhi_p_tsd_p;

    double op = 0.;
    double delta = 1.;
    
    if(Dp(0) > iota)
        op = 1.;
    
    double Bbar = op*op*Bhat + delta*(1-op*op);

    double invBbar = 0.;
    double invBhat = 0.;
    invBbar = 1./Bbar;
    invBhat = op*op*invBbar;
    
    vec Pjay0 = zeros(6);
    Pjay0 = L*(invBhat*dPhi_p_tsd_sigma);
    
    Lt = L_tilde + L_tilde*(kappamat0*trans(Pjay0));
    
/*    if (Y_t > Y_22_u) {
        Lt = L_iso(1, 0.3, "Enu");
        sigma = Lt*Eel;
    }*/
    
    statev(0) = Tinit;
    statev(1) = d_12;
    statev(2) = d_22;
    
    statev(3) = p_ts;
               
    statev(4) = EP(0);
    statev(5) = EP(1);
    statev(6) = EP(2);
    statev(7) = EP(3);
    statev(8) = EP(4);
    statev(9) = EP(5);
    
    //Returning the energy
    vec DEel = Eel - Eel_start;
    vec Dsigma = sigma - sigma_start;
    
    double Dtde = 0.5*sum((sigma_start+sigma)%DEtot);
    double Dsse = sum(sigma_start%DEel) + 0.5*sum(Dsigma%DEel);
    
    sse += Dsse;
    spd += Dtde - Dsse;
}
    
} //namespace smart
