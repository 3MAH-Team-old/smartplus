/* This file is part of SMART+ private.
 
 Only part of SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file num_solve.cpp
///@brief random number generators
///@author Chemisky & Anagnostou
///@version 1.0
///@date 10/23/2014

#include <iostream>
#include <assert.h>
#include <armadillo>
#include <smartplus/parameter.hpp>

using namespace std;
using namespace arma;

/* Notice*/
// Recall that these methods are intended to solve an problem by means of the Kuhn-Tucker conditions
//This means : Phi <= 0 ; Dp >= 0; Phi*Dp = 0, with Dp = pdot * Dt, a non-zero time state

namespace smart{

void Newton_Raphon(const vec &Phi, const vec &Y_crit, const mat &denom, vec &Dp, vec &dp, double &error)
{
    int n=Phi.n_elem;
    
    for (int i=0; i<n; i++) {
        assert(fabs(Y_crit(i)) > 0.);
    }
    
    if (fabs(det(denom)) > limit) {
        dp = -1.*solve(denom,Phi);
    }
    else {
        dp = zeros(n);
    }
    
    //New update of the vector of the variable increment throughout the step
    Dp = Dp + dp;
    
    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(Phi(i))/fabs(Y_crit(i));
    }
    
}
    
void Fischer_Burmeister(const vec &Phi, const vec &Y_crit, const mat &denom, vec &Dp, vec &dp, double &error)
{
    
    int n=Phi.n_elem;
    
    for (int i=0; i<n; i++) {
        assert(fabs(Y_crit(n)) > 0);
    }
    
    vec FB = zeros(n);
    mat denomFB = zeros(n,n);
    mat delta = eye(n,n);
    
    //Normalized Fischer-Burmeister set of equations
    for (int i=0; i<n; i++) {
        
        if ((fabs(Phi(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
            FB(i) = sqrt(pow(Phi(i),2.) + pow(Dp(i),2.)) + Phi(i) - Dp(i);
        }
        else if(fabs(Phi(i)) > 0.) {
            FB(i) = sqrt(pow(Phi(i),2.)) + Phi(i);
        }
        else if(fabs(Dp(i)) > 0.) {
            FB(i) = sqrt(pow(Dp(i),2.)) - Dp(i);
        }
        else {
            FB(i) = 0.;
        }
        
        for (int j=0; j<n; j++) {
            if ((fabs(Phi(i)) > 0.)&&(fabs(Dp(i)) > 0.)) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.) + pow(Dp(i),2.)))+1.)*denom(i,j) + delta(i,j)*(Dp(i)/(sqrt(pow(Phi(i),2.) + pow(Dp(i),2.))) - 1.);
            }
            else if(fabs(Phi(i)) > 0.) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.)))+1.)*denom(i,j) - delta(i,j);
            }
            else if(fabs(Dp(i)) > 0.) {
                denomFB(i,j) = denom(i,j) + delta(i,j)*(Dp(i)/(sqrt(pow(Dp(i),2.))) - 1.);
            }
            else {
                denomFB(i,j) = 1.E12;
            }
        }
        
    }
    
    if (fabs(det(denomFB)) > limit) {
        dp = -1.*solve(denomFB,FB);
    }
    else {
        dp = zeros(n);
    }
    
    //New update of the vector of the variable increment throughout the step
    Dp = Dp + dp;
    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(FB(i))/fabs(Y_crit(i));
    }
    
}
    
void Fischer_Burmeister_m(const vec &Phi, const vec &Y_crit, const mat &denom, vec &Dp, vec &dp, double &error)
{
    int n=Phi.n_elem;
        
    for (int i=0; i<n; i++) {
        assert(fabs(Y_crit(i)) > 0);
    }
    
    vec FB = zeros(n);
    mat denomFB = zeros(n,n);
    mat delta = eye(n,n);
    
    //Determine the eigenvalues of denom
    vec factor_denom = zeros(n);
    for (int i=0; i<n; i++) {
        factor_denom(i) = fabs(denom(i,i));
    }
    
    vec Dpstar = Dp%factor_denom;
    vec dDpstar = factor_denom;
    
    //Normalized Fischer-Burmeister set of equations
    for (int i=0; i<n; i++) {
        
        if ((fabs(Phi(i)) > 0.)&&(fabs(Dpstar(i)) > 0.)) {
            FB(i) = sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.)) + Phi(i) - Dpstar(i);
        }
        else if(fabs(Phi(i)) > 0.) {
            FB(i) = sqrt(pow(Phi(i),2.)) + Phi(i);
        }
        else if(fabs(Dpstar(i)) > 0.) {
            FB(i) = sqrt(pow(Dpstar(i),2.)) - Dpstar(i);
        }
        else {
            FB(i) = 0.;
        }
        
        for (int j=0; j<n; j++) {
            if ((fabs(Phi(i)) > 0.)&&(fabs(Dpstar(i)) > 0.)) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.)))+1.)*denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.))) - 1.);
            }
            else if(fabs(Phi(i)) > 0.) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.)))+1.)*denom(i,j) - delta(i,j)*dDpstar(i);
            }
            else if(fabs(Dpstar(i)) > 0.) {
                denomFB(i,j) = denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Dpstar(i),2.))) - 1.);
            }
            else {
                if (i==j)
                    denomFB(i,j) = 1.E12;
                else
                    denomFB(i,j) = 0.;
            }
        }
        
    }
    
    if (fabs(det(denomFB)) > limit) {
        dp = -1.*solve(denomFB, FB);
    }
    else {
        dp = zeros(n);
    }
    
    //New update of the transformation/orientation multipliers
    Dp = Dp + dp;

    error = 0.;
    for (int i=0; i<n; i++) {
        error+= fabs(FB(i))/fabs(Y_crit(i));
    }
    
}
    
mat denom_FB_m(const vec &Phi, const mat &denom, const vec &Dp)
{
    
    int n=Phi.n_elem;
    
    vec FB = zeros(n);
    mat denomFB = zeros(n,n);
    mat delta = eye(n,n);
    
    //Determine the eigenvalues of denom
    vec factor_denom = zeros(n);
    for (int i=0; i<n; i++) {
        factor_denom(i) = fabs(denom(i,i));
    }
    
    vec Dpstar = Dp%factor_denom;
    vec dDpstar = factor_denom;
    
    //Normalized Fischer-Burmeister set of equations
    for (int i=0; i<n; i++) {
        
        for (int j=0; j<n; j++) {
            if ((fabs(Phi(i)) > 0.)&&(fabs(Dpstar(i)) > 0.)) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.)))+1.)*denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Phi(i),2.) + pow(Dpstar(i),2.))) - 1.);
            }
            else if(fabs(Phi(i)) > 0.) {
                denomFB(i,j) = (Phi(i)/(sqrt(pow(Phi(i),2.)))+1.)*denom(i,j) - delta(i,j)*dDpstar(i);
            }
            else if(fabs(Dpstar(i)) > 0.) {
                denomFB(i,j) = denom(i,j) + delta(i,j)*dDpstar(i)*(Dpstar(i)/(sqrt(pow(Dpstar(i),2.))) - 1.);
            }
            else {
                denomFB(i,j) = 1.E12;
            }
        }
        
    }
    return denomFB;
}
    
} //namespace smart
