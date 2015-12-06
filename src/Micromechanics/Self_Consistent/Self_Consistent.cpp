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

///@file Self-Consistent.hpp
///@brief User subroutine for non-linear N-phases heterogeneous materials using
///@brief Self-Consistent scheme
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <smartplus/Micromechanics/Self_Consistent/Self_Consistent.hpp>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Homogenization/eshelby.hpp>
#include <smartplus/Umat/umat_smart.hpp>
#include <smartplus/Libraries/Homogenization/ellipsoid_characteristics.hpp>

using namespace std;
using namespace arma;

namespace smart{

///@brief props[0] : Number of phases
///@brief props[1] : Number of the file NPhase[i].dat utilized
///@brief props[2] : Number of integration points in the 1 direction
///@brief props[3] : Number of integration points in the 2 direction

///@brief The table phases.dat will store the necessary informations about the geometry of the phases and the material properties

void umat_SC_N(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &sse, double &spd, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{	
	//This will create the N phase object, that stores the props, stran, Dstran, stress and stetev for each phase:
	
    UNUSED(Etot);
    UNUSED(nprops);
    UNUSED(nstatev);
    
	//This will create the N phase object, that stores the props, stran, Dstran, stress and stetev for each phase:
	int nphases = props(0); // Number of phases
	int filenumber = props(1); //file # that stores the microstructure properties
	
	int nItg1 = props(2);
	int nItg2 = props(3);
	
	std::vector<ellipsoid_characteristics> rvesvs(nphases); // vector of ellipsoid phases
    
	///@brief Properties of the phases reading, use "phases.dat" to specify the parameters of each phase
    
	string buffer;
	stringstream sstm;
	const char *car;
	string sar;
	
	sstm << "data/Nphases" << filenumber << ".dat";
	sar = sstm.str();
	sstm.str(std::string());
	car = sar.c_str();
    
	///@brief Properties of the phases reading, use "phases.dat" to specify the parameters of each phase
	ifstream paramphases;
	paramphases.open(car, ios::in);
	if(paramphases) {
		string chaine1;
		paramphases >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1;
		for(int i=0; i<nphases; i++) {
			int nprops, nstatev;
			paramphases >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> nprops >> nstatev;
			rvesvs[i].resize(nprops,nstatev);
			for(int j=0; j<nprops; j++) {
				paramphases >> chaine1;
			}
		}
	}
	else {
		cout << "Error: cannot open phases.dat file \n";
	}
	
	paramphases.close();
	
	paramphases.open(car, ios::in);
	if(paramphases) {
		string chaine1;
		paramphases >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1;
		for(int i=0; i<nphases; i++) {
			paramphases >> rvesvs[i].number >> rvesvs[i].coatingof >> rvesvs[i].umat_name >> rvesvs[i].concentration >> rvesvs[i].psi_mat >> rvesvs[i].theta_mat >> rvesvs[i].phi_mat >> rvesvs[i].a1 >> rvesvs[i].a2 >> rvesvs[i].a3 >> rvesvs[i].psi_geom >> rvesvs[i].theta_geom >> rvesvs[i].phi_geom >> chaine1 >> chaine1;
			for(int j=0; j<rvesvs[i].dimprops(); j++) {
				paramphases >> rvesvs[i].props[j];
			}
			rvesvs[i].psi_mat*=(pi/180.);
			rvesvs[i].theta_mat*=(pi/180.);
			rvesvs[i].phi_mat*=(pi/180.);
            
			rvesvs[i].psi_geom*=(pi/180.);
			rvesvs[i].theta_geom*=(pi/180.);
			rvesvs[i].phi_geom*=(pi/180.);
		}
	}
	paramphases.close();
	
	//Statev of all the phases
	int Sum_nstatev_i = 0;
	for (int i=0; i<nphases; i++) {

		for (int l=0; l<6; l++) {
			rvesvs[i].global.Etot(l) = statev(l + Sum_nstatev_i);
			rvesvs[i].global.DEtot(l) = statev(l + 6 + Sum_nstatev_i);
			rvesvs[i].global.sigma(l) = statev(l + 12 + Sum_nstatev_i);
		}

		for(int m=0; m<6; m++) {
			for(int n=0; n<6; n++) {
				rvesvs[i].global.L(m,n) = statev(18 + m*6 + n + Sum_nstatev_i);
				rvesvs[i].global.Lt(m,n) = statev(54 + m*6 + n + Sum_nstatev_i);
			}
		}
		
		for (int k=0; k<rvesvs[i].nstatev; k++) {
			rvesvs[i].statev(k) = statev(90 + Sum_nstatev_i + k);
 		}	
		Sum_nstatev_i += rvesvs[i].nstatev + 90;
	}

	//Initial rotation of the state variables (from global to local)
	for (int i=0; i<nphases; i++) {
		rvesvs[i].global2local();
	}

	//Initialization
	if (start) {
		for (int i=0; i<nphases; i++) {
		
			//Run the appropriate constitutive model		
			select_umat(rvesvs[i].umat_name, rvesvs[i].local.Etot, rvesvs[i].local.DEtot, rvesvs[i].local.sigma, rvesvs[i].local.Lt, DR, rvesvs[i].nprops, rvesvs[i].props, rvesvs[i].nstatev, rvesvs[i].statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
			//Initialization of the elastic stiffness tensor
			rvesvs[i].local.L = rvesvs[i].local.Lt;
            
		}
        
        // Compute the effective tangent modulus, necessary in the self_consistent scheme
        Lt = zeros(6,6);
        for(int i=0; i<nphases; i++) {
            rvesvs[i].local2global();
            rvesvs[i].A = eye(6,6);
            Lt += rvesvs[i].concentration*(rvesvs[i].global.Lt*rvesvs[i].A);
        }
        
	}
    
	std::vector<ellipsoid_characteristics> rvesvs_start; //Table that stores all initial states variables of the RVE
	rvesvs_start = rvesvs;

	// ###################### Preliminaries of the convergence loop #################################
	int nbiter = 0;
	double error = 1.;
	std::vector<vec> DEtot_N(nphases); //Table that stores all the previous increments of strain
        
	//################################################
	// ********** Convergence Loop ******************/
	//################################################
	
	//Definition of the Gauss integration points for the numerical evaluation of the Eshelby tensor
	int mp=nItg1;
	int np=nItg2;
	vec x(mp);
	vec wx(mp);
	vec y(np);
	vec wy(np);
	points(x, wx, y, wy, mp, np);
    mat Lt_local = zeros(6,6);
    
	while ((error > precision_micro)&&(nbiter <= maxiter_micro)) {
        
		for(int i=0; i<nphases; i++) {
            DEtot_N[i] = rvesvs[i].global.DEtot;
        }
        
		//Compute the Eshelby tensor and the interaction tensor for each phase
        #pragma omp parallel for
		for(int i=0; i<nphases; i++) {
			//Compute the local matrix tensor
            rvesvs[i].fillT(Lt, x, wx, y, wy, mp, np);
        }

		#pragma omp parallel for       
		for(int i=0; i<nphases; i++) {
		//Compute the Self-COnsistent tensor for each phase
		
			rvesvs[i].A = rvesvs[i].T;
			rvesvs[i].global.DEtot = rvesvs[i].A*DEtot;
            rvesvs[i].statev = rvesvs_start[i].statev;
            
			//Rotation From global to local (in fact, just for DEtotV)
			rvesvs[i].global2local();
			rvesvs[i].local.sigma = rvesvs_start[i].local.sigma;
			
			//Run the appropriate constitutive model
			select_umat(rvesvs[i].umat_name, rvesvs[i].local.Etot, rvesvs[i].local.DEtot, rvesvs[i].local.sigma, rvesvs[i].local.Lt, DR, rvesvs[i].nprops, rvesvs[i].props, rvesvs[i].nstatev, rvesvs[i].statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
			
			//This is the semi-implicit method with a alph=2./3.
			rvesvs[i].local.Lt = (1 - (2./3.))*rvesvs_start[i].local.Lt + (2./3.)*rvesvs[i].local.Lt;
			//rvesvs[i].local.Lt = rvesvs[i].local.Lt;			
			//From local to global
			rvesvs[i].local2global();			
		}
        
        Lt = zeros(6,6);
        // Compute the effective tangent modulus
        for(int i=0; i<nphases; i++) {
            Lt += rvesvs[i].concentration*(rvesvs[i].global.Lt*rvesvs[i].A);
        }
        
        error = 0.;
		for(int i=0; i<nphases; i++) {
            error += norm(DEtot_N[i] - rvesvs[i].global.DEtot,2)*(1./nphases);
        }
		error = sqrt(error);
        
		nbiter++;
	}

//	######################  UMAT -> ABAQUS #################################
    
	//Return the updated strain for each phase
	for(int i=0; i<nphases; i++) {
		rvesvs[i].global.Etot += rvesvs[i].global.DEtot;
	}

	//Compute the effective stress
	sigma=0.*sigma;
	for(int i=0; i<nphases; i++) {
		sigma += rvesvs[i].concentration*rvesvs[i].global.sigma;
	}

	// Compute the effective tangent modulus, and the effective stress
    Lt = zeros(6,6);
	for(int i=0; i<nphases; i++) {
		Lt += rvesvs[i].concentration*(rvesvs[i].global.Lt*rvesvs[i].A);
	}
	
	Sum_nstatev_i=0;

	for (int i=0; i<nphases; i++) {

		for (int l=0; l<6; l++) {
			statev(l + Sum_nstatev_i) = rvesvs[i].global.Etot(l);
			statev(l + 6 + Sum_nstatev_i) = rvesvs[i].global.DEtot(l);
			statev(l + 12 + Sum_nstatev_i) = rvesvs[i].global.sigma(l);
		}

		for(int m=0; m<6; m++) {
			for(int n=0; n<6; n++) {
				statev(18 + m*6 + n + Sum_nstatev_i) = rvesvs[i].global.L(m,n);
				statev(54 + m*6 + n + Sum_nstatev_i) = rvesvs[i].global.Lt(m,n);
			}
		}
		
		for (int k=0; k<rvesvs[i].nstatev; k++) {
			statev(90 + Sum_nstatev_i + k) = rvesvs[i].statev(k);
 		}	
		Sum_nstatev_i += rvesvs[i].nstatev + 90;
	}

}

} //namespace smart
