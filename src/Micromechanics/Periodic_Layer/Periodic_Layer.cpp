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

///@file periodic_layer.cpp
///@brief User subroutine for non-linear periodic layers
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <smartplus/Micromechanics/Periodic_Layer/Periodic_Layer.hpp>
#include <smartplus/parameter.hpp>
#include <smartplus/Umat/umat_smart.hpp>
#include <smartplus/Libraries/Homogenization/layer_characteristics.hpp>

using namespace std;
using namespace arma;

namespace smart{

///@brief props[0] : Number of phases
///@brief props[1] : Number of the file NPhase[i].dat utilized
///@brief props[2] : Number of integration points in the 1 direction
///@brief props[3] : Number of integration points in the 2 direction

///@brief The table phases.dat will store the necessary informations about the geometry of the phases and the material properties

void umat_PL_N(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &sse, double &spd, const int &ndi, const int &nshr, const bool &start)
{
    
    double tnew_dt = 1.;
    
	//This will create the N phase object, that stores the props, stran, Dstran, stress and stetev for each phase:
	int nlayers = props(0); // Number of phases
	int filenumber = props(1); //file # that stores the microstructure properties
		
	std::vector<layer_characteristics> layersvs(nlayers); // This is the pointor on the state_variable object
    
	///@brief Properties of the phases reading, use "phases.dat" to specify the parameters of each phase
    
	string buffer;
	stringstream sstm;
	const char *car;
	string sar;
	
	sstm << "data/Nlayers" << filenumber << ".dat";
	sar = sstm.str();
	sstm.str(std::string());
	car = sar.c_str();

	///@brief Properties of the phases reading, use "phases.dat" to specify the parameters of each phase
	ifstream paramlayers;
	paramlayers.open(car, ios::in);
	if(paramlayers) {
		string chaine1;
		paramlayers >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1;
		for(int i=0; i<nlayers; i++) {
			int nprops, nstatev;
			paramlayers >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> nprops >> nstatev;
			layersvs[i].resize(nprops,nstatev);
			for(int j=0; j<nprops; j++) {
				paramlayers >> chaine1;
			}			
		}		
	}
	else {
		cout << "Error: cannot open phases.dat file \n";
	}
	
	paramlayers.close();
	
	paramlayers.open(car, ios::in);
	if(paramlayers) {
		string chaine1;		
		paramlayers >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1;
		for(int i=0; i<nlayers; i++) {  
			paramlayers >> layersvs[i].number >> layersvs[i].umat_name >> layersvs[i].concentration >> layersvs[i].psi_mat >> layersvs[i].theta_mat >> layersvs[i].phi_mat >> chaine1 >> chaine1;
			for(int j=0; j<layersvs[i].dimprops(); j++) {
				paramlayers >> layersvs[i].props[j];
			}
			layersvs[i].psi_mat*=(pi/180.);
			layersvs[i].theta_mat*=(pi/180.);
			layersvs[i].phi_mat*=(pi/180.);
		}
	}
		
	paramlayers.close();
	
	//Statev of all the layers
	int Sum_nstatev_i = 0;
	for (int i=0; i<nlayers; i++) {

		for (int l=0; l<6; l++) {
			layersvs[i].global.Etot(l) = statev(l + Sum_nstatev_i);
			layersvs[i].global.DEtot(l) = statev(l + 6 + Sum_nstatev_i);
			layersvs[i].global.sigma(l) = statev(l + 12 + Sum_nstatev_i);
		}

		for(int m=0; m<6; m++) {
			for(int n=0; n<6; n++) {
				layersvs[i].global.L(m,n) = statev(18 + m*6 + n + Sum_nstatev_i);
				layersvs[i].global.Lt(m,n) = statev(54 + m*6 + n + Sum_nstatev_i);
			}
		}
		
		for (int k=0; k<layersvs[i].nstatev; k++) {
			layersvs[i].statev(k) = statev(90 + Sum_nstatev_i + k);
 		}	
		Sum_nstatev_i += layersvs[i].nstatev + 90;
	}

	std::vector<mat> Dnn(nlayers);
	std::vector<mat> Dtn(nlayers);
	std::vector<mat> Dtt(nlayers);		
	std::vector<vec> sigma_hat(nlayers);
	std::vector<vec> dzdx1(nlayers);		
	
	for (int i=0; i<nlayers; i++) {
		Dnn[i] = zeros(3,3);
		Dtn[i] = zeros(3,3);
		Dtt[i] = zeros(3,3);
		sigma_hat[i] = zeros(3);
		dzdx1[i] = zeros(3);		
		layersvs[i].global.DEtot = DEtot;
	}			
	
	mat Dnn_eff = zeros(3,3);
	mat Dtn_eff = zeros(3,3);
	mat	Dtt_eff = zeros(3,3);
	
	mat sumDnn = zeros(3,3);
	mat sumDtn = zeros(3,3);
	mat	sumDtt = zeros(3,3);	
	
	vec m = zeros(3);	
	vec sumcDsig = zeros(3);		

	//Initial rotation of the state variables (from global to local)
	for (int i=0; i<nlayers; i++) {
		layersvs[i].global2local();
	}

	std::vector<layer_characteristics> layersvs_start; //Table that stores all initial states variables of the RVE
	layersvs_start = layersvs;
    
	//Initialization
    for (int i=0; i<nlayers; i++) {
    
        //Run the appropriate constitutive model		
        select_umat(layersvs[i].umat_name, layersvs[i].local.Etot, layersvs[i].local.DEtot, layersvs[i].local.sigma, layersvs[i].local.Lt, DR, layersvs[i].nprops, layersvs[i].props, layersvs[i].nstatev, layersvs[i].statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
        //Initialization of the elastic stiffness tensor
        layersvs[i].local.L = layersvs[i].local.Lt;
    }

	// ###################### Preliminaries of the convergence loop #################################
	int nbiter = 0;
	double error = 1.;
            
	//################################################
	// ********** Convergence Loop ******************/
	//################################################

	while ((error > precision_micro)&&(nbiter <= maxiter_micro)) {

		for(int i=0; i<nlayers; i++) {

			//Compute the derivative of z with respect to x1 in the global
			layersvs[i].local2global();
			Dnn[i](0,0) = layersvs[i].global.Lt(0,0); 
			Dnn[i](0,1) = layersvs[i].global.Lt(0,3);	
			Dnn[i](0,2) = layersvs[i].global.Lt(0,4); 
			Dnn[i](1,0) = layersvs[i].global.Lt(3,0);
			Dnn[i](1,1) = layersvs[i].global.Lt(3,3);			
			Dnn[i](1,2) = layersvs[i].global.Lt(3,4);			
			Dnn[i](2,0) = layersvs[i].global.Lt(4,0);
			Dnn[i](2,1) = layersvs[i].global.Lt(4,3);
			Dnn[i](2,2) = layersvs[i].global.Lt(4,4);			
		
			sigma_hat[i](0) = layersvs[i].global.sigma(0);
			sigma_hat[i](1) = layersvs[i].global.sigma(3);
			sigma_hat[i](2) = layersvs[i].global.sigma(4);
		}
		
		sumDnn = zeros(3,3);
		for(int i=0; i<nlayers; i++) {
			sumDnn += layersvs[i].concentration*inv(Dnn[i]);
		}
		
		sumcDsig = zeros(3);
		for(int i=0; i<nlayers; i++) {
			sumcDsig += layersvs[i].concentration*inv(Dnn[i])*sigma_hat[i];
		}		
		
		m = inv(sumDnn)*sumcDsig;

		#pragma omp parallel for        
		for(int i=0; i<nlayers; i++) {

			dzdx1[i] = inv(Dnn[i])*(m-sigma_hat[i]);
			
			layersvs[i].global.DEtot(0) += dzdx1[i](0);
			layersvs[i].global.DEtot(3) += dzdx1[i](1);
			layersvs[i].global.DEtot(4) += dzdx1[i](2);	
            
            layersvs[i].statev = layersvs_start[i].statev;
            
			//Rotation From global to local (in fact, just for DEtotV)
			layersvs[i].global2local();
			layersvs[i].local.sigma = layersvs_start[i].local.sigma;
			
			//Run the appropriate constitutive model
			select_umat(layersvs[i].umat_name, layersvs[i].local.Etot, layersvs[i].local.DEtot, layersvs[i].local.sigma, layersvs[i].local.Lt, DR, layersvs[i].nprops, layersvs[i].props, layersvs[i].nstatev, layersvs[i].statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
			
			//From local to global
			layersvs[i].local2global();								
		}
			
		//	Definition of the error
		error = 0.;
		for(int i=0; i<nlayers; i++) {		
			error += norm(dzdx1[i],2);
		}
		nbiter++;
	}

//	######################  UMAT -> ABAQUS #################################

	//Return the updated strain for each phase
	for(int i=0; i<nlayers; i++) {
		layersvs[i].global.Etot += layersvs[i].global.DEtot;
	}

	//Compute the effective stress
	sigma=0.*sigma;
	for(int i=0; i<nlayers; i++) {
		sigma += layersvs[i].concentration*layersvs[i].global.sigma;
	}

	// Compute the effective tangent modulus

	for(int i=0; i<nlayers; i++) {
		Dnn[i](0,0) = layersvs[i].global.Lt(0,0);
		Dnn[i](0,1) = layersvs[i].global.Lt(0,3);	
		Dnn[i](0,2) = layersvs[i].global.Lt(0,4);
		Dnn[i](1,0) = layersvs[i].global.Lt(3,0);
		Dnn[i](1,1) = layersvs[i].global.Lt(3,3);			
		Dnn[i](1,2) = layersvs[i].global.Lt(3,4);			
		Dnn[i](2,0) = layersvs[i].global.Lt(4,0);
		Dnn[i](2,1) = layersvs[i].global.Lt(4,3);
		Dnn[i](2,2) = layersvs[i].global.Lt(4,4);			

		Dtn[i](0,0) = layersvs[i].global.Lt(1,0);
		Dtn[i](0,1) = layersvs[i].global.Lt(1,3);	
		Dtn[i](0,2) = layersvs[i].global.Lt(1,4);
		Dtn[i](1,0) = layersvs[i].global.Lt(2,0);
		Dtn[i](1,1) = layersvs[i].global.Lt(2,3);			
		Dtn[i](1,2) = layersvs[i].global.Lt(2,4);			
		Dtn[i](2,0) = layersvs[i].global.Lt(5,0);
		Dtn[i](2,1) = layersvs[i].global.Lt(5,3);
		Dtn[i](2,2) = layersvs[i].global.Lt(5,4);
		
		Dtt[i](0,0) = layersvs[i].global.Lt(1,1);
		Dtt[i](0,1) = layersvs[i].global.Lt(1,2);	
		Dtt[i](0,2) = layersvs[i].global.Lt(1,5);
		Dtt[i](1,0) = layersvs[i].global.Lt(2,1);
		Dtt[i](1,1) = layersvs[i].global.Lt(2,2);			
		Dtt[i](1,2) = layersvs[i].global.Lt(2,5);			
		Dtt[i](2,0) = layersvs[i].global.Lt(5,1);
		Dtt[i](2,1) = layersvs[i].global.Lt(5,2);
		Dtt[i](2,2) = layersvs[i].global.Lt(5,5);	
	}
	
	sumDnn = zeros(3,3);
	for(int i=0; i<nlayers; i++) {
		sumDnn += layersvs[i].concentration*inv(Dnn[i]);
	}
	Dnn_eff = inv(sumDnn);		
	
	sumDtn = zeros(3,3);
	for(int i=0; i<nlayers; i++) {
		sumDtn += layersvs[i].concentration*Dtn[i]*inv(Dnn[i]);
	}
	Dtn_eff = sumDtn*Dnn_eff;

	sumDtt = zeros(3,3);
	for(int i=0; i<nlayers; i++) {
		sumDtt += layersvs[i].concentration*(Dtn[i]*inv(Dnn[i])*(trans(Dtn_eff)-trans(Dtn[i]))+Dtt[i]);
	}
	Dtt_eff = sumDtt;	
	
	Lt(0,0) = Dnn_eff(0,0);
	Lt(0,1) = Dtn_eff(0,0);		
	Lt(0,2) = Dtn_eff(1,0);		
	Lt(0,3) = Dnn_eff(0,1);
	Lt(0,4) = Dnn_eff(0,2);	
	Lt(0,5) = Dtn_eff(2,0);
	
	Lt(1,0) = Dtn_eff(0,0);
	Lt(1,1) = Dtt_eff(0,0);		
	Lt(1,2) = Dtt_eff(0,1);		
	Lt(1,3) = Dtn_eff(0,1);
	Lt(1,4) = Dtn_eff(0,2);	
	Lt(1,5) = Dtt_eff(0,2);
	
	Lt(2,0) = Dtn_eff(1,0);
	Lt(2,1) = Dtt_eff(1,0);		
	Lt(2,2) = Dtt_eff(1,1);		
	Lt(2,3) = Dtn_eff(1,1);
	Lt(2,4) = Dtn_eff(1,2);	
	Lt(2,5) = Dtt_eff(1,2);
	
	Lt(3,0) = Dnn_eff(1,0);
	Lt(3,1) = Dtn_eff(0,1);		
	Lt(3,2) = Dtn_eff(1,1);		
	Lt(3,3) = Dnn_eff(1,1);
	Lt(3,4) = Dnn_eff(1,2);	
	Lt(3,5) = Dtn_eff(2,1);
	
	Lt(4,0) = Dnn_eff(2,0);
	Lt(4,1) = Dtn_eff(0,2);		
	Lt(4,2) = Dtn_eff(1,2);		
	Lt(4,3) = Dnn_eff(2,1);
	Lt(4,4) = Dnn_eff(2,2);	
	Lt(4,5) = Dtn_eff(2,2);
	
	Lt(5,0) = Dtn_eff(2,0);
	Lt(5,1) = Dtt_eff(2,0);		
	Lt(5,2) = Dtt_eff(2,1);		
	Lt(5,3) = Dtn_eff(2,1);
	Lt(5,4) = Dtn_eff(2,2);	
	Lt(5,5) = Dtt_eff(2,2);
	
	//Return the statev
	Sum_nstatev_i=0;

	for (int i=0; i<nlayers; i++) {

		for (int l=0; l<6; l++) {
			statev(l + Sum_nstatev_i) = layersvs[i].global.Etot(l);
			statev(l + 6 + Sum_nstatev_i) = layersvs[i].global.DEtot(l);
			statev(l + 12 + Sum_nstatev_i) = layersvs[i].global.sigma(l);
		}

		for(int m=0; m<6; m++) {
			for(int n=0; n<6; n++) {
				statev(18 + m*6 + n + Sum_nstatev_i) = layersvs[i].global.L(m,n);
				statev(54 + m*6 + n + Sum_nstatev_i) = layersvs[i].global.Lt(m,n);
			}
		}
		
		for (int k=0; k<layersvs[i].nstatev; k++) {
			statev(90 + Sum_nstatev_i + k) = layersvs[i].statev(k);
 		}	
		Sum_nstatev_i += layersvs[i].nstatev + 90;
	}
    
}

} //namespace smart
