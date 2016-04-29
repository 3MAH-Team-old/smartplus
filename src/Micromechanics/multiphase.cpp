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

///@file multiphase.cpp
///@brief User subroutine for non-linear N-phases heterogeneous materials
///@version 1.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <smartplus/Libraries/Continuum_Mechanics/constitutive.hpp>
#include <smartplus/Micromechanics/multiphase.hpp>
#include <smartplus/parameter.hpp>
#include <smartplus/Umat/umat_smart.hpp>
#include <smartplus/Libraries/Phase/state_variables_M.hpp>
#include <smartplus/Libraries/Phase/read.hpp>
#include <smartplus/Libraries/Homogenization/ellipsoid_multi.hpp>
#include <smartplus/Libraries/Homogenization/eshelby.hpp>
#include <smartplus/Micromechanics/schemes.hpp>
#include <smartplus/Umat/umat_smart.hpp>

using namespace std;
using namespace arma;

namespace smart{
    
///@brief The elastic isotropic UMAT requires 4 constants, plus the number of constants for each materials:
///@brief props[0] : Number of phases
///@brief props[1] : File # that stores the microstructure properties
///@brief props[2] : Number of integration points in the 1 direction
///@brief props[3] : Number of integration points in the 2 direction

///@brief The table Nphases.dat will store the necessary informations about the geometry of the phases and the material properties

void umat_multi(phase_characteristics &phase, const mat &DR, const double &Time, const double &DTime, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &method)
{

    int nphases = phase.sptr_matprops->props(0); // Number of phases
    int filenumber = phase.sptr_matprops->props(1); //file # that stores the microstructure properties
    
    shared_ptr<state_variables_M> umat_phase_M = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local); //shared_ptr on state variables of the rve
    shared_ptr<state_variables_M> umat_sub_phases_M; //shared_ptr on state variables
    
    //1 - We need to figure out the type of geometry and read the phase
    if(start) {
        switch (method) {
                
            case 100: case 101: case 102: case 103: {
                //Definition of the static vectors x,wx,y,wy
                ellipsoid_multi::mp = phase.sptr_matprops->props(2);
                ellipsoid_multi::np = phase.sptr_matprops->props(3);
                ellipsoid_multi::x.set_size(ellipsoid_multi::mp);
                ellipsoid_multi::wx.set_size(ellipsoid_multi::mp);
                ellipsoid_multi::y.set_size(ellipsoid_multi::np);
                ellipsoid_multi::wy.set_size(ellipsoid_multi::np);
                points(ellipsoid_multi::x, ellipsoid_multi::wx, ellipsoid_multi::y, ellipsoid_multi::wy,ellipsoid_multi::mp, ellipsoid_multi::np);
                
                read_ellipsoid(phase, filenumber);
                break;
            }
            case 104: {
                read_layer(phase, filenumber);
                break;
            }
        }
    }
    
	//Initialization
	if (start) {
        
        for (auto r: phase.sub_phases) {
            //Run the appropriate constitutive model
            select_umat_M(r, DR, Time, DTime, ndi, nshr, start, tnew_dt);
        }
    }

    // Preliminaries of the convergence loop
    int nbiter = 0;
    double error = 1.;
    std::vector<vec> DEtot_N(nphases); //Table that stores all the previous increments of strain
    
	//Convergence loop, localization
	while ((error > precision_micro)&&(nbiter <= maxiter_micro)) {
	  
        for(int i=0; i<nphases; i++) {
            DEtot_N[i] = phase.sub_phases[i].sptr_sv_global->DEtot;
        }

        //Compute the strain concentration tensor for each phase:
        
        switch (method) {
                
            case 100: {
                Homogeneous_E(phase);
                break;
            }
            case 101: {
                Mori_Tanaka(phase);
                break;
            }
            case 102: {
                Self_Consistent(phase, start, 0);
                break;
            }
            case 104: {
                Periodic_Layer(phase);
                break;
            }
        
        }
    
        for (auto r : phase.sub_phases) {
            r.sptr_sv_global->DEtot = r.sptr_multi->A*umat_phase_M->DEtot; //Recall that the global coordinates of subphases is the local coordinates of the generic phase
            r.sptr_sv_global->to_start();

            //Theta method for the tangent modulus
            //mat Lt_start = umat_sub_phases_M->Lt
            select_umat_M(r, DR, Time, DTime, ndi, nshr, start, tnew_dt);

            //Theta method for the tangent modulus
            //umat_sub_phases_M = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
            //Lt* = (1 - (2./3.))*Lt_start + 2./3.*Lt;
        }
        
        error = 0.;
        for(int i=0; i<nphases; i++) {
            error += norm(DEtot_N[i] - phase.sub_phases[i].sptr_sv_global->DEtot,2);
        }
        error*=(1./nphases);
    
        nbiter++;
	}
    
    //	Homogenization
	//Compute the effective stress
	umat_phase_M->sigma = zeros(6);
    for (auto r : phase.sub_phases) {
		umat_phase_M->sigma += r.sptr_shape->concentration*r.sptr_sv_global->sigma;
	}
    
    umat_phase_M->Lt = zeros(6,6);
	// Compute the effective tangent modulus, and the effective stress
    for (auto r : phase.sub_phases) {
        umat_sub_phases_M = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
		umat_phase_M->Lt += r.sptr_shape->concentration*(umat_sub_phases_M->Lt*r.sptr_multi->A);
	}
    
}

} //namespace smart
