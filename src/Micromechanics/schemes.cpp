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

///@file schemes.cpp
///@brief micromechanical schemes for non-linear N-phases heterogeneous materials:
///@brief Mori-Tanaka scheme
///@version 1.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <smartplus/Libraries/Maths/rotation.hpp>
#include <smartplus/Libraries/Phase/phase_characteristics.hpp>
#include <smartplus/Libraries/Phase/state_variables.hpp>
#include <smartplus/Libraries/Phase/state_variables_M.hpp>
#include <smartplus/Libraries/Geometry/geometry.hpp>
#include <smartplus/Libraries/Geometry/layer.hpp>
#include <smartplus/Libraries/Geometry/ellipsoid.hpp>
#include <smartplus/Libraries/Homogenization/phase_multi.hpp>
#include <smartplus/Libraries/Homogenization/layer_multi.hpp>
#include <smartplus/Libraries/Homogenization/ellipsoid_multi.hpp>
#include <smartplus/Libraries/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;

namespace smart{
  
void Homogeneous_E(phase_characteristics &phase) {
    
    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        r.sptr_multi->A = eye(6,6);
    }
}
    
void Mori_Tanaka(phase_characteristics &phase) {
    
    mat sumT = zeros(6,6);
    mat inv_sumT = zeros(6,6);

    std::shared_ptr<ellipsoid_multi> elli_multi;
    std::shared_ptr<ellipsoid> elli;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_0 = std::dynamic_pointer_cast<state_variables_M>(phase.sub_phases[0].sptr_sv_global);
    std::shared_ptr<state_variables_M> sv_r;
    
    //Compute the Eshelby tensor and the interaction tensor for each phase
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(r.sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        elli_multi->fillT(sv_0->Lt, sv_r->Lt, *elli);
        //Compute the normalization interaction tensir sumT
		sumT += elli->concentration*elli_multi->T;
    }
				
    inv_sumT = inv(sumT);
    
    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        elli_multi->A = elli_multi->T*inv_sumT;
    }
}
    
void Self_Consistent(phase_characteristics &phase, const bool &start, const int &option_start = 0) {
    
    std::shared_ptr<ellipsoid_multi> elli_multi;
    std::shared_ptr<ellipsoid> elli;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_r;
    std::shared_ptr<state_variables_M> sv_eff = std::dynamic_pointer_cast<state_variables_M>(phase.sptr_sv_local);

    //In the self_consistent scheme we need to have the effective tangent modulus first, based on some guessed initial concentration tensor.
    if(start) {
        //ensure that the Lt of phase is actually 0
        if(option_start == 0)
            Homogeneous_E(phase);
        else if(option_start == 1)
            Mori_Tanaka(phase);
        else {
            cout << "error , option is not valid for the start option of Self-Consistent scheme (0 : MT, 1 : h_E)";
        }
        
        //Compute the effective tensor from the previous strain localization tensors
        mat Lt_eff = zeros(6,6);
        for(auto r : phase.sub_phases) {
            sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
            Lt_eff += r.sptr_shape->concentration*r.sptr_multi->A*sv_r->Lt;
        }
        sv_eff->Lt = Lt_eff;
    }
    
    //Compute the Eshelby tensor and the interaction tensor for each phase
    for(auto r : phase.sub_phases) {
        elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(r.sptr_multi);
        elli = std::dynamic_pointer_cast<ellipsoid>(r.sptr_shape);
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        
        //Note The tangent modulus are turned in the coordinate system of the ellipspoid in the fillT function
        elli_multi->fillT(sv_eff->Lt, sv_r->Lt, *elli);
        
        //Compute the strain concentration tensor A
        elli_multi->A = elli_multi->T;
    }

}
    
void Periodic_Layer(phase_characteristics &phase) {
    
    std::shared_ptr<layer_multi> lay_multi;
    std::shared_ptr<layer> lay;
    //ptr on the matrix properties
    std::shared_ptr<state_variables_M> sv_r;
    
    //We need to get the axis
    
    //Compute the strain concentration tensor A
    for(auto r : phase.sub_phases) {
        
        sv_r = std::dynamic_pointer_cast<state_variables_M>(r.sptr_sv_global);
        lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);
        
        lay_multi->Dnn(0,0) = sv_r->Lt(0,0);
        lay_multi->Dnn(0,1) = sv_r->Lt(0,3);
        lay_multi->Dnn(0,2) = sv_r->Lt(0,4);
        lay_multi->Dnn(1,0) = sv_r->Lt(3,0);
        lay_multi->Dnn(1,1) = sv_r->Lt(3,3);
        lay_multi->Dnn(1,2) = sv_r->Lt(3,4);
        lay_multi->Dnn(2,0) = sv_r->Lt(4,0);
        lay_multi->Dnn(2,1) = sv_r->Lt(4,3);
        lay_multi->Dnn(2,2) = sv_r->Lt(4,4);
        
        lay_multi->Dnt(0,0) = sv_r->Lt(0,1);
        lay_multi->Dnt(0,1) = sv_r->Lt(0,2);
        lay_multi->Dnt(0,2) = sv_r->Lt(0,5);
        lay_multi->Dnt(1,0) = sv_r->Lt(3,1);
        lay_multi->Dnt(1,1) = sv_r->Lt(3,2);
        lay_multi->Dnt(1,2) = sv_r->Lt(3,5);
        lay_multi->Dnt(2,0) = sv_r->Lt(4,1);
        lay_multi->Dnt(2,1) = sv_r->Lt(4,2);
        lay_multi->Dnt(2,2) = sv_r->Lt(4,5);
    }
    
    mat sumDnn = zeros(3,3);
    mat sumDnt = zeros(3,3);
    for(auto r : phase.sub_phases) {
        lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);
        lay = std::dynamic_pointer_cast<layer>(r.sptr_shape);
        sumDnn += lay->concentration*inv(lay_multi->Dnn);
        sumDnt += lay->concentration*inv(lay_multi->Dnn)*lay_multi->Dnt;
    }
    mat m_n = inv(sumDnn);
    mat m_t = m_n*sumDnt;

    for(auto r : phase.sub_phases) {
        lay_multi = std::dynamic_pointer_cast<layer_multi>(r.sptr_multi);
        lay_multi->dXn = inv(lay_multi->Dnn)*(m_n-lay_multi->Dnn);
        lay_multi->dXt = inv(lay_multi->Dnn)*(m_t-lay_multi->Dnt);

        lay_multi->A = eye(6,6);
        
        lay_multi->A(0,0) += lay_multi->dXn(0,0);
        lay_multi->A(0,1) += lay_multi->dXt(0,0);
        lay_multi->A(0,2) += lay_multi->dXt(0,1);
        lay_multi->A(0,3) += lay_multi->dXn(0,1);
        lay_multi->A(0,4) += lay_multi->dXn(0,2);
        lay_multi->A(0,5) += lay_multi->dXt(0,2);
        
        lay_multi->A(3,0) += lay_multi->dXn(1,0);
        lay_multi->A(3,1) += lay_multi->dXt(1,0);
        lay_multi->A(3,2) += lay_multi->dXt(1,1);
        lay_multi->A(3,3) += lay_multi->dXn(1,1);
        lay_multi->A(3,4) += lay_multi->dXn(1,2);
        lay_multi->A(3,5) += lay_multi->dXt(1,2);
        
        lay_multi->A(4,0) += lay_multi->dXn(2,0);
        lay_multi->A(4,1) += lay_multi->dXt(2,0);
        lay_multi->A(4,2) += lay_multi->dXt(2,1);
        lay_multi->A(4,3) += lay_multi->dXn(2,1);
        lay_multi->A(4,4) += lay_multi->dXn(2,2);
        lay_multi->A(4,5) += lay_multi->dXt(2,2);
    }
    
}
    
    
} //namespace smart
