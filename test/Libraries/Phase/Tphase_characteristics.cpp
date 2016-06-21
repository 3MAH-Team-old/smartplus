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

///@file Tconstitutive.cpp
///@brief Test for Constitutive tensors in Voigt notation
///@version 1.0

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "rotation"
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iterator>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Phase/phase_characteristics.hpp>
#include <smartplus/Libraries/Phase/read.hpp>
#include <smartplus/Libraries/Phase/write.hpp>

using namespace std;
using namespace arma;
using namespace smart;

BOOST_AUTO_TEST_CASE( read_write )
{

    
    string umat_name;
//    int nprops = 2;
//    int nstatev = 0;
    vec props = {2,0};
    
    double rho = 1.12;
    double c_p = 1.68;
    
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;

    //Phases
    phase_characteristics rve_phase;
    rve_phase.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props, rho, c_p);
    rve_phase.construct(0,1); //The rve is supposed to be mechanical only here
    read_phase(rve_phase, 0);
    write_phase(rve_phase, 1);
    
    std::ifstream ifs1_phase("data/Nphases0.dat");
    std::ifstream ifs2_phase("data/Nphases1.dat");
    
    std::istream_iterator<char> b1_phase(ifs1_phase), e1_phase;
    std::istream_iterator<char> b2_phase(ifs2_phase), e2_phase;
        
    BOOST_CHECK_EQUAL_COLLECTIONS(b1_phase, e1_phase, b2_phase, e2_phase);
    
    //Layers
    phase_characteristics rve_layer;
    rve_layer.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props, rho, c_p);
    rve_layer.construct(1,1); //The rve is supposed to be mechanical only here
    read_layer(rve_layer, 0);
    write_layer(rve_layer, 1);
    
    std::ifstream ifs1_layer("data/Nlayers0.dat");
    std::ifstream ifs2_layer("data/Nlayers1.dat");
    
    std::istream_iterator<char> b1_layer(ifs1_layer), e1_layer;
    std::istream_iterator<char> b2_layer(ifs2_layer), e2_layer;
    
    BOOST_CHECK_EQUAL_COLLECTIONS(b1_layer, e1_layer, b2_layer, e2_layer);
    
    //Ellipsoid
    phase_characteristics rve_ellipsoid;
    rve_ellipsoid.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props, rho, c_p);
    rve_ellipsoid.construct(2,1); //The rve is supposed to be mechanical only here
    read_ellipsoid(rve_ellipsoid, 0);
    write_ellipsoid(rve_ellipsoid, 1);
    
    std::ifstream ifs1_ellipsoid("data/Nellipsoids0.dat");
    std::ifstream ifs2_ellipsoid("data/Nellipsoids1.dat");
    
    std::istream_iterator<char> b1_ellipsoid(ifs1_ellipsoid), e1_ellipsoid;
    std::istream_iterator<char> b2_ellipsoid(ifs2_ellipsoid), e2_ellipsoid;
    
    BOOST_CHECK_EQUAL_COLLECTIONS(b1_ellipsoid, e1_ellipsoid, b2_ellipsoid, e2_ellipsoid);
    
    //Cylinder
    phase_characteristics rve_cylinder;
    rve_cylinder.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props, rho, c_p);
    rve_cylinder.construct(3,1); //The rve is supposed to be mechanical only here
    read_cylinder(rve_cylinder, 0);
    write_cylinder(rve_cylinder, 1);

    cout << "La tete a Toto!" << endl;
    
    std::ifstream ifs1_cylinder("data/Ncylinders0.dat");
    std::ifstream ifs2_cylinder("data/Ncylinders1.dat");
    
    std::istream_iterator<char> b1_cylinder(ifs1_cylinder), e1_cylinder;
    std::istream_iterator<char> b2_cylinder(ifs2_cylinder), e2_cylinder;
    
    BOOST_CHECK_EQUAL_COLLECTIONS(b1_cylinder, e1_cylinder, b2_cylinder, e2_cylinder);
    
}
