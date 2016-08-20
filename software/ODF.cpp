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

///@file ODF
///@brief ODF: get		//
//	for a homogeneous loading path, allowing repeatable steps
///@version 1.9

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Solver/read.hpp>
#include <smartplus/Libraries/Phase/phase_characteristics.hpp>
#include <smartplus/Libraries/Phase/read.hpp>
#include <smartplus/Libraries/Phase/write.hpp>
#include <smartplus/Libraries/Material/peak.hpp>
#include <smartplus/Libraries/Material/ODF.hpp>
#include <smartplus/Libraries/Material/ODF2Nphases.hpp>
#include <smartplus/Libraries/Material/read.hpp>

using namespace std;
using namespace arma;
using namespace smart;

int main() {

    int nphases_rve = 36;
    
    phase_characteristics rve_init;
    
    string umat_name = "MIPLN";
    vec props = {2,0};
    
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    
    rve_init.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    
    string path_data = "data";
    string inputfile;
    string outputfile;
    
    std::map<string, int> list_umat;
    list_umat = {{"MIHEN",100},{"MIMTN",101},{"MISCN",102},{"MIPCW",103},{"MIPLN",104}};
    
    //Here we read everything about the initial rve
    switch (list_umat[rve_init.sptr_matprops->umat_name]) {
            
        case 100: case 101: case 102: case 103: {
            rve_init.construct(2,1); //The rve is supposed to be mechanical only here
            
            inputfile = "Nellipsoids" + to_string(int(rve_init.sptr_matprops->props(1))) + ".dat";
            read_ellipsoid(rve_init, path_data, inputfile);
            break;
        }
        case 104: {
            rve_init.construct(1,1); //The rve is supposed to be mechanical only here
            
            inputfile = "Nlayers" + to_string(int(rve_init.sptr_matprops->props(1))) + ".dat";
            read_layer(rve_init, path_data, inputfile);
            break;
        }
    }
    
    double angle_min = 0.;
    double angle_max = 180.;
    
    ODF odf_rve(0, false, angle_min, angle_max);
    read_peak(odf_rve, path_data, "Npeaks0.dat");
    
    vec x = linspace<vec>(0,  178,  90);
    cout << "x = " << x.t() << endl;
    
    vec y = get_densities(x, path_data, "Npeaks0.dat", false);
    cout << "y = " << y.t() << endl;
    
    phase_characteristics rve = discretize_ODF(rve_init, odf_rve, 1, nphases_rve,0);

    if(rve.shape_type == 0) {
        outputfile = "Nphases1.dat";
        write_phase(rve, path_data, outputfile);
    }
    if(rve.shape_type == 1) {
        outputfile = "Nlayers1.dat";
        write_layer(rve, path_data, outputfile);
    }
    else if(rve.shape_type == 2) {
        outputfile = "Nellipsopids1.dat";
        write_ellipsoid(rve, path_data, outputfile);
    }
    else if(rve.shape_type == 3) {
        outputfile = "Ncylinders1.dat";
        write_cylinder(rve, path_data, outputfile);
    }

	return 0;
}
