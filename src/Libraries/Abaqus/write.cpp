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

///@file write.cpp
///@brief To write NphasesX.dat and NlayerX.dat files
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <smartplus/Libraries/Phase/phase_characteristics.hpp>
#include <smartplus/Libraries/Abaqus/section_characteristics.hpp>
#include <smartplus/Libraries/Abaqus/materials.hpp>
#include <smartplus/Libraries/Abaqus/steps.hpp>
#include <smartplus/Libraries/Abaqus/write.hpp>
#include <smartplus/Libraries/Solver/block.hpp>
#include <smartplus/Libraries/Phase/read.hpp>

using namespace std;
using namespace arma;

namespace smart{

void update_sections(section_characteristics &section_rve, const int &nsections, const int &nscale, const int &geom_type, const string &path_data)
{
    
    string umat_name_macro = "ABAPC";
    phase_characteristics rve;
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    vec props_rve = {(double)nsections,(double)nscale};
    
    rve.sptr_matprops->update(0, umat_name_macro, 1, psi_rve, theta_rve, phi_rve, props_rve.n_elem, props_rve);
    rve.construct(geom_type,1); //The rve is supposed to be mechanical only here
    string inputfile; //file # that stores the microstructure properties
    
    switch (geom_type) {
            
        case 0 : {
            //Definition from Nphases.dat
            inputfile = "Nphases" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_phase(rve, path_data, inputfile);
            break;
        }
        case 1: {
            //Definition from Nlayers.dat
            inputfile = "Nlayers" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_layer(rve, path_data, inputfile);
            break;
        }
        case 2: {
            //Definition from Nellipsoids.dat
            inputfile = "Nellipsoids" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_ellipsoid(rve, path_data, inputfile);
            break;
        }
        case 3: {
            //Definition from Ncylinders.dat
            inputfile = "Ncylinders" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_cylinder(rve, path_data, inputfile);
            break;
        }
    }
    
    int id = 99999;
    section_rve.update(umat_name_macro, id, rve);
}

void write_sections(section_characteristics &section_rve, const string &path_data, const string &outputfile) {
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream param_mats;
    
    param_mats.open(filename, ios::out);
    param_mats << "** ==================\n";
    param_mats << "** MATERIALS\n";
    param_mats << "** ==================\n";
    param_mats << "**" << endl;
    param_mats.close();
    
    for(unsigned int i=0; i<section_rve.sub_sections.size(); i++) {
        section_rve.sub_sections[i].abamat.write(path_data,outputfile,true);
    }
    
    param_mats.open(filename, ios::app);
    
    for(unsigned int i=0; i<section_rve.sub_sections.size(); i++) {
        param_mats << "*Solid Section, ElSet=" << section_rve.sub_sections[i].elset_name << ", Material=" << section_rve.sub_sections[i].abamat.umat_name << "-" << section_rve.sub_sections[i].abamat.id << endl;
    }
    param_mats.close();
}
    
void update_materials(std::vector<aba_material> &aba_mats, const int &nphases, const int &nscale, const int &geom_type, const string &path_data) {

    aba_mats.clear();
    phase_characteristics rve;
    string umat_name_macro = "ABAPC";
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    vec props_rve = {(double)nphases,(double)nscale};
    
    rve.sptr_matprops->update(0, umat_name_macro, 1, psi_rve, theta_rve, phi_rve, props_rve.n_elem, props_rve);
    rve.construct(geom_type,1); //The rve is supposed to be mechanical only here
    string inputfile; //file # that stores the microstructure properties
    
    switch (geom_type) {
            
        case 0 : {
            //Definition from Nphases.dat
            inputfile = "Nphases" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_phase(rve, path_data, inputfile);
            break;
        }
        case 1: {
            //Definition from Nlayers.dat
            inputfile = "Nlayers" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_layer(rve, path_data, inputfile);
            break;
        }
        case 2: {
            //Definition from Nellipsoids.dat
            inputfile = "Nellipsoids" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_ellipsoid(rve, path_data, inputfile);
            break;
        }
        case 3: {
            //Definition from Ncylinders.dat
            inputfile = "Ncylinders" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_cylinder(rve, path_data, inputfile);
            break;
        }
    }
    
    int id_mat = 0;
    
    for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
        aba_material temp;
        id_mat = 100000 + rve.sptr_matprops->props(1)*1000 + rve.sub_phases[i].sptr_matprops->number;
        temp.update(*rve.sub_phases[i].sptr_matprops, *rve.sub_phases[i].sptr_sv_global, id_mat);
        aba_mats.push_back(temp);
    }
    
}
    
void write_materials(std::vector<aba_material> &aba_mats, const string &path_data, const string &outputfile) {
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream param_mats;
    
    param_mats.open(filename, ios::out);
    param_mats << "** ==================\n";
    param_mats << "** MATERIALS\n";
    param_mats << "** ==================\n";
    param_mats << "**" << endl;
        
    for(unsigned int i=0; i<aba_mats.size(); i++) {
        aba_mats[i].write(path_data,outputfile,true);
    }
    param_mats.close();
}
    
void update_steps(std::vector<aba_step_meca> &aba_steps, const std::vector<block> &blocks, const bool &nlgeom) {

    string name_step_ini = "Step";
    string name_step = name_step_ini;
    int type = 0;
    
    for (unsigned int i=0; i<blocks.size(); i++) {
        for (int j=0; j<blocks[i].ncycle; j++) {
            for (int k=0; k<blocks[i].nstep; k++) {
                aba_step_meca temp;
                name_step = name_step_ini + to_string(i+1) + to_string(j+1) + to_string(k+1);
                shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[k]);
                if ((sptr_meca->Dn_init == 1.)&&(sptr_meca->Dn_mini == 1.)) {
                    type = 1;
                }
                else {
                    type = 0;
                }
                temp.update(*sptr_meca, name_step, nlgeom, type);
                aba_steps.push_back(temp);
            }
        }
    }
    
}
    
void write_steps(std::vector<aba_step_meca> &aba_steps, const double &temp_init, const string &path_data, const string &outputfile) {
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream param_mats;
    
    param_mats.open(filename, ios::out);

    param_mats << "*Boundary\n";
    param_mats << "CentreNode, 1, 1\n";
    param_mats << "CentreNode, 2, 2\n";
    param_mats << "CentreNode, 3, 3\n";
    
    param_mats << "**\n";
    param_mats << "**Initial Conditions, type=TEMPERATURE\n";
    param_mats << "**AllNodes, " << temp_init << "\n";
    param_mats << "**\n";
    
    param_mats << "** ==================\n";
    param_mats << "** STEPS\n";
    param_mats << "** ==================\n";
    param_mats << "**" << endl;
    
    for(unsigned int i=0; i<aba_steps.size(); i++) {
        aba_steps[i].write(path_data,outputfile,true);
    }
    param_mats.close();
}
    
} //namespace smart
