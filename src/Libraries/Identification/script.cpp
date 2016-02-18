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

///@file script.cpp
///@brief Scripts that allows to run identification algorithms based on Smart+ Control functions
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <smartplus/Libraries/Identification/parameters.hpp>
#include <smartplus/Libraries/Identification/constants.hpp>
#include <smartplus/Libraries/Identification/generation.hpp>
#include <smartplus/Libraries/Identification/read.hpp>
#include <smartplus/Libraries/Identification/script.hpp>
#include <smartplus/Libraries/Solver/read.hpp>
#include <smartplus/Libraries/Solver/solver.hpp>

using namespace std;
using namespace arma;
using namespace arma;

namespace smart{

//This function will copy the parameters files
void copy_parameters(const vector<parameters> &params, const string &src_path, const string &dst_path) {

    string src_files;
    string dst_files;
    
    for (auto pa : params) {
        for(auto ifiles : pa.input_files) {
            src_files = src_path + ifiles;
            dst_files = dst_path + ifiles;
            boost::filesystem::copy_file(src_files,dst_files,boost::filesystem::copy_option::overwrite_if_exists);
        }
    }
}

//This function will copy the parameters files
void copy_constants(const vector<constants> &consts, const string &src_path, const string &dst_path) {
    
    string src_files;
    string dst_files;
    
    for (auto co : consts) {
        for(auto ifiles : co.input_files) {
            src_files = src_path + ifiles;
            dst_files = dst_path + ifiles;
            boost::filesystem::copy_file(src_files,dst_files,boost::filesystem::copy_option::overwrite_if_exists);
        }
    }
}
    
//This function will replace the keys by the parameters
void apply_parameters(const vector<parameters> &params, const string &dst_path) {
  
    string mod_files;
    string buffer;
    
    ifstream in_files;
    ofstream ou_files;
    
    for (auto pa : params) {
        for(auto ifiles : pa.input_files) {
            mod_files = dst_path + ifiles;

            in_files.open(mod_files, ios::in);
            
            std::vector<string> str;
            while (!in_files.eof())
            {
                getline(in_files,buffer);
                str.push_back(buffer);
            }
            in_files.close();
            
            ou_files.open(mod_files);
            for (auto s : str) {
                boost::replace_all(s, pa.key, to_string(pa.value));
                ou_files << s << "\n";
            }
            ou_files.close();
        }
    }

}

//This function will replace the keys by the parameters
void apply_constants(const vector<constants> &consts, const string &dst_path) {
    
    string mod_files;
    string buffer;
    
    ifstream in_files;
    ofstream ou_files;
    
    for (auto co : consts) {
        for(auto ifiles : co.input_files) {
            mod_files = dst_path + ifiles;
            
            in_files.open(mod_files, ios::in);
            
            std::vector<string> str;
            while (!in_files.eof())
            {
                getline(in_files,buffer);
                str.push_back(buffer);
            }
            in_files.close();
            
            ou_files.open(mod_files);
            for (auto s : str) {
                boost::replace_all(s, co.key, to_string(co.value));
                ou_files << s << "\n";
            }
            ou_files.close();
        }
    }
    
}
    
void launch_solver(const generation &genrun, const int &nfiles, vector<parameters> &params, vector<constants> &consts, const string &folder, const string &name, const string &ext)
{
	stringstream outputfile;
	stringstream pathfile;
    
	//#pragma omp parallel for private(sstm, path)
	for (int j = 0; j <genrun.nindividuals; j++) {
		for (int i = 0; i<nfiles; i++) {
			///Creating the right path & output filenames
			outputfile << folder << "/" << name << i+1 << "_" << j+1 << "." << ext;
			pathfile << "path_id_" << i+1 << ".txt";
            
            string umat_name;
            int nprops = 0;
            int nstatev = 0;
            vec props;
            
            double rho = 0.;
            double c_p = 0.;
            
            //Replace the constants
            for (unsigned int k=0; k<consts.size(); k++) {
                consts[k].value = consts[k].input_values(i);
            }
                        //Replace the parameters
            for (unsigned int k=0; k<params.size(); k++) {
                params[k].value = genrun.pop[j].p(k);
            }
            string data = "data/";
            string key = "data/key/";

            copy_constants(consts, key, data);
            copy_parameters(params, key, data);
            
            apply_constants(consts, data);
            apply_parameters(params, data);
            
            //Then read the material properties
            read_matprops(umat_name, nprops, props, nstatev, rho, c_p);
            
			///Launching the solver with relevant parameters
            solver(umat_name, props, nstatev, rho, c_p, pathfile.str(), outputfile.str());
            
            outputfile.str(std::string());
            pathfile.str(std::string());
		}
	}
}

    
} //namespace smart
