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

///@file read.cpp
///@brief To read ODF from Npeak.dat
///@version 1.0

#include <assert.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Material/ODF.hpp>

using namespace std;
using namespace arma;

namespace smart{

void read_peak(ODF &odf_rve, const string &path_data, const string &inputfile) {
    
    unsigned int npeaks = 0;
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream paramodf;
    
    paramodf.open(path_inputfile, ios::in);
    if(paramodf) {
        while (!paramodf.eof())
        {
            getline (paramodf,buffer);
            if (buffer != "") {
                npeaks++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the peak characteristics in the folder: " << path_data << endl;
        return;
    }
    paramodf.close();
    npeaks--;
    
    //Generate the sub_phase vector and har-create the objects pointed buy the shared_ptrs
    odf_rve.construct(npeaks);
    
    int nparams = 0;
    
    paramodf.open(path_inputfile, ios::in);
    paramodf >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(unsigned int i=0; i<npeaks; i++) {
        
        paramodf >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> nparams;
        
        odf_rve.peaks[i].params = zeros(nparams);
        for(unsigned int j=0; j<odf_rve.peaks[i].params.n_elem; j++) {
            paramodf >> buffer;
        }
    }
    paramodf.close();
    
    paramodf.open(path_inputfile, ios::in);
    paramodf >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(unsigned int i=0; i<npeaks; i++) {
        
        paramodf >> odf_rve.peaks[i].number >> odf_rve.peaks[i].method >> odf_rve.peaks[i].mean >> odf_rve.peaks[i].s_dev >> odf_rve.peaks[i].width >> odf_rve.peaks[i].ampl >> buffer;
        
        for(unsigned int j=0; j<odf_rve.peaks[i].params.n_elem; j++) {
            paramodf >> odf_rve.peaks[i].params(j);
        }
        
        if(odf_rve.radian == false) {
            odf_rve.peaks[i].mean *=(pi/180.);
            odf_rve.peaks[i].s_dev *=(pi/180.);
            odf_rve.peaks[i].width *=(pi/180.);
        }
    }
    
    paramodf.close();
}

} //namespace smart
