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
///@brief To read from material.dat and path.dat
///@version 1.0

#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Solver/block.hpp>
#include <smartplus/Libraries/Solver/step.hpp>
#include <smartplus/Libraries/Solver/step_meca.hpp>
#include <smartplus/Libraries/Solver/step_thermomeca.hpp>
#include <smartplus/Libraries/Solver/output.hpp>

using namespace std;
using namespace arma;

namespace smart{

Col<int> subdiag2vec() {
    
    Col<int> temp;
    temp = zeros<Col<int> >(6);
    temp(0) = 0;
	temp(1) = 3;
	temp(2) = 1;
	temp(3) = 4;
	temp(4) = 5;
	temp(5) = 2;
    
    return temp;
}

/// Function that fills the matrix Tdsde for mix strain/stress conditions
void Lt_2_K(const mat &Lt, mat &K, const Col<int> &cBC_meca, const double &lambda)
{
	K = zeros(6,6);
    
    for (int i=0; i<6; i++) {
        if (cBC_meca(i)) {
            K.row(i) = Lt.row(i);
        }
        else
            K(i,i) = lambda;
    }
}

/// Function that fills the matrix Tdsde for mix strain/stress conditions
void Lth_2_K(const mat &dSdE, mat &dSdT, mat &dQdE, mat &dQdT, mat &K, const Col<int> &cBC_meca, const int &cBC_T, const double &lambda)
{
	K = zeros(7,7);
    
    K.submat(0, 0, 5, 5) = dSdE;
    K.submat(0, 6, 5, 6) = dSdT;
    K.submat(6, 0, 6, 5) = dQdE;
    K.submat(6, 6, 6, 6) = dQdT;
    
    for (int i=0; i<6; i++) {
        if (cBC_meca(i) == 0) {
            K.row(i) = 0.*K.row(i);
            K(i,i) = lambda;
        }
    }
    if (cBC_T == 0) {
            K.row(6) = 0.*K.row(6);
            K(6,6) = lambda;
    }
}

/*void Lt_2_K(const mat &Lt, mat &K, const Col<int> &cBC_meca)
{
	int kT;
	int lT = 0;
	for (int j = 0; j < 6; j++){
		if (cBC_meca(j)){
			kT = 0;
			for (int i = 0; i < 6; i++){
				if (i == j){
					K(kT,lT) = Lt(i,j);
					kT++;
				}
				else if (cBC_meca(i))
				{
					K(kT,lT) = Lt(i,j);
					kT++;
				}
			}
			lT++;
		}
	}
}*/

void read_matprops(string &umat_name, int &nprops, vec &props, int &nstatev, vec &statev, double &rho, double &c_p) {

    ///Material properties reading, use "material.dat" to specify parameters values
	string buffer;
	ifstream propsmat;
	propsmat.open("data/material.dat", ios::in);
	if(propsmat) {
        
		string buffer;
		propsmat >> buffer >> buffer >> umat_name >> buffer >> nprops >> buffer >> nstatev;
	}
	else {
		cout << "Error: cannot open material.dat file \n";
	}
	
	char *cmname = new char [umat_name.length()];
	strcpy (cmname, umat_name.c_str());
    
	propsmat.close();
    	
	props = zeros(nprops);
	statev = zeros(nstatev);
    
	propsmat.open("data/material.dat", ios::in);
	if(propsmat) {
		string buffer;
		propsmat >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> rho >> buffer >> c_p >> buffer;
        
		for(int i=0;i<nprops;i++)
			propsmat >> buffer >> props(i);
	}
	else {
		cout << "Error: cannot open the file material.dat \n";
	}
    
	propsmat.close();
    
}

/*void read_path_init(int &type, int &nblock) {

    string buffer;
    ifstream path;
    path.open("path.txt", ios::in);
    if(!path) {
        cout << "Error: cannot open the file path.txt \n";
	}
    
	///temperature is initialized
	path >> buffer >> type >> buffer >> buffer >> buffer >> nblock;
    path.close();
}*/

void read_output(solver_output &so, const int &nblock, const int &nstatev) {
        
    string buffer;
    string file;
    
    ifstream cyclic_output;
	cyclic_output.open("data/output.dat", ios::in);
	if(!cyclic_output)
	{
		cout << "Error: cannot open the file data/output.dat \n";
	}
    
	cyclic_output >> buffer;
    cyclic_output >> buffer >> so.o_nb_meca;
    so.o_meca.zeros(so.o_nb_meca);
    for (int i=0; i<so.o_nb_meca; i++) {
        cyclic_output >> so.o_meca(i);
    }
    cyclic_output >> buffer >> so.o_nb_T;
    
    ///Selection of the wanted umat statev, use "cyclic.dat" to specify wanted internal variables
    cyclic_output >> buffer >> buffer;
    if ((buffer == "all") || (buffer == "All") || (buffer == "ALL")){
        so.o_wanted_statev.zeros(1);
        so.o_wanted_statev(0) = -1;
    }
    else if(atoi(buffer.c_str()) != 0){
        so.o_nw_statev = atoi(buffer.c_str());
        so.o_wanted_statev.zeros(so.o_nw_statev);
        so.o_range_statev.zeros(so.o_nw_statev);
        for (int i = 0; i < so.o_nw_statev; i++){
            cyclic_output >> buffer >> buffer;
            if ((buffer == "from") || (buffer == "From") || (buffer == "FROM")){
                cyclic_output >> so.o_wanted_statev(i) >> buffer >> so.o_range_statev(i);
            }
            else{
                so.o_wanted_statev(i) = atoi(buffer.c_str());
                so.o_range_statev(i) = so.o_wanted_statev(i);
            }
            
            if(so.o_range_statev(i) > nstatev -1) {
                cout << "Error : The range of outputed statev is greater than the actual number of statev!\n";
                cout << "Check output.dat and/or material.dat\n\n";
                
                exit(0);
            }
        }
    }
    else {
        so.o_nw_statev = 0;
    }
    
    cyclic_output >> buffer >> buffer >> buffer;
    for(int i = 0 ; i < nblock ; i++){
        cyclic_output >> buffer >> so.o_type(i);
        if(so.o_type(i) == 1)
            cyclic_output >> so.o_nfreq(i);
        else if(so.o_type(i) == 2)
            cyclic_output >> so.o_tfreq(i);
        else
            cyclic_output >> buffer;
    }

	cyclic_output.close();
    
}


void read_path(std::vector<block> &blocks, double &T, const string &pathfile) {
    
	/// Reading the loading path file, Path.txt
    string buffer;
    int conver;
    char bufferchar;
    int nblock;
    Col<int> Equiv = subdiag2vec();
    
	ifstream path;
	path.open(pathfile, ios::in);
	if(!path)
	{
		cout << "Error: cannot open the file path.txt \n";
	}

	///temperature is initialized
	path >> buffer >> T >> buffer >> nblock;
    blocks.resize(nblock);
    
    /// Reading path_file
    for(int i = 0 ; i < nblock ; i++){
        
        path >> buffer >> blocks[i].number >> buffer >> blocks[i].type >> buffer >> blocks[i].ncycle >> buffer >> blocks[i].nstep;

        if (blocks[i].number != i+1) {
            cout << "could not find the UMAT: This is due to a polymorphism exception mode 8875 (Vérifier Numero de block...)";
        }
        
        blocks[i].generate();
        
        switch(blocks[i].type) {
            case 1: {

                for(int j = 0; j < blocks[i].nstep; j++){
                    
                    path >> buffer >> blocks[i].steps[j]->mode;
                    blocks[i].steps[j]->number = j+1;
                    
                    if ((blocks[i].steps[j]->mode == 1)||(blocks[i].steps[j]->mode == 2)) {
                        
                        shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[j]);
                        
                        path >> buffer >> sptr_meca->Dn_init >> buffer >> sptr_meca->Dn_mini >> buffer >> sptr_meca->Dn_maxi >> buffer >> sptr_meca->BC_Time >> buffer;
                        for(int k = 0 ; k < 6 ; k++) {
                            path >> bufferchar;
                            conver = bufferchar;
                            if (conver == 83){
                                sptr_meca->cBC_meca(Equiv(k)) = 1;
                                path >> sptr_meca->BC_meca(Equiv(k));
                            }
                            else if (conver == 69){
                                sptr_meca->cBC_meca(Equiv(k)) = 0;
                                path >> sptr_meca->BC_meca(Equiv(k));
                            }
                        }
                        path >> buffer >> bufferchar;
                        conver = bufferchar;
                        if (conver == 84){
                            path >> sptr_meca->BC_T;
                        }
                        else
                            cout << "Error, This is a mechanical step, only temperature boundary condition is allowed here\n";

                    }
                    else if (blocks[i].steps[j]->mode == 3) {
                        
                        shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[j]);                        
                        
                        path >> buffer >> sptr_meca->file >> buffer;
                        
                        for(int k = 0 ; k < 6 ; k++) {
                            path >> bufferchar;
                            conver = bufferchar;
                            if (conver == 83){
                                sptr_meca->cBC_meca(Equiv(k)) = 1;
                            }
                            else if (conver == 69) {
                                sptr_meca->cBC_meca(Equiv(k)) = 0;
                            }
                            else if (conver == 48) {
                                sptr_meca->cBC_meca(Equiv(k)) = 2;      // this is a special stress-controlled one (it is different since it does not read data from the path_inc tabular file)
                            }
                                
                        }
                        path >> buffer >> bufferchar;
                        conver = bufferchar;
                        if (conver == 84){
                            sptr_meca->cBC_T = 0;                       //This is the classical temperature imposed in the file
                        }
                        else if (conver == 48) {                        //This is a special case where the temperature is constant
                            sptr_meca->cBC_T = 2;
                        }
                    }
                    else {
                        cout << "The SMART+ is not compatible with your OS. Please think about upgrading to Mac OS 10.9 or later, it is much better than your shitty Windows 8";
                    }
                }
                break;
            }
            case 2: {
                
                for(int j = 0; j < blocks[i].nstep; j++){
                    
                    path >> buffer >> blocks[i].steps[j]->mode;
                    blocks[i].steps[j]->number = j+1;
                    
                    if ((blocks[i].steps[j]->mode == 1)||(blocks[i].steps[j]->mode == 2)) {
                        
                        shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[j]);
                        
                        path >> buffer >> sptr_thermomeca->ninc >> buffer >> sptr_thermomeca->BC_Time >> buffer;
                        for(int k = 0 ; k < 6 ; k++) {
                            path >> bufferchar;
                            conver = bufferchar;
                            if (conver == 83){
                                sptr_thermomeca->cBC_meca(Equiv(k)) = 1;
                                path >> sptr_thermomeca->BC_meca(Equiv(k));
                            }
                            else if (conver == 69){
                                sptr_thermomeca->cBC_meca(Equiv(k)) = 0;
                                path >> sptr_thermomeca->BC_meca(Equiv(k));
                            }
                        }
                        path >> buffer >> bufferchar;
                        conver = bufferchar;
                        if (conver == 81){
                            sptr_thermomeca->cBC_T = 1;
                            path >> sptr_thermomeca->BC_T;
                        }
                        else if (conver == 84){
                            sptr_thermomeca->cBC_T = 0;
                            path >> sptr_thermomeca->BC_T;
                        }
                        
                    }
                    else if (blocks[i].steps[j]->mode == 3) {
                        
                        shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[j]);
                        
                        path >> buffer >> sptr_thermomeca->file >> buffer;
                        
                        for(int k = 0 ; k < 6 ; k++) {
                            path >> bufferchar;
                            conver = bufferchar;
                            if (conver == 83){
                                sptr_thermomeca->cBC_meca(Equiv(k)) = 1;
                            }
                            else if (conver == 69) {
                                sptr_thermomeca->cBC_meca(Equiv(k)) = 0;
                            }
                            else if (conver == 48) {
                                sptr_thermomeca->cBC_meca(Equiv(k)) = 2;      // this is a special stress-controlled one (it is different since it does not read data from the path_inc tabular file)
                            }
                            
                        }
                        path >> buffer >> bufferchar;
                        conver = bufferchar;
                        if (conver == 84){
                            sptr_thermomeca->cBC_T = 0;                       //This is the classical temperature imposed in the file
                        }
                        else if (conver == 81){
                            sptr_thermomeca->cBC_T = 1;                       //This is the classical heat flux quantitt imposed in the file
                        }
                        else if (conver == 48) {                        //This is a special case where the temperature is constant
                            sptr_thermomeca->cBC_T = 2;
                        }
                        
                    }
                    else {
                        cout << "The SMART+ is not compatible with your OS. Please think about upgrading to Mac OS 10.9 or later, it is much better than your shitty Windows 8";
                    }
                    
                }
                break;
            }
            default: {
                cout << "Seriously?\n";
                break;
            }
            
        }
    }
    path.close();
    
}

} //namespace smart
