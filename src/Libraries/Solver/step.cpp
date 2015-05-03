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

///@file step.cpp
///@brief object that defines an step
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <smartplus/Libraries/Solver/step.hpp>
#include <smartplus/Libraries/Solver/output.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for ellipsoid_characteristics===================================

//=====Public methods for ellipsoid_characteristics============================================

//@brief default constructor
//-------------------------------------------------------------
step::step()
//-------------------------------------------------------------
{
	number=0;
	ninc=0;
	mode=0;
    
    Time = 0.;
    BC_Time = 0.;
    
    file = "";
    
}

/*!
 \brief Constructor with parameters
 \param mnumber : number of the step
 \param mninc : number of increments
 \param mmode : mode : type of loading
 */

//-------------------------------------------------------------
step::step(int mnumber, int mninc, int mmode)
//-------------------------------------------------------------
{
	assert(mnumber>=0);
	assert(mninc>0);
	assert(mmode>0);
    
	number = mnumber;
	ninc = mninc;
	mode = mmode;
    
    Time = 0.;
    times = zeros(ninc);
    BC_Time = 0.;
    
    file = "";
}

/*!
 \brief Copy constructor
 \param st step object to duplicate
 */

//------------------------------------------------------
step::step(const step& st)
//------------------------------------------------------
{    
	number = st.number;
	ninc = st.ninc;
	mode = st.mode;
    
    Time = st.Time;
    times = st.times;
    BC_Time = st.BC_Time;
    
    file = st.file;
    
}

/*!
 \brief destructor
 */

step::~step() {}

//-------------------------------------------------------------
void step::generate(const double &mTime, const vec &msigma, const vec &mEtot, const double &mT)
//-------------------------------------------------------------
{
	assert(number>=0);
	assert(ninc>0);
	assert(mode>0);

    times = zeros(ninc);
}

/*!
 \brief Standard operator = for block
 */

//----------------------------------------------------------------------
step& step::operator = (const step& st)
//----------------------------------------------------------------------
{
//	assert(st.ninc>0);
//	assert(st.mode>0);
    
	number = st.number;
	ninc = st.ninc;
	mode = st.mode;
        
    Time = st.Time;
    times = st.times;
    BC_Time = st.BC_Time;
    
    file = st.file;
        
	return *this;
}

void step::output(ostream& output, const solver_output &so, const int &kblock, const int&kcycle, const int &kinc, const vec &statev) {
    
    output << kblock << "\t";
    output << kcycle << "\t";
    output << number << "\t";
    output << kinc << "\t";
    output << times(kinc) << "\t\t";
    
    if (so.o_nb_T) {
        output << 0  << "\t";
        output << 0 << "\t";                //This is for the flux
    }
    if (so.o_nb_meca) {
        for (int z=0; z<so.o_nb_meca; z++) {
            output << 0 << "\t";
        }
        for (int z=0; z<so.o_nb_meca; z++) {
            output << 0 << "\t";
        }
    }
    output << "\t";
    if(so.o_nw_statev != 0){
        if (so.o_wanted_statev(0) < 0) {
            for(int k = 0 ; k < statev.n_elem ; k++)
                output << statev(k) << "\t";
        }
        else{
            for(int k = 0 ; k < so.o_nw_statev ; k++){
                for (int l = so.o_wanted_statev(k); l < (so.o_range_statev(k)+1); l++){
                    output << statev(l) << "\t";
                }
            }
        }
    }
    output << endl;
    
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const step& st)
//--------------------------------------------------------------------------
{    
	s << "Display info on the step " << st.number << "\n";
	s << "Number of increments: " << st.ninc << "\twithin " << st.BC_Time << " s\n";
	s << "Loading mode: " << st.mode << "\n";
	   
	return s;
}

} //namespace smart
