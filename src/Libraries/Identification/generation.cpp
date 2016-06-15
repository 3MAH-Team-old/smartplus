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

///@file generation.cpp
///@brief generation for genetic algorithm (among others)
///@author Chemisky & Despringre
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <smartplus/Libraries/Identification/individual.hpp>
#include <smartplus/Libraries/Identification/generation.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for phases_characteristics===================================

//=====Public methods for phases_characteristics============================================

///@brief default constructor
//----------------------------------------------------------------------
generation::generation()
//----------------------------------------------------------------------
{

}

///@brief Constructor
///@param n : number of individuals
///@param init boolean that indicates if the constructor has to initialize  (default value is true)
//----------------------------------------------------------------------
generation::generation(const int &n, const int &m, int &idnumber)
//----------------------------------------------------------------------
{
	assert(n>0);
	
    for (int i=0; i<n; i++) {
        pop.push_back(individual(m, idnumber));
        idnumber++;
    }
}

///@brief Copy constructor
///@param gp generation object to duplicate
//----------------------------------------------------------------------
generation::generation(const generation& gp)
//----------------------------------------------------------------------
{
    pop = gp.pop;
}

///@brief destructor
//----------------------------------------------------------------------
generation::~generation() {}
//----------------------------------------------------------------------

///@brief Construct : A method to construct the generation after its initial construction
//----------------------------------------------------------------------
void generation::construct(const int &n, const int &m, int &idnumber)
//----------------------------------------------------------------------
{
	assert(n>0);

    for (int i=0; i<n; i++) {
        pop.push_back(individual(m, idnumber));
        idnumber++;
    }
}

///@brief classify : A method to classify the individuals in a generation
//----------------------------------------------------------------------
void generation::classify()
//----------------------------------------------------------------------
{
    
	assert(pop.size()>0);

	double mini = -1.;
	int posmini = 0;
	individual temp;
	
	for(unsigned int i=0; i < pop.size()-1; i++) {
		mini=pop[i].cout;
		posmini=i;
		
		for(unsigned int j=i; j < pop.size(); j++) {
			if(pop[j].cout < mini) {
				mini=pop[j].cout;
				posmini=j;
			}
		}
		temp=pop[posmini];
		pop[posmini]=pop[i];
		pop[i]=temp;
	} 	

	for(unsigned int i=0; i < pop.size(); i++) {
		pop[i].rank=i+1;
	}

}

///@brief newid : A method to assign new id's to a generation
//----------------------------------------------------------------------
void generation::newid(int &idnumber)
//----------------------------------------------------------------------
{
    
    assert(pop.size()>0);

    for (unsigned int i=0; i<pop.size(); i++) {
		pop[i].id = idnumber;
        idnumber++;
    }
}

///@brief newid : A method to destruct and free memory of a generation
//----------------------------------------------------------------------
void generation::destruct() {}
//----------------------------------------------------------------------

///@brief Standard operator = for generation
//----------------------------------------------------------------------
generation& generation::operator = (const generation& gp)
//----------------------------------------------------------------------
{
    pop = gp.pop;    

	return *this;
}

///@brief Standard operator << for generation
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const generation& gp)
//--------------------------------------------------------------------------
{
	s << "Display info on the generation\n";
	s << "Number of individuals: " << gp.pop.size() << "\n";
	
	s << "Characteristics of each individual: \n";
    for(auto ind : gp.pop) {
		s << ind;
	}
	s << "\n\n";

	return s;
}

} //namespace smart