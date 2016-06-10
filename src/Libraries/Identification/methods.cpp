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

///@file methods.cpp
///@brief methods for genetic algorithm (among others)
///@author Chemisky & Despringre
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>

#include <smartplus/Libraries/Maths/random.hpp>
#include <smartplus/Libraries/Maths/lagrange.hpp>

#include <smartplus/Libraries/Identification/parameters.hpp>
#include <smartplus/Libraries/Identification/methods.hpp>
#include <smartplus/Libraries/Identification/generation.hpp>
#include <smartplus/Libraries/Identification/optimize.hpp>

using namespace std;
using namespace arma;

namespace smart{
    
//Genetic method
void genetic(generation &gen_g, generation &gensons, int &idnumber, const double &probaMut, const double &pertu, const vector<parameters> &params){
    
    int n_param = params.size();
    int maxpop = gensons.nindividuals;
    
    //Generate two genitoers
	individual dad(n_param, 0);
	individual mom(n_param, 0);

    int chromosome = 0;
    //Very small pertubation
    
    gensons.newid(idnumber);
    for(int i=0; i<gensons.nindividuals; i++) {
        /// Random determination of "father" and "mother"
        dad = gen_g.pop[alea(maxpop-1)];
        mom = gen_g.pop[alea(maxpop-1)];
        while(dad.id==mom.id)
            mom = gen_g.pop[alea(maxpop-1)];
            
        for(int j=0; j<n_param; j++) {
            chromosome = alea(1);
            if(chromosome==0) {
                gensons.pop[i].p(j)=dad.p(j)*alead(1.-pertu,1.+pertu);
            }
            else
                gensons.pop[i].p(j)=mom.p(j)*alead(1.-pertu,1.+pertu);
                
                if (gensons.pop[i].p(j) > params[j].max_value)
                    gensons.pop[i].p(j) = params[j].max_value;
                    if (gensons.pop[i].p(j) < params[j].min_value)
                        gensons.pop[i].p(j) = params[j].min_value;
                        
                        ///Apply a mutation
                        if (alea(99)<probaMut)
                            gensons.pop[i].p(j) = alead(params[j].min_value, params[j].max_value);
        }
    }
    
}


///Genrun creation
void to_run(generation &gensons, generation &gboys, generation &genrun, const double &delta, const vector<parameters> &params) {

    int z = 0;
    int maxpop = gensons.nindividuals;
    int n_param = params.size();
    
    //genrun part of the genetic
    if (maxpop > 1) {
        for(int i=0; i<gensons.nindividuals; i++) {
            genrun.pop[z] = gensons.pop[i];
            z++;
        }
    }

    //genrun part of the gradient
    for(int k=0; k<gboys.nindividuals ;k++) {
        genrun.pop[z].p = gboys.pop[k].p;
        z++;
        for(int j=0; j<n_param; j++) {
            ///set up the new set of parameters
            genrun.pop[z].p = gboys.pop[k].p;
            genrun.pop[z].p(j) = gboys.pop[k].p(j)*(1. + delta);
            z++;
        }
    }

}

} //namespace smart
