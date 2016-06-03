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

///@file identification.cpp
///@brief The main identification function
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <algorithm>

#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Maths/random.hpp>

#include <smartplus/Libraries/Identification/parameters.hpp>
#include <smartplus/Libraries/Identification/constants.hpp>
#include <smartplus/Libraries/Identification/optimize.hpp>
#include <smartplus/Libraries/Identification/generation.hpp>
#include <smartplus/Libraries/Identification/methods.hpp>
#include <smartplus/Libraries/Identification/opti_data.hpp>
#include <smartplus/Libraries/Identification/doe.hpp>
#include <smartplus/Libraries/Identification/read.hpp>
#include <smartplus/Libraries/Identification/script.hpp>

using namespace std;
using namespace arma;

namespace smart{
    
void run_identification_solver(const int &n_param, const int &n_consts, const int &nfiles, const int &ngen, const int &aleaspace, int &apop, int &spop, const int &ngboys, const int &maxpop, const std::string &path_data, const std::string &path_keys, const std::string &outputfile, const std::string &data_num_name, const double &probaMut, const double &pertu, const double &c, const double &p0, const double &lambdaLM) {

    std::string data_num_ext = data_num_name.substr(data_num_name.length()-4,data_num_name.length());
    std::string data_num_name_root = data_num_name.substr(0,data_num_name.length()-4); //to remove the extension
    
    ///Allow non-repetitive pseudo-random number generation
    srand(time(0));
    int TOOL = 1;   ///Which code is going to compute numerical files
    ofstream result;    ///Output stream, with parameters values and cost function

    //Define the parameters
    vector<parameters> params(n_param);  //vector of parameters
    vector<constants> consts(n_consts);  //vector of parameters
    vec Dp = zeros(n_param);
    vec p = zeros(n_param);
    //Read the parameters

    read_parameters(n_param, params);
    read_constants(n_consts, consts, nfiles);
    
    //string path_data = "data/";
    //string path_keys = "data/key/";
    //Copy of the parameters to the keys folder

    //copy_parameters(params, path_data, path_keys);
    //copy_constants(consts, path_data, path_keys);

    int idnumber=1;
    int id0=0;

    //Define the necessary vectors to compute the cost functions and constrain optimization
    vec GBcost = zeros(ngboys);
    vec lambdaLMV = zeros(ngboys);
    for(int i=0; i<ngboys; i++) {
        lambdaLMV(i) = lambdaLM;
    }

    //Get the data structures
    vector<opti_data> data_exp(nfiles);
    vector<opti_data> data_weight(nfiles);
    vector<opti_data> data_num(nfiles);

    Col<int> weight_types(3);
    vec weight_files = zeros(nfiles);
    vector<vec> weight_cols(nfiles);

    read_data_exp(nfiles, data_exp);
    read_data_weights(nfiles, weight_types, weight_files, weight_cols, data_weight, data_exp);
    read_data_num(nfiles, data_exp, data_num);

    /// Get the data vectors
    ///Import of the experimental data
    string data_exp_folder="exp_data";
    int sizev = 0;
    for(int i=0; i<nfiles;i++) {
        data_exp[i].import(data_exp_folder);
        data_weight[i].import(data_exp_folder);
        sizev += data_exp[i].ndata * data_exp[i].ninfo;
    }

    //Data structure has been created. Next is the generation of structures to compute cost function and associated derivatives
    vec vexp = calcV(data_exp, nfiles, sizev);
    
    vec W = calcW(sizev, nfiles, weight_types, weight_files, weight_cols, data_weight, data_exp);
    vec vnum = zeros(sizev);   //weight vector
    mat S(sizev,n_param);
    Col<int> pb_col;
    pb_col.zeros(n_param);

    result.open(outputfile,  ios::out);
    result << "g" << "\t";
    result << "nindividual" << "\t";
    result << "cost" << "\t";
    for(int i=0; i<n_param; i++) {
        result << "p(" << i << ")" << "\t";
    }
    result << "\n";
    result.close();

    ///Optimization process
    ///Creation of the generation table
    vector<generation> gen(ngen+1);
    vector<generation> gboys(ngen+1);

    generation geninit;
    int g=0;

    gen[g].nindividuals = maxpop;
    gen[g].construct(n_param, id0);

    if(ngboys) {
        gboys[g].nindividuals = ngboys;
        gboys[g].construct(n_param, id0);
    }

    gen_initialize(geninit, spop, apop, idnumber, aleaspace, n_param, params);

    string data_num_folder = "num_data";

    /// Run the simulations corresponding to each individual
    /// The simulation input files should be ready!
    /// Tool defines the numerical algorithm that will compute num files.
    switch (TOOL) {
        case 1: { ///Use of the solver
            
            launch_solver(geninit, nfiles, params, consts, data_num_folder, data_num_name, path_data, path_keys);
            break;
        }
            /*		case 2: {///Use of the specific ODF code
             launch_ODF(geninit, nfiles, p_umat, data_num_folder, data_num_name, data_num_ext, Angle, Nb_ligne);
             break;
             }	*/
        default: {
            cout << "\n\nCRITICAL ERROR : The specified optimization Tool (" << TOOL << ") does not exist.\n";
            exit(0);
        }
    }
    
    ///determination of the numerical data to use
    for(int j=0; j<geninit.nindividuals; j++) {
        for(int i=0; i<nfiles; i++) {

            data_num[i].name = data_num_name_root + to_string(i+1) + "_" + to_string(j+1) + "_global-0" + data_num_ext;
            data_num[i].import(data_num_folder);
            assert(data_exp[i].ndata==data_num[i].ndata);
        }
        
        ///Computation of the cost function
        vnum = calcV(data_num, nfiles, sizev);
        geninit.pop[j].cout = calcC(vexp, vnum, W);
    }

    //	lambdaLM = gen[g].pop[0].cout*lambdaLM;
    for(int i=0; i<maxpop; i++) {
        gen[0].pop[i]=geninit.pop[i];
    }
    gen[0].classify();

    result.open(outputfile,  ios::out | ios::app);
    for(int i=0; i<maxpop; i++) {
        result << 0 << "\t" << gen[0].pop[i].id << "\t" << gen[0].pop[i].cout << "\t";
        for(int j=0; j<n_param;j++) {
            result << gen[0].pop[i].p(j);
            if(j==n_param-1)
                result << "\n";
            else
                result << "\t";
        }
    }
    result.close();

    cout << "Cout[Meilleur Individu] = " << gen[0].pop[0].cout << "\n";

    ///Next step : Classify the best one, and compute the next generation!
    // Here we can choose the type of optimizer we want
    ///Creation of the generation table

    generation gensons(maxpop, n_param, id0);
    generation genall((maxpop > 1) ?  2*maxpop : maxpop, n_param, id0);
    generation genrun((maxpop > 1) ?  (maxpop + ngboys*(n_param+1)) : (ngboys*(n_param+1)),n_param, id0);
    generation gboys_n(ngboys, n_param, id0);

    vec temp_p;
    double temp_cout;

    for(int i=0; i<ngboys; i++) {
        gboys[0].pop[i] = gen[0].pop[i];
    }

    double delta = 0.01;
    bool mauvaise_descendance;
    
    double costnm1 = 0.;
    double stationnarity = 1.E-12;		/// Stationnary stopping criteria (no more evolution of the cost function)
        
    while((g<ngen)&&(fabs(gen[g].pop[0].cout - costnm1) > stationnarity)) {
        cout << "Generation " << g+1 << " sur " << ngen << "...\t";
        
        costnm1 = gen[g].pop[0].cout;
        
        if (maxpop > 1) {
            genetic(gen[g], gensons, idnumber, probaMut, pertu, params);
        }
        
        ///prepare the individuals to run
        to_run(gensons, gboys[g], genrun, delta, params);
        
        ///Run the simulations for the new individuals : Sons generation & gboys + delta
        switch (TOOL) {
            case 1: {
                launch_solver(genrun, nfiles, params, consts, data_num_folder, data_num_name, path_data, path_keys);
                break;
            }
                /*			case 2: {
                 launch_ODF(genrun, nfiles, p_umat, data_num_folder, data_num_name, data_num_ext, Angle, Nb_ligne);
                 break;
                 }		*/
            default: {
                exit(0);
            }
        }
        
        ///Import result for the sons generation and compute their cost
        int z = 0;
        if (maxpop > 1) {
            for(int j=0; j<gensons.nindividuals; j++) {
                for(int i=0; i<nfiles; i++) {
                    data_num[i].name = data_num_name_root + to_string(i+1) + "_" + to_string(z+1) + "_global-0" + data_num_ext;
                    data_num[i].import(data_num_folder);
                    assert(data_exp[i].ndata == data_num[i].ndata);
                }
                ///Computation of the cost function
                vnum = calcV(data_num, nfiles, sizev);
                gensons.pop[j].cout = calcC(vexp, vnum, W);
                z++;
            }
        }
        
        ///Import result for the gboys + delta and compute their cost
        ///temp gboys generation
        gboys_n = gboys[g];
        ///set up the S and V for the gradient based
        for(int k=0; k<ngboys; k++) {
            ///Get the result for the unmodified Gboys
            for(int i=0; i<nfiles; i++) {
                data_num[i].name = data_num_name_root + to_string(i+1) + "_" + to_string(z+1) + "_global-0" + data_num_ext;
                data_num[i].import(data_num_folder);
                assert(data_exp[i].ndata == data_num[i].ndata);
            }
            vnum = calcV(data_num, nfiles, sizev);
            gboys[g].pop[k].cout = calcC(vexp, vnum, W);
            z++;
            
            ///Get the result for +delta specimens
            S.set_size(sizev,n_param);
            for(int j=0; j<n_param; j++) {
                for(int i=0; i<nfiles; i++) {
                    data_num[i].name = data_num_name_root + to_string(i+1) + "_" + to_string(z+1) + "_global-0" + data_num_ext;
                    data_num[i].import(data_num_folder);
                    assert(data_exp[i].ndata == data_num[i].ndata);
                }
                calcS(gboys[g].pop[k], S, vnum, data_num, nfiles, j, delta);
                z++;
            }
            
            ///Check the consistency of the sensitivity matrix
            pb_col.zeros(n_param);
            pb_col = checkS(S);
            
            p = gboys[g].pop[k].p;
            
            ///Compute the parameters increment
            Dp = calcDp(S, vexp, vnum, W, p, params, lambdaLMV(k), c, p0, n_param, pb_col);
            
            p += Dp;
            
            for(int j=0; j < n_param; j++) {
                if(p(j) > params[j].max_value)
                    p(j) = params[j].max_value;
                if(p(j) < params[j].min_value)
                    p(j) = params[j].min_value;
            }
            
            gboys_n.pop[k].p = p;
        }
        
        switch (TOOL) {
            case 1: {
                launch_solver(gboys_n, nfiles, params, consts, data_num_folder, data_num_name, path_data, path_keys);
                break;
            }
                /*			case 2: {
                 launch_ODF(gboys_n, nfiles, p_umat, data_num_folder, data_num_name, data_num_ext, Angle, Nb_ligne);
                 break;
                 }		*/
            default: {
                exit(0);
            }
        }
        
        z=0;
        mauvaise_descendance = true;
        temp_p = gen[g].pop[0].p;
        temp_cout = gen[g].pop[0].cout;
        
        for(int k=0; k<ngboys; k++) {
            for(int i=0; i<nfiles; i++) {
                data_num[i].name = data_num_name_root + to_string(i+1) + "_" + to_string(z+1) + "_global-0" + data_num_ext;
                data_num[i].import(data_num_folder);
                assert(data_exp[i].ndata == data_num[i].ndata);
            }
            vnum = calcV(data_num, nfiles, sizev);
            
            gboys_n.pop[k].cout = calcC(vexp, vnum, W);
            z++;
            if (gboys_n.pop[k].cout < gen[g].pop[0].cout) {
                mauvaise_descendance = false;
            }
            gboys[g].pop[k].p = gboys_n.pop[k].p;
            gboys[g].pop[k].cout = gboys_n.pop[k].cout;
        }
        if (mauvaise_descendance) {
            gboys[g].pop[0].p = temp_p;
            gboys[g].pop[0].cout = temp_cout;
            delta = delta/2.;
        }
        else {
            delta = 0.01;
        }
        if (delta < 0.002) {
            delta = alead(0.03, 0.005);
        }
        
        ///Find the bests
        for(int i=0; i<ngboys; i++) {
            genall.pop[i] = gboys[g].pop[i];
        }
        for(int i=ngboys; i<maxpop; i++) {
            genall.pop[i]=gen[g].pop[i];
        }
        
        if (maxpop > 1) {
            for(int i=0; i<gensons.nindividuals; i++) {
                genall.pop[i+maxpop]=gensons.pop[i];	  
            }			
            genall.classify();
        }
        
        g++;
        gen[g].nindividuals=maxpop;
        gen[g].construct(n_param, id0);		  
        gboys[g].nindividuals = ngboys;
        
        if(ngboys) {
            gboys[g].construct(n_param, id0);
        }
        
        for(int i=0; i<maxpop; i++) {
            gen[g].pop[i] = genall.pop[i];
        }
        
        for(int i=0; i<ngboys; i++) {
            gboys[g].pop[i] = genall.pop[i];
        }
        
        result.open(outputfile,  ios::out | ios::app);	
        for(int i=0; i<maxpop; i++) {
            
            result << g << "\t" << gen[g].pop[i].id << "\t" << gen[g].pop[i].cout << "\t";
            for(int j=0; j<n_param;j++) {		
                result << gen[g].pop[i].p(j);
                if(j==n_param-1)
                    result << "\n";
                else
                    result << "\t";
            }
        }
        result.close();
        cout << "Cout[Meilleur Individu] = " << gen[g].pop[0].cout << "\n";		
    }

//    copy_parameters(params, path_keys, path_data);
//    copy_constants(consts, path_keys, path_data);

    gensons.destruct();
    genall.destruct();
    genrun.destruct();
    gboys_n.destruct();
    
}

} //namespace smart
