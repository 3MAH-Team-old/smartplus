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

///@file optimize.cpp
///@brief functions for optimization
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>

#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Maths/lagrange.hpp>

#include <smartplus/Libraries/Identification/parameters.hpp>
#include <smartplus/Libraries/Identification/opti_data.hpp>
#include <smartplus/Libraries/Identification/individual.hpp>
#include <smartplus/Libraries/Identification/optimize.hpp>

using namespace std;
using namespace arma;

namespace smart{

///This function constructs the vector of exp/num
vec calcV(vector<opti_data> &exp_data, const int &nfiles, const int &sizev) {
    vec v = zeros(sizev);
    int z=0;
    
	for(int i=0; i<nfiles; i++) {
		for(int l=0; l<exp_data[i].ninfo; l++) {
			for(int k=0; k<exp_data[i].ndata; k++) {
				v(z) = exp_data[i].data(k,l);
				z++;
			}
		}
	}
    return v;
}

///This function constructs the sensitivity matrix
void calcS(const individual &ind, mat &S, const vec &vnum0, const vector<opti_data> &num_data, const int& nfiles, const int& j, const double &delta) {
	double vnum_pertu=0.;
	
	int z=0;
	for(int i=0; i<nfiles;i++) {
		for(int l=0; l<num_data[i].ninfo;l++) {
			for(int k=0; k<num_data[i].ndata;k++) {
				vnum_pertu = num_data[i].data(k,l);
				S(z,j) = (vnum_pertu-vnum0(z))*(1./(ind.p(j)*delta));
				z++;
			}
		}
	}
}

///This function checks the sensitivity matrix.
///This ensures that if a parameter didn't modify at all the result, the sensibility matrix doesn't have a column of "0" (inversion) issues
Col<int> checkS(mat &S) {
	double somme;
	int problem = 0;
	Col<int> pb_col;
	pb_col.zeros(S.n_cols + 1);
	
	for (int j = 0; j<fabs(S.n_cols); j++) {
		somme = 0.;
		for(int i=0; i<fabs(S.n_rows);i++) {
			somme += fabs(S(i,j));
		}
		if (somme < limit) {
			problem++;
			pb_col(j)++;
			pb_col(S.n_cols)++;
		}
	}
	if (problem > 0) {
		/*if (problem == fabs(S.n_cols)){
			cout << "\nERROR -- PLEASE HANDLE THIS SPECIFIC CASE IN FUNCTION checkS FROM optimize.hpp\n";
		}*/
		mat newS = zeros(S.n_rows,S.n_cols - problem);
		int k = 0;
		
		for (int j = 0; j < (fabs(S.n_cols)); j++) {
			if (pb_col(j) == 0) {
				for(int i=0; i<fabs(S.n_rows);i++) {
					newS(i,k) = S(i,j);
				}
				k++;
			}
		}
		S.set_size(S.n_rows,S.n_cols - problem);
		S = newS;
	}
	return pb_col;
}

///This function computes the Cost function (Square differnces) from the components of experimental values and numerically evaluated values 
/*double calcC(vec vector<opti_data> &data_exp, vector<opti_data> &data_num, const int &nfiles) {
	double Cout = 0.;
	vector<mat> mdiff(nfiles);
	    
	for(int z=0; z<nfiles; z++) {
		assert(data_exp[z].data.n_elem == data_num[z].data.n_elem);
        mdiff[z] = zeros(data_exp[z].ndata, data_exp[z].ninfo);
	}
    
	double denom = 0.;
	double CoutCol = 0.;	
	
	for(int z=0; z<nfiles; z++) {
		mdiff[z] = data_num[z].data - data_exp[z].data;	
		
		for(int j=0; j<data_exp[z].ninfo; j++) {
			for(int i=0; i<data_exp[z].ndata; i++) {
				CoutCol += pow(mdiff[z](i,j), 2.);
				denom += pow(data_exp[z].data(i,j), 2.);				
			}
			Cout += CoutCol/denom;
			CoutCol = 0.;
			denom = 0.;
		}
	}
	return Cout;
}*/
   
double calcC(const vec &vexp, const vec &vnum, const vec &W) {
    double Cout = 0.;
    for(unsigned int z=0; z<vexp.n_elem; z++) {
        if (W(z) > iota) {
            Cout += pow((vexp(z)-vnum(z)), 2.)*W(z);
        }
    }
	return Cout;    
}
    
//Hessian matrix
mat Hessian(const int &size, const mat &S, const vec &W) {
    ///Hessian matrix
    mat H = zeros(size, size);
    for(int i=0; i<(size); i++) {
        for(unsigned int j=0; j<W.n_elem; j++) {
            for(int l=0; l<(size); l++) {
                H(i,l) += S(j,i)*W(j)*S(j,l);
            }
        }
    }
    return H;
}
    
///Gradient matrix
mat diagJtJ(const mat &H){
    return diagmat(H);
}
    
//Minimal bound Lagrange multiplier vector
vec bound_min(const int &size, const vec &p, const vector<parameters> &params, const double &c, const double &p0) {
    vec L_min = zeros(size);
    for(int k=0; k<(size); k++) {
        L_min(k) = -1.*lagrange_exp(params[k].min_value*(1.-p(k)/params[k].min_value), c*params[k].min_value, p0);
    }
    return L_min;
}
    
//Minimal bound Lagrange multiplier vector derivative
vec dbound_min(const int &size, const vec &p, const vector<parameters> &params, const double &c, const double &p0) {
    vec dL_min = zeros(size);
    for(int k=0; k<(size); k++) {
		dL_min(k) = dlagrange_exp(params[k].min_value*(1.-p(k)/params[k].min_value), c*params[k].min_value, p0);
    }
    return dL_min;
}

    
//Maximal bound Lagrange multiplier vector
vec bound_max(const int &size, const vec &p, const vector<parameters> &params, const double &c, const double &p0) {
    vec L_max = zeros(size);
    for(int k=0; k<(size); k++) {
        L_max(k) = -1.*lagrange_exp(p(k)*(p(k)/params[k].max_value-1.), c*p(k), p0);
    }
    return L_max;
}

    
//Maximal bound Lagrange multiplier vector derivative
vec dbound_max(const int &size, const vec &p, const vector<parameters> &params, const double &c, const double &p0) {
    vec dL_max = zeros(size);
    for(int k=0; k<(size); k++) {
		dL_max(k) = -1.*dlagrange_exp(params[k].max_value*(p(k)/params[k].max_value-1.), c*params[k].max_value, p0);
    }
    return dL_max;
}

    
vec calcW(const int &sizev, const int &nfiles, const Col<int> &weight_types, const vec &weight_files, const vector<vec> &weight_cols, const vector<opti_data> &weight, const vector<opti_data> &data_exp) {
    
    vec W = ones(sizev);
    double denom = 0.;
    int z=0;
    
    //Load info for the weight type 1 : Weight for each data file
    //if (weight_types(0) == 0) : Nothing to do
    if(weight_types(0) == 1) {      //Add the weight per file
        for(int i=0; i<nfiles; i++) {
            for(int k=0; k<data_exp[i].ninfo; k++) {
                for(int j=0; j<data_exp[i].ndata; j++) {
                    W(z) *= weight_files(i);
                    z++;
                }
            }
        }
    }

    //Load info for the weight type 2 : Weight for each data columns
    //if (weight_types(1) == 0) : Nothing to do
    z=0;
    if(weight_types(1) == 1) {      //Add the weight per columns corresponding to the sum
        for(int i=0; i<nfiles; i++) {
            for(int k=0; k<data_exp[i].ninfo; k++) {
                for(int j=0; j<data_exp[i].ndata; j++) {
                    denom += pow(data_exp[i].data(j,k),2.);
                }
                for(int j=0; j<data_exp[i].ndata; j++) {
                    W(z) *= (1./denom);
                    z++;
                }
                denom = 0.;
            }
        }
    }
    else if(weight_types(1) == 2) {      //Add the weight per columns corresponding to the sum + weight
        for(int i=0; i<nfiles; i++) {
            for(int k=0; k<data_exp[i].ninfo; k++) {
                for(int j=0; j<data_exp[i].ndata; j++) {
                    denom += pow(data_exp[i].data(j,k),2);
                }
                for(int j=0; j<data_exp[i].ndata; j++) {
                    W(z) *= weight_cols[i](k)/denom;
                    z++;
                }
                denom = 0.;
            }
        }
    }
    else if(weight_types(1) == 3) {      //Add the weight per columns corresponding to weight only
        for(int i=0; i<nfiles; i++) {
            for(int k=0; k<data_exp[i].ninfo; k++) {
                for(int j=0; j<data_exp[i].ndata; j++) {
                    W(z) *= weight_cols[i](k);
                    z++;
                }
            }
        }
    }
        
    //Load info for the weight type 3 : Weight for each data
    //if (weight_types(2) == 0) : Nothing to do
    z=0;
    if(weight_types(2) == 1) {
        for(int i=0; i<nfiles; i++) {
            for(int k=0; k<data_exp[i].ninfo; k++) {
                for(int j=0; j<data_exp[i].ndata; j++) {
                    W(z) *= fabs(weight[i].data(j,k));
                    z++;
                }
            }
        }
    }
    return W;
}
    
vec G_cost(const int &size, const mat &S, const vec &W, const vec &Dv, const vec &L_min, const vec &L_max) {
    vec G = zeros(size);
    for(int i=0; i<(size); i++) {
        for(unsigned int j=0; j<W.n_elem; j++) {
            G(i) += S(j,i)*W(j)*Dv(j);
        }
    }
    
    //Integrate the limits
    for(int k=0; k<(size); k++) {
        G(k) += -1.*fabs(G(k))*(L_min(k) + L_max(k));
    }
    return G;
}

///Levenberg-Marquardt matrix, with bounds
mat LevMarq(const mat &H, const double &lambdaLM, const vec &dL_min, const vec &dL_max){

    int size=H.n_cols;
    mat LM = H + lambdaLM*diagmat(H);

    for(int k=0; k<(size); k++) {
        LM(k,k) += fabs(LM(k,k))*(dL_min(k) + dL_max(k));
    }
    return LM;
}
    
    
vec calcDp(const mat &S, const vec &vexp, const vec &vnum, const vec &W, const vec &p, const vector<parameters> &params, const double &lambdaLM, const double &c, const double &p0, const int &nprops, const Col<int>& pb_col) {
    
    vec FullDp = zeros(nprops);
    int problem = pb_col(nprops);
    int sizepb = nprops-problem;
    vec Dp = zeros(nprops-problem);
    vec Dv = (vexp-vnum);
    mat H = Hessian(sizepb, S, W);
    
    ///Constrain optimization
    vec L_min = bound_min(sizepb, p, params, c, p0);
    vec L_max = dbound_max(sizepb, p, params, c, p0);
    vec dL_min = bound_min(sizepb, p, params, c, p0);
    vec dL_max = dbound_max(sizepb, p, params, c, p0);

    vec G = G_cost(sizepb, S, W, Dv, L_min, L_max);
    mat LM = LevMarq(H, lambdaLM, dL_min, dL_max);
    
    Dp = inv(LM)*G;
    
    int z = 0;
    for(int i=0; i<nprops; i++) {
        if (pb_col(i)==0) {
            FullDp(i) = Dp(z);
            z++;
        }
        else {
            FullDp(i) = 0.;
        }
    }
    
    ///Safety to avoid any big change in the values of p(i) due to over-sensibility
    for(int i=0; i<nprops; i++) {
        if(fabs(FullDp(i)) > 0.1*p(i)) {
            if(FullDp(i) > 0) {
                FullDp(i) = 0.1*fabs(p(i));
            }
            else  {
                FullDp(i) = -1.*0.1*fabs(p(i));
            }
        }
    }
    return FullDp;
}
    
} //namespace smart