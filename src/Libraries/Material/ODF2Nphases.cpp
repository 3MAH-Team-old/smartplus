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

///@file ODF2Nphases.cpp
///@brief ODF2Nphases discretization of ODFs
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Material/ODF2Nphases.hpp>
#include <smartplus/Libraries/Homogenization/ellipsoid_characteristics.hpp>
#include <smartplus/Libraries/Maths/random.hpp>
#include <smartplus/Libraries/Maths/stats.hpp>

using namespace std;
using namespace arma;

namespace smart{

double ODF(const double& theta, const int& method, const vec& param, const bool& radian, const double& dec){
	
	double Theta = theta + dec;
  
	switch (method) {
		case 1: {
			double Mean;
			if (radian) {
				Mean = param(0);
			}
			else {
				Mean = param(0)*pi/180.;
			}
			return ODF_sd(Theta, Mean, param(1), param(2), param(3), param(4));
		}
		case 2: {
			double Mean;
			double Sd;
			if (radian) {
				Mean = param(0);
				Sd = param(2);
			}
			else {
				Mean = param(0)*pi/180.;
				Sd = param(2)*pi/180.;
			}
			return ODF_hard(Theta, Mean, param(1), Sd) + ODF_hard(Theta - pi, Mean, param(1), Sd) + ODF_hard(Theta + pi, Mean, param(1), Sd);
		}
		case 3: {
			double Mean;
			double SD;
			double Ampl = param(2);
			if (radian) {
				Mean = param(0);
				SD = param(1);
			}
			else {
				Mean = param(0)*pi/180.;
				SD = param(1)*pi/180.;
			}
			if (fabs(Ampl) < precision) {
				Ampl = 1.;
			}
			return Gaussian(Theta, Mean, SD, Ampl) + Gaussian(Theta - pi, Mean, SD, Ampl) + Gaussian(Theta + pi, Mean, SD, Ampl);
		}
		case 4: {
			double Mean;
			double Width;
			double Ampl = param(2);
			if (radian) {
				Mean = param(0);
				Width = param(1);
			}
			else {
				Mean = param(0)*pi/180.;
				Width = param(1)*pi/180.;
			}
			if (fabs(Ampl) < precision) {
				Ampl = 1.;
			}
			return Lorentzian(Theta, Mean, Width, Ampl) + Lorentzian(Theta - pi, Mean, Width, Ampl) + Lorentzian(Theta + pi, Mean, Width, Ampl);
		}
		case 5: {
			double Mean;
			double Width;
			double SD;
			double Ampl = param(4);
			if (radian) {
				Mean = param(1);
				Width = param(2);
				SD = param(3);
			}
			else {
				Mean = param(1)*pi/180.;
				Width = param(2)*pi/180.;
				SD = param(3)*pi/180.;
			}
			if (fabs(Ampl) < precision) {
				Ampl = 1.;
			}
			return PseudoVoigt(Theta, param(0), Mean, Width, SD, Ampl) + PseudoVoigt(Theta - pi, param(0), Mean, Width, SD, Ampl) + PseudoVoigt(Theta + pi, param(0), Mean, Width, SD, Ampl);
		}
		case 6: {
			double Mean;
			double Inv_Width;
			double Max = param(3);
			if (radian) {
				Mean = param(0);
				Inv_Width = param(1);
			}
			else {
				Mean = param(0)*pi/180.;
				Inv_Width = param(1)*180./pi;
			}
			if (fabs(Max) < precision) {
				Max = 1.;
			}
			return Pearson7(Theta, Mean, Inv_Width, param(2), Max) + Pearson7(Theta - pi, Mean, Inv_Width, param(2), Max) + Pearson7(Theta + pi, Mean, Inv_Width, param(2), Max);
		}
		case 7: {
			return 1.;
		}
		case 30: {
			vec Mean = zeros(param(0));
			vec SD = ones(param(0));
			vec Ampl = zeros(param(0));
			for (int i = 0; i < param(0); i++) {
				if (radian) {
					Mean(i) = param(3*i + 1);
					SD(i) = param(3*i + 2);				
				}
				else {
					Mean(i) = param(3*i + 1)*pi/180.;
					SD(i) = param(3*i + 2)*pi/180.;	
				}
				Ampl(i) = param(3*i + 3);	
			}
			return Mult_Gaussian(Theta, param(0), Mean, SD, Ampl) + Mult_Gaussian(Theta - pi, param(0), Mean, SD, Ampl) + Mult_Gaussian(Theta + pi, param(0), Mean, SD, Ampl);
		}
		case 40: {
			vec Mean = zeros(param(0));
			vec Width = ones(param(0));
			vec Ampl = zeros(param(0));
			for (int i = 0; i < param(0); i++) {
				if (radian) {
					Mean(i) = param(3*i + 1);
					Width(i) = param(3*i + 2);	
				}
				else {
					Mean(i) = param(3*i + 1)*pi/180.;
					Width(i) = param(3*i + 2)*pi/180.;
				}
				Ampl(i) = param(3*i + 3);
			}
			return Mult_Lorentzian(Theta, param(0), Mean, Width, Ampl) + Mult_Lorentzian(Theta - pi, param(0), Mean, Width, Ampl) + Mult_Lorentzian(Theta + pi, param(0), Mean, Width, Ampl);
		}
		case 50: {
			vec Eta = zeros(param(0));
			vec Mean = zeros(param(0));
			vec Width_Lor = ones(param(0));
			vec SD_Gau = ones(param(0));
			vec Ampl = zeros(param(0));
			for (int i = 0; i < param(0); i++) {
				Eta(i) = param(5*i + 1);
				if (radian) {		
					Mean(i) = param(5*i + 2);
					Width_Lor(i) = param(5*i + 3);
					SD_Gau(i) = param(5*i + 4);
				}
				else {
					Mean(i) = param(5*i + 2)*pi/180.;
					Width_Lor(i) = param(5*i + 3)*pi/180.;
					SD_Gau(i) = param(5*i + 4)*pi/180.;
				}
				Ampl(i) = param(5*i + 5);
			}
			return Mult_PseudoVoigt(Theta, param(0), Eta, Mean, Width_Lor, SD_Gau, Ampl) + Mult_PseudoVoigt(Theta - pi, param(0), Eta, Mean, Width_Lor, SD_Gau, Ampl) + Mult_PseudoVoigt(Theta + pi, param(0), Eta, Mean, Width_Lor, SD_Gau, Ampl);
		}
		case 60: {
			vec Mean = zeros(param(0));
			vec Inv_Width = zeros(param(0));
			vec Shape = ones(param(0));
			vec Max = zeros(param(0));
			for (int i = 0; i < param(0); i++) {
				if (radian) {		
					Mean(i) = param(4*i + 1);
					Inv_Width(i) = param(4*i + 2);						
				}
				else {
					Mean(i) = param(4*i + 1)*pi/180.;
					Inv_Width(i) = param(4*i + 2)*180./pi;					
				}
				Shape(i) = param(4*i + 3);
				Max(i) = param(4*i + 4);
			}
			return Mult_Pearson7(Theta, param(0), Mean, Inv_Width, Shape, Max) + Mult_Pearson7(Theta - pi, param(0), Mean, Inv_Width, Shape, Max) + Mult_Pearson7(Theta + pi, param(0), Mean, Inv_Width, Shape, Max);
		}
		default: {
			cout << "ERROR : cannot find the function labeled " << method << ".\n";
			return 0.;			
		}
	}
}

void ODF2Nphases(const Col<int> &Nphases, const Col<int> &Angle, const Col<int> &method, const vector<string> &filename, const mat &paramODF, const bool &radian, const double& dec) {
/// WARNING: The paramODF is not passed as a parameter to the function ODF for the moment!
		
	cout << "Initializing...\n";
    ///Determination of the number of phases
    int phases_init = Nphases.n_elem;
    int Number_phase = 0.;
	string chaine1;

	///Creation of each object
	std::vector<ellipsoid_characteristics> rvesvs_init(phases_init);
    std::vector < std::vector<ellipsoid_characteristics> > rvesvs(phases_init);
    for (int i=0; i<phases_init; i++) {
        rvesvs[i].resize(Nphases(i));
    }
    
	cout << "Getting phases.dat file...\n";
	///@brief Properties of the phases reading, use "phases.dat" to specify the parameters of each phase
	ifstream phases;
	phases.open("data/phases.dat", ios::in);
	if(phases) {
		phases >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1;
		for(int i=0; i < phases_init; i++) {
			int nprops, nstatev;
			phases >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> nprops >> nstatev;
			rvesvs_init[i].resize(nprops,nstatev);
			for(int j=0; j<nprops; j++) {
				phases >> chaine1;
			}
		}
	}
	else {
		cout << "Error: cannot read  phases.dat file using ODF2phases.hpp\n";
	}
	phases.close();
    
	phases.open("data/phases.dat", ios::in);
	if(phases) {
		phases >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1 >> chaine1;
		for(int i=0; i < phases_init; i++) {
			phases >> rvesvs_init[i].number >> rvesvs_init[i].coatingof >> rvesvs_init[i].umat_name >> rvesvs_init[i].concentration >> rvesvs_init[i].psi_geom >> rvesvs_init[i].theta_geom >> rvesvs_init[i].phi_geom >> rvesvs_init[i].a1 >> rvesvs_init[i].a2 >> rvesvs_init[i].a3 >> rvesvs_init[i].nprops >> rvesvs_init[i].nstatev;
			for(int j=0; j<rvesvs_init[i].dimprops(); j++) {
				phases >> rvesvs_init[i].props[j];
			}
		}
	}
	phases.close();
    
	double normODF = 0.;
	double dalpha = 0.;
	double alpha = 0.;
	cout << "Coupling with ODF and orientation families...\n\n";
	stringstream path;
	mat full_ODF;
	int limiteinf = 0;
	int limitesup;
	int compteur;
	double corr;
	bool flag = false;
	char answer;
	int conver;
	int ndata = 0;

	for(int i = 0; i < phases_init; i++) {
        normODF = 0.;
        alpha = 0.;
        if(Nphases(i) == 1){
            rvesvs[i][0].resize(rvesvs_init[i].nprops,rvesvs_init[i].nstatev);
            rvesvs[i][0] = rvesvs_init[i];
            rvesvs[i][0].number = Number_phase;
            Number_phase++;
        }
        else{
			///Print chosen method
			switch (method(i)) {
				case 0: {
					break;
				}
				case 1: {
					cout << "Phase " << i << ": use of the standard ODF.\n";
					break;
				}
				case 2: {
					cout << "Phase " << i << ": use of the hardening-like ODF.\n";
					break;
				}
				case 3: {
					cout << "Phase " << i << ": use of a Gaussian ODF.\n";
					break;
				}
				case 4: {
					cout << "Phase " << i << ": use of a Lorentzian ODF.\n";
					break;
				}
				case 5: {
					cout << "Phase " << i << ": use of a Pseudo-Voigt ODF.\n";
					break;
				}
				case 6: {
					cout << "Phase " << i << ": use of a Pearson VII ODF.\n";
					break;
				}
				case 7: {
					cout << "Phase " << i << ": use of an uniform ODF.\n";
					break;
				}
				case 30: {
					cout << "Phase " << i << ": use of a multi-peak (" << paramODF(i,0) << ") Gaussian ODF.\n";
					break;
				}
				case 40: {
					cout << "Phase " << i << ": use of a multi-peak (" << paramODF(i,0) << ") Lorentzian ODF.\n";
					break;
				}
				case 50: {
					cout << "Phase " << i << ": use of a multi-peak (" << paramODF(i,0) << ") Pseudo-Voigt ODF.\n";
					break;
				}
				case 60: {
					cout << "Phase " << i << ": use of a multi-peak (" << paramODF(i,0) << ") Pearson VII  ODF.\n";
					break;
				}
				default: {
					cout << "ERROR : cannot find the function labeled " << method << " (phase " << i << ").\n";
					break;			
				}
			}
			///Chargement du fichier externe .dat dans le cas d'une ODF Tabulé
			if (method(i) == 0) {
				///Composition du chemin d'accès
				path.str(std::string());
				path << "data/" << filename[i] << ".dat";
				
				phases.open(path.str().c_str(), ios::in);
				///Décompte du nombre de données
				if(phases) {
					ndata = 0;
					while (!phases.eof())
					{
						getline (phases,chaine1);
						if (chaine1 != "") {
							ndata++;
						}
					}
				}
				else {
					cout << "ERROR: cannot read the file related to the ODF of phase " << i << " using " << path.str() << "\n";
				}
				phases.close();
				
				///Récupération des données
				phases.open(path.str().c_str(), ios::in);
				if(phases) {
					full_ODF.set_size(ndata,2);
					full_ODF = zeros(ndata,2);
					for (int k = 0; k < ndata; k++) {
						phases >> full_ODF(k,0) >> full_ODF(k,1);
					}
				}
				else {
					cout << "ERROR: cannot read the file related to the ODF of phase " << i << " using " << path.str() << "\n";
				}
				phases.close();
				cout << "Using a dedicated ODF file for phase " << i << ": " << path.str() << "...\n";
				
				///Vérification des données
				if (ndata < 2) {
					cout << "ERROR: a file containing less than 2 points can't represent an ODF.\n";
				}
								
				for (int k = 0; k < ndata; k++) {
					if ((full_ODF(k,0) < 0.)||(full_ODF(k,0)>=180.)) {
					cout << "ERROR: Angle must belong to [0 ; 180[. Error line " << k+1 << ".\n";
					}
					full_ODF(k,0) = full_ODF(k,0) * pi/180.;
				}				
				for (int k = 0; k < (ndata-1); k++) {
					if (full_ODF(k,0) == full_ODF(k+1,0)) {
					cout << "ERROR: Only one value for one angle (" << full_ODF(k,0) << "). Error line " << k+2 << ".\n";
					}
					if (full_ODF(k,0) > full_ODF(k+1,0)) {
					cout << "ERROR: Angles must have an ascending order. Error line " << k+2 << ".\n";
					}
				}
				///Utilisation éventuel du fichier ODF en tant que concentrations directes
				flag = false;
				if (ndata == Nphases(i)) {
					cout << "\nThe number of phases is equivalent to the number of points in the ODF. Is the input file considering directly the volume fraction for each orientation family (if Y Nphases.dat will consider directly the input values, if N volume fraction of the family will be built based on the ODF result)? Y or N (Default answer is N)";
					cin >> answer;
					conver = answer;
					cout << "\n";
					if ((conver == 121)||(conver == 89)) {
						flag = true;
						cout << "Warning this mode require the correct construction of the input file : angles with equal distance, starting at 0.\n";
					}
					else if ((conver == 110)||(conver == 78)) {
						flag = false;
					}
					else {
						flag = false;
						cout << "Neither Y nor N was answered. Default answer is N\n";
					}				
				}
			}
			
            dalpha = pi/Nphases(i);			
            for(int j = 0; j < Nphases(i); j++) {
                                
                rvesvs[i][j].resize(rvesvs_init[i].nprops,rvesvs_init[i].nstatev);
                rvesvs[i][j] = rvesvs_init[i];                
                rvesvs[i][j].number = Number_phase;
                
                ///Respect of the chosen angle
                switch (Angle(i)) {
					case 0: {
						rvesvs[i][j].psi_geom = alpha*180./pi;
						break;
					}
					case 1: {
						rvesvs[i][j].theta_geom = alpha*180./pi;
						break;
					}
					case 2: {
						rvesvs[i][j].phi_geom = alpha*180./pi;
						break;
					}
					default: {
						rvesvs[i][j].psi_geom = alpha*180./pi;
						break;					
					}
				}
                
                ///Cas d'une ODF tabulé
                if (method(i) == 0) {
					if (flag) {
						if (fabs((full_ODF(j,0) - alpha)) > precision) {
							cout << "ERROR: At line " << j+1 << ", the file angle " << full_ODF(j,0)*180./pi << " does not match its orientation family, centered in " << alpha*180./pi << ".\n";
						}
						rvesvs[i][j].concentration = full_ODF(j,1);
					}
					else {
						if (j==0) {
							///Détermination de la borne inférieure de la première famille d'orientation (toujours supérieure ou égale à la borne théorique)
							compteur = 0;
							if (full_ODF(ndata-1,0) + precision < (pi - dalpha/2.)) {
								///Cas particulier où le fichier ne contient aucune donnée pour des angles approchant 180° par la limite inférieure, pour la première famille d'orientation
								limiteinf = 0;
							}
							else {							
								while (full_ODF(compteur,0) + precision < (pi - dalpha/2.)) {
										compteur++;
								}
							}
							limiteinf = compteur;
							
							///Détermination de la borne supérieure de la première famille d'orientation  (toujours supérieure ou égale à la borne théorique)
							compteur = 0;
							while (full_ODF(compteur,0) + precision < dalpha/2.) {
									compteur++;
							}
							limitesup = compteur;			
							
							///Trapezoidal integration
							rvesvs[i][j].concentration = 0.;
							if (limiteinf > 0) {
								for (int k = limiteinf; k < (ndata-1); k++) {
									rvesvs[i][j].concentration += (full_ODF(k+1,0)-full_ODF(k,0)) * (full_ODF(k+1,1)+full_ODF(k,1)) / 2.;
								}
								rvesvs[i][j].concentration += ((full_ODF(0,0)+pi)-full_ODF(ndata-1,0)) * (full_ODF(0,1)+full_ODF(ndata-1,1)) / 2.;
							}
							for (int k = 0; k < limitesup; k++) {
								rvesvs[i][j].concentration += (full_ODF(k+1,0)-full_ODF(k,0)) * (full_ODF(k+1,1)+full_ODF(k,1)) / 2.;
							}	
					//cout << (pi - dalpha/2.)*180./pi << "\t" << full_ODF(limiteinf,0)*180./pi << "\t" ;
					//cout << (dalpha/2.*180./pi) << "\t" << full_ODF(limitesup,0)*180./pi << "\t\t" ;		
					//cout << rvesvs[i][j].concentration *180./pi<< "\n";				
							
							///Correction for intervals between 2 orientations family
							///Correction regarding lower bound, always something to be add
							if (limiteinf > 0) {
								if (fabs(full_ODF(limiteinf,0) - (pi - dalpha/2.)) > precision) {
									corr = (full_ODF(limiteinf,0) - (pi - dalpha/2.))/(full_ODF(limiteinf,0) - full_ODF(limiteinf-1,0));
									rvesvs[i][j].concentration += (full_ODF(limiteinf,0)-full_ODF(limiteinf-1,0)) * corr * (full_ODF(limiteinf,1)*(2.-corr)+full_ODF(limiteinf-1,1)*corr) / 2.;
					//cout << " Corr Inf: " << corr << "\t+ " << ((full_ODF(limiteinf,0)-full_ODF(limiteinf-1,0)) * corr * 180./pi * (full_ODF(limiteinf,1)*(2.-corr)+full_ODF(limiteinf-1,1)*corr) / 2.) << "\t\t" << rvesvs[i][j].concentration *180./pi<< "\n";	
								}
							}
							else {
								corr = (full_ODF(0,0) + dalpha/2.)/(full_ODF(0,0)+pi - full_ODF(ndata-1,0));
								rvesvs[i][j].concentration += (full_ODF(0,0)+pi-full_ODF(ndata-1,0)) * corr * (full_ODF(0,1)*(2.-corr)+full_ODF(ndata-1,1)*corr) / 2.;
					//cout << "***SPECIAL MISS_0*** Corr Inf: " << corr << "\t+ " << ((full_ODF(0,0)+pi-full_ODF(ndata-1,0)) * 180./pi * corr * (full_ODF(0,1)*(2.-corr)+full_ODF(ndata-1,1)*corr) / 2.) << "\t\t" << rvesvs[i][j].concentration *180./pi<< "\n";	
							}
							///Correction regarding higher bound, always something to be soustracted
							if (fabs(full_ODF(limitesup,0) - (dalpha/2.)) > precision) {
								if (limitesup == 0) {
									///Cas particulier où le fichier ne contient aucune donnée pour des angles approchant 0° par la limite supérieure, pour la première famille d'orientation
									corr = (full_ODF(0,0) - dalpha/2.)/(full_ODF(0,0)+pi - full_ODF(ndata-1,0));
									rvesvs[i][j].concentration -= (full_ODF(0,0)+pi-full_ODF(ndata-1,0)) * corr * (full_ODF(0,1)*(2.-corr)+full_ODF(ndata-1,1)*corr) / 2.;
					//cout << "***SPECIAL MISS_0*** Corr Sup: " << corr << "\t- " << ((full_ODF(0,0)+pi-full_ODF(ndata-1,0)) * corr * 180./pi * (full_ODF(0,1)*(2.-corr)+full_ODF(ndata-1,1)*corr) / 2.) << "\t\t" << rvesvs[i][j].concentration *180./pi<< "\n";
								}
								else {
									corr = (full_ODF(limitesup,0) - dalpha/2.)/(full_ODF(limitesup,0) - full_ODF(limitesup-1,0));
									rvesvs[i][j].concentration -= (full_ODF(limitesup,0)-full_ODF(limitesup-1,0)) * corr * (full_ODF(limitesup,1)*(2.-corr)+full_ODF(limitesup-1,1)*corr) / 2.;
					//cout << " Corr Sup: " << corr << "\t- " << ((full_ODF(limitesup,0)-full_ODF(limitesup-1,0)) * corr * 180./pi * (full_ODF(limitesup,1)*(2.-corr)+full_ODF(limitesup-1,1)*corr) / 2.) << "\t\t" << rvesvs[i][j].concentration *180./pi<< "\n";
								}
					
							}
							limiteinf = 0;
						}
						
						else {			
							///Détermination de la borne inférieure de la i+1 ème famille d'orientation	(toujours supérieure ou égale à la borne théorique)
							compteur = limiteinf;
							if (full_ODF(ndata-1,0) + precision < (alpha - dalpha/2.)) {
								///Cas particulier où le fichier ne contient aucune donnée pour des angles supérieurs à la borne inférieure de la i+1 ème famille d'orientation
								limiteinf = 0;
							}
							else {
								while (full_ODF(compteur,0) + precision < (alpha - dalpha/2.)) {
										compteur++;
								}
								limiteinf = compteur;
							}
							///Détermination de la borne supérieure de la i+1 ème famille d'orientation (toujours supérieure ou égale à la borne théorique)
							if (full_ODF(ndata-1,0) + precision < (alpha + dalpha/2.)) {
								///Cas particulier où le fichier ne contient aucune donnée pour des angles supérieurs à la borne supérieure de la i+1 ème famille d'orientation
								limitesup = 0;
							}
							else {
								while (full_ODF(compteur,0) + precision < (alpha + dalpha/2.)) {
										compteur++;
								}
								limitesup = compteur;
							}
							
							///Trapezoidal integration
							rvesvs[i][j].concentration = 0.;
							if (limitesup ==0) {
								///Cas particulier où la borne supérieure de la i+1 ème famille d'orientation est supérieure au dernier point
								if (limiteinf > 0) {
									for (int k = limiteinf; k < ndata-1; k++) {
										rvesvs[i][j].concentration += (full_ODF(k+1,0)-full_ODF(k,0)) * (full_ODF(k+1,1)+full_ODF(k,1)) / 2.;
									}
									rvesvs[i][j].concentration += (full_ODF(0,0)+pi-full_ODF(ndata-1,0)) * (full_ODF(0,1)+full_ODF(ndata-1,1)) / 2.;
								}
							}
							else {
								for (int k = limiteinf; k < limitesup; k++) {
									rvesvs[i][j].concentration += (full_ODF(k+1,0)-full_ODF(k,0)) * (full_ODF(k+1,1)+full_ODF(k,1)) / 2.;
								}
							}
					//cout << (alpha - dalpha/2.)*180./pi << "\t" <<  full_ODF(limiteinf,0)*180./pi<< "\t" ;
					//cout << (alpha + dalpha/2.)*180./pi << "\t" <<  full_ODF(limitesup,0)*180./pi<< "\t\t" ;	
					//cout << rvesvs[i][j].concentration *180./pi<< "\n";				
							
							///Correction for intervals between 2 orientations family
							///Correction regarding lower bound, always something to be add
							if (fabs(full_ODF(limiteinf,0) - (alpha - dalpha/2.)) > precision) {
								if (limiteinf == 0) {
									///Cas particulier où la borne inférieure de la i+1 ème famille d'orientation est le premier point
									corr = (full_ODF(0,0)+pi - (alpha - dalpha/2.))/(full_ODF(0,0)+pi - full_ODF(ndata-1,0));
									if (fabs(corr) > 1.) {
										///Cas spécial pour accorder le "modulo pi" avec les corrections à faire, afin d'intégrer le plus petit intervalle
										corr = (full_ODF(0,0) - (alpha - dalpha/2.))/(full_ODF(0,0)+pi - full_ODF(ndata-1,0));
									}
									rvesvs[i][j].concentration += (full_ODF(0,0)+pi-full_ODF(ndata-1,0)) * corr * (full_ODF(0,1)*(2.-corr)+full_ODF(ndata-1,1)*corr) / 2.;
					//cout << "***SPECIAL MISS*** Corr Inf: " << corr << "\t+ " << ((full_ODF(0,0)+pi-full_ODF(ndata-1,0)) * 180./pi * corr * (full_ODF(0,1)*(2.-corr)+full_ODF(ndata-1,1)*corr) / 2.) << "\t\t" << rvesvs[i][j].concentration *180./pi<< "\n";
								}
								else {
									corr = (full_ODF(limiteinf,0) - (alpha - dalpha/2.))/(full_ODF(limiteinf,0) - full_ODF(limiteinf-1,0));
									rvesvs[i][j].concentration += (full_ODF(limiteinf,0)-full_ODF(limiteinf-1,0)) * corr * (full_ODF(limiteinf,1)*(2.-corr)+full_ODF(limiteinf-1,1)*corr) / 2.;
					//cout << " Corr Inf: " << corr << "\t+ " << ((full_ODF(limiteinf,0)-full_ODF(limiteinf-1,0)) * corr * 180./pi * (full_ODF(limiteinf,1)*(2.-corr)+full_ODF(limiteinf-1,1)*corr) / 2.) << "\t\t" << rvesvs[i][j].concentration *180./pi<< "\n";	
								}
							}
							///Correction regarding higher bound, always something to be soustracted		
							if (limitesup == 0) {
								///Cas particulier où la borne supérieure de la i+1 ème famille d'orientation est le premier point
								corr = (full_ODF(0,0)+pi -(alpha + dalpha/2.))/(full_ODF(0,0)+pi - full_ODF(ndata-1,0));
								if (fabs(corr) > 1.) {
									///Cas spécial pour accorder le "modulo pi" avec les corrections à faire, afin d'intégrer le plus petit intervalle
									corr = (full_ODF(0,0) -(alpha + dalpha/2.))/(full_ODF(0,0)+pi - full_ODF(ndata-1,0));
								}
								rvesvs[i][j].concentration -= (full_ODF(0,0)+pi-full_ODF(ndata-1,0)) * corr * (full_ODF(0,1)*(2.-corr)+full_ODF(ndata-1,1)*corr) / 2.;
					//cout << "***SPECIAL MISS*** Corr Sup: " << corr << "\t- " << ((full_ODF(0,0)+pi-full_ODF(ndata-1,0)) * corr * 180./pi * (full_ODF(0,1)*(2.-corr)+full_ODF(ndata-1,1)*corr) / 2.) << "\t\t" << rvesvs[i][j].concentration *180./pi<< "\n";
							}
							else {
								if (fabs(full_ODF(limitesup,0) - (alpha + dalpha/2.)) > precision) {
									corr = (full_ODF(limitesup,0) - (alpha + dalpha/2.))/(full_ODF(limitesup,0) - full_ODF(limitesup-1,0));
									rvesvs[i][j].concentration -= (full_ODF(limitesup,0)-full_ODF(limitesup-1,0)) * corr * (full_ODF(limitesup,1)*(2.-corr)+full_ODF(limitesup-1,1)*corr) / 2.;
					//cout << " Corr Sup: " << corr << "\t- " << ((full_ODF(limitesup,0)-full_ODF(limitesup-1,0)) * corr  * 180./pi* (full_ODF(limitesup,1)*(2.-corr)+full_ODF(limitesup-1,1)*corr) / 2.) << "\t\t" << rvesvs[i][j].concentration *180./pi<< "\n";
								}
							}
						}
					}
				}
				
                ///Cas d'une ODF basée sur une fonction prédéfinie
				else {
                ///Simpson integration
					if (j==0) {
						rvesvs[i][j].concentration = dalpha/6. * (ODF(pi-dalpha/2, method(i), paramODF.row(i).t(), radian, dec)+4.*ODF(alpha, method(i), paramODF.row(i).t(), radian, dec)+ODF(alpha+dalpha/2, method(i), paramODF.row(i).t(), radian, dec));
					}
					else {
						rvesvs[i][j].concentration = dalpha/6. * (ODF(alpha-dalpha/2, method(i), paramODF.row(i).t(), radian, dec)+4.*ODF(alpha, method(i), paramODF.row(i).t(), radian, dec)+ODF(alpha+dalpha/2, method(i), paramODF.row(i).t(), radian, dec));
					}
				}
				
                ///Normalization
                normODF += rvesvs[i][j].concentration;
                alpha += dalpha;
                Number_phase++;
            }
            for(int j = 0; j < Nphases(i); j++) {	
                rvesvs[i][j].concentration *= (rvesvs_init[i].concentration / normODF);
            }
        }
    }

	for(int i = 0; i < phases_init; i++) {
        if (rvesvs_init[i].coatingof != 0) {
            for(int j = 0; j < Nphases(i); j++) {
                rvesvs[i][j].coatingof = rvesvs[rvesvs_init[i].coatingof][j].number;
            }
        }
    }
	cout << "Writing Nphases.dat file.\n";
    
	ofstream Nphases_file("data/Nphases.dat", ios::out);
	if(Nphases_file) {
        Nphases_file << "Number\t" << "Coatingof\t" << "umat\t" << "c\t" << "psi\t" << "theta\t" << "phi\t" << "a1\t" << "a2\t" << "a3\t" << "nprops\t" << "nstatev\t" << "props\n";        

        for(int i = 0; i < phases_init; i++) {
            for(int j = 0; j < Nphases(i); j++) {
                Nphases_file << rvesvs[i][j].number << "\t" << rvesvs[i][j].coatingof << "\t\t" << rvesvs[i][j].umat_name << "\t" << rvesvs[i][j].concentration << "\t" << rvesvs[i][j].psi_geom << "\t" << rvesvs[i][j].theta_geom << "\t" << rvesvs[i][j].phi_geom << "\t" << rvesvs[i][j].a1 << "\t" << rvesvs[i][j].a2 << "\t" << rvesvs[i][j].a3 << "\t" << rvesvs[i][j].nprops << "\t" << rvesvs[i][j].nstatev;
                    for(int k=0; k<rvesvs[i][j].dimprops(); k++) {
                        Nphases_file << "\t" << rvesvs[i][j].props[k];
                    }
                Nphases_file << "\n";
            }
        }
    }
	else
		cout << "Error: cannot write in Nphases.dat file using ODF2phases.hpp\n";
    
	Nphases_file.close();
}

} //namespace smart
