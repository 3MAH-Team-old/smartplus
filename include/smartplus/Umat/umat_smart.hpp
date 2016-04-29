///@file abaqus2smart.hpp
///@brief Procedure that transfer the abaqus format to a SMART+ format:
///@brief Implemented in 1D-2D-3D
///@author Chemisky
///@version 1.0
///@date 12/03/2013

#pragma once
#include <armadillo>
#include "../Libraries/Phase/phase_characteristics.hpp"

namespace smart{

///@param stress array containing the components of the stress tensor (dimension ntens)
///@param stran array containing total strain component (dimension ntens) at the beginning of increment
///@param dstran array containing the component of total strain increment (dimension ntens)
///@param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
///@param dtime time increment
///@param temperature temperature avlue at the beginning of increment
///@param Dtemperature temperature increment
///@param ndi number of direct stress components
///@param nshr number of shear stress components
///@param drot rotation increment matrix (dimension 3*3)
///@param pnewdt ratio of suggested new time increment
///@param celent characteristic element length
///@param dfgrd0 array containing the deformation gradient at the beginning of increment (dimension 3*3)
///@param dfgrd1 array containing the deformation gradient at the end of increment (dimension 3*3)
///@param noel element number
///@param npt integration point number
///@param layer layer number - not used
///@param kspt section point number within the current layer - not used
///@param kstep step number
///@param kinc increment number    
    
void abaqus2smart(double *, double *, const double *, const double *, const double *, const double &, const double &, const double &, const int &,const double *, const int &, double *, const double &, const int &, const int &, const double *, arma::vec &, arma::mat &, arma::vec &, arma::vec &, double &, double &, double &, double &, arma::vec &, arma::vec &, double &, arma::mat &, bool &);

void select_umat_T(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, const bool &, double &);
    
void select_umat_M(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, const bool &, double &);

void run_umat_T(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, double &);

void run_umat_M(phase_characteristics &, const arma::mat &, const double &, const double &, const int &, const int &, bool &, double &);

void smart2abaqus(double *, double *, double *, const int &, const int &, const arma::vec &, const arma::mat &, const arma::vec &, double &, const double &);
    
} //namespace smart
