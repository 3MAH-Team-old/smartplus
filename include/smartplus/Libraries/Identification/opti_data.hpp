/* This file is part of SMART+ private.
 
 Only part of SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file opti_data.hpp
///@brief Handle of data from optimization
///@author Chemisky & Despringre
///@version 1.0

#pragma once
#include <iostream>
#include <armadillo>

namespace smart{

//======================================
class opti_data
//======================================
{
	private:

	protected:

	public :
        std::string name;
		int number;
		int ndata;
		int ninfo;
		int ncolumns;
        arma::Col<int> c_data;
        arma::mat data;
		
		opti_data(); 	//default constructor
		opti_data(int, int);	//constructor - allocates memory for statev
        opti_data(std::string, int, int, int, int); //Constructor with parameters
		opti_data(const opti_data &);	//Copy constructor
		~opti_data();
		
		int dimdata () const {return ndata;}       // returns the number of data points
		int diminfo () const {return ninfo;}       // returns the number of informations at each datapoint

		void constructc_data();
		void constructdata();

        void import(std::string, std::string, int=0);
				
		virtual opti_data& operator = (const opti_data&);
		
        friend  std::ostream& operator << (std::ostream&, const opti_data&);
};

} //namespace smart
