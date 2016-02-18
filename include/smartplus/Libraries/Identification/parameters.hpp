/* This file is part of SMART+ private.
 
 Only part of SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file parameters.hpp
///@brief Handle of input parameters
///@author Chemisky
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

namespace smart{

//======================================
class parameters
//======================================
{
	private:

	protected:

	public :
    
		int number;     //number of the parameters
        double value;       //Value of the parameter
    
        double min_value;   //Minimum value of the parameter
        double max_value;   //Maximum value of the parameter
    
        std::string key;    //A unique key utilized to replace the parameters in file(s)
        int ninput_files;
        std::vector<std::string> input_files;   //vector of files impacted (automaticaly filed for some parameter types)
    
		parameters(); 	//default constructor
		parameters(const int&, const double&, const double&);	//constructor - number, min and max values
        parameters(const int&, const double&, const double&, const std::string&, const int &, const std::vector<std::string>&); //Constructor with parameters
		parameters(const parameters &);	//Copy constructor
		~parameters();
		
		int dimfiles () const {return ninput_files;}  // returns the number of files associated to this parameter
    
        void update(const double &);
        void resize(const int&);
				
		virtual parameters& operator = (const parameters&);
		
        friend  std::ostream& operator << (std::ostream&, const parameters&);
};

} //namespace smart
