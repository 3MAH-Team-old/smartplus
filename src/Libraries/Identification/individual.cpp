/* This file is part of SMART+ private.
 
 Only part of SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
*/

///@file individual.hpp
///@brief individual for genetic algorithm (among others)
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

///@brief default constructor
//----------------------------------------------------------------------
individual::individual()
//----------------------------------------------------------------------
{
	np=0;
	cout = 0.;
	id = 0;
	rank=0;
}

/*!
  \brief Constructor
  \param m : number of parameters
  \param init boolean that indicates if the constructor has to initialize  (default value is true) \n
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/
//-------------------------------------------------------------
individual::individual(const int &n, const int &idnumber)
//-------------------------------------------------------------
{
	np=n;
	cout = 0.;
	id = idnumber;
	rank = 0;
	
    if (n>0) {
        p = zeros(n);
    }

}

/*!
  \brief Copy constructor
  \param gp individual object to duplicate
*/
//------------------------------------------------------
individual::individual(const individual& gp)
//------------------------------------------------------
{
	np=gp.np;
	cout=gp.cout;
	id=gp.id;
	rank=gp.rank;
	p=gp.p;
}

/*!
  \brief destructor
*/
individual::~individual() {}

///@brief Construct method to construct an element that had size zero
//-------------------------------------------------------------
void individual::construct()
//-------------------------------------------------------------
{
	assert(np>0);
	p = zeros(np);
}

/*!
  \brief Standard operator = for individual
*/
//----------------------------------------------------------------------
individual& individual::operator = (const individual& gp)
//----------------------------------------------------------------------
{
    np=gp.np;
    cout=gp.cout;
    id=gp.id;
    rank=gp.rank;
    p=gp.p;
    
    return *this;
}

/*!
\brief Ostream operator << for individual
*/
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const individual& gp)
//--------------------------------------------------------------------------
{
	assert(gp.np>0);	
	
	s << "Parameters of individual\n";
	s << "id = " << gp.id << "\n" ;
	s << "rank = " << gp.id << "\n" ;	
	s << "cost = " << gp.cout << "\n" ;
	s << "p = " << gp.p << "\n" ;

	return s;
}
    
} //namespace smart