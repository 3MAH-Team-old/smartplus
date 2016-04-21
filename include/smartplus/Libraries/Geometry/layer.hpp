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

///@file layer.hpp
///@brief Characteristics of an layer geometry, which hereditates from:
///-geometry
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include "geometry.hpp"

namespace smart{

//======================================
class layer : public geometry
//======================================
{
	private:

	protected:

	public :

        int layerup;
        int layerdown;
    
		layer(); 	//default constructor
        layer(const double &, const int &, const int &); //Constructor with parameters

		layer(const layer&);	//Copy constructor
        virtual ~layer();
    
		virtual layer& operator = (const layer&);
		
        friend std::ostream& operator << (std::ostream&, const layer&);
};

} //namespace smart
