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

///@file peak.cpp
///@brief peak charcteristics from a density function:
///@version 0.9

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <smartplus/parameter.hpp>
#include <smartplus/Libraries/Material/peak.hpp>
#include <smartplus/Libraries/Maths/stats.hpp>

using namespace std;
using namespace arma;

namespace smart{
    
//=====Private methods for peak===================================

//=====Public methods for peak============================================

/*!
 \brief default constructor
 */
    
//-------------------------------------------------------------
peak::peak()
//-------------------------------------------------------------
{
    method = 0;
    mean = 0.;
    s_dev = 0.;
    width = 0.;
    ampl = 0.;
}
    
    /*!
     \brief Constructor with parameters
     \f$ \textbf{Examples :} \f$ \n
     */
    
//-------------------------------------------------------------
peak::peak(const int &mmethod, const double &mmean, const double &ms_dev, const double &mwidth, const double &mampl, const vec &mparam, const vec &mdensities)
//-------------------------------------------------------------
{
    method = mmethod;
    mean = mmean;
    s_dev = ms_dev;
    width = mwidth;
    ampl = mampl;
    param = mparam;
    densities = mdensities;
}
    
/*!
 \brief Copy constructor
 \param s phase_characteristics object to duplicate
 */

//------------------------------------------------------
peak::peak(const peak &pc)
//------------------------------------------------------
{
    method = pc.method;
    mean = pc.mean;
    s_dev = pc.s_dev;
    width = pc.width;
    ampl = pc.ampl;
    param = pc.param;
    densities = pc.densities;
}

/*!
 \brief Destructor
 
 Deletes peak
 */
    
//-------------------------------------
peak::~peak() {}
//-------------------------------------
    
/*!
 \brief Standard operator = for peak
 */

//-------------------------------------------------------------
void peak::construct(const double &mmean, const double &ms_dev, const double &mwidth, const double &mampl, const double &ampl_limit, const bool &radian)
//-------------------------------------------------------------
{

    if (radian) {
        mean = mmean;
        s_dev = ms_dev;
        width = mwidth;
    }
    else {
        mean = mmean*pi/180.;
        s_dev = ms_dev*pi/180.;
        width = mwidth*pi/180.;
    }
    ampl = mampl;
    if (fabs(ampl) < ampl_limit) {
        ampl = 1.;
    }
}
    
//-------------------------------------------------------------
double peak::get_density(const double &theta, const double &ampl_limit)
//-------------------------------------------------------------
{
    
/*    switch (method) {
        case 1: {
            return ODF_sd(theta, mean, param);
        }
        case 2: {
            return ODF_hard(theta, mean, s_dev, param) + ODF_hard(theta - pi, mean, s_dev, param) + ODF_hard(theta + pi, mean, s_dev, param);
        }
        case 3: {
            return Gaussian(theta, mean, s_dev, ampl) + Gaussian(theta - pi, mean, s_dev, ampl) + Gaussian(theta + pi, mean, s_dev, ampl);
        }
        case 4: {
            return Lorentzian(theta, mean, width, ampl) + Lorentzian(theta - pi, mean, width, ampl) + Lorentzian(theta + pi, mean, width, ampl);
        }
        case 5: {
            return PseudoVoigt(theta, param(0), mean, width, s_dev, ampl) + PseudoVoigt(theta - pi, param(0), mean, width, s_dev, ampl) + PseudoVoigt(theta + pi, param(0), mean, width, s_dev, ampl);
        }
        case 6: {
            assert(Width > 0.);
            double inv_width = 1./width;
            double max = param(1);
            if (fabs(max) < ampl_limit) {
                max = 1.;
            }
            return Pearson7(theta, mean, inv_width, param(0), max) + Pearson7(theta - pi, mean, inv_width, param(0), max) + Pearson7(theta + pi, mean, inv_width, param(0), max);
        }
        case 7: {
            return 1.;
        }
    }*/
}
    
//----------------------------------------------------------------------
peak& peak::operator = (const peak& pc)
//----------------------------------------------------------------------
{
    method = pc.method;
    mean = pc.mean;
    s_dev = pc.s_dev;
    width = pc.width;
    ampl = pc.ampl;
    param = pc.param;
    densities = pc.densities;
    
    return *this;
}
    
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const peak& pc)
//--------------------------------------------------------------------------
{
    s << "Display peak characteristics:\n";
    s << "method: \t" << pc.method << "\n";
    s << "mean = \t" << pc.mean << "\n";
    s << "standard deviation = \t" << pc.s_dev << "\n";
    s << "width = \t" << pc.width << "\n";
    
    s << "\n\n";
    
    return s;
}
    
} //namespace smart
