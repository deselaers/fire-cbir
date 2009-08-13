/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

FIRE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __gaussiandensity_hpp__
#define __gaussiandensity_hpp__

#include <vector>
#include "diag.hpp"

typedef ::std::vector<double> DoubleVector;

class GaussianDensity {
public:
  uint dim;
  uint elements;
  DoubleVector mean;
  DoubleVector sigma;
  
  GaussianDensity(uint d=0) : dim(d), elements(0) , mean(d,0.0), sigma(d,0.0) {}

};

#endif
 
