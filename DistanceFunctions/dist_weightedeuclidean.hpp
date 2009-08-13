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
#ifndef __dist_weightedeuclidean_hpp__
#define __dist_weightedeuclidean_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include <iostream>
#include <sstream>

class WeightedEuclideanDistance : public BaseDistance {

private:
  double alpha_;
  uint numItemsFirstWeight;

public:

  WeightedEuclideanDistance(double alpha) {
    alpha_ = alpha;
    numItemsFirstWeight = 2;
  }

  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
    double result1 = 0.0, result2 = 0.0;
    
    const VectorFeature* db=dynamic_cast<const VectorFeature*>(databaseFeature);
    const VectorFeature* query=dynamic_cast<const VectorFeature*>(queryFeature);

    if(db && query && (db->size() > 2) && (query->size() > 2)) {
    
      double tmp;
      uint db_size = db->size();
      for(uint i = 0; i < db_size - 2; ++i) {
        tmp = (db->operator[](i)) - (query->operator[](i));
        result1 += tmp * tmp;
      }
      for(uint i = db_size - 2; i < db_size; ++i) {
        tmp = (db->operator[](i)) - (query->operator[](i));
        result2 += tmp * tmp;
      }

      return result1 * 2.0 * (1.0 - alpha_) + result2 * 2.0 * alpha_;
      
    } else {
      ERR << "Features not comparable" << ::std::endl;
      return -1.0;
    }
  }

  virtual ::std::string name() { 
    ::std::ostringstream strstream(::std::ostringstream::out);
    strstream << "weightedeuclidean, alpha=" << alpha_;
    return strstream.str();
  }
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
};

#endif
