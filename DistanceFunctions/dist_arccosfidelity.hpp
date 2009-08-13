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
#ifndef __dist_arccosfidelity_hpp__
#define __dist_arccosfidelity_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include <iostream>

class ArccosFidelityDistance : public BaseDistance {
public:
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
    double result=0.0;
    
    const VectorFeature* db=dynamic_cast<const VectorFeature*>(databaseFeature);
    const VectorFeature* query=dynamic_cast<const VectorFeature*>(queryFeature);

    if(db && query) {
    
      double n1, n2;

      for(uint i=0;i<db->size();++i) {
        n1=(db->operator[](i));
        n2=(query->operator[](i));
        result+=sqrt(n1)*sqrt(n2);
      }
      
      return 2.0/M_PI*acos(result);
      
    } else {
      ERR << "Features not comparable" << ::std::endl;
      return -1.0;
    }
  }

  virtual ::std::string name() {return "arccosfidelity";}
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
};

#endif
