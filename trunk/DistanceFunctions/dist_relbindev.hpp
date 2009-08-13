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
#ifndef __dist_relativebindeviation_hpp__
#define __dist_relativebindeviation_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include <iostream>

class RelativeBinDeviationDistance : public BaseDistance {
public:
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
    double result=0.0;
    
    const VectorFeature* db=dynamic_cast<const VectorFeature*>(databaseFeature);
    const VectorFeature* query=dynamic_cast<const VectorFeature*>(queryFeature);

    if(db && query) {
      double tmp;
      double n1, n2, div;
      for(uint i=0;i<db->size();++i) {
        n1=(db->operator[](i));
        n2=(query->operator[](i));

        div=n1*n1+n2*n2;
        div*=0.5;
        
        tmp=n1-n2;
        tmp*=tmp;
        tmp=sqrt(tmp);
        
        result+=tmp/div;
      }
      return result/db->size();
    } else {
      ERR << "Features not comparable" << ::std::endl;
      return -1.0;
    }
  }

  virtual ::std::string name() {return "relativebindeviation";}
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
};

#endif
