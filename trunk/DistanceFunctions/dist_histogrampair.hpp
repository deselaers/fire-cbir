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
#ifndef __dist_histogrampair_hpp__
#define __dist_histogrampair_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include "histogrampairfeature.hpp"
#include <iostream>
#include <limits>

class HistogramPairDistance : public BaseDistance {
private:
  BaseDistance * groundDist_;
  double centerWeight_;
public:
  HistogramPairDistance(BaseDistance *gdist, double cweight) : groundDist_(gdist), centerWeight_(cweight) {
    DBG(10) << "groundist=" << groundDist_->name() << " centerWeight=" << centerWeight_ << ::std::endl;
  }
  
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
    double result=0.0;
    
    const HistogramPairFeature* db=dynamic_cast<const HistogramPairFeature*>(databaseFeature);
    const HistogramPairFeature* query=dynamic_cast<const HistogramPairFeature*>(queryFeature);
    
    if(db && query) {
      double cdist, bgdist;
      double minDist=::std::numeric_limits<double>::max();
      double aktDist;

      for(uint i=0;i<db->nOfHistograms();++i) {
        for(uint j=0;j<query->nOfHistograms();++j) {

          //distance between foreground histograms
          HistogramFeature qhisto=query->histo(j);
          HistogramFeature dbhisto=db->histo(i);
          cdist=groundDist_->distance(&qhisto,&dbhisto);
          
          //distance between background histograms
          qhisto=query->bghisto(j);
          dbhisto=db->bghisto(i);
          bgdist=groundDist_->distance(&qhisto,&dbhisto);
          
          // combine and find smallest
          aktDist=centerWeight_*cdist+(1.0-centerWeight_)*bgdist;
          if(aktDist<minDist) {
            minDist=aktDist;
          }
        }
      }
      result=minDist;
      return result;
    } else {
      ERR << "Features not comparable" << ::std::endl;
      return -1.0;
    }
  }

  virtual ::std::string name() {return "histogrampair";}
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
};

#endif
