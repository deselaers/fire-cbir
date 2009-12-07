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

#include <limits>
#include "dist_distfile.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include "basescoring.hpp"
#include "distancefilefeature.hpp"
#include <iostream>
#include "getscoring.hpp"

using namespace std;

DistanceFileDistance::DistanceFileDistance(const ::std::string& sname, bool clearing) :  currentLine_(0), clearing_(clearing) {
    scoring_=getScoring(sname, 0);
}


 double DistanceFileDistance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
    
    const DistanceFileFeature* db=dynamic_cast<const DistanceFileFeature*>(databaseFeature);
    const DistanceFileFeature* query=dynamic_cast<const DistanceFileFeature*>(queryFeature);
    
    if (db && query) {
        //find the right line in the queryFeature
        
        double score=scoring_->getScore((*query)[currentLine_]);
        ++currentLine_;
        return -log(score); // here we have to retransform the scoring
        // into the -log(score) because this is put
        // into exp(- X) in the linear scoring where
        // it is treated as a distance.
        
    } else {
        ERR << "Features not comparable: need DistanceFileFeatures" << ::std::endl;
        return -1.0;
    }
}
