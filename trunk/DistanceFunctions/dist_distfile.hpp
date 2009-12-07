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
#ifndef __dist_distancefile_hpp__
#define __dist_distancefile_hpp__

#include "diag.hpp"
#include "basedistance.hpp"
#include "basescoring.hpp"
#include "distancefilefeature.hpp"

class DistanceFileDistance : public BaseDistance {
private:
    BaseScoring *scoring_;
    DistanceFileFeature *queryFeat_;

    uint currentLine_;
    bool clearing_;

public:

    DistanceFileDistance(const ::std::string& sname="linear", bool clearing=true);

    virtual ~DistanceFileDistance() {
        delete scoring_;
    }
    
    virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature);
    
    virtual ::std::string name() {
        return "distfile";
    }
    
    virtual void start(const BaseFeature *q) {
        currentLine_=0;
        queryFeat_=dynamic_cast<DistanceFileFeature*>(const_cast<BaseFeature *>(q));
        if(not queryFeat_->loaded()) {
            queryFeat_->loadYourself();
        }
    }

    virtual void stop() {
        if(clearing_) {
            queryFeat_->clear();
        }
    }
};

#endif
