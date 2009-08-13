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
#ifndef __dist_rast_hpp__
#define __dist_rast_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include <iostream>
#include "lfposclusteridfeature.hpp"

class RASTDistance : public BaseDistance {
private:
  int eps_;
  double tolerance_;
  int mindx_,maxdx_,mindy_,maxdy_;
  int minq_;
  double amin_,amax_;
  double minscale_,maxscale_;

  int maximumMatches_;

public:

  RASTDistance(int eps,double tolerance,int mindx,int maxdx,int mindy,int maxdy,int minq, double amin,double amax,double minscale,double maxscale):
    eps_(eps),tolerance_(tolerance),mindx_(mindx),maxdx_(maxdx),
    mindy_(mindy),maxdy_(maxdy),minq_(minq),
    amin_(amin),amax_(amax),minscale_(minscale),maxscale_(maxscale) {
  }
  
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature);
  
  virtual ::std::string name() {return "RAST";}
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
};

#endif
