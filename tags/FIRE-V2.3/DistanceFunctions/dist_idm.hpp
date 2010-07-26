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
#ifndef __dist_imagedistortionmodel_hpp__
#define __dist_imagedistortionmodel_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include "imagefeature.hpp"
#include <iostream>

class ImageDistortionModelDistance : public BaseDistance {
private: 
  uint wr1_, wr2_;
  double threshold_;
  bool sobel_;
  bool gray_;

  double pixelDist(const ImageFeature& query, 
                   const ImageFeature& db, 
                   const uint x, const uint y, const uint xx, const uint yy) const;
  
public:
  ImageDistortionModelDistance(uint wr1=2, uint wr2=1, double threshold=-2.0, bool sobel=true, bool gray=true) :
    wr1_(wr1), wr2_(wr2), threshold_(threshold), sobel_(sobel) , gray_(gray){
    DBG(35) << "WR1=" << wr1_ 
            << " WR2=" << wr2_ 
            << " TH=" << threshold_ 
            << " SOBEL=" << sobel_ 
            << " GRAY="<< gray_ << ::std::endl;
  }
  
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature);
  
  virtual ::std::string name() {return "idm";}
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
};

#endif
