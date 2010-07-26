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
#ifndef __dist_weighted_l1_hpp__
#define __dist_weighted_l1_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include <iostream>

class WeightedL1Distance : public BaseDistance {
public:
  WeightedL1Distance(int maxiter, double stepwidth, double regularisationweight) : FLT_MAX(100000000.0), FLT_MIN(-100000000.0), maxIterations_(maxiter), stepWidth_(stepwidth), regularisationWeight_(regularisationweight), gradientWeight_(1.0-regularisationweight) , rampa(10){
    DBG(10) << VAR(maxIterations_) << " " << VAR(stepWidth_) << " " << VAR(regularisationweight) << std::endl;
  }
  
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature);
  virtual ::std::string name() {return "weightedl1";}
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
  
  virtual void tune(const std::vector<const BaseFeature*>& posFeat, const std::vector<const BaseFeature*>& negFeat);
  
private:
  float wdisl1(float * i, float *j, int dim); 
  float sigmoide(float x);

  std::vector<double> weight_;
  float FLT_MAX, FLT_MIN;
  int maxIterations_;
  double  stepWidth_,regularisationWeight_,gradientWeight_;
  float rampa;
};

#endif
