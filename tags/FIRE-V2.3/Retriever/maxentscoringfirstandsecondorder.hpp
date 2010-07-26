/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

FIRE is distributed in the hope that it will be useful, but WITHOyT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software Foundation, Inc.,
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#ifndef __maxentsecondandfirstorderscoring_hpp__
#define __maxentsecondandfirstorderscoring_hpp__
#include "diag.hpp"
#include "basescoring.hpp"


class MaxEntFirstAndSecondOrderScoring : public BaseScoring {
private:
  ::std::vector<double> lambdas_;
  double factor_, offset_;
  uint numLambdasCls_, numLambdasCom_, numLambdasPool_, numCls_;  
public:
  MaxEntFirstAndSecondOrderScoring(uint no=0) : lambdas_(no) {type_="maxent1st2nd";}
  MaxEntFirstAndSecondOrderScoring(const ::std::string &filename)  {type_="maxent1st2nd";load(filename);}
  MaxEntFirstAndSecondOrderScoring(const ::std::string &filename, const uint&)  {type_="maxent1st2nd";if(filename!="") {load(filename);} else {DBG(10) << "need configfile" << ::std::endl;}}

  /** read a lambda file trained with ME-first for using in this
      Scoring function */
  void load(const ::std::string &filename);

  /** 
   * This function is strongly inspired by calcPC from the MaxEnt
   * Toolkit by Daniel Keysers. We assume that class 1 is relevant and
   * return the probability for an image to be from class 1 as score.
   * the parameter is a distance vector of an image and the score is
   * returned. For this the distance vector has to be of the same sort
   * as the distance vectors that are saved by save distances.
   */
  virtual double getScore(const ::std::vector<double>& dists);

  virtual const ::std::string settings();
};

#endif
