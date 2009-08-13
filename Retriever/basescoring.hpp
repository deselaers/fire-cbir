/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

FIRE is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software Foundation, Inc.,
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#ifndef __basescoring_hpp__
#define __basescoring_hpp__
#include "diag.hpp"
#include "imagecomparator.hpp"

/** 
 * a class to calculate scores from distance vectors.  so far this
 * class is implemented in LinearSoring and MaxEntScoring.
 * 
 * a scorer is expected to return high values for relevant images and
 * low values for irrelevant images. Basically a scorer can be seen as
 * a classfifier classifying images into relevant and irrelevant
 * ones. And it returns the confidence for an image to be relevant. 
 */

class BaseScoring {
protected:
  ::std::string type_;
public:

  
  BaseScoring() {};
  BaseScoring(const ::std::string&, const uint&) {};
  
  virtual void load(const ::std::string& filename)=0;
  
  /// destructor :-)
  virtual ~BaseScoring() {}
  /// given a normalized distance vector return the score
  virtual double getScore(const ::std::vector<double>& dists)=0;
  
  /// give the name of the scoring as reference
  virtual ::std::string& type() {return type_;}
  
  /// give the name of the scoring as const reference
  virtual const ::std::string& type() const {return type_;}
  
  /// give settings for this scoring algorithm
  virtual const ::std::string settings() {return "";}

  // this is necessary for da factory
  template<class T>
  static BaseScoring* create(const ::std::string& configfile, const uint& numberOfDistances) {
    return new T(configfile, numberOfDistances);
  }

};



#endif
