/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

FIRE is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software Foundation, Inc.,
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#ifndef __vectorvectorfeature_hpp__
#define __vectorvectorfeature_hpp__

#include <vector>
#include "basefeature.hpp"

/** A vector of vectors as feature. E.g. can be used for gabor
    features, local features, image sets, ... */

class VectorVectorFeature : public BaseFeature { 
protected:

  /// the data, outer vector: vector of vectors
  ::std::vector< ::std::vector<double> > data_;
public:
  
  /// destructor
  virtual ~VectorVectorFeature() { data_.clear();}
  
//  virtual VectorVectorFeature* clone() const {return new VectorVectorFeature(*this);}
  
  
  /// access the idx-th vector 
  virtual const ::std::vector<double>& operator[](uint idx) const { return data_[idx];}
  
  /// access the idx-th vector 
  virtual ::std::vector<double>& operator[](uint idx) { return data_[idx];}
  
  /// how many vectors are there?
  virtual const uint numberOfVectors() const {return data_.size();}
  
  /// calculate the total size
  virtual const uint size() const {
    uint res=0;
    for(uint i=0;i<data_.size();++i) {
      res+=data_[i].size();
    }
    return res;
  }

  virtual void push_back(const ::std::vector<double>& d) {
    data_.push_back(d);
  }
  
  
};

#endif
