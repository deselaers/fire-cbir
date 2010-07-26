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

#ifndef __lfposclusteridfeature__hpp__
#define __lfposclusteridfeature__hpp__

#include <vector>
#include <string>
#include <iostream>
#include "basefeature.hpp"
#include "localfeatures.hpp"

/** cannot be derived from vectorvectorfeature as we do not want to
    store doubles, but ints */
class LFPositionClusterIdFeature : public BaseFeature {
private:
  ::std::vector< PositionPair > positions_;
  uint dim_;
  ::std::string comment_;
  ::std::vector< ::std::vector <int> > data_;
  
public:
  
  LFPositionClusterIdFeature() :  positions_(::std::vector< PositionPair >()), dim_(0){
    type_=FT_LFPOSCLSIDFEAT;
  }

  virtual ~LFPositionClusterIdFeature() {}
  
  virtual LFPositionClusterIdFeature* clone() const {return new LFPositionClusterIdFeature(*this);}
  
  /// read this feature from the given stream
   virtual bool read(::std:: istream & is);

  /// write this feature to the given stream
  virtual void write(::std::ostream & os); 
  
  ///  return the number stored features
  const uint size() const {return data_.size();}

  /// position of the idx-th patch
  PositionPair& position(uint idx) {return positions_[idx];}
  const PositionPair& position(uint idx) const {return positions_[idx];}
  
  // return number of features in feature vectors
  uint& dim() {return dim_;}
  const uint& dim() const {return dim_;}

  const ::std::string& comment() const {return comment_;}
  ::std::string& comment() {return comment_;}
   
  /// access the idx-th vector 
  virtual const ::std::vector<int>& operator[](uint idx) const { return data_[idx];}
  virtual ::std::vector<int>& operator[](uint idx) { return data_[idx];}

  virtual void push_back(const PositionPair& pos, const ::std::vector<int>& feat);

};

#endif
