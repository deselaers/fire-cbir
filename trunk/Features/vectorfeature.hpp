
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
#ifndef __vectorfeature_hpp__
#define __vectorfeature_hpp__

#include <vector>
#include "basefeature.hpp"
#include <string>


class VectorFeature : public BaseFeature {
protected:
  ::std::vector<double> data_;
public:
  
  /// initialize a vector feature, prepared to take no entries
  VectorFeature() : data_() {type_=FT_VEC;}
  
  /// initialize a vector feature, prepared to take "size" entries
  VectorFeature(uint size) : data_(size) {type_=FT_VEC;}
  
  /// take a vector and make this vector feature contain this data
  VectorFeature(const ::std::vector<double> &in) : data_(in) {type_=FT_VEC;}

  //  VectorFeature::VectorFeature(const VectorFeature &src) : BaseFeature(), data_(src.data()) {type_=FT_VEC;}
  
  virtual VectorFeature* clone() const {
    return new VectorFeature(*this);
  }
  
  virtual ~VectorFeature() {}
  
  //inherited from BaseFeature
  //virtual void load(const ::std::string &filename);
  //virtual void save(const ::std::string &filename);

  VectorFeature operator-(const VectorFeature &v ) const;
  VectorFeature & operator-=(const VectorFeature &v );

  ///inherited from BaseFeature
  virtual bool read(::std:: istream & is);
  ///inherited from BaseFeature
  virtual bool readBinary(::std:: istream & is);
  ///inherited from BaseFeature
  virtual void write(::std::ostream & os);  
  ///inherited from BaseFeature
  virtual void writeBinary(::std::ostream & os);  
  
  
  ///access to data from vector, get entry number idx, this allows for modification
  virtual double operator[](uint idx) const {  return data_[idx];}
  ///access to data from vector, get entry number idx, this does not allow for modification
  virtual double& operator[](uint idx) {  return data_[idx];}
  
  ///return the number of entries of this vector feature
  virtual const uint size() const;
  /// inhertied form basefeature, calculate the binary size of an class instance
  virtual const unsigned long int calcBinarySize() const;
  
  
  ///return the internal data structure (::std::vector) of this
  ///feature. This allows for various modifications.
  virtual ::std::vector<double> & data() {return data_;}
  virtual const ::std::vector<double> & data() const {return data_;}
  
};

#endif
