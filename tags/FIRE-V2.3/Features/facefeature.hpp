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

#ifndef __facefeature_hpp__
#define __facefeature_hpp__

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "basefeature.hpp"
#include "vectorfeature.hpp"
#include "vectorvectorfeature.hpp"


class Face :public VectorFeature {
protected: 
  uint posx_,posy_,width_,height_;
  uint size_;
  
public:

  /// default constructor
  Face(uint size=1024) : VectorFeature() ,posx_(0), posy_(0), width_(0), height_(0), size_(size){
  }
  

  virtual bool read(::std::istream &is) {
    ::std::string line,tmp;
    uint no;
    
    is >> tmp  >> no >> posx_ >> posy_ >> width_ >> height_ ;
    if(tmp != "face") {
      ERR << "This is not a face feature. Expected \"face\" got \"" << tmp << "\"." << ::std::endl;
      return false;
    }
    double v;
    for(uint i=0;i<size_;++i) {
      is >> v;
      data_.push_back(v);
    }

    DBG(20) << "Read face " << no << " for position ("<< posx_ << "," <<posy_ <<") of size (" << width_ << "x"<< height_ << ") with " << data_.size() << " coefficients." << ::std::endl;
    return true;
  }
  
  virtual void write(::std:: ostream &os) {
    os << posx_ << " " << posy_ << " " << width_ << " " << height_ ;
    for(uint i=0;i<data_.size();++i) {
      os << " " << data_[i];
    }
  }
    
  
  uint &posx() {return posx_;}
  uint &posy() {return posy_;}
  uint &width() {return width_;}
  uint &height() {return height_;}
  const uint &posx() const {return posx_;}
  const uint &posy() const {return posy_;}
  const uint &width() const {return width_;}
  const uint &height() const {return height_;}
};

class FaceFeature : public VectorVectorFeature {
protected:
  uint noffaces_, nofcoefficients_;
  ::std::vector< Face > data_;

public:
  FaceFeature() {
    type_=FT_FACEFEAT;
  }

  virtual FaceFeature* clone() const {return new FaceFeature(*this);}
  virtual const ::std::vector<double>& operator[](uint idx) const { return data_[idx].data();}
  virtual ::std::vector<double>& operator[](uint idx) { return data_[idx].data();}
  virtual const uint numberOfVectors() const {return data_.size();}
  
  virtual bool read(::std::istream& is) {
    ::std::string tmp;

    is>>tmp;
    if(tmp!="FIRE_facefeature") {
      ERR << "This is not a face feature file?!?. Expected \"FIRE_facefeature \", got \""<< tmp<<"\"."<< ::std::endl;
      return false;
    }
    is >> tmp >> noffaces_;
    is >> tmp >> nofcoefficients_;

    DBG(35) << tmp << " " << noffaces_ << ::std::endl;
    DBG(35) << tmp << " " << nofcoefficients_ << ::std::endl;
    
    Face ftmp;
    for(uint i=0;i<noffaces_;++i) {
      ftmp=Face(nofcoefficients_);
      ftmp.read(is);
      data_.push_back(ftmp);
      DBG(20) << "Got face no " << i << ::std::endl;
    }
    return true;
  }

  virtual void write(::std::ostream &os) {
    os << "FIRE_facefeature" << ::std::endl
       << "noffaces " << noffaces_ << ::std::endl
       << "nofcoefficients " << nofcoefficients_ << ::std::endl;
    for(uint i=0;i<noffaces_;++i) {
      os << "face " << i << " ";
      data_[i].write(os);
      os << ::std::endl;
    }
  }

  const Face& face(const uint i) const {return data_[i];}
  Face& face(const uint i) {return data_[i];}

  
  const uint nOfFaces() {return noffaces_;}
  const uint nOfCoefficients() const {return nofcoefficients_;}
  uint& nOfCoefficients() {return nofcoefficients_;}
  
};

#endif
