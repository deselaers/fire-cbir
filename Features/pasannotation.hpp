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
#ifndef __pascalannotationfeature_hpp__
#define __vectorfeature_hpp__

#include <string>
#include <sstream>
#include <iostream>
#include "basefeature.hpp"


struct PASObject {
  ::std::string label_;
  uint Xmin_,Ymin_,Xmax_,Ymax_;

  void resize(double fac) {
    Xmin_=uint(double(Xmin_)*fac);
    Xmax_=uint(double(Xmax_)*fac);
    Ymin_=uint(double(Ymin_)*fac);
    Ymax_=uint(double(Ymax_)*fac);
  }

  bool inside(int x, int y) const {
    return ( (x>=int(Xmin_)) && (x<=int(Xmax_)) && (y>=int(Ymin_)) && (y<=int(Ymax_)) );
  }

};

class PascalAnnotationFeature : public BaseFeature {
protected:
  /// image size
  uint X_,Y_,C_;
  
  /// number of object annotated in image
  uint nOfObjects_;

  /// image file name
  ::std::string filename_;

  ///object annotations
  ::std::vector<PASObject> objects_;
  
public:
 
  /// constructor, if no parameter (or 0) given -> defaults to false,
  /// else to true
  PascalAnnotationFeature()  {
    type_=FT_PASCALANNOTATION;
  }
  
  virtual PascalAnnotationFeature *clone() const {return new PascalAnnotationFeature(*this);}

  uint& X() {return X_;}
  uint& Y() {return Y_;}
  uint& C() {return C_;}
  uint& nOfObjects() {return nOfObjects_;}
  ::std::string filename() {return filename_;}
  PASObject& object(const uint i) {return objects_[i];}

  const uint& X() const {return X_;}
  const uint& Y() const {return Y_;}
  const uint& C() const {return C_;}
  const uint& nOfObjects() const {return nOfObjects_;}
  const ::std::string filename() const {return filename_;}
  const PASObject & object(const uint i) const {return objects_[i];}
  

  void resize(double fac) {
    X_=uint(double(X_)*fac); Y_=uint(double(Y_)*fac);
    for(std::vector<PASObject>::iterator i=objects_.begin();i!=objects_.end();++i) {
      i->resize(fac);
    }
  }
  
  /// derived from base feature, read from stream
  bool read(::std::istream &is) {
    ::std::string line,token;
    ::std::istringstream iss;

    getline(is,line);
    if(line!="PASAnnotation") {
      ERR << "reading. Expected 'PASAnnotation', got '" << token << "'.";
      return false;
    }
    
    token=""; getline(is,line); iss.clear(); iss.str(line); iss >> token;
    if(token=="filename") {
      iss >> filename_;
    } else {
      ERR << "reading. Expected 'filename', got '" << token << "'.";
      return false;
    }

    token=""; getline(is,line); iss.clear(); iss.str(line); iss >> token;
    if(token=="size") {
      iss >> X_ >> Y_ >> C_;
    } else {
      ERR << "reading. Expected 'size', got '" << token << "'."<< ::std::endl;
      return false;
    }
    
    token=""; getline(is,line); iss.clear(); iss.str(line); iss >> token;
    if(token=="objects") {
      iss >> nOfObjects_;
    } else {
      ERR << "reading. Expected 'size', got '" << token << "'."<< ::std::endl;
      return false;
    }

    for(uint i=0;i<nOfObjects_;++i) {
      PASObject tmp;

      token=""; getline(is,line); iss.clear(); iss.str(line); iss >> token;
      if(token=="PASObject") {
        iss >> tmp.label_ >> tmp.Xmin_ >> tmp.Ymin_ >> tmp.Xmax_ >> tmp.Ymax_ ;
      } else {
        ERR << "reading. Expected 'size', got '" << token << "'."<< ::std::endl;
        return false;
      }
      objects_.push_back(tmp);
    }
    return true;
  }
  
  /// derived from base feature, write to stream
  void write(::std::ostream &os) {
    os << "PASAnnotation" << ::std::endl
       << "filename " << filename_ << ::std::endl 
       << "size " << X_ << " " << Y_ << " " << C_ << ::std::endl
       << "objects " << nOfObjects_ << ::std::endl;
    for(uint i=0;i<nOfObjects_;++i) {
      os << "PASObject " << objects_[i].label_ << " " << objects_[i].Xmin_ << " " <<objects_[i].Ymin_ << " " << objects_[i].Xmax_ << " " <<objects_[i].Ymax_ << ::std::endl ;
    }
  }
  
};

#endif
