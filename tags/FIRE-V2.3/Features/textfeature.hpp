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
#ifndef __textfeature_hpp__
#define __textfeature_hpp__
#include "diag.hpp"
#include <string>
#include <fstream>

class TextFeature: public BaseFeature {
private:
  ::std::string textfilename_;

public:
  TextFeature() : textfilename_(""){
    type_ = FT_TEXT;
  }

  TextFeature(::std::string val) : textfilename_(val) {
    type_ = FT_TEXT;
  }
    
  virtual ~TextFeature() {}

  virtual TextFeature* clone() const {return new TextFeature(*this);}
  bool read(::std::istream &is) {
    ::std::string line;
    ::std::getline(is,line);
    textfilename_ = line;
    if(textfilename_=="") return false;
    return true;
  }
 
  void write(::std::ostream &os) {
    os << textfilename_ << ::std::endl;
  }
  
  const ::std::string& value() const {return textfilename_;}

};

#endif
