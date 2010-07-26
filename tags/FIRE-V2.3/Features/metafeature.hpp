// 14.09.2005 :: CT :: modifications for Sun C++ compiler

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
#ifndef __metafeature_hpp__
#define __metafeature_hpp__
#include "diag.hpp"
#include <string>
#include <map>
#include <fstream>
#include <string.h>
#include <strings.h>

class MetaFeature: public BaseFeature {
private:
  std::string value_;
  std::map<std::string,std::string> values_;

public:
  MetaFeature() {
    type_ = FT_META;
  }

  MetaFeature(std::map<std::string,std::string> val):values_(val) {
    type_ = FT_META;
  }
  

    
  virtual ~MetaFeature() {}
  // This generates a map<str,str> of the form: map[key] == value from a text
  // file that has the form
  // key1: value1
  // key2: value2
  // ...
  // Bugs: No ':'s are allowed in key or value as they are the delimiters!
  
  virtual MetaFeature* clone() const {return new MetaFeature(*this);}
  
  bool read(::std::istream &is) {
	::std::string line;
	::std::getline(is,line);
    bool returnvalue=false;
    while(!is.eof()) {
		char linec[160];
		char *namepos, *valpos;
		
		strcpy(linec, line.c_str());
		
		if(( namepos = strtok(linec, ":\n") )) {
          valpos = strtok(NULL, ":\n");
          valpos++; // Kill the space
          values_[namepos] = valpos;
          //std::cout << "name: '" << namepos << "'" << std::endl << "value: '" << valpos << "'" << std::endl;
          returnvalue=true;
		} else {
          DBG(15) << "Invalid format in metafeature description!" << std::endl;
          return false;
		}
        value_ += line;
        if(!is.eof()) {
            value_ += "\n";
        }
		::std::getline(is,line);
    }
	//std::cout << std::endl;
    DBG(100) << value_ << std::endl;
    return returnvalue;
  }
  
  void write(::std::ostream &os) {
    os << value_;
  }
  
  const std::string& value() const {return value_;}
  const std::map<std::string,std::string>& values() const {return values_;}
};


#endif
