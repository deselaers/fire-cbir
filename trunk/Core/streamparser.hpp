#ifndef __streamparser_hpp__
#define __streamparser_hpp__

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "streamparser.hpp"
#include "gzstream.hpp"
#include "diag.hpp"

class StreamParser {

public:
  StreamParser(::std::istream &is) : is_(is) {}
  
  bool getLine(::std::string &line) {
    ::std::getline(is_,line);
    DBG(50) << "read: " <<  VAR(line) << ::std::endl;
    return true;
  }

  template<class T>
  bool getFromLine(const ::std::string& keyword, T& b) {
    ::std::string line, key;
    ::std::getline(is_,line);
    ::std::istringstream iss(line);
    iss >> key;
    if(keyword!=key) {
      ERR << "Unable to parse stream. Expecting keyword '" << keyword <<"' but got '" << key <<"'." << ::std::endl;
      return false;
    } else {
      iss >> b;
      DBG(50) << VAR(line) << " " <<VAR(keyword) << " " << VAR(key) << " " << VAR(b) << ::std::endl;
      return true;
    }
  }

  template<class T>
  bool getVectorFromLine(const ::std::string& keyword, const int dim, ::std::vector<T>& v) {
    ::std::string line, key;
    ::std::getline(is_,line);
    ::std::istringstream iss(line);
    iss >> key;
    if(keyword!=key) {
      ERR << "Unable to parse stream. Expecting keyword '" << keyword <<"' but got '" << key <<"'." << ::std::endl;
      return false;
    } else {
      v.resize(dim);
      for(int i=0;i<dim;++i) {
        iss >> v[i];
      }
      DBG(50) << VAR(line) << " " <<VAR(keyword) << " " << VAR(key) << " " << VAR(v) << ::std::endl;
      return true;
    }
  }
  
private:
  ::std::istream &is_;
  
};

#endif
