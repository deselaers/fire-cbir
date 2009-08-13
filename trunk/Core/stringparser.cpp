#include "stringparser.hpp"
#include "diag.hpp"
#include <sstream>
using namespace std;

/// function returning the integer which follows the string t in string s.
/// examples:
///    getIntAfter("egal=7asdf","egal=",1) will return 7
///    getIntAfter("egal=7asdf","doof=",1) will return 1, as the string doof= is not found

const int getIntAfter(const string& s, const string& t, const int def) {
  int val=def;

  uint pos=s.find(t);
  string INTstring="";
  if(pos<s.size()) {
    for (uint i=pos+t.size(); s[i]>='0' && s[i]<='9' && i<s.size() ; ++i) {
      INTstring+=s[i];
    }
    istringstream istr(INTstring);
    istr >> val;
  }
  return val;
}

/// function returning the double which follows the string t in string s.
/// examples:
///    getDoubleAfter("egal=7asdf","egal=",1) will return 7
///    getDoubleAfter("egal=7asdf","doof=",1) will return 1, as the string doof= is not found

const double getDoubleAfter(const string& s, const string& t, const double def) {
  double val=def;

  uint pos=s.find(t);
  string DOUBLEstring="";
  if (pos<s.size()) {
    for (uint i=pos+t.size(); (s[i]>='0' && s[i]<='9' || s[i]=='.') && i<s.size(); ++i) {
      DOUBLEstring+=s[i];
    }
    istringstream istr(DOUBLEstring);
    istr >> val;
  }
  return val;
}

/// function returning the string which follows the string t in string s.
/// delimiter can be set
/// examples:
///    getDoubleAfter("egal=7asdf","egal=",1) will return 7
///    getDoubleAfter("egal=7asdf","doof=",1) will return 1, as the string doof= is not found

const string getStringAfter(const string& s, const string& t, const string& def, const char delimiter) {
  string val=def;

  uint pos=s.find(t);
  string STRINGstring="";
  if (pos<s.size()) {
    for (uint i=pos+t.size(); s[i]!=delimiter && i<s.size() ; ++i) {
      STRINGstring+=s[i];
    }
    val=STRINGstring;
  }
  return val;
}

///function returning whether a certain string is part of another string
const bool getBooleanString(const string& s, const string &t) {
  return s.find(t)<s.size();
}
