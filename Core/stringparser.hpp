#ifndef __stringparser_hpp__
#define __stringparser_hpp__
#include <string>

const int getIntAfter(const std::string& s, const std::string& t, const int def);
const double getDoubleAfter(const std::string& s, const std::string& t, const double def);
const std::string getStringAfter(const std::string& s, const std::string& t, const std::string& def, const char delimiter=':');
const bool getBooleanString(const std::string& s, const std::string &t); 

#endif
