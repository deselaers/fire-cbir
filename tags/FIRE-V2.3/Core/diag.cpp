#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <string>
#include "diag.hpp"
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>

#ifdef OS_linux
#include <execinfo.h>
#endif
#include <unistd.h>


using namespace std;

// common function, where else to place ??
bool CheckRuntimeDebugLevel(int level, bool set) {
  static int s_Level = -1;
  if (set)  s_Level=level;
  if (s_Level == -1) {
    char * p = getenv("DEBUG_LEVEL");
    if (p) {
      string s(p);
      cout << "DEBUG_LEVEL is " << s << endl;
      istringstream is(s);
      is >> s_Level;
      cout << "DEBUG_LEVEL set to " << s_Level << endl;
    }
  }
  if (s_Level == -1) s_Level=99;
  return level<=s_Level;
}

string GetCurrentWorkingDirectory() {
  char cwd[1024];
  getcwd(cwd, 1024);
  return ::std::string(cwd);
}

void printCmdline(uint argc, char **argv) {
  cout << argv[0];
  for(uint i=1;i<argc;++i) {
    cout << " "  << argv[i];
  }
  cout << endl;
}


inline char my_tolower(char c) { return tolower(c); }
inline char my_toupper(char c) { return toupper(c); }

::std::string& string_tolower(::std::string &s) {
  ::std::transform(s.begin(),s.end(),s.begin(),my_tolower);
    return s;
}

::std::string& string_toupper(::std::string &s) {
  ::std::transform(s.begin(),s.end(),s.begin(),my_toupper);
  return s;
}


void stackTrace(int cutof) {
#ifdef OS_linux
  ERR << "stack trace (innermost first):" << std::endl;
  static const size_t maxTraces = 100;
  void *array[maxTraces];
  size_t nTraces = backtrace(array, maxTraces);
  ERR << nTraces << endl;
  char **strings = backtrace_symbols(array, nTraces);
  for (size_t i = cutoff+1; i < nTraces; i++)
	ERR << '#' << i << "  " << strings[i] << std::endl;
  free(strings);
#endif
}


