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
#ifndef __logfile_hpp__
#define __logfile_hpp__

#include <fstream>
#include <iostream>
#include <time.h>

class LogFile {
private:
  ::std::string logfile_;

public:
  LogFile(::std::string filename="") {
    logfile_=filename;
  }

  void log(::std::string logstring) {
    if(logfile_ != "") {
      ::std::ofstream os(logfile_.c_str(),::std::ios::out|::std::ios::app);
      os << getTime() << " " << logstring << ::std::endl;
      os.close();
    }
  }
  
  ::std::string getTime() {
    char buffer[40];
    struct tm *TM;
    time_t t;
    time(&t);
    TM=localtime(&t);
    strftime(buffer,40,"%Y-%m-%d %H:%M:%S",TM);
    ::std::string result(buffer);
    return result;
  }
};


#endif
