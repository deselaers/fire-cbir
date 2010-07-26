// 13.09.2005 :: CT :: modifications for Sun C++ compiler

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
#ifndef __diag_hpp
#define __diag_hpp
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdlib>


typedef ::std::vector<double> DoubleVector;
typedef ::std::vector< DoubleVector* > DoubleVectorVector;
typedef ::std::vector<int> IntVector;
typedef unsigned int uint;

//#include "colors.hpp"

/** Library for diagnostic output.
    
    @author Nils Springob <nils.springob@crazy-idea.de>, 
            Thomas Deselaers <thomas@deselaers.de>
    
    This is a quite efficient and easy way to have some debug messages
    in your code. It is possible to disable this while run or while
    compile time in as fine steps as you like.

    -# while compiletime: to leave all the messages with a debug-level
       < x out of the programm compile with -DDEBUG_LEVEL=x

    -# while runtime: set the environment variable to the value you
       prefer.

    -# also it is possible to specify this directly in your code using
       the CheckRuntimeDebugLevel function

 */
#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif

#include <iostream>

///simple debug-output macro

#define DBGI(level,instruction) \
 if(DEBUG_LEVEL>=level && CheckRuntimeDebugLevel(level)) instruction
#ifdef _SUNPRO_CC
#define DBG(level) \
 if(DEBUG_LEVEL>=level && CheckRuntimeDebugLevel(level)) ::std::cout << "(" << ((int)(level)) << ") [" << __FILE__<<":"<<__LINE__<<"] "
#else
#define DBG(level) \
 if(DEBUG_LEVEL>=level && CheckRuntimeDebugLevel(level)) ::std::cout << "(" << ((int)(level)) << ") [" << __FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<"] "
#endif
#define BLINK(level) \
 if(DEBUG_LEVEL>=level && CheckRuntimeDebugLevel(level)) ::std::cout 

///even simpler error-message macro
#ifdef _SUNPRO_CC
#define ERR ::std::cerr << "[" << __FILE__<<":"<<__LINE__<<"] ERROR: "
#else
#define ERR ::std::cerr << "[" << __FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<"] ERROR: "
#endif

// print variable name and value for use with DBG() etc
#define VAR(x)  #x " = " << x 

///function to specify and check the desired debug level
bool CheckRuntimeDebugLevel(int level, bool set=false);

::std::string GetCurrentWorkingDirectory();


void printCmdline(uint argc, char **argv);

///print a stack-trace. works in linux only
void stackTrace(int cutoff=0);


::std::string& string_tolower(::std::string &s);
::std::string& string_toupper(::std::string &s);

#endif
