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
  59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#include <stdio.h>
#include "dist_mpeg7.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <sstream>
#include "mpeg7feature.hpp"
#include "diag.hpp"

#ifndef _SUNPRO_CC
#include "pstream.hpp"
using namespace redi;
#endif
using namespace std;

double MPEG7Distance::distance(const BaseFeature*, const BaseFeature* databaseFeature) {
  double result=1.0e20; // return a high variance if image does
                            // not occur (probably something has gone
                            // wrong)
  string bn=dynamic_cast<const MPEG7Feature*>(databaseFeature)->basename();
  if(distances_.find(bn)!=distances_.end()) {
    result=distances_[bn];
  } else {
    ERR << "Probably something has gone wrong. Image " << bn << " has no associated MPEG7Distance." << endl;
  }
  return result;
}


void MPEG7Distance::start(const BaseFeature *queryFeature) {
  distances_.clear();
  vector<string> output;
  
  const MPEG7Feature *q=dynamic_cast< const MPEG7Feature * >(queryFeature); 
  if(!q) {
    ERR << "This is not a valid MPEG7Feature:" << queryFeature->type() << endl;
    return;
  }

  string mpegtype=q->getMPEG7Type();

  ostringstream cmdline;
  cmdline << xmmain_ 
          << " -p" << mpeg7data_<<".par"     //parfile
          << " -a" << mpegtype << "Client"   //application
          << " -l" << mpeg7data_ <<".lst"    //list
          << " -b" << mpeg7data_<<".mp7"     //bitstream
          << " -q" << q->basename()          //name of the query image
          << " -n 1000000000";               //number of results is sufficiently large to get all results :-)
  
  // does not work, is not ansi compliant
  //  FILE *file=popen(cmdline.str().c_str(),"r"); // run the program
  //  ifstream is(fileno(file));
  
#ifdef _SUNPRO_CC 
  FILE *file=popen(cmdline.str().c_str(),"r"); // run the program
  ifstream is(fileno(file));
#else
  ipstream is(cmdline.str());
#endif

  string line;
  /** currently the output of the XMMain looks like:
   *
   * MPEG-7 XM Version 6.1_alpha
   * Initializing Xerces Parser
   * Loading : /u/deselaers/tmp/555.jpg
   * /work/deselaers/imageretrieval/wang/images/555.jpg	0.000000
   * /work/deselaers/imageretrieval/wang/images/504.jpg	7.111817
   * /work/deselaers/imageretrieval/wang/images/506.jpg	7.232543
   * /work/deselaers/imageretrieval/wang/images/355.jpg	8.012952
   * /work/deselaers/imageretrieval/wang/images/953.jpg	8.014173
   * ...
   *
   * that is, we can parse the input from the 3rd line.
   */
  

  uint lines=0;
  getline(is,line);
  getline(is,line);
  getline(is,line);
  string filename; double dist;
  while(!is.eof()) {
    getline(is,line);
    istringstream iss(line);
    iss >> filename >> dist;

    distances_[filename]=dist;
    DBG(20) << VAR(filename) << " " << VAR(dist) << " " << VAR(distances_[filename]) << endl;
    ++lines;
  }
  
  DBG(10) << lines << " lines from output of " << cmdline.str() << " read." << endl;
  is.close();
}

void MPEG7Distance::stop() {
  distances_.clear();
}


