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
#include <string>
#include <vector>
#include <limits>
#include "getpot.hpp"
#include "gzstream.hpp"
#include "diag.hpp"
#include "imagefeature.hpp"
#include "vectorfeature.hpp"
using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "extractaspectratio [options] (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "    -h, --help   show this help" << endl
       << "    -s, --suffix <suffix> default: ar.vec.gz" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  //command line parsing
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();
  
  string suffix=cl.follow(".ar.vec.gz","--suffix");

  //get list of files to be processed
  vector<string> infiles;
  if(cl.search("--image")) {
    string filename=cl.next(" ");;
    while(filename!=" ") {
      infiles.push_back(filename);
      filename=cl.next(" ");
    }
  } else if (cl.search("--filelist")) {
    string filename="test";
    igzstream ifs; ifs.open(cl.follow("list","--filelist"));
    if(!ifs.good() || !ifs) {
      ERR << "Cannot open filelist " <<cl.follow("list","--filelist")  << ". Aborting." << endl;
      exit(20);
    }
    while(!ifs.eof() && filename!="") {
      getline(ifs,filename);
      if(filename!="") {
        infiles.push_back(filename);
      }
    }
    ifs.close();
  } else {
    USAGE();
    exit(20);
  }
  

  // processing the files
  ImageFeature img;
  VectorFeature out(2);
  VectorFeature min(2), max(2), avg(2);
  min[0]=numeric_limits<double>::max();   min[1]=numeric_limits<double>::max(); 
    
  for(uint i=0;i<infiles.size();++i) {
    string filename=infiles[i];
    DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << "): ";
    img.load(filename);
    out[0]=img.xsize();
    out[1]=img.ysize();
    out.save(filename+"."+suffix);
    
    BLINK(10) << out[0] << " " << out[1] << endl;

    if(min[0]>out[0]) min[0]=out[0];
    if(min[1]>out[1]) min[1]=out[1];
    
    if(max[0]<out[0]) max[0]=out[0];
    if(max[1]<out[1]) max[1]=out[1];

    avg[0]+=out[0]; avg[1]+=out[1];
    DBG(20) << "Finished with '" << filename << "'." << endl;
  }
  
  avg[0]/=infiles.size(); avg[1]/=infiles.size();
  DBG(10) << "min: " << min[0] << " " << min[1] << endl;
  DBG(10) << "max: " << max[0] << " " << max[1] << endl;
  DBG(10) << "avg: " << avg[0] << " " << avg[1] << endl;

  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
