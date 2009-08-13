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
#include "imagefeature.hpp"
#include "sparsehistogramfeature.hpp"
#include "getpot.hpp"
#include "gzstream.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "extractsparesecolorhistogram [options] (--gray|--color) (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "      --suffix <suffix>  set the suffix of the output file" << endl
       << "      --steps <steps>    set the steps per dimension of the histogram" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv) {
  
  GetPot cl(argc,argv);
  
  if(cl.search(2,"--help","-h")) USAGE();
  vector<string> infiles;

  string suffix=cl.follow("sparsehisto.gz","--suffix");

  uint steps=cl.follow(8,"--steps");

  if(cl.search("--images")) {
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
  
  if(cl.search("--gray")) {

    string filename;
    for (uint i = 0; i < infiles.size(); i++) {
      ImageFeature im;
      filename = infiles[i];
      im.load(filename);

      SparseHistogramFeature result(1, steps);
      result.min()=vector<double>(1,0.0);
      result.max()=vector<double>(1,1.0);
      
      result.initStepsize();
      for(uint x=0;x<im.xsize();++x) {
	for(uint y=0;y<im.ysize();++y) {
	  result.feed(vector<double>(1,im(x,y,0)));
	}
      }
    
      result.save(filename+"."+suffix);
    }
    
  } else if(cl.search("--color")) {

    string filename;
    for (uint i = 0; i < infiles.size(); i++) {
      ImageFeature im;
      filename = infiles[i];
      im.load(filename);
      SparseHistogramFeature result(vector<uint>(3,steps));
      result.min()=vector<double>(3,0.0);
      result.max()=vector<double>(3,1.0);
      result.initStepsize();
      
      vector<double> tofeed(3);

      for(uint x=0;x<im.xsize();++x) {
	for(uint y=0;y<im.ysize();++y) {
	  tofeed[0]=im(x,y,0);
	  tofeed[1]=im(x,y,1);
	  tofeed[2]=im(x,y,2);
	  result.feed(tofeed);
	}
      }
      result.save(filename+"."+suffix);
    }
    
  } else {
    USAGE();
  }
  

}
