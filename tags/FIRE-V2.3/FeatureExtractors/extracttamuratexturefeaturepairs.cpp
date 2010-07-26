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
#include "getpot.hpp"
#include "gzstream.hpp"
#include "diag.hpp"
#include "imagefeature.hpp"
#include "histogrampairfeature.hpp"
#include "histogramfeature.hpp"
#include "imagelib.hpp"
#include "tamurafeature.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "extracttemplate [options] (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "    -h, --help   show this help" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  //command line parsing
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();
  string suffix=cl.follow("tamura.histopairs.gz","--suffix");
  
  vector<double> sizes;
  double s=cl.follow(1.0,"--sizes");
  while(s>0) {
    sizes.push_back(s);
    s=cl.next(-10.0);
  }
  
  //get list of files to be processed
  vector<string> infiles;
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
  

  // processing the files
  for(uint i=0;i<infiles.size();++i) {
    string filename=infiles[i];
    DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;

    
    ImageFeature im;
    im.load(filename);

    uint ximgsize=im.xsize();
    uint yimgsize=im.ysize();
    
    ImageFeature tamuraImage=calculate(im);
    normalize(tamuraImage);
    
    vector<uint> binImage=histogramImageBinImage(tamuraImage);
    HistogramFeature completeHisto=histogramize(tamuraImage,binImage,0,0,ximgsize,yimgsize);

    HistogramPairFeature histopair(completeHisto);
    
    uint centerx=ximgsize/2; uint centery=yimgsize/2;
    for(uint i=0;i<sizes.size();++i) {
      double s=sizes[i];
      uint startx, stopx, starty, stopy;

      if(s<1.0) {
        // get half window size
        uint sizex2=uint(double(ximgsize)*s/2.0); uint sizey2=uint(double(yimgsize)*s/2.0);
        
        // get the start and stop position for histogramization
        startx=centerx-sizex2; starty=centery-sizey2;
        stopx=centerx+sizex2; stopy=centery+sizey2; 
      } else {
        startx=0;starty=0; 
        stopx=ximgsize; stopy=yimgsize;
      }
        
      HistogramFeature partHisto=histogramize(tamuraImage,binImage,startx,starty,stopx,stopy);
      histopair.push_back(partHisto);
    }

    histopair.save(filename+"."+suffix);

    DBG(20) << "Finished with '" << filename << "'." << endl;
  }
  
  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
