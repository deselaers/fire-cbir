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
#include "invariantfeaturehistogram.hpp"

using namespace std;

void USAGE() {
  cout  << "USAGE:" << endl
        << "extracttemplate [options] (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
        << "   Options:" << endl
        << "    -h, --help   show this help" << endl
        << "    --suffix <suffix>" << endl
        << "          default: pairs.inv.histo.gz" << endl
        << "    --sizes <list of sizes>" << endl
        << "          default: 1.0" << endl
        << "    --kernelfunction <description> " << endl
        << "          description of the kernelfunction" << endl
        << "          default: mon2:x1=4:y1=0:x2=0:y2=8" << endl
        << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  //command line parsing
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();
  string suffix=cl.follow("inv.histopair.gz","--suffix");
  
  bool gray=cl.search("--gray");
  
  vector<double> sizes;
  double s=cl.follow(1.0,"--sizes");
  while(s>0) {
    DBG(30) << "reading size: " << s << endl;
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
  
  // reading the R
  int R = 36;
  if (cl.search("--R")) {
    R = cl.follow(R,"-R");
  }

  // reading of the number of steps
  int steps = 8;
  if (cl.search("--steps")) {
    steps = cl.follow(steps,"--steps");
  }		
	
  // reading of the number of samples
  int samples = -1;
  if (cl.search("--samples")) {
    samples = cl.follow(samples, "--samples");	
  }
	
  // reading of the kernelfunction's description
  string kernelfunction = "mon2:x1=4:y1=0:x2=0:y2=8";
  if (cl.search("--kernelfunction")) {
    kernelfunction = cl.next(" ");
  }
  
  
  
  // processing the files
  for(uint i=0;i<infiles.size();++i) {
    string filename=infiles[i];
    DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;

    
    ImageFeature img;
    img.load(filename);

    uint ximgsize=img.xsize();
    uint yimgsize=img.ysize();
    
    if(gray) {
      img=makeGray(img);
    }
    
    // Histogram berechnen
    HistogramFeature completeHisto = getInvariantFeatureHistogram(img, kernelfunction, samples, steps, R);
 
    HistogramPairFeature histopair(completeHisto);
    
    // iterate over all given sizes
    uint centerx=ximgsize/2; uint centery=yimgsize/2;
    
    DBG(30) << "sizes.size():" << sizes.size() << endl;
    
    for(uint i=0;i<sizes.size();++i) {
      double s=sizes[i];
      uint startx, stopx, starty, stopy;

      DBG(30) << "size:" << s << endl;
      
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
        
      DBG(30) << "startx: " << startx << " stopx: " << stopx << " starty: " << starty << " stopy: " << stopy << endl;
      
      // calculate histogram of the part of the image
      HistogramFeature partHisto= getInvariantFeatureHistogramPart(img, kernelfunction, samples, steps, R, startx, stopx, starty, stopy);
      histopair.push_back(partHisto);
    }

    histopair.save(filename+"."+suffix);

    DBG(20) << "Finished with '" << filename << "'." << endl;
  }
  
  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
