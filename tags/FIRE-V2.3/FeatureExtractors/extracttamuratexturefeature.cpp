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
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/
#include <string>
#include <vector>
#include "getpot.hpp"
#include "gzstream.hpp"
#include "diag.hpp"

#include "imagefeature.hpp"
#include "histogramfeature.hpp"
#include "tamurafeature.hpp"
#include "imagelib.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "extracttamurafeature [options] (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "      -h, --help   show this help" << endl
       << "      --saveHistogram [<suffix>]  to save the texture histogram, suffix default to tamura.histo.gz" << endl
       << "      --saveTextureImage [<suffix>]  to save the image, suffix defaults to tamura.png" << endl
       << "      --savePartialTextureImages to save the single layers of the texture image" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  //command line parsing
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();

  bool saveHisto=false, saveImage=false, savePartialImage=false;
  string histosuffix, imgsuffix;
  if(cl.search("--saveHistogram")) {
    saveHisto=true;
    histosuffix=cl.follow("tamura.histo.gz","--saveHistogram");
  }
  
  if(cl.search("--saveTextureImage")) {
    saveImage=true;
    imgsuffix=cl.follow("tamura.png","--saveTextureImage");
  }
  
  if(cl.search("--savePartialTextureImages")) {
    savePartialImage=true;
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
  ImageFeature im, tamuraImage;
  for(uint i=0;i<infiles.size();++i) {
    string filename=infiles[i];
    DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;
    im.load(filename);
    ImageFeature tamuraImage=calculate(im);
    normalize(tamuraImage);

    if(saveHisto) {
      HistogramFeature h=histogramize(tamuraImage);
      h.save(filename+"."+histosuffix);
    }
    
    if(saveImage) {
      normalize(tamuraImage);
      tamuraImage.save(filename+"."+imgsuffix);
    }
    
    if(savePartialImage) {
      normalize(tamuraImage);
      saveLayer(0,tamuraImage,filename+".directionality.png");
      saveLayer(1,tamuraImage,filename+".coarseness.png");
      saveLayer(2,tamuraImage,filename+".contrast.png");
    }
    
    DBG(20) << "Finished with '" << filename << "'." << endl;
  }

  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
