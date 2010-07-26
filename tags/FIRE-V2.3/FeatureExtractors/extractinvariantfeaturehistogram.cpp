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
#include "histogramfeature.hpp"
#include "invariantfeaturehistogram.hpp"
#include "imagelib.hpp"

using namespace std;

void USAGE() {
   cout << "USAGE:" << endl
        << "extractinvariantfeaturehistogram [options] (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
        << "   Options:" << endl
        << "    -h, --help   show this help" << endl
        << "    --gray   treat all images as gray images" << endl
        << "    --suffix <suffix>" << endl
        << "          default: inv.histo.gz" << endl
        << "    --steps <steps> " << endl
        << "          steps in histogram per dimension" << endl
        << "          default: 8" << endl
        << "    --samples <noofsamples> " << endl
        << "          number of samples used in the monte carlo method" << endl
        << "          default -1" << endl
        << "    --kernelfunction <description> " << endl
        << "          description of the kernelfunction" << endl
        << "          default: mon2:x1=4:y1=0:x2=0:y2=8" << endl
        << "    --R <noofsteps> " << endl
        << "          number of steps used for integration over the rotation group " << endl
        << "          default: 36" << endl
        << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  string suffix=cl.follow("inv.histo.gz","--suffix");

  bool gray=cl.search("--gray");

  //command line parsing
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();
  
  //get list of files to be processed
  vector<string> infiles;
  if(cl.search("--images")) {
    string filename=cl.next(" ");
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
  

  // initializing image and histogram
  ImageFeature img;
  HistogramFeature histo;
  
  // processing the files
  for(uint i=0;i<infiles.size();++i) {
    string filename=infiles[i];
    DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;

    // load Image
    img.load(filename);
    
    
    if(gray) {
      img=makeGray(img);
    }
    
    // Histogram berechnen
    histo = getInvariantFeatureHistogram(img, kernelfunction, samples, steps, R);
 
    // Histogram speichern
    histo.save(filename+"."+suffix);
    
    DBG(20) << "Finished with '" << filename << "'." << endl;
  }
  
  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
