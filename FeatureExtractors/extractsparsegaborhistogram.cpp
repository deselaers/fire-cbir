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
#include "createsparsehisto.hpp"
#include "getpot.hpp"
#include "gabor.hpp"
#include "gzstream.hpp"
#include "imagelib.hpp"
#include "jflib.hpp"
#include <vector>
#include <iostream>

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "createsparsegaborhisto [options] (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "    -s  --suffix           suffix of sparse histogram files" << endl      
       << "    -ph, --phases          number of phases for the gabor feature extraction" << endl
       << "                            default value is 2" << endl
       << "    -fr, --frequencies     number of frequences for the gabor feature extraction" << endl
       << "                            default value is 2" << endl
       << "    -st, --steps           number of steps for each dimension in the histogram" << endl
       << "                           default value is 4" << endl
       << "    -h, --help   show this help" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  const int STEP_SIZE_DEFAULT = 4;
  const int NUM_PHASES_DEFAULT  = 2;
  const int NUM_FREQUENCIES_DEFAULT = 2;
  int numPhases = NUM_PHASES_DEFAULT, numFrequencies = NUM_FREQUENCIES_DEFAULT, steps = STEP_SIZE_DEFAULT;

  //command line parsing
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();
  if (cl.search(2, "-ph", "--phases")) {
    numPhases = cl.next(NUM_PHASES_DEFAULT);
  }
  if (cl.search(2, "-fr", "--frequencies")) {
    numFrequencies = cl.next(NUM_FREQUENCIES_DEFAULT);
  }
  if (cl.search(2, "-st", "--stepsize")) {
    steps = cl.next(STEP_SIZE_DEFAULT);
  }
    
  string suffix = ".gabor.sparsehisto.gz";
  if (cl.search(2, "--suffix", "-s")) {
    suffix = cl.next(" ");
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

    CreateSparseHisto csh(steps);

    ImageFeature img;
    img.load(filename);
    
    Gabor gabor(img);
    DBG(10) << "calculating gabor features..." << endl;
    gabor.calculate(numPhases, numFrequencies);
    
    DBG(10) << "adding gabor features to histogram..." << endl;
    normalize(gabor);
    for (int i = 0; i < numPhases * numFrequencies; i++) {
      ImageFeature gaborImage = gabor.getImage(i);
      normalize(gaborImage);
      csh.addLayers(gaborImage);
    }
    
    string output = filename + suffix;
    DBG(10) << "saving histogram to file '" << output  << "'..." << endl;
    csh.write(output);
    
    DBG(20) << "Finished with '" << filename << "'." << endl;
  }

  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
