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
#include "gabor.hpp"
#include "imagelib.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "extractgaborfeaturevector [options] (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "    -h, --help      show this help" << endl
       << "    --numPhases,-p  " << endl
       << "    --numFrequencies,-f" << endl
       << "    --suffix <suffix>" << endl
       << "    --hMargin" << endl
       << "    --vMargin" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  //command line parsing
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();

  uint numPhases=cl.follow(5,2,"--numPhases","-p");
  uint numFrequencies=cl.follow(3,2,"--numFrequencies","-f");
  string suffix=cl.follow("gabor.vec.gz","--suffix");
  int hMargin=cl.follow(32,"--hMargin");
  int vMargin=cl.follow(32,"--vMargin");
  bool saveImgs=cl.search("--saveImages");

  
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
  ImageFeature img;
  VectorFeature vecfeat;
  for(uint i=0;i<infiles.size();++i) {
    string filename=infiles[i];
    DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;
    img.load(filename);

    vecfeat=VectorFeature(numPhases*numFrequencies*2);

    DBG(20) << "Calculating Gabor features" << endl;
    Gabor gabor(img);
    gabor.calculate(numPhases, numFrequencies, hMargin, vMargin);
    DBG(20) << "Gabor feature calculated" << endl;

    DBG(20) << "Making GaborImage to vector" << endl;
    double mean, variance;
    for(uint i=0;i<numPhases*numFrequencies;++i) {
      ImageFeature gaborImage=gabor.getImage(i);
      normalize(gaborImage);
      if(saveImgs) {
        ostringstream iss; iss << filename << "-" << i << ".png";
        gaborImage.save(iss.str());
      }
      meanandvariance(gaborImage, mean, variance);
      vecfeat[2*i]=mean; vecfeat[2*i+1]=sqrt(variance);
    }
    
    vecfeat.save(filename+"."+suffix);
    DBG(20) << "Finished with '" << filename << "'." << endl;
  }
  

  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
