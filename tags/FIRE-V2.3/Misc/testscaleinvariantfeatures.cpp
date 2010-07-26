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
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include "getpot.hpp"
#include "imagelib.hpp"
#include "gzstream.hpp"
#include "differenceofgaussian.hpp"


using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "testscaleinvariantfeatures [options] (--color|--gray) (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "    -h, --help         show this help" << endl
       << "    -n, --num 'num'    extracts 'num' patches (default 200)" << endl;
}



int main(int argc, char** argv) {

  GetPot cl(argc,argv);
  if (cl.search(2,"--help","-h")) {
    USAGE();
  }

  vector<string> images;
  if (cl.search("--images")) {
    string filename =  cl.next(" ");
    while (filename != " ") {
      images.push_back(filename);
      filename = cl.next(" ");
    }
  } else if (cl.search("--filelist")) {
    string filename = "test";
    igzstream ifs; ifs.open(cl.follow("list","--filelist"));
    if(!ifs.good() || !ifs) {
      ERR << "Cannot open filelist " <<cl.follow("list","--filelist")  << ". Aborting." << endl;
      exit(20);
    }
    while(!ifs.eof() && filename!="") {
      getline(ifs,filename);
      if(filename!="") {
        images.push_back(filename);
      }
    }
    ifs.close();
  } else {
    USAGE();
    exit(20);
  }

  bool forceGray;
  if (cl.search(1, "--gray")) {
    forceGray = true;
  } else if (cl.search(1, "--color")) {
    forceGray = false;
  } else {
    USAGE();
    exit(20);
  }

  int numPatches = 200;
  if (cl.search(2, "-n", "--num")) {
    numPatches = cl.next(200);
  }

  for(uint i = 0; i < images.size(); i++) {
    string filename=images[i];
    DBG(10) << "Processing '"<< filename << "'.(" << i << "/" << images.size()<< ")" <<endl;
    ImageFeature img; img.load(filename,forceGray);
    DifferenceOfGaussian sift(img);
    vector<InterestPoint> interestPoints = sift.getInterestPoints(numPatches);
    
    ImageFeature padded = img;
    vector<double> color(3,1.0);
    
    for (int i = 0; i < (int) interestPoints.size(); i++) {
      InterestPoint ip = interestPoints[i];
      box(padded, ip.x, ip.y, color, ip.scale / 2);
    }
    padded.display();

    break;
  }


  DBG(10) << "cmdline was: "; printCmdline(argc,argv);

}


