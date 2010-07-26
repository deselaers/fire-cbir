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
#include "relationalfeaturehistogram.hpp"


using namespace std;

void USAGE() {
   cout << "USAGE:" << endl
        << "extractrelationalfeaturehistogram [options] (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
        << "   Options:" << endl
 	      << "    -h, --help   show this help" << endl
        << "    --suffix <suffix>" << endl
        << "          default: rel.histo.gz" << endl
        << "    --kernelfunction <description> " << endl
        << "          description of the kernelfunction" << endl
        << "          default: rel:x1=0:y1=0:x2=0:y2=4" << endl

        << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  string suffix=cl.follow("rel.histo.gz","--suffix");

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
    igzstream ifs(cl.follow("list","--filelist"));
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

	
  // reading of the kernelfunction's description
  string kernelfunction = "rel:x1=0:y1=0:x2=0:y2=4";
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
    
    // Histogram berechnen
    histo = getRelationalFeatureHistogram(img, kernelfunction);
 
    // Histogram speichern
    histo.save(filename+"."+suffix);

    
    DBG(20) << "Finished with '" << filename << "'." << endl;
  }
  
  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
