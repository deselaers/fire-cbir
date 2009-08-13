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
#include "lfsignaturefeature.hpp"
#include "em.hpp"
#include "distancemaker.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "lf2lfsignature [options] (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "    -h, --help   show this help" << endl
       << "    -s, --suffix suffix of output files" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  //command line parsing
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();
   
  //get list of files to be processed
  vector<string> infiles;
  
  string suffix="lfsig.gz"; 
  if (cl.search(2, "--suffix", "-s")) {
    suffix = cl.next(" ");
  }
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


  string splitMode = cl.follow("allSplit","-splitMode");
  string disturbMode = cl.follow("meanDisturb","-disturbMode");
  string poolMode = cl.follow("bothPooling","-poolMode");
  uint maxSplit = cl.follow(5,"-maxSplit");
  uint stopWithNClusters = cl.follow(128,"-stopWithNClusters");
  uint dontSplitBelow = cl.follow(8,"-dontSplitBelow");
  uint iterationsBetweenSplits = cl.follow(10,"-iterationsBetweenSplits");
  uint minObservationsPerCluster = cl.follow(5,"-minObservationsPerCluster");
  double epsilon = cl.follow(0.1,"-epsilon");
  string distanceMaker = cl.follow("","-dist");
  //setting up clustering
  EM em;
  em.splitMode(splitMode);
  em.disturbMode(disturbMode);
  em.poolMode(poolMode);
  em.maxSplit()=maxSplit;
  em.stopWithNClusters()=stopWithNClusters;
  em.dontSplitBelow()=dontSplitBelow;
  em.iterationsBetweenSplits()=iterationsBetweenSplits;
  em.minObservationsPerCluster()=minObservationsPerCluster;
  em.epsilon()=epsilon;
  if(distanceMaker!="") {
    DistanceMaker dm;
    em.setDist(dm.makeDistance(distanceMaker));
    DBG(10) << "making distance " << distanceMaker << endl;
  }  else {
    DBG(10) << "not making any distance" << endl;
  }
  em.printConfig();

  

  // processing the files
  for(uint i=0;i<infiles.size();++i) {
    string filename=infiles[i];
    DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;
    
    LFSignatureFeature lfsig;
    LocalFeatures lf; lf.load(filename);
    vector< GaussianDensity > & model=lfsig.model();
    vector< uint > & histo=lfsig.signature();
    ResultVector clusterinformation;
    
    DoubleVectorVector vec(0);
    for(uint l=0;l<lf.size();++l) {
      vector<double> * tmpvec=new vector<double>(lf[l]);
      vec.push_back(tmpvec);
    }
    em.clear();
    em.run(vec,clusterinformation);

    //    for(uint l=0;l<lf.size();++l) {
    //  delete vec[i];
    //}

    for(int c=0;c<em.numberOfClusters();++c) {
      model.push_back(em.center(c));
    } 
    histo.resize(model.size());
    for(uint l=0;l<clusterinformation.size();++l) {
      if(clusterinformation[l]>0)  // there might be some observations that could not be classified in the last reestimation iteration
        histo[clusterinformation[l]]++;
    }
    BLINK(10) << endl;
    
    lfsig.save(filename+"."+suffix);
    
    DBG(20) << "Finished with '" << filename << "'." << endl;
  }
 

  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
