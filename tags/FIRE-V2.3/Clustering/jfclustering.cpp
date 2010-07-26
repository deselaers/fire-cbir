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
  
#include <iostream>
#include <string>
#include "baseclusterer.hpp"
#include "em.hpp"
#include "dbscan.hpp"
#include "getpot.hpp"
#include "jflib.hpp"
#include "distancemaker.hpp"

using namespace std;

void USAGE() {
  
  cout << "USAGE: jfclustering -jf <joergfile to be clustered>" << endl
       << " General parameters:" << endl
       << "    -jf <filename of joerfile to be clustered>" << endl
       << "    -c <clusterer> <clusterer in dbscan, lbg" << endl
       << "    -dist <distance name, if not specified mahalanobis>" << endl
       << "    -addCluster <offset> to specify an offset for the output" << endl
       << " EM/LBG specific parameters:" << endl
       << "    -splitMode <splitmodename>" <<endl
       << "    -poolMode <poolmodename>" << endl
       << "    -maxSplit <maximal number of splits>" << endl
       << "    -stopWithNClusters <maximal number of clusters>" << endl
       << "    -dontSplitBelow <minimal number of observations in a cluster to be split>" << endl
       << "    -iterationisBetweenSplits <number of reestimate iterations after each split>" << endl
       << "    -minObservationPerCluster <minimal number of observation in a cluster not to be purged>" << endl
       << "    -epsilon <epsilon to disturb means>" << endl
       << " DBSCAN specific parameters:" << endl
       << "    -minPts <MinPts>" << endl
       << "    -epsilon <epsilon>" << endl
       << endl;
}

int main(int argc, char **argv) {
  GetPot cl(argc, argv);
  
  if(cl.search("-h") || !cl.search("-jf")) USAGE();
  
  string joergfilename=cl.follow("data.jf","-jf");
  int offset=cl.follow(0,"-addCluster");
  
  string clusterername=cl.follow("lbg","-c");
  BaseClusterer * clusterer=NULL;
  if("lbg"==clusterername) {
    clusterer=new EM();
    dynamic_cast<EM*>(clusterer)->splitMode(cl.follow("allSplit","-splitMode"));
    dynamic_cast<EM*>(clusterer)->disturbMode(cl.follow("meanDisturb","-disturbMode"));
    dynamic_cast<EM*>(clusterer)->poolMode(cl.follow("bothPooling","-poolMode"));
    dynamic_cast<EM*>(clusterer)->maxSplit()=cl.follow(10,"-maxSplit");
    dynamic_cast<EM*>(clusterer)->stopWithNClusters()=cl.follow(128,"-stopWithNClusters");
    dynamic_cast<EM*>(clusterer)->dontSplitBelow()=cl.follow(8,"-dontSplitBelow");
    dynamic_cast<EM*>(clusterer)->iterationsBetweenSplits()=cl.follow(10,"-iterationsBetweenSplits");
    dynamic_cast<EM*>(clusterer)->minObservationsPerCluster()=cl.follow(5,"-minObservationsPerCluster");
    dynamic_cast<EM*>(clusterer)->epsilon()=cl.follow(0.1,"-epsilon");
    if(cl.search("-dist")) {
      DistanceMaker dm;
      dynamic_cast<EM*>(clusterer)->setDist(dm.makeDistance(cl.follow("l1","-dist")));
    }
    dynamic_cast<EM*>(clusterer)->printConfig();
  } else if (clusterername=="dbscan") {
    clusterer=new DBSCAN();
    dynamic_cast<DBSCAN*>(clusterer)->epsilon()=cl.follow(0.1,"-epsilon");
    dynamic_cast<DBSCAN*>(clusterer)->minPts()=cl.follow(3,"-minPts");
    if(!cl.search("-dist")) {
      DistanceMaker dm;
      dynamic_cast<DBSCAN*>(clusterer)->setDist(dm.makeDistance("euclidean"));
    } else {
      DistanceMaker dm;
      dynamic_cast<DBSCAN*>(clusterer)->setDist(dm.makeDistance(cl.follow("l1","-dist")));
    }
  } else {
    ERR << "Unknown clusterer: " << cl.follow("lbg","-c") << endl;
    USAGE();
    exit(20);
  }
  
  pair<DoubleVectorVector,IntVector> JF=readJF(joergfilename);
  cout << "Read " << JF.first.size() << " entries from JF " << joergfilename << "." << endl;
  
  ResultVector clusterinformation;
  


  clusterer->run(JF.first,clusterinformation);
  for(uint i=0;i< clusterinformation.size();++i) {
    cout << clusterinformation[i] + offset;
    for(uint j=0;j<JF.first[i]->size();++j) {
      cout << " " << JF.first[i]->operator[](j);
    }
    cout << endl;
  }
}
