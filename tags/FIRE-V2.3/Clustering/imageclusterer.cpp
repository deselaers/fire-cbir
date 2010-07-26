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
#include "imagelib.hpp"
#include "database.hpp"
#include "imagecontainer.hpp"
#include "em.hpp"
#include "baseclusterer.hpp"
#include "distancemaker.hpp"
#include "getscoring.hpp"

using namespace std;

void USAGE() {
   cout << "USAGE: imageclusterer - cluster images" << endl
        << "imageclusterer <options> --filelist <list of localfeature files>" << endl
       << "  Available options are: " << endl
       << "    -splitMode <splitmodename> in {allSplit,largestSplit,varianceSplit}" <<endl
       << "    -poolMode <poolmodename> in {noPooling, clusterPooling, dimensionPooling, bothPooling}" << endl
       << "    -disturbMode <distmodename> in {constDisturb, meanDisturb, meanDisturb2, varianceDisturb}" << endl
       << "    -maxSplit <maximal number of splits> " << endl
       << "    -stopWithNClusters <maximal number of clusters>" << endl
       << "    -dontSplitBelow <minimal number of observations in a cluster to be split>" << endl
       << "    -iterationsBetweenSplits <number of reestimate iterations after each split>" << endl
       << "    -minObservationsPerCluster <minimal number of observation in a cluster not to be purged>" << endl
       << "    -epsilon <epsilon to disturb means>" << endl
       << "    -dist <distance name, if not specified mahalanobis>" << endl
       << "    --saveModel <filename> to save the model of the clustering to the specified file after clustering" << endl
       << endl;

  exit(20);
}

void setupClusterer(EM &clusterer, const string& splitMode, const uint maxSplit, const uint stopWithNClusters, const string& disturbMode, const string& poolMode,
                    const uint dontSplitBelow, const uint iterationsBetweenSplits, const uint minObservationsPerCluster, const double epsilon,
                    const string& distanceMaker) {
  clusterer.splitMode(splitMode);
  clusterer.disturbMode(disturbMode);
  clusterer.poolMode(poolMode);
  clusterer.maxSplit() = maxSplit;
  clusterer.stopWithNClusters() = stopWithNClusters;
  clusterer.dontSplitBelow() = dontSplitBelow;
  clusterer.iterationsBetweenSplits() = iterationsBetweenSplits;
  clusterer.minObservationsPerCluster() = minObservationsPerCluster;
  clusterer.epsilon() = epsilon;

  if(distanceMaker!="") {
    DistanceMaker dm;
    clusterer.setDist(dm.makeDistance(distanceMaker));
    DBG(10) << "making distance " << distanceMaker << endl;
  }

  else {
    DBG(10) << "not making any distance" << endl;
  }

  clusterer.printConfig();
}

int main(int argc, char** argv) {
  GetPot cl(argc, argv);
  
  if(!cl.search("--filelist") || cl.search("-h")) { USAGE(); exit(20);}
  
  string filename;
  if (cl.search("--filelist"))  filename=cl.follow("list","--filelist");
  cout << "filelist: " << filename << endl;
 
  string splitMode = cl.follow("allSplit","-splitMode");
  cout << "splitMode: " << splitMode << endl;
  string disturbMode = cl.follow("meanDisturb","-disturbMode");
  cout << "disturbMode: " << disturbMode << endl;
  string poolMode = cl.follow("bothPooling","-poolMode");
  cout << "poolMode: " << poolMode << endl;
  uint maxSplit = cl.follow(3,"-maxSplit");
  cout << "maxSplit: " << maxSplit << endl;
  uint stopWithNClusters = cl.follow(8,"-stopWithNClusters");
  cout << "stopWithNClusters: " << stopWithNClusters << endl;
  uint dontSplitBelow = cl.follow(8,"-dontSplitBelow");
  cout << "dontSplitBelow: " << dontSplitBelow << endl;
  uint iterationsBetweenSplits = cl.follow(10,"-iterationsBetweenSplits");
  cout << "iterationsBetweenSplits: " << iterationsBetweenSplits << endl;
  uint minObservationsPerCluster = cl.follow(5,"-minObservationsPerCluster");
  cout << "minObservationsPerCluster: " << minObservationsPerCluster << endl;
  double epsilon = cl.follow(0.1,"-epsilon");
  cout << "epsilon: " << epsilon << endl;
  string distanceMaker="";
  if(cl.search("-dist")) {
    distanceMaker = cl.follow("euclidean","-dist");
    cout << "distanceMaker: " << distanceMaker << endl;
  }
  string modelfile=cl.follow("","--saveModel");
  cout << "modelfile: " << modelfile << endl;
  
  // setting up clusterer
  ResultVector clusterinformation;
  EM em;      
  setupClusterer(em, splitMode, maxSplit, stopWithNClusters, disturbMode, poolMode, dontSplitBelow, 
                 iterationsBetweenSplits, minObservationsPerCluster, epsilon, distanceMaker);

  // reading in data                 
  Database db;
  db.loadFileList(filename);
  db.loadFeatures();
  DoubleVectorVector imageData;
  for (uint i = 0; i < db.size(); i++) {
    DoubleVector *vec=new DoubleVector();
    *vec = db[i]->asVector();
    imageData.push_back(vec);
  }
  
  DBG(10)  <<imageData.size() << "images represented by " << imageData[0]->size() << " dimensional vectors loaded." << endl;


  

  // let clusterer work
  if(modelfile!="") {
    em.run(imageData, clusterinformation, modelfile);
    em.saveModel(modelfile);
  } else {
    em.run(imageData, clusterinformation);
  }
  
  // printing result
  for (uint i = 0; i < clusterinformation.size(); i++) {
    cout << "clusterresult: " << clusterinformation[i] << " " << db[i]->basename() <<endl;
  }
  
  
  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
