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


#include "clusterlocalfeatures.hpp"
#include "em.hpp"
#include "baseclusterer.hpp"
#include "distancemaker.hpp"

using namespace std;

void setupClusterer(EM* clusterer, string splitMode, uint maxSplit, uint stopWithNClusters, string disturbMode, string poolMode,
                    uint dontSplitBelow, uint iterationsBetweenSplits, uint minObservationsPerCluster, double epsilon,
                    string distanceMaker, bool saveBeforeSplits)
{
  clusterer->splitMode(splitMode);
  clusterer->disturbMode(disturbMode);
  clusterer->poolMode(poolMode);
  clusterer->maxSplit() = maxSplit;
  clusterer->stopWithNClusters() = stopWithNClusters;
  clusterer->dontSplitBelow() = dontSplitBelow;
  clusterer->iterationsBetweenSplits() = iterationsBetweenSplits;
  clusterer->minObservationsPerCluster() = minObservationsPerCluster;
  clusterer->epsilon() = epsilon;
  clusterer->saveBeforeSplits()=saveBeforeSplits;

  if(distanceMaker!="")
  {
    DistanceMaker dm;
    clusterer->setDist(dm.makeDistance(distanceMaker));
    DBG(10) << "making distance " << distanceMaker << endl;
  }

  else
  {
    DBG(10) << "not making any distance" << endl;
  }

  clusterer->printConfig();
}

ResultVector cluster(vector<LocalFeatures*> &locfeat, string splitMode, uint maxSplit, uint stopWithNClusters, string disturbMode,
                     string poolMode, uint dontSplitBelow, uint iterationsBetweenSplits, uint minObservationsPerCluster,
                     double epsilon, string distanceMaker, string saveModelTo, bool saveBeforeSplits,string loadModelFrom)
{
  return cluster(locfeat, false, splitMode, maxSplit, stopWithNClusters, disturbMode, poolMode, dontSplitBelow, iterationsBetweenSplits,
                 minObservationsPerCluster, epsilon, distanceMaker, saveModelTo, saveBeforeSplits,loadModelFrom);
}



ResultVector cluster(vector<LocalFeatures*> &locfeat, bool clusterWithPositions, string splitMode, uint maxSplit, uint stopWithNClusters, string disturbMode, string poolMode, uint dontSplitBelow, uint iterationsBetweenSplits, uint minObservationsPerCluster, double epsilon, string distanceMaker, string saveModelTo, bool saveBeforeSplits, string loadModelFrom)
{

  BaseClusterer *clusterer = new EM();

  if(loadModelFrom!="")
  {
    clusterer->loadModel(loadModelFrom);
  }
  setupClusterer(dynamic_cast<EM*>(clusterer), splitMode, maxSplit, stopWithNClusters, disturbMode, poolMode, dontSplitBelow,
                 iterationsBetweenSplits, minObservationsPerCluster, epsilon, distanceMaker, saveBeforeSplits);


  ResultVector clusterinformation;

  if (clusterWithPositions)
  {
    DBG(10) << "clustering with position information included" << endl;
  }

  DoubleVectorVector localFeaturesData;
  for (uint j = 0; j < locfeat.size(); j++)
  {
    LocalFeatures* localFeatures = locfeat[j];
    int xSize = localFeatures->imageSizeX();
    int ySize = localFeatures->imageSizeY();
    for (uint i = 0; i < localFeatures->numberOfFeatures(); i++)
    {
      DoubleVector* lfvector = &((*localFeatures)[i]);
      if (clusterWithPositions)
      {
        pair<double, double> pos = localFeatures->relativePosition(i);
        lfvector->push_back((double) localFeatures->position(i).x / (double) xSize);
        lfvector->push_back((double) localFeatures->position(i).y / (double) ySize);
      }
      localFeaturesData.push_back( lfvector );
    }
  }


  if(saveModelTo!="")
  {
    dynamic_cast<EM*>(clusterer)->run(localFeaturesData, clusterinformation, saveModelTo);
    clusterer->saveModel(saveModelTo);
  }
  else
  {
    clusterer->run(localFeaturesData, clusterinformation);
  }

  delete clusterer;
  return clusterinformation;
}
