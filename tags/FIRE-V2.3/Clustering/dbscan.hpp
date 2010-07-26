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
#ifndef __dbscan_hpp__
#define __dbscan_hpp__

#include "baseclusterer.hpp"
#include "diag.hpp"
#include <vector>
#include <string>
#include <set>
#include "basedistance.hpp"

class DBSCAN : public BaseClusterer {
private:
  BaseDistance *dist_;
  DoubleVectorVector distanceTable_;
  
  double epsilon_;
  uint minPts_;
  
  // optimized distance function with lookup table
  double distance(const DoubleVectorVector &inputdata, const uint i, const uint j);
  

  // get neighborhood of a point (the indices in inputdata are returned in the set)
  ::std::set<uint> neighborhood(const DoubleVectorVector& inputdata, uint startObjectID);

  // expand a cluster with given startObject and given clusterId
  bool expandClusters(const DoubleVectorVector& inputdata, ResultVector& clusterInformation, const uint startObjectId, const uint clusterId);

public:
  static const int NOISE=-2;
  static const int UNCLASSIFIED;

  //methods derived from BaseClusterer
  DBSCAN();
  virtual void run(const DoubleVectorVector& inputdata, ResultVector &clusterInformation);
  virtual ~DBSCAN();
  virtual int classify(const DoubleVector& toClassify);
  virtual int classify(const DoubleVector& toClassify, DoubleVector& dists) {
    ERR << "NYI" << ::std::endl;
    return 0;
  }

  void printConfig();

  double& epsilon() { return epsilon_;}
  uint& minPts() {return minPts_;}
  const double& epsilon() const  { return epsilon_;}
  const uint& minPts() const {return minPts_;}
  virtual int numberOfClusters() {return -1;}
  void setDist(BaseDistance* d) {dist_=d;}
  BaseDistance* getDist() {return dist_;}

  
  virtual void saveModel(const ::std::string filename);
  virtual void loadModel(const ::std::string filename);
};

#endif
