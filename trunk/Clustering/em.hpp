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
#ifndef __em_hpp__
#define __em_hpp__

#include <vector>
#include <string>
#include "diag.hpp"
#include "baseclusterer.hpp"
#include "gaussiandensity.hpp"
#include "basedistance.hpp"

/** set poolMode to noPooling, to have variances unpooled. */
static const uint noPooling=0;
/** set poolMode to ClusterPmooling, to have variances for all
    clusters equal, but different in every dimension. */
static const uint clusterPooling=1;
/** set poolMode to dimensionPooling, to have variances equal for
    all dimensions, but different for every cluster. */
static const uint dimensionPooling=2;
/** set poolMode to bothPooling, to have variances pooled over
    cluster and dimensions */
static const uint bothPooling=3;
  
/** set disturbMode to constDisturb to have
    mean_new1=mean_old+epsilon, mean_new2=mean_old-epsilon after
    split.*/ 
static const uint constDisturb=0; 
/** set disturbMode to meanDisturb to have
    mean_new1=mean_old+epsilon*mean_old, mean_new2=mean_old-epsilon*mean_old after
    split.*/ 
static const uint meanDisturb=1;
/** set disturbMode to meanDisturb to have
    mean_new1=mean_old+epsilon*sigma, mean_new2=mean_old-epsilon*sigma after
    split.*/ 
static const uint varianceDisturb=2;
/** set disturbMode to meanDisturb2 to have meandisturb, but switch
    between + and - after half of the entries */
static const uint meanDisturb2=3;



/** set splitMode to allSplit to split all densities in each split
    iteration */
static const uint allSplit=0;
/** set splitMode to largestSplit to split only the largest density
    in each split iteration */
static const uint largestSplit=1;
/** set splitMode to varianceSplit to split only density with
    hightes variance in each split iteration */
static const uint varianceSplit=2;

class EM : public BaseClusterer {
private:
  uint splitMode_;
  uint maxSplits_;
  uint stopWithNClusters_;
  uint disturbMode_;
  uint poolMode_;
  uint dontSplitBelow_;
  uint iterationsBetweenSplits_;
  uint minObservationsPerCluster_;
  double epsilon_;
  bool saveBeforeSplits_;
  
  BaseDistance *dist_;
  
  ::std::vector<GaussianDensity> clusters_;
  // the following vector stores values which are precomputed to make classification faster
  // they are only needed if the build-in mahalanobis distance is used and no pooling
  // over dimensions is applied
  ::std::vector<double> clusterMahalanobisNormalizationValues_;

  void poolVariances(::std::vector<GaussianDensity> &clusters);
  const uint classify(const DoubleVector *observation,const ::std::vector<GaussianDensity> &clusters,DoubleVector& dists) const;
  void unzeroSigma(GaussianDensity &gd, double minAllowed=1E-8);
  GaussianDensity getDensity(const DoubleVectorVector &inputdata, ResultVector clusterinformation, int cluster);
  GaussianDensity splitCluster(GaussianDensity& cluster);
  void deleteToSmallClusters(::std::vector<GaussianDensity>& clusters,ResultVector &clusterInformation);
  // calculates the normalization values for mahalanobis distance, see comment above
  double calcClusterMahalanobisNormalizationValue(const GaussianDensity& cluster);

public:
  //methods derived from BaseClusterer
  EM();
  virtual void run(const DoubleVectorVector& inputdata, ResultVector &clusterInformation, ::std::string filename);
  virtual void run(const DoubleVectorVector& inputdata, ResultVector &clusterInformation);
  
  // clear a model to allow reuse of this clusterer
  virtual void clear();


  virtual ~EM();
  
  virtual int classify(const DoubleVector& toClassify,DoubleVector& dists);

  virtual int classify(const DoubleVector& toClassify) {
    DoubleVector dists;
    return classify(toClassify,dists);
  }
  
  virtual int numberOfClusters() { return clusters_.size();}
  
  void printConfig();
  
  //methods to set and get options
  uint &splitMode() {return splitMode_;}
  uint &maxSplit() {return maxSplits_;}
  uint &stopWithNClusters() {return stopWithNClusters_;}
  uint &disturbMode() {return disturbMode_;}
  uint &poolMode() {return poolMode_;}
  uint &dontSplitBelow() {return dontSplitBelow_;}
  uint &iterationsBetweenSplits() {return iterationsBetweenSplits_;}
  uint &minObservationsPerCluster() {return minObservationsPerCluster_;}
  double &epsilon() {return epsilon_;}

  bool& saveBeforeSplits() {return saveBeforeSplits_;}
  const bool saveBeforeSplits() const {return saveBeforeSplits_;}
  
  void setDist(BaseDistance* d) {dist_=d;}
  BaseDistance* getDist() {return dist_;}
  
  const uint &poolMode(const ::std::string poolModeString);
  const uint &disturbMode(const ::std::string disturbModeString);
  const uint &splitMode(const ::std::string splitModeString);

  virtual void saveModel(const ::std::string filename);
  virtual void loadModel(const ::std::string filename);
  virtual const GaussianDensity& center(const uint i) { return clusters_[i];}
};



#endif
