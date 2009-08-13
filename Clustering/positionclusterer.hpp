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
#ifndef __positionclusterer_hpp__
#define __positionclusterer_hpp__

#include <vector>
#include <string>
#include "math.h"
#include "diag.hpp"
#include "baseclusterer.hpp"

// data structure to represent and manipulate 2x2 covariance matrices
struct PositionCovarianceMatrix {
  double xx;
  double xy;
  double yx;
  double yy;

  double normValue;
  double inverted_xx;
  double inverted_xy;
  double inverted_yx;
  double inverted_yy;

  // function to invert a covariance matrix
  // be assured, it works ;-)
  PositionCovarianceMatrix invert() const {
    PositionCovarianceMatrix inverse;
    double tmp1 = 1.0 / xx;
    double tmp2 = yx / xx;
    double tmp3 = xy / xx;
    double tmp4 = yy - yx * tmp3;
    inverse.xx = tmp1 + (tmp2 * tmp3 / tmp4);
    inverse.xy = -tmp3 / tmp4;
    inverse.yx = -tmp2 / tmp4;
    inverse.yy = 1.0 / tmp4;
    return inverse;
  }

  // multiplies this covariance matrix with another one
  // and returns the result
  PositionCovarianceMatrix multiply(const PositionCovarianceMatrix& b) const {
    PositionCovarianceMatrix c;
    c.xx = xx * b.xx + xy * b.yx;
    c.xy = xx * b.xy + xy * b.yy;
    c.yx = yx * b.xx + yy * b.yx;
    c.yy = yx * b.xy + yy * b.yy;
    return c;
  }

  // add another covariance matrix to this one
  // and returns the result
  PositionCovarianceMatrix add(const PositionCovarianceMatrix& b) const {
    PositionCovarianceMatrix c;
    c.xx = xx + b.xx;
    c.xy = xy + b.xy;
    c.yx = yx + b.yx;
    c.yy = yy + b.yy;
    return c;
  }

  // subtracts another covariance matrix from this one
  // and returns the result
  PositionCovarianceMatrix subtract(const PositionCovarianceMatrix& b) const {
    PositionCovarianceMatrix c;
    c.xx = xx - b.xx;
    c.xy = xy - b.xy;
    c.yx = yx - b.yx;
    c.yy = yy - b.yy;
    return c;
  }

  double getDeterminant() const {
    return xx * yy - yx * xy;
  }

  double getTheta() const {
    if (fabs(xx - yy) < 1E-10) {
      return M_PI / 2.0;
    } else {
      return 0.5 * atan((2.0 * xy) / (xx - yy));
    }
  }

  void rotate(double theta) {
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double new_xx = (xx * cos_theta - yx * sin_theta) * cos_theta + (xy * cos_theta - yy * sin_theta) * -sin_theta;
    double new_xy = (xx * cos_theta - yx * sin_theta) * sin_theta + (xy * cos_theta - yy * sin_theta) * cos_theta;
    double new_yx = (xx * sin_theta + yx * cos_theta) * cos_theta + (xy * sin_theta + yy * cos_theta) * -sin_theta;
    double new_yy = (xx * sin_theta + yx * cos_theta) * sin_theta + (xy * sin_theta + yy * cos_theta) * cos_theta;
    xx = new_xx;
    xy = new_xy;
    yx = new_yx;
    yy = new_yy;
  }

  void calcCachedData() {
    normValue = log(1.0 / (2.0 * M_PI * sqrt(getDeterminant())));
    PositionCovarianceMatrix inv = invert();
    inverted_xx = inv.xx;
    inverted_xy = inv.xy;
    inverted_yx = inv.yx;
    inverted_yy = inv.yy;
  }

};

// data structure to represent the position of a cluster, given by its mean value and covariance matrix
struct ClusterPosition {
  ::std::pair<double, double> mean;
  PositionCovarianceMatrix covariance;
  uint elements;
  double elements_log;
  ClusterPosition() {
    mean.first = 0.0; mean.second = 0.0;
    covariance.xx = 0.0; covariance.xy = 0.0; covariance.yx = 0.0; covariance.yy = 0.0;
  }
};

class PositionClusterer : public BaseClusterer {
private:
  uint maxSplits_;
  uint stopWithNClusters_;
  uint dontSplitBelow_;
  uint iterationsBetweenSplits_;
  uint minObservationsPerCluster_;
  double epsilon_;
  bool diagonalize_;
  double minVariance_;

  ::std::vector<ClusterPosition> clusters_;

  const uint classify(const DoubleVector& observation, const ::std::vector<ClusterPosition> &clusters, DoubleVector& dists) const;
  ClusterPosition getDensity(const ::std::vector<DoubleVector>& inputdata, ResultVector clusterinformation, int cluster);
  ClusterPosition splitCluster(ClusterPosition& cluster);
  void deleteTooSmallClusters(::std::vector<ClusterPosition>& clusters, ResultVector &clusterInformation);

public:

  virtual void run(const ::std::vector<DoubleVector>& inputdata, ResultVector& clusterInformation, ::std::vector< ::std::vector< ::std::vector<ClusterPosition> > >& allClusters, bool saveAllClusters);
  virtual void run(const ::std::vector<DoubleVector>& inputdata, ::std::vector<ClusterPosition>& clusters);
  virtual void run(const ::std::vector<DoubleVector>& inputdata, ::std::vector<ClusterPosition>& clusters, ::std::vector< ::std::vector< ::std::vector<ClusterPosition> > >& allClusters, bool saveAllClusters);

  //methods derived from BaseClusterer
  PositionClusterer();
  virtual void run(const DoubleVectorVector& inputdata, ResultVector &clusterInformation);
  virtual ~PositionClusterer();
  virtual int classify(const DoubleVector& toClassify, DoubleVector& dists);
  virtual int classify(const DoubleVector& toClassify) {
    DoubleVector dists;
    return classify(toClassify,dists);
  }
  virtual int numberOfClusters() { return clusters_.size();}
  void printConfig();
  
  //methods to set and get options
  uint &maxSplit() {return maxSplits_;}
  uint &stopWithNClusters() {return stopWithNClusters_;}
  uint &dontSplitBelow() {return dontSplitBelow_;}
  uint &iterationsBetweenSplits() {return iterationsBetweenSplits_;}
  uint &minObservationsPerCluster() {return minObservationsPerCluster_;}
  double &epsilon() {return epsilon_;}
  double &minVariance() {return minVariance_;}
  bool &diagonalize() { return diagonalize_; }

  // the following methods are derived from BaseClusterer and do nothing
  virtual void saveModel(const ::std::string filename) {
    DBG(25) << "not writing model " << filename << ::std::endl; // just to avoid compiler warning because of unused variable
  }
  virtual void loadModel(const ::std::string filename) { 
    DBG(25) << "not loading model " << filename << ::std::endl; // just to avoid compiler warning because of unused variable
  }
};



#endif
