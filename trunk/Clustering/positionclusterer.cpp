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
#include <sstream>
#include <limits>
#include "em.hpp"
#include "diag.hpp"
#include "baseclusterer.hpp"
#include "vectorfeature.hpp"
#include "basefeature.hpp"
#include "positionclusterer.hpp"
#include "gzstream.hpp"

using namespace std;

PositionClusterer::PositionClusterer() {
  stopWithNClusters_=0;
}


void PositionClusterer::printConfig() {
  cout << "LBG clustering settings" << endl
       << "  maxSplits               = " << maxSplits_ << endl
       << "  stopWithNClusters       = " << stopWithNClusters_ << endl
       << "  dontSplitBelow          = " << dontSplitBelow_ << endl
       << "  iterationsBetweenSplits = " << iterationsBetweenSplits_ << endl
       << "  minObsPerCluster        = " << minObservationsPerCluster_ << endl
       << "  epsilon                 = " << epsilon_ << endl
       << endl;
}

PositionClusterer::~PositionClusterer() {
}

ClusterPosition PositionClusterer::splitCluster(ClusterPosition& cluster) {

  ClusterPosition tmpCluster;

  double angle = cluster.covariance.getTheta();
  if (cluster.covariance.xx - cluster.covariance.yy < 0.0) {
    angle += M_PI / 2;
  }
  DBG(25) << "mean is: " << cluster.mean.first << " " << cluster.mean.second << endl;
  DBG(25) << "covariance matrix is:" << endl;
  DBG(25) << cluster.covariance.xx << " " << cluster.covariance.xy << endl;
  DBG(25) << cluster.covariance.yx << " " << cluster.covariance.yy << endl;
  DBG(25) << "angle: " << angle << " (" << (angle * 180 / M_PI) << ")" << endl;
  double disturbX = cos(angle);
  double disturbY = sin(angle);
  DBG(25)  << "disturbing by " << disturbX << " " << disturbY << endl;


  // calculate the two eigenvalues of the covariance matrix
  double eigen1 = 0.5 * ((cluster.covariance.xx + cluster.covariance.yy) + sqrt(4 * cluster.covariance.xy * cluster.covariance.yx + ((cluster.covariance.xx - cluster.covariance.yy) * (cluster.covariance.xx - cluster.covariance.yy))));
  double eigen2 = 0.5 * ((cluster.covariance.xx + cluster.covariance.yy) - sqrt(4 * cluster.covariance.xy * cluster.covariance.yx + ((cluster.covariance.xx - cluster.covariance.yy) * (cluster.covariance.xx - cluster.covariance.yy))));
  DBG(25) << "eigenvalues: " << eigen1 << ", " << eigen2 << endl;
  if (isnan(eigen1) || isnan(eigen2)) {
    ERR << "eigenvalues: " << eigen1 << ", " << eigen2 << endl;
    ERR << "covariance matrix is:" << endl;
    ERR << cluster.covariance.xx << " " << cluster.covariance.xy << endl;
    ERR << cluster.covariance.yx << " " << cluster.covariance.yy << endl;
    exit(1);
  }

  double maxEigen;
  if (abs(eigen1) > abs(eigen2)) {
    maxEigen = eigen1;
  } else {
    maxEigen = eigen2;
  }
  if (fabs(cluster.covariance.xy) < 1E-10) {
    if (cluster.covariance.xx > cluster.covariance.yy) {
      disturbX = cluster.covariance.xx;
      disturbY = 0.0;
    } else {
      disturbX = 0.0;
      disturbY = cluster.covariance.yy;
    }
    DBG(25) << "trivial eigenvector for eigenvalue " << maxEigen << ": " << disturbX << "," << disturbY << endl;
  } else {
    disturbX = 1.0;
    disturbY = ((maxEigen - cluster.covariance.xx) * disturbX) / cluster.covariance.xy;
    DBG(25) << "eigenvector for eigenvalue " << maxEigen << ": " << disturbX << "," << disturbY << endl;

    double lengtheigenvector = sqrt(disturbX * disturbX + disturbY * disturbY);
    disturbX /= lengtheigenvector;
    disturbY /= lengtheigenvector;
    DBG(25) << "normalized eigenvector for eigenvalue " << maxEigen << ": " << disturbX << "," << disturbY << endl;
  }

  if (isnan(disturbX) || isnan(disturbY)) {
    ERR << "disturbX = " << disturbX << ", disturbY = " << disturbY << endl;
    exit(1);
  }

  tmpCluster.mean.first = cluster.mean.first + epsilon_ * disturbX;
  tmpCluster.mean.second = cluster.mean.second + epsilon_ * disturbY;
  tmpCluster.covariance = cluster.covariance;
  cluster.mean.first -= epsilon_ * disturbX;
  cluster.mean.second -= epsilon_ * disturbY;

  return tmpCluster;
}

void PositionClusterer::run(const vector<DoubleVector>& inputdata, vector<ClusterPosition>& clusters,
                            vector< vector< vector<ClusterPosition> > >& allClusters, bool saveAllClusters) {
  clusters_.clear();
  ResultVector clusterInformation;
  vector< vector< vector<ClusterPosition> > > dummy;
  run(inputdata, clusterInformation, allClusters, saveAllClusters);
  clusters = clusters_;
}

void PositionClusterer::run(const vector<DoubleVector>& inputdata, vector<ClusterPosition>& clusters) {
  clusters_.clear();
  ResultVector clusterInformation;
  vector< vector< vector<ClusterPosition> > > dummy;
  run(inputdata, clusterInformation, dummy, false);
  clusters = clusters_;
}


void PositionClusterer::run(const DoubleVectorVector& inputdata, ResultVector &clusterInformation) {
  vector<DoubleVector> input(inputdata.size());
  for (uint i = 0; i < inputdata.size(); i++) {
    input[i] = *(inputdata[i]);
  }
  vector< vector< vector<ClusterPosition> > > dummy;
  run(input, clusterInformation, dummy, false);
}

void PositionClusterer::run(const vector<DoubleVector>& inputdata, ResultVector &clusterInformation,
                            vector< vector< vector<ClusterPosition> > >& allClusters, bool saveAllClusters) {
  DBG(20) << "Starting to cluster with " << inputdata.size() << " features" << endl;

  //base variables
  uint nOfObservations=inputdata.size();
  
  ClusterPosition tmpCluster;
  clusterInformation = vector<int>(nOfObservations, 0);

//  if(inputdata.size()<dontSplitBelow_) {
//    DBG(20) << "Too few observations to estimate clustermodel" << endl;
//    tmpCluster=getDensity(inputdata,clusterInformation,0);
//    return;
//  }
  
  
  if (clusters_.size() ==0 ) {
    //estimate inital gaussian
    tmpCluster = getDensity(inputdata, clusterInformation, 0);
    clusters_.push_back(tmpCluster);
    DBG(20) << "Init distribution estimated" << endl;
  } else {
    DBG(20) << "Skipping initialization as we already have a clustermodel" << endl;
  }
 
  if (saveAllClusters) {
    allClusters.resize(maxSplits_ + 1);
    allClusters[0].resize(1);
    allClusters[0][0] = clusters_;
  }
  // split iteration
  for(uint split = 0; split < maxSplits_; split++) {
    DBG(20) << "Splitting " <<split+1 << "/" << maxSplits_<<": "<< clusters_.size() << " clusters." <<endl;
    uint nOfClusters = clusters_.size();
    
    if ((nOfClusters < stopWithNClusters_) || (stopWithNClusters_ == 0)) {
      for (int i = nOfClusters - 1; i >= 0; --i) {
      	if (clusters_[i].elements > dontSplitBelow_) {
          DBG(20) << "Splitting cluster " << i << " with " << clusters_[i].elements << " elements."<< endl;
          tmpCluster = splitCluster(clusters_[i]);
          clusters_[i].elements = 0;
          tmpCluster.elements = 0;
          clusters_.push_back(tmpCluster);
        } else {
      	  DBG(20) << "Not splitting. Cluster has only " << clusters_[i].elements << " elements." << endl;
      	}
      }

      if (saveAllClusters) {
        allClusters[split + 1].resize(iterationsBetweenSplits_);
      }

      // reestimate iterations
      nOfClusters = clusters_.size();
      uint reestimation = 0; 
      bool deleted = false;
      while ((reestimation < iterationsBetweenSplits_) || deleted) {
        DBG(20) << "Reestimation Iteration " << (reestimation + 1) << "/" << iterationsBetweenSplits_ << " " << clusters_.size() << " clusters: ";
        for(uint i = 0; i < clusters_.size(); ++i) {
          BLINK(20) << clusters_[i].elements << " ";
          clusters_[i].elements=0;
        }
        BLINK(20) << endl;

        DBG(25) << "Reassigning" <<endl;
        //reassigning the features to the clusters
        DoubleVector dists;
        for(uint i = 0; i < inputdata.size(); ++i) {
          int bestCluster = classify(inputdata[i], clusters_, dists);
          clusterInformation[i] = bestCluster;
          ++clusters_[bestCluster].elements;
        }
        
        DBG(25)  <<"#(Clusters):" << clusters_.size() << " after: "<< (split + 1) << "/" << maxSplits_ << " splits, " << (reestimation + 1) << "/" << iterationsBetweenSplits_ << " reestimations." << endl;
        for (vector<ClusterPosition>::iterator aktClust  =clusters_.begin(); aktClust <clusters_.end(); ++aktClust) {
          DBG(25) << "   Reassigned: " << aktClust->elements << " entries." << endl;
        }
        
        //check whether there are any clusters which are too small
        DBG(25) << "Deleting too small clusters" << endl;
        uint nOfClustersBeforeDeletion = clusters_.size();
        deleteTooSmallClusters(clusters_, clusterInformation);
        nOfClusters = clusters_.size();
        deleted = (nOfClusters != nOfClustersBeforeDeletion);
	
        DBG(30) << "Now having " << nOfClusters << "clusters." << endl;
        if (nOfClusters == 0) {
          // all clusters where deleted...
          clusterInformation = vector<int>(nOfObservations, 0);
          tmpCluster = getDensity(inputdata, clusterInformation, 0);
          clusters_.push_back(tmpCluster);
        } else {
          // Reestimation
          DBG(30) << "Reestimating" << endl;
          for(uint i = 0; i < nOfClusters; ++i) {
            clusters_[i] = getDensity(inputdata, clusterInformation, i);
          }
        }
        if (saveAllClusters) {
          allClusters[split + 1][reestimation] = clusters_;
        }
        reestimation++;
      }
      
    } else {
      DBG(15) << "stopWithNClusters=" << stopWithNClusters_ << " reached." << endl;
    }
  } 

  DBG(30) << " ending" << endl;
}


void PositionClusterer::deleteTooSmallClusters(vector<ClusterPosition>&clusters, ResultVector &clusterInformation) {
  vector<ClusterPosition>::iterator aktClust = clusters_.begin();
  int aktClustNo = 0;
  while (aktClust != clusters.end()) {
    DBG(30) << "cluster" << aktClustNo << " has " << aktClust->elements << " entries." << endl;
    if (aktClust->elements < minObservationsPerCluster_) {
      DBG(20) << "deleting cluster with " << aktClust->elements << " elements." << endl;
      aktClust = clusters.erase(aktClust);
      for(uint i = 0; i < clusterInformation.size(); ++i) {
        if (clusterInformation[i]==aktClustNo) {
          clusterInformation[i] = -1;
        } else if (clusterInformation[i] > aktClustNo) {
          --clusterInformation[i];
        }
      }
    } else {
      ++aktClust;
      ++aktClustNo;
    }
  }
}

ClusterPosition PositionClusterer::getDensity(const vector<DoubleVector>& inputdata, ResultVector clusterInformation, int cluster) {
  uint nOfelements=0;
  ClusterPosition result;
  DBG(25) << "getDensity(...) entered!" << endl;

  for(uint i = 0; i < clusterInformation.size(); ++i) {
    if (clusterInformation[i] == cluster) {
      ++nOfelements;
      double xValue = inputdata[i][0];
      double yValue = inputdata[i][1];
      result.mean.first += xValue;
      result.mean.second += yValue;
    }
  }
  if (nOfelements == 0) {
    ERR << "must not be 0 here!" << endl;
    exit(1);
  }
  result.elements = nOfelements;
  result.mean.first /= (double) nOfelements;
  result.mean.second /= (double) nOfelements;
  for(uint i = 0; i < clusterInformation.size(); ++i) {
    if (clusterInformation[i] == cluster) {
      double xValue = inputdata[i][0];
      double yValue = inputdata[i][1];
      result.covariance.xx += (xValue - result.mean.first) * (xValue - result.mean.first);
      result.covariance.xy += (xValue - result.mean.first) * (yValue - result.mean.second);
      result.covariance.yx += (yValue - result.mean.second) * (xValue - result.mean.first);
      result.covariance.yy += (yValue - result.mean.second) * (yValue - result.mean.second);
    }
  }
  result.covariance.xx /= nOfelements;
  result.covariance.xy /= nOfelements;
  result.covariance.yx /= nOfelements;
  result.covariance.yy /= nOfelements;

  /*

  // this is a quicker, but numerically critical way to estimate the density

  for(uint i = 0; i < clusterInformation.size(); ++i) {
  if (clusterInformation[i] == cluster) {
  ++nOfelements;
  double xValue = inputdata[i][0];
  double yValue = inputdata[i][1];
  result.mean.first += xValue;
  result.mean.second += yValue;
  result.covariance.xx += xValue * xValue;
  result.covariance.xy += xValue * yValue;
  result.covariance.yx += yValue * xValue;
  result.covariance.yy += yValue * yValue;
  }
  }

  DBG(25) << "calculating density from " << nOfelements << " observations" << endl;

  if (nOfelements == 0) {
  ERR << "must not be 0 here!" << endl;
  exit(1);
  }
  result.elements = nOfelements;

  result.mean.first /= nOfelements;
  result.mean.second /= nOfelements;
  result.covariance.xx /= nOfelements;
  result.covariance.xy /= nOfelements;
  result.covariance.yx /= nOfelements;
  result.covariance.yy /= nOfelements;
  result.covariance.xx -= result.mean.first * result.mean.first;
  result.covariance.xy -= result.mean.first * result.mean.second;
  result.covariance.yx -= result.mean.second * result.mean.first;
  result.covariance.yy -= result.mean.second * result.mean.second;
  */

  DBG(25) << "mean is: " << result.mean.first << " " << result.mean.second << endl;
  DBG(25) << "cov. matrix is: " << result.covariance.xx << " " << result.covariance.xy << " " << result.covariance.yx << " " << result.covariance.yy << endl;

  if (diagonalize_) {
    result.covariance.xy = 0.0;
    result.covariance.yx = 0.0;
  }

  if (isnan(result.covariance.xx) || isnan(result.covariance.xy) || isnan(result.covariance.yx) || isnan(result.covariance.yy)) {
    ERR << "invalid covariance matrix" << endl;
    cout << result.covariance.xx << " " << result.covariance.xy << endl;
    cout << result.covariance.yx << " " << result.covariance.yy << endl;
    cout << result.elements << endl;
    exit(1);
  }

  double theta = result.covariance.getTheta();
  bool adjusted = false;
  PositionCovarianceMatrix savedCovarianceMatrix = result.covariance;
  result.covariance.rotate(-theta);
  // the non-diagonal elements should be 0.0 after rotation, but in practice they might have 
  // values which slightly differ from that
  result.covariance.xy = 0.0;
  result.covariance.yx = 0.0;
  if (result.covariance.xx < minVariance_) {
    result.covariance.xx = minVariance_;
    adjusted = true;
  }
  if (result.covariance.yy < minVariance_) {
    result.covariance.yy = minVariance_;
    adjusted = true;
  }
  PositionCovarianceMatrix adjustedCovarianceMatrix = result.covariance;
  // back-rotation only if the covariance matrix has been adjusted, otherwise use saved copy
  if (adjusted) {
    result.covariance.rotate(theta);
  } else {
    result.covariance = savedCovarianceMatrix;
  }

  if (result.covariance.getDeterminant() - minVariance_ * minVariance_ < -1E-15) {
    cout << result.covariance.xx << " " << result.covariance.xy << endl;
    cout << result.covariance.yx << " " << result.covariance.yy << endl;
    cout << result.covariance.getDeterminant() << endl;
    ERR << "determinant should now be > " << minVariance_ * minVariance_ << " !" << endl;
    cout << "after setting minimal variance, covariance matrix was:" << endl;
    cout << adjustedCovarianceMatrix.xx << " " << adjustedCovarianceMatrix.xy << endl;
    cout << adjustedCovarianceMatrix.yx << " " << adjustedCovarianceMatrix.yy << endl;
    cout << adjustedCovarianceMatrix.getDeterminant() << endl;
    cout << "before setting minimal variance, covariance matrix was:" << endl;
    cout << savedCovarianceMatrix.xx << " " << savedCovarianceMatrix.xy << endl;
    cout << savedCovarianceMatrix.yx << " " << savedCovarianceMatrix.yy << endl;
    cout << savedCovarianceMatrix.getDeterminant() << endl;
    
    exit(1);
  }
  if (isnan(result.covariance.xx) || isnan(result.covariance.xy) || isnan(result.covariance.yx) || isnan(result.covariance.yy)) {
    ERR << "invalid covariance matrix after rotation, theta was " << theta << endl;
    cout << result.covariance.xx << " " << result.covariance.xy << endl;
    cout << result.covariance.yx << " " << result.covariance.yy << endl;
    cout << result.elements << endl;
    exit(1);
  }
  DBG(25) << "result cov. matrix: " << result.covariance.xx << " " << result.covariance.xy << " " << result.covariance.yx << " " << result.covariance.yy << endl;

  if (fabs(result.covariance.xy - result.covariance.yx) > 1E-10) {
    ERR << "covariance matrix is not symmetric !" << endl;
    ERR << result.covariance.xx << " " << result.covariance.xy << " " << result.covariance.yx << " " << result.covariance.yy << endl;
    exit(1);
  }
  if (result.covariance.xy * result.covariance.yx < 0.0) {
    result.covariance.xy = 0.0;
    result.covariance.yx = 0.0;
  }
  
  return result;
}

const uint PositionClusterer::classify(const DoubleVector& aktEl, const vector<ClusterPosition> &clusters, DoubleVector &dists) const {
  double maxProb = -1.0; // any value < 0.0 does the job here
  int bestCluster=-1;

  dists = DoubleVector(clusters.size());

  for (uint j = 0; j < clusters.size(); ++j) {
    double xDiff = aktEl[0] - clusters[j].mean.first;
    double yDiff = aktEl[1] - clusters[j].mean.second;
    PositionCovarianceMatrix inverted = clusters[j].covariance.invert();
    double aktProb = (1.0 / (2.0 * M_PI * sqrt(clusters[j].covariance.getDeterminant()))) *
      exp(-0.5 * (xDiff * xDiff * inverted.xx + xDiff * yDiff * (inverted.xy + inverted.yx) + yDiff * yDiff * inverted.yy));

    if(isnan(aktProb)) {
      ERR << "isnan(aktProb) "<< endl;
      ERR << "x=" << clusters[j].mean.first << endl;
      ERR << "y=" << clusters[j].mean.second << endl;
      ERR << "cov_xx=" << clusters[j].covariance.xx << endl;
      ERR << "cov_xy=" << clusters[j].covariance.xy << endl;
      ERR << "cov_yx=" << clusters[j].covariance.yx << endl;
      ERR << "cov_yy=" << clusters[j].covariance.yy << endl;
      ERR << "inverted covariance matrix:" << endl;
      ERR << "cov_xx=" << inverted.xx << endl;
      ERR << "cov_xy=" << inverted.xy << endl;
      ERR << "cov_yx=" << inverted.yx << endl;
      ERR << "cov_yy=" << inverted.yy << endl;
      ERR << "determinant: " << clusters[j].covariance.getDeterminant() << endl;
      aktProb = 0.0;
      exit(1);
    }

    dists[j] = aktProb;
    if (aktProb > maxProb) {
      maxProb = aktProb;
      bestCluster = j;
    }
  }

  if (bestCluster == -1) {
    ERR << "x=" << aktEl[0] << endl;
    ERR << "y=" << aktEl[1] << endl;
    ERR << "bestcluster == -1, maxprob==" << maxProb << endl;
    ERR << "numClusters = " << clusters.size() << endl;
    for (uint j = 0; j < clusters.size(); j++) {
      ERR << "cluster " << j << ":" << endl;
      ERR << "x=" << clusters[j].mean.first << endl;
      ERR << "y=" << clusters[j].mean.second << endl;
      ERR << "cov_xx=" << clusters[j].covariance.xx << endl;
      ERR << "cov_xy=" << clusters[j].covariance.xy << endl;
      ERR << "cov_yx=" << clusters[j].covariance.yx << endl;
      ERR << "cov_yy=" << clusters[j].covariance.yy << endl;
      ERR << "elements=" << clusters[j].elements << endl;
    }
    exit(1);
  } else {
    BLINK(30) << bestCluster << " ";
  }
  DBG(30) << "dist=" << maxProb << " bestCluster=" << bestCluster<< endl;
  return bestCluster;
}

int PositionClusterer::classify(const DoubleVector& toClassify,DoubleVector &dists) {
  return classify(toClassify, clusters_,dists);
}
