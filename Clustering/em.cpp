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
#include "gaussiandensity.hpp"
#include "vectorfeature.hpp"
#include "basefeature.hpp"
#include "gzstream.hpp"
#include "distancemaker.hpp"

using namespace std;

#include "ScopeTimer.h"

EM::EM() {
  poolMode_=noPooling;
  disturbMode_=varianceDisturb;
  splitMode_=allSplit;
  stopWithNClusters_=0;
  dist_=NULL;
}


void EM::clear() {
  clusters_.resize(0);
}

void EM::printConfig() {
  string distname="built-in";
  if(dist_) distname=getDist()->name();
  cout << "LBG clustering settings" << endl
       << "  splitMode               = " << splitMode_ << endl
       << "  maxSplits               = " << maxSplits_ << endl
       << "  stopWithNClusters       = " << stopWithNClusters_ << endl
       << "  disturbMode             = " << disturbMode_ << endl
       << "  poolMode                = " << poolMode_ << endl
       << "  dontSplitBelow          = " << dontSplitBelow_ << endl
       << "  iterationsBetweenSplits = " << iterationsBetweenSplits_ << endl
       << "  minObsPerCluster        = " << minObservationsPerCluster_ << endl
       << "  epsilon                 = " << epsilon_ << endl
       << "  dist                    = " << distname << endl
       << endl;
}

EM::~EM() {
}

GaussianDensity EM::splitCluster(GaussianDensity& cluster) {
  uint dim=cluster.dim;
  GaussianDensity tmpCluster(dim);
  
  DBG(150) << tmpCluster.mean.size() << " " <<tmpCluster.sigma.size() << endl;

  for(unsigned int j=0;j<dim;++j) {
    switch(disturbMode_) {
    case varianceDisturb:
      tmpCluster.mean[j]=cluster.mean[j]+epsilon_*cluster.sigma[j];
      cluster.mean[j]-=epsilon_*cluster.sigma[j];
      tmpCluster.sigma[j]=cluster.sigma[j];
      break;
    case meanDisturb:
      tmpCluster.mean[j]=cluster.mean[j]+epsilon_*cluster.mean[j];
      cluster.mean[j]-=epsilon_*cluster.mean[j];
      tmpCluster.sigma[j]=cluster.sigma[j];
      break;
    case meanDisturb2:
      if(2*j<dim) {
        tmpCluster.mean[j]=cluster.mean[j]+epsilon_*cluster.mean[j];
        cluster.mean[j]-=epsilon_*cluster.mean[j];
        tmpCluster.sigma[j]=cluster.sigma[j];
      } else {
        tmpCluster.mean[j]=cluster.mean[j]-epsilon_*cluster.mean[j];
        cluster.mean[j]+=epsilon_*cluster.mean[j];
        tmpCluster.sigma[j]=cluster.sigma[j];
      }
      break;
    case constDisturb:
      tmpCluster.mean[j]=cluster.mean[j]+epsilon_;
      cluster.mean[j]-=epsilon_;
      tmpCluster.sigma[j]=cluster.sigma[j];
      break;
    default:
      ERR << "Unknown disturbMode: " << disturbMode_ << endl;
      break;
    }
  }
  return tmpCluster;
}

void EM::run(const DoubleVectorVector& inputdata, ResultVector &clusterInformation) {
  run(inputdata, clusterInformation, "");
}

void EM::run(const DoubleVectorVector& inputdata, ResultVector &clusterInformation, string filename) {
#ifdef SCOPETIMER
  ScopeTimer st1("EM::run");
#endif

  DBG(20) << " Starting " << endl;
  
  if(inputdata.size()==0) return;

  //base variables
  uint dim=inputdata[0]->size();
  DBG(10) << "Clustering " << dim << " dimensional data." << endl;
  uint nOfObservations=inputdata.size();
  
  GaussianDensity tmpCluster;
  clusterInformation=vector<int>(nOfObservations,0);
  clusterMahalanobisNormalizationValues_.clear();
  if(clusters_.size()==0) {
    clusters_=vector<GaussianDensity>();
    
    //estimate inital gaussian
    tmpCluster=getDensity(inputdata,clusterInformation,0);
    unzeroSigma(tmpCluster);
    clusters_.push_back(tmpCluster);
    poolVariances(clusters_);
    if ((poolMode_ == noPooling) || (poolMode_ == clusterPooling)) {
      clusterMahalanobisNormalizationValues_.push_back(calcClusterMahalanobisNormalizationValue(tmpCluster));
    }
    DBG(20) << " Init distribution estimated" << endl;
  } else {
    if ((poolMode_ == noPooling) || (poolMode_ == clusterPooling)) {
      for (int c = 0; c < (int) clusters_.size(); c++) {
	clusterMahalanobisNormalizationValues_.push_back(calcClusterMahalanobisNormalizationValue(clusters_[c]));
      }
    }
    DBG(10) << "Skipping initialization as we already have a clustermodel" << endl;
  }
  

  
#ifdef SCOPETIMER
{
   ScopeTimer st2("EM::run -> split iteration");
#endif
  //--------------------------------------------------------------------
  // split iteration
  
   for(uint split=0;split<maxSplits_;++split) {
     if(saveBeforeSplits_) {
       ostringstream filenamebeforesplit;
      filenamebeforesplit << filename << "-beforeSplit" << split+1 << ".clustermodel.gz";

      DBG(10) << "Saving model to '" << filenamebeforesplit.str() << "'." <<endl;
      saveModel(filenamebeforesplit.str());
    }



    DBG(5) << "Splitting " <<split+1 << "/" << maxSplits_<<": "<< clusters_.size() << " clusters." <<endl;
    //-----------------------------------------------    
    //split
    uint nOfClusters=clusters_.size();
    
    if (nOfClusters < stopWithNClusters_ || stopWithNClusters_ == 0) {
      if(splitMode_==allSplit) {      
        DBG(15) << "allSplit Mode" << endl;
        for(int i=nOfClusters-1;i>=0;--i) {
          if(clusters_[i].elements>dontSplitBelow_) {
            DBG(15) << "Splitting cluster with " << clusters_[i].elements << " elements."<< endl;
            tmpCluster=splitCluster(clusters_[i]);
            clusters_[i].elements=0;
            tmpCluster.elements=0;
            clusters_.push_back(tmpCluster);
	    if ((poolMode_ == noPooling) || (poolMode_ == clusterPooling)) {
	      clusterMahalanobisNormalizationValues_.push_back(clusterMahalanobisNormalizationValues_[i]);
	    }
          } else {
            DBG(10) << "Not splitting. Cluster has only " << clusters_[i].elements << " elements." << endl;
          }
        }
      } else if(splitMode_==largestSplit) {
        DBG(15) << "largestSplit Mode" << endl; 
        int maxCl=-1,maxSize=-1;
        for(unsigned int i=0;i<nOfClusters;++i) {
          if(int(clusters_[i].elements)>maxSize) {
            maxSize=clusters_[i].elements;
            maxCl=i;
          }
        }
        if(maxCl==-1) {
          ERR << "STRANGE: No largest cluster found" << endl;
        } else {
          tmpCluster=splitCluster(clusters_[maxCl]);
          clusters_[maxCl].elements=0;
          tmpCluster.elements=0;
          clusters_.push_back(tmpCluster);
	  if ((poolMode_ == noPooling) || (poolMode_ == clusterPooling)) {
	    clusterMahalanobisNormalizationValues_.push_back(clusterMahalanobisNormalizationValues_[maxCl]);
	  }
        }
      } else if(splitMode_==varianceSplit) {
        DBG(15) << "varianceSplit Mode" << endl; 
        int maxCl=-1;
        double maxVar=-1;
        for(unsigned int i=0;i<nOfClusters;++i) {
          double tmpVar=0.0;
          for(unsigned int j=0;j<dim;++j) {
            tmpVar+=clusters_[i].sigma[j];
          }
          if(tmpVar>maxVar) {
            maxVar=tmpVar;
            maxCl=i;
          }
        }
        if(maxCl==-1) {
          ERR << "STRANGE: No highest variance cluster found" << endl;
        } else {
          tmpCluster=splitCluster(clusters_[maxCl]);
          clusters_[maxCl].elements=0;
          tmpCluster.elements=0;
          clusters_.push_back(tmpCluster);
	  if ((poolMode_ == noPooling) || (poolMode_ == clusterPooling)) {
	    clusterMahalanobisNormalizationValues_.push_back(clusterMahalanobisNormalizationValues_[maxCl]);
	  }
        }
      }
      

      DBGI(200,{ostringstream filenames(""); filenames << "aftersplit-" << split<<".clustermodel";saveModel(filenames.str());});

#ifdef SCOPETIMER
{
   ScopeTimer st3("EM::run -> reestimate iterations");
#endif
      //-----------------------------------------------
      //reestimate iterations
      nOfClusters=clusters_.size();
      for(unsigned int reestimation=0;reestimation<iterationsBetweenSplits_;++reestimation) {
        DBG(10) << "Reestimation Iteration " << reestimation+1 << "/" << iterationsBetweenSplits_ << " " << clusters_.size() << " clusters: ";
        for(unsigned int i=0;i<clusters_.size();++i) {
          BLINK(15) << clusters_[i].elements << " ";
          clusters_[i].elements=0;
        }
        BLINK(10) << endl;
	
        DBG(25) << "Reassigning" <<endl;

#ifdef SCOPETIMER
{
   ScopeTimer st4("EM::run -> reassigning the features");
#endif
        //-----------------------------------------------
        //reassigning the features to the clusters
        

#ifdef _OPENMP
        vector<int> bestClusterResult(inputdata.size());
#endif

	int ids = inputdata.size();
#pragma omp parallel
{

#pragma omp for
        for(int i=0;i<ids;++i)
        {
#ifdef SCOPETIMER
ScopeTimer st6("EM::run -> reassigning the features - for");
#endif

          DoubleVector dists;
          int bestCluster=classify(inputdata[i],clusters_,dists);
          clusterInformation[i]=bestCluster;

#ifdef _OPENMP
          bestClusterResult[i] = bestCluster;
#else
          ++clusters_[bestCluster].elements;
#endif
        }

} // end omp parallel

#ifdef _OPENMP
        for(long i=0;i<inputdata.size();++i) {
          ++clusters_[bestClusterResult[i]].elements;
        }
#endif

#ifdef SCOPETIMER
}
#endif

        DBG(15)  <<"#(Clusters):" << clusters_.size() << " after: "<< split+1 << "/" << maxSplits_ << " splits, " << reestimation+1 << "/" <<iterationsBetweenSplits_ << " reestimations." << endl;
        for(vector<GaussianDensity>::iterator aktClust=clusters_.begin();aktClust<clusters_.end();++aktClust) {
          DBG(15) << "   Reassigned: " << aktClust->elements << " entries." << endl;
        }
        

#ifdef SCOPETIMER
{
   ScopeTimer st5("EM::run -> reestimation");
#endif
        //check whether there are any clusters which are to small
        DBG(25) << "Deleting too small clusters" << endl;
        deleteToSmallClusters(clusters_,clusterInformation);
        nOfClusters=clusters_.size();
        DBG(25) << "Now having " << nOfClusters << "clusters." << endl;
        
        
        //------------------------------------------------------------
        // Reestimation
        DBG(25) << "Reestimating" << endl;

        for(unsigned int i=0;i<nOfClusters;++i) {
          clusters_[i]=getDensity(inputdata,clusterInformation,i);
          //unzeroSigma(clusters_[i], 0.5 / M_PI);
          unzeroSigma(clusters_[i]);
        }
        poolVariances(clusters_);
        // calculate new normalization values, if necessary
        if ((poolMode_ == noPooling) || (poolMode_ == clusterPooling)) {
          clusterMahalanobisNormalizationValues_.resize(nOfClusters);
          for(unsigned int i=0;i<nOfClusters;++i) {
            clusterMahalanobisNormalizationValues_[i] = calcClusterMahalanobisNormalizationValue(clusters_[i]);
          }
        }
        

	DBGI(200,{ostringstream filenames("");filenames << filename <<  "-aftersplit-" << split << ".clustermodel"; saveModel(filenames.str());});

#ifdef SCOPETIMER
}
#endif

      }

#ifdef SCOPETIMER
}
#endif

    } else {
      DBG(15) << "stopWithNClusters=" << stopWithNClusters_ << " reached." << endl;
    }

  }

#ifdef SCOPETIMER
}
#endif
 
  DBG(25) << " ending" << endl;
}


void EM::deleteToSmallClusters(vector<GaussianDensity>&clusters,ResultVector &clusterInformation) {
  
  vector<GaussianDensity>::iterator aktClust=clusters_.begin();
  
  int aktClustNo=0;
  while(aktClust!=clusters.end()) {
    DBG(25) << "cluster" << aktClustNo<< " has " << aktClust->elements << " entries." << endl;
    if (aktClust->elements< minObservationsPerCluster_) {
      DBG(15) << "deleting cluster with "<<aktClust->elements<< " elements."<< endl;
      aktClust=clusters.erase(aktClust);
      for(unsigned int i=0;i<clusterInformation.size();++i) {
        if(clusterInformation[i]==aktClustNo) {
          clusterInformation[i]=-1;
        } else if (clusterInformation[i]>aktClustNo) {
          --clusterInformation[i];
        }
      }
    } else {
      ++aktClust;
      ++aktClustNo;
    }
  }
}

const uint & EM::splitMode(const ::std::string splitModeString) {
  if(splitModeString=="allSplit") {
    splitMode_=allSplit;
  } else if(splitModeString =="largestSplit") {
    splitMode_=largestSplit;
  } else if(splitModeString =="varianceSplit") {
    splitMode_=varianceSplit;
  } else {
    ERR << "cannot understand splitModeString: '" << splitModeString << "'." << endl;
  }
  return splitMode_;
}

const uint & EM::poolMode(const ::std::string poolModeString) {
  if(poolModeString=="noPooling") {
    poolMode_=noPooling;
  } else if( poolModeString=="clusterPooling") {
    poolMode_=clusterPooling;
  } else if( poolModeString=="dimensionPooling") {
    poolMode_=dimensionPooling;
  } else if( poolModeString=="bothPooling") {
    poolMode_=bothPooling;
  } else {
    ERR << "Cannot understand poolModeString: " << poolModeString << endl;
  }
  return poolMode_;
}

const uint &EM::disturbMode(const ::std::string disturbModeString) {
  if(disturbModeString=="constDisturb") {
    disturbMode_=constDisturb;
  } else if(disturbModeString=="meanDisturb") {
    disturbMode_=meanDisturb;
  } else if(disturbModeString=="meanDisturb2") {
    disturbMode_=meanDisturb2;
  } else if(disturbModeString=="varianceDisturb") {
    disturbMode_=varianceDisturb;
  } else {
    ERR << "Cannot understand disturbModeString: " << disturbModeString << endl;
  }
  return disturbMode_;
}

GaussianDensity EM::getDensity(const DoubleVectorVector &inputdata, ResultVector clusterInformation, int cluster) {
  uint nOfelements=0;
  uint dim=inputdata[0]->size();
    
  GaussianDensity result(dim);
  
  DoubleVector *aktEl;
  for(uint i=0;i<clusterInformation.size();++i) {
    if(clusterInformation[i]==cluster) {
      ++nOfelements;
      aktEl=inputdata[i];
      
      for(unsigned int j=0;j<dim;++j) {
        result.mean[j]+=(*aktEl)[j];
        result.sigma[j]+=(*aktEl)[j]*(*aktEl)[j];
      }
      
    }
  }

  if (nOfelements == 0) {
    ERR << "must not be 0 here!" << endl;
  }

  for(unsigned int j=0;j<dim;++j) {
    result.mean[j]/=nOfelements;
    result.sigma[j]/=nOfelements;
  }
  
  result.elements=nOfelements;
  for(unsigned int j=0;j<dim;++j) {
    result.sigma[j]-=(result.mean[j]*result.mean[j]);
  }
    
  return result;
}

double EM::calcClusterMahalanobisNormalizationValue(const GaussianDensity& cluster) {
  double normValue = 0.0;
  for (int d = 0; d < (int) cluster.dim; d++) {
    normValue += log(2.0 * M_PI * cluster.sigma[d]);
  }
  return normValue;
}

void EM::unzeroSigma(GaussianDensity &gd, double minAllowedSigma) {
	/*
  double minSigma=numeric_limits<double>::max();
  int dim=gd.sigma.size();
  for(int i=0;i<dim;++i) {
    if(gd.sigma[i]>minAllowedSigma) {
      if(gd.sigma[i]<minSigma) {
        minSigma=gd.sigma[i];
      }
    }
  }

  minSigma=max(minAllowedSigma,minSigma*0.0001);
  */
  int dim=gd.sigma.size();
  for(int i=0;i<dim;++i) {
    if(gd.sigma[i]<=minAllowedSigma) {
      //gd.sigma[i]=minSigma;
      gd.sigma[i]=minAllowedSigma;
    }
  }
}

const uint EM::classify(const DoubleVector *aktEl, const vector<GaussianDensity> &clusters,DoubleVector &dists) const {
  double minDist=numeric_limits<double>::max();
  int bestCluster=-1;
  double aktDist;
  // these are never freed, but the memory leak is neglectible, as they are allocated only once
  static BaseFeature *t1=new VectorFeature();
  static BaseFeature *t2=new VectorFeature();

  dists=DoubleVector(clusters.size());

  for(uint j=0;j<clusters.size();++j) {
    if(dist_) { // a distance was specified, use it
      dynamic_cast<VectorFeature*>(t1)->data()=*aktEl;
      dynamic_cast<VectorFeature*>(t2)->data()=clusters[j].mean;
      aktDist=dist_->distance(t1, t2);
    } else { // no distance was specified, use mahalanobis with
             // diagonal covariance matrix
      uint dim=aktEl->size();
      aktDist=0.0;
      double tmp;
      double norm=0.0;

      for(uint i=0;i<dim;++i) {
        tmp=(*aktEl)[i]-clusters[j].mean[i];
        tmp*=tmp;
        if(clusters[j].sigma[i]>=1E-8) {
          tmp/=clusters[j].sigma[i];
        } else {
          ERR << "Sigma[" << i << "] too small: " << clusters[j].sigma[i] << endl;
          exit(1);
        }
        aktDist+=tmp;
      }
      
      if ((poolMode_ == dimensionPooling) || (poolMode_ == bothPooling)) {
        norm += dim * log(2*M_PI*clusters[j].sigma[0]);
      } else {
        norm += clusterMahalanobisNormalizationValues_[j];
      }
      
      aktDist+=norm;
      
      DBG(25) << VAR(aktDist) << " " << VAR(norm) << endl;
    }
    
    if(isnan(aktDist)) {ERR << "isnan(aktDist) "<< endl;}
    
    dists[j]=aktDist;
    if (aktDist<minDist) {
      minDist=aktDist;
      bestCluster=j;
    }
  }
  
  if(bestCluster == -1) {
    ERR << "bestcluster == -1, mindist=="<< minDist << endl;
  } else {
    BLINK(25) << bestCluster << " ";
  }
  DBG(25) << "dist=" << minDist << " bestCluster=" << bestCluster<< endl;
  return bestCluster;
}

void EM::poolVariances(vector<GaussianDensity> &clusters) {
  if(poolMode_==noPooling) {
    //no pooling, thus ready
  } else if(poolMode_==clusterPooling) {
    for(unsigned int j=0;j<clusters[0].sigma.size();++j) {
      double pool=0.0;
      int cnt=0;
      for(unsigned int i=0;i<clusters.size();++i) {
        pool+=clusters[i].elements*clusters[i].sigma[j];
        cnt+=clusters[i].elements;
      }
      pool/=double(cnt);
      for(unsigned int i=0;i<clusters.size();++i) {
        clusters[i].sigma[j]=pool;
      }
    }
  } else if(poolMode_==dimensionPooling) {
    for(unsigned int i=0;i<clusters.size();++i) {
      double pool=0.0;
      int cnt=0;
      for(unsigned int j=0;j<clusters[i].sigma.size();++j) {
        pool+=clusters[i].sigma[j];
        ++cnt;
      }
      pool/=double(cnt);
      for(unsigned int j=0;j<clusters[i].sigma.size();++j) {
        clusters[i].sigma[j]=pool;
      }
    }
  } else if(poolMode_==bothPooling) {
    double pool=0.0;
    int cnt=0;
    for(unsigned int i=0;i<clusters.size();++i) {
      for(unsigned int j=0;j<clusters[i].sigma.size();++j) {
        pool+=clusters[i].elements*clusters[i].sigma[j];
        cnt+=clusters[i].elements;
      }
    }
    pool/=double(cnt);
    for(unsigned int i=0;i<clusters.size();++i) {
      for(unsigned int j=0;j<clusters[i].sigma.size();++j) {
        clusters[i].sigma[j]=pool;
      }
    }
  } else {
    ERR << "Unknown Pooling Mode: " << poolMode_ << endl;
  }
}

int EM::classify(const DoubleVector& toClassify,DoubleVector &dists) {
  return classify(&toClassify, clusters_,dists);
}

void EM::saveModel(const ::std::string filename) {
  ogzstream ofs; ofs.open(filename.c_str());
  if(!ofs.good()) {
    ERR << "Unable to open file " << filename << " to save the model." << endl;
    return;
  } else {
    ofs << "# EM/LBG clustering model file" << endl
        << "splitmode " << splitMode_ << endl
        << "maxsplits " << maxSplits_ << endl
        << "stopWithNClusters " <<  stopWithNClusters_ << endl
        << "disturbMode " << disturbMode_<< endl
        << "poolMode " << poolMode_<< endl
        << "dontSplitBelow " << dontSplitBelow_<< endl
        << "iterationsBetweenSplits " <<iterationsBetweenSplits_<< endl
        << "minObservationsPerCluster " << minObservationsPerCluster_<< endl
	<< "distance " << (dist_? dist_->name(): "built-in") << endl
        << "epsilon " <<  epsilon_ << endl
        << "gaussians " << clusters_.size() << endl;
    
    for(uint i=0;i<clusters_.size();++i) {
      GaussianDensity& c=clusters_[i];
      ofs << "density " << i << " elements " << c.elements << endl;
      ofs << "density " << i << " dim " << c.dim << endl;
      ofs << "density " << i << " mean";
      for(uint k=0;k<c.dim;++k) { ofs << " " << c.mean[k]; } ofs << endl;
      
      ofs << "density " << i << " variance";
      for(uint k=0;k<c.dim;++k) { ofs << " " << c.sigma[k]; } ofs << endl;
    }

    ofs.close();
  }
}
void EM::loadModel(const ::std::string filename ) {
  igzstream is; is.open(filename.c_str());
  if(!is||!is.good()) {
    ERR << "Unable to read file '" << filename << "'." << endl;
    return;
  } else {
    string line;
    getline(is,line);
    if(!is.good()) {
      ERR << "Error reading file '" << filename << "'." << endl;
      return;
    } else {
      if (line!="# EM/LBG clustering model file") {
        ERR << "This is probably not an EM-Model file" << endl
            << "Continuing anyway" << endl;
      }
      while(!is.eof()) {
        getline(is,line);
        if(!is.eof()) {
          istringstream iss(line);
          string keyword;
          iss >> keyword;
          if(keyword=="splitmode" ) {
            iss >> splitMode_;
          } else if(keyword=="maxsplits" ) {
            iss >> maxSplits_;
          } else if(keyword=="stopWithNClusters" ) {
            iss >> stopWithNClusters_;
          } else if(keyword=="disturbMode" ) {
            iss >> disturbMode_;
          } else if(keyword=="poolMode" ) {
            iss >> poolMode_;
          } else if(keyword=="dontSplitBelow" ) {
            iss >> dontSplitBelow_;
          } else if(keyword=="iterationsBetweenSplits" ) {
            iss >> iterationsBetweenSplits_;
          } else if(keyword=="minObservationsPerCluster" ) {
            iss >> minObservationsPerCluster_;
          } else if(keyword=="distance" ) {
            string distname;
            iss >> distname;
            if(distname!="built-in") {
              DistanceMaker dm;
              dist_=dm.makeDistance(distname);
            }
          } else if(keyword=="epsilon" ) {
            iss >> epsilon_;
          } else if(keyword=="gaussians" ) {
            uint size;
            iss >> size;
            clusters_.resize(size);
          } else if(keyword=="density" ) {
            uint no;
            iss >> no;
            iss >> keyword;
            if(keyword=="elements") {
              iss >> clusters_[no].elements;
            } else if(keyword=="dim") {
              iss >> clusters_[no].dim;
            } else if(keyword=="mean") {
              GaussianDensity& c=clusters_[no];
              c.mean=vector<double>(c.dim);
              for(uint i=0;i<c.dim;++i) {
                iss >> c.mean[i];
              }
            } else if(keyword=="variance") {
              GaussianDensity& c=clusters_[no];
              c.sigma=vector<double>(c.dim);
              for(uint i=0;i<c.dim;++i) {
                iss >> c.sigma[i];
              }
            } else {
              ERR << "Reading density received unknown keyword in position 3: '" << keyword << "'." << endl;
            }
          } else {
            ERR << "Unknown keyword '" << keyword << "'." << endl;
          }
        }
      } //while
    } //file is ok
    is.close();

	if ((poolMode_ == noPooling) || (poolMode_ == clusterPooling)) {
	  clusterMahalanobisNormalizationValues_.resize(clusters_.size());
	  for(unsigned int i=0;i<clusters_.size();++i) {
	    clusterMahalanobisNormalizationValues_[i] = calcClusterMahalanobisNormalizationValue(clusters_[i]);
	  }
	}


  }
}
