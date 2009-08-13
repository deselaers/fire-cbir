#include "maxentscoringfirstandsecondorder.hpp"
#include "diag.hpp"
#include "gzstream.hpp"
#include <limits>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

/** read a lambda file trained with ME-first for using in this
    Scoring function */
void MaxEntFirstAndSecondOrderScoring::load(const ::std::string &filename) {
  igzstream is; is.open(filename.c_str());
  if(is.good()) {
    is >> numCls_ >> factor_ >> offset_ 
       >> numLambdasCls_ >> numLambdasCom_ >> numLambdasPool_;
    
    if(numCls_ != 2) {
      ERR << "Expecting two classes only, got " << numCls_ << endl;
    }
    
    if(numLambdasCls_ !=0) {
      ERR << "Expecting no class specific lambdas, got " << numLambdasCls_ << endl;
    }
    
    if(numLambdasPool_ !=0) {
      ERR << "Expecting no pooled lambdas, got " << numLambdasPool_ << endl;
    }
    

    DBG(10) << "Reading " << numCls_<< " * "<< numLambdasCom_ << " lambdas from " << filename << endl;
    lambdas_=vector<double>(numCls_*numLambdasCom_);
    for(uint c=0;c<numCls_;++c) {
      for(uint i=0;i<numLambdasCom_;++i) {
        is>> lambdas_[c*numLambdasCom_+i];
      }
    }
    is.close();
  } else {
    ERR << "Error opening lambda file: " << filename << endl;
  }
}


/** 
 * This function is strongly inspired by calcPC from the MaxEnt
 * Toolkit by Daniel Keysers. We assume that class 1 is relevant and
 * return the probability for an image to be from class 1 as score.
 * the parameter is a distance vector of an image and the score is
 * returned. For this the distance vector has to be of the same sort
 * as the distance vectors that are saved by save distances.
 */
double MaxEntFirstAndSecondOrderScoring::getScore(const ::std::vector<double>& dists) {
  double Z=0;
  double pc;
  double max=-numeric_limits<double>::max(); uint best=0;
  uint D=dists.size();

  vector<double> dist2nd(D+D*(D+1)/2+1);
  uint cnt=0;

  for(uint i=0;i<D;++i) {
    dist2nd[cnt++]=dists[i];
  }
  for(uint i=0;i<D;++i) {
    for(uint j=i;j<D;++j) {
      dist2nd[cnt++]=dists[i]*dists[j];
    }
  }
  dist2nd[cnt++]=1.0;
  
  DBG(50) << VAR(cnt) << " features generated" << endl;
  
  vector<double> P(numCls_);
  for(uint c=0;c<numCls_;++c) {
    pc=0.0;
    for(uint i=0;i<numLambdasCom_;++i) {
      pc+=lambdas_[c*numLambdasCom_+i]*(dist2nd[i]*factor_+offset_);
    }
    P[c]=pc;
    if(max<pc) {max=pc; best=c;}
  }
  
  DBG(50) << VAR(P[0]) << " " VAR(P[1]) << " "<< VAR(max)<<  endl;
  
  for(uint c=0;c<numCls_;++c) {
    P[c]-=max; 
    if(P[c]>-700.0) {
      P[c]=exp(P[c]);
    } else {
      P[c]=0.0;
    }
    Z+=P[c];
  }

  for(uint c=0;c<numCls_;++c) {
    P[c]/=Z;
  }
  DBG(30) << VAR(P[0]) << " " VAR(P[1]) << endl;
  return P[1];
}

const ::std::string MaxEntFirstAndSecondOrderScoring::settings() {
  ostringstream oss;
  oss << "size " << lambdas_.size() 
      << " factor " << factor_
      << " offset " << offset_;
  for(uint i=0;i<lambdas_.size();++i) {
    oss << " lambda " << i << " " << lambdas_[i];
  }
  return oss.str();
}
