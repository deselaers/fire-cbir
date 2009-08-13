#include "dbscan.hpp"
#include "diag.hpp"
#include "basefeature.hpp"
#include "vectorfeature.hpp"
using namespace std;


const int DBSCAN::UNCLASSIFIED=-1;

DBSCAN::DBSCAN() {
}

void DBSCAN::run(const DoubleVectorVector& inputdata, ResultVector &clusterInformation){

  DBG(25) << "Init ...";
  uint nOfObservations=inputdata.size();
  distanceTable_=DoubleVectorVector(nOfObservations);
  for(uint i=0;i<nOfObservations;++i) {
    distanceTable_[i]=new DoubleVector(nOfObservations,-1);
  }
  clusterInformation=vector<int>(nOfObservations,UNCLASSIFIED); // all observations are unclassified
  uint tenpercent=nOfObservations/10;
  uint percent=0;
  BLINK(25) << "done" << endl;
  
  uint clusterId=0;
  for(uint i=0;i<nOfObservations;++i) {
    if(i%tenpercent==0) {DBG(10) << percent << "% done" << endl; percent+=10;}
    if(clusterInformation[i]==UNCLASSIFIED) {
      if(expandClusters(inputdata,clusterInformation,i,clusterId)) {
        ++clusterId;
        DBG(10) << "ClusterId=" << clusterId << endl;
      }
    }
  }
}

double DBSCAN::distance(const DoubleVectorVector &inputdata, const uint i, const uint j) {
  double result;
  if((*distanceTable_[i])[j]>0) {
    result=(*distanceTable_[i])[j];
  } else {
    static BaseFeature *t1=new VectorFeature();
    static BaseFeature *t2=new VectorFeature();
    dynamic_cast<VectorFeature*>(t1)->data()=*(inputdata[i]);
    dynamic_cast<VectorFeature*>(t2)->data()=*(inputdata[j]);
    result=dist_->distance(t1,t2);
    (*distanceTable_[i])[j]=result;
  }
  return result;
}

set<uint> DBSCAN::neighborhood(const DoubleVectorVector& inputdata, uint startObjectID) {
  set<uint> result;
  uint nOfObjects=inputdata.size();
  for(uint i=0;i<nOfObjects;++i) {
    if(distance(inputdata,i,startObjectID)<=epsilon_) {
      result.insert(i);
    }
  }
  return result;
}

bool DBSCAN::expandClusters(const DoubleVectorVector& inputdata, ResultVector& clusterInformation, const uint startObjectId, const uint clusterId) {
  set<uint> seeds=neighborhood(inputdata, startObjectId);
  if(seeds.size()<minPts_) { // this object is noise
    clusterInformation[startObjectId]=NOISE;
    return false; 
  } else { // this object is not noise
    for(set<uint>::const_iterator i=seeds.begin();i!=seeds.end();++i) {
      clusterInformation[*i]=clusterId;
    }
    while(!seeds.empty()) {
      set<uint>::iterator iterator_o=seeds.begin();
      uint o=*iterator_o;
      seeds.erase(iterator_o);
      
      set<uint> neighborHood=neighborhood(inputdata,o);
      if(neighborHood.size() > minPts_) {
        for(set<uint>::iterator i=neighborHood.begin();i!=neighborHood.end();++i) {
          if(clusterInformation[*i]<0) { //UNCLASSIFIED or NOISE
            if(clusterInformation[*i]==UNCLASSIFIED) {
              seeds.insert(*i);
            } // if UNCLASSIFIED
            clusterInformation[*i]=clusterId;
          } // if UNCLASSIFIED or NOISE
        } // for i in neighborHood
      } // if neighborHood > minPts
    } // while (!seeds.empty)
    return true;
  }
}

DBSCAN::~DBSCAN(){
  delete dist_;
  for(uint i=0;i<distanceTable_.size();++i) {
    delete distanceTable_[i];
  }
}

int DBSCAN::classify(const DoubleVector& ){
  ERR << "NYI" << endl;
  return UNCLASSIFIED;
}

void DBSCAN::printConfig(){
  DBG(2) << "Epsilon: " << epsilon_ << endl
         << "MinPts:  " << minPts_ << endl;
}

void DBSCAN::saveModel(const ::std::string ){
  ERR << "NYI" << endl;
}
void DBSCAN::loadModel(const ::std::string){
  ERR << "NYI" << endl;
}
