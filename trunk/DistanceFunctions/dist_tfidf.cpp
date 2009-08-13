#include "dist_tfidf.hpp"
#include <map>
#include <math.h>
using namespace std;


//this initializes the term frequencies of the features contained in the query image
void TFIDFDistance::start(const BaseFeature * queryFeature) {
  
  queryFeature_=dynamic_cast<const SparseHistogramFeature*>(queryFeature);
  
  //get feature map to avoid casting everytime some distance is computed
  queryMap_=queryFeature_->getMap();

  //precompute
  queryLength_=queryFeature_->length();
  //compute score for query Image for GIFT-TFIDF normalization
  queryScore_=scoreFeatureSet(queryFeature_);

 
  
}

//score a featureset according to some scoring function
double TFIDFDistance::scoreFeatureSet(const SparseHistogramFeature * featureset){

  MapTypeDouble tmpMap=featureset->getMap();
  double tmpscore=0;
  for (MapTypeDouble::iterator i=tmpMap.begin();i!=tmpMap.end();i++){
    tmpscore+=scoreFeature(featureset,i);
  }
  return tmpscore;
}


//score single feature
double TFIDFDistance::scoreFeature(const SparseHistogramFeature * featureset, MapTypeDouble::iterator F){

  //gift tf/idf
  return (F->second)*collectionFrequencies_[F->first];
}

//Dummy for the Moment
MapTypeDouble TFIDFDistance::getTermFrequencies(MapTypeDouble inQuery){
  return inQuery;
}

//euclidean length of feature vector
double TFIDFDistance::getLength(MapTypeDouble inMap){

  double length=0;
  for (MapTypeDouble::iterator i=inMap.begin();i!=inMap.end();i++){
    double tmp=i->second;//*collectionFrequencies_[i->first];
    length+=tmp*tmp;
  }
  return sqrt(length);

}


//distance computation with tf/idf scoring
double TFIDFDistance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
  
  double result=0.0;
  const SparseHistogramFeature * data=dynamic_cast<const SparseHistogramFeature*>(databaseFeature);
  const SparseHistogramFeature* query=dynamic_cast<const SparseHistogramFeature*>(queryFeature);
  //get map
  MapTypeDouble databaseMap=data->getMap();
  MapTypeDouble queryMap= query->getMap();

  //iterate over query features
  
  for (MapTypeDouble::iterator i=queryMap.begin();i!=queryMap.end();i++){

    //look for that feature in document
    MapTypeDouble::iterator j=databaseMap.find(i->first); 

    if (j!=databaseMap.end()){
  
      
	//Machery tf idf
	//adding up with tf factors
      //      double wtq=scoreFeature(queryFeature_,i);
      //double wtd=scoreFeature(data,j);
      //result+=wtq*wtd;
      
      //fidelity with idf
      result+=sqrt(i->second)*sqrt(j->second);
      //gift tfidf
      //result+=scoreFeature(query,i);
    }
  }

//  double databaseFeatureLength;
//  databaseFeatureLength=getLength(databaseMap);
//  result/=(queryLength_*databaseFeatureLength);
//  result/=fabs(queryScore_);

  return 1-result;
  //  return result/sqrt(qL*dL);
//  return sqrt(result/count);//result;
}


//clears query-relevant information
void TFIDFDistance::stop(){
  
  //data
  queryMap_.erase(queryMap_.begin(),queryMap_.end());
  //term frequencies
  queryTermFrequencies_.erase(queryTermFrequencies_.begin(),queryTermFrequencies_.end());
  //length and score for convenience
  queryLength_=0;
  queryScore_=0;
}


//this initializes the collection frequencies of all features by iterating over the database
void TFIDFDistance::initialize(Database &db, uint distanceIndex){

  
  //iterate over database
  for (uint i=0;i<db.size();i++){

    //get ImageContainer
    ImageContainer * ic=db[i];
    //get correct Featuremap
    const SparseHistogramFeature * shf=dynamic_cast<const SparseHistogramFeature*>((*ic)[distanceIndex]);
    MapTypeDouble documentMap=shf->getMap(); 
    //compute collectionfrequencies
    for(MapTypeDouble::iterator i=documentMap.begin();i!=documentMap.end();i++){
     
      //search whether Feature is already registered in CollectionFrequencies
      MapTypeDouble::iterator tmp=collectionFrequencies_.find(i->first);
     
      if (tmp==collectionFrequencies_.end()){
	//if not, insert
	collectionFrequencies_[i->first]=1;
      }
      else{
	//if, increase
	collectionFrequencies_[i->first]+=1;
      }
     
    } 
  }
  

  for (MapTypeDouble::iterator i=collectionFrequencies_.begin();i!=collectionFrequencies_.end();i++){
  
    i->second=log(db.size()/i->second);
    i->second*=i->second;
  }
  //save value for speed
  dataBaseSize_=db.size();

}



