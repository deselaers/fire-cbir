#include "dist_smart2.hpp"
#include "dist_tfidf.hpp"
#include <map>
#include <math.h>

using namespace std;

void SMART2Distance::start(const BaseFeature * queryFeature) {
  
  queryFeature_=dynamic_cast<const SparseHistogramFeature*>(queryFeature);
  
  //get feature map to avoid casting everytime some distance is computed
  queryMap_=queryFeature_->getMap();

  //precompute
  queryLength_=queryFeature_->length();
  //compute score for query Image for GIFT-TFIDF normalization
  //  queryScore_=scoreFeatureSet(queryFeature_);

}

//score a featureset according to some scoring function
double SMART2Distance::scoreFeatureSet(const SparseHistogramFeature * featureset){

  MapTypeDouble tmpMap=featureset->getMap();
  double tmpscore=0;
  for (MapTypeDouble::iterator i=tmpMap.begin();i!=tmpMap.end();i++){
    tmpscore+=scoreFeature(i);
  }
  return tmpscore;
}

double SMART2Distance::scoreFeature(MapTypeDouble::iterator F){

  return 1+log(F->second);

}
    
double SMART2Distance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature){
  double result=0.0;

  const SparseHistogramFeature * data=dynamic_cast<const SparseHistogramFeature*>(databaseFeature);
  double idf, wtq, gtd, wtd;
  //get map
  MapTypeDouble databaseMap=data->getMap();
  queryScore_=0;
  //iterate over query features
  for (MapTypeDouble::iterator i=queryMap_.begin();i!=queryMap_.end();i++){

    //look for that feature in document
    MapTypeDouble::iterator j=databaseMap.find(i->first); 
    if (j!=databaseMap.end()){

      idf=collectionFrequencies_[i->first];//log(dataBaseSize_/collectionFrequencies_[i->first]);
      wtq=scoreFeature(i)*idf;
      gtd=(1+scoreFeature(j))/(1+ADI_[data].meanTF);
      wtd=gtd/((0.8)*pivot_+0.3*ADI_[data].numSingletons);
      

      result+=wtd;///wtq;    
	queryScore_+=wtq;
     

//      DBG(10)<<idf<<" "<<wtq<<" "<<gtd<<" "<<wtd<<" "<<pivot_<<" "<<ADI_[data].numSingletons<<" "<<ADI_[data].meanTF<<endl;
//      DBG(10)<<result<<endl;
//      counter++;
    }

  }
  

  return 1-result/queryScore_;

}

  




  
void SMART2Distance::initialize(Database &db,uint distanceIndex){


  pivot_=0;
  DBG(10)<<"INITIALIZING"<<endl;
  //initialize collection frequencies
  for (uint i=0;i<db.size();i++){

    //get ImageContainer
    ImageContainer * ic=db[i];
    //get correct Featuremap
    const SparseHistogramFeature * shf=dynamic_cast<const SparseHistogramFeature*>((*ic)[distanceIndex]);
    MapTypeDouble documentMap=shf->getMap(); 
    //compute collectionfrequencies, mean TF, singletons and pivot
    tfsin container;
    container.meanTF=0;
    container.numSingletons=0;
    uint counter=0;
    for(MapTypeDouble::iterator i=documentMap.begin();i!=documentMap.end();i++){
      
      //search whether Feature is already registered in CollectionFrequencies
      MapTypeDouble::iterator tmp=collectionFrequencies_.find(i->first);
     
      if (tmp==collectionFrequencies_.end()){
	//if not, insert
	collectionFrequencies_[i->first]=1;
      }
      else{
	//if, increase
	collectionFrequencies_[i->first]+=1;//i->second;
      }
      //mean term frequency
      container.meanTF+=i->second;
      counter++;

      //count singletons
      if (i->second==1){
	container.numSingletons++;
      }     
    } 
    //pivot=number of singletons in database
    //container.numSingletons/=counter;


    pivot_+=container.numSingletons;
    DBG(10)<<" NS "<<container.numSingletons<<" pivot:"<<pivot_<<endl;
    container.meanTF/=counter;
    ADI_.insert(make_pair(shf,container));
    

  }
  //normalized
  pivot_/=db.size();
  DBG(10)<<pivot_<<endl;
  dataBaseSize_=db.size();
}
