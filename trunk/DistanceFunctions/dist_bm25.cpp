#include "dist_bm25.hpp"
#include "dist_tfidf.hpp"
#include <map>
#include <math.h>

using namespace std;

void BM25Distance::start(const BaseFeature * queryFeature) {

	queryFeature_=dynamic_cast<const SparseHistogramFeature*>(queryFeature);

	//get feature map to avoid casting everytime some distance is computed
	queryMap_=queryFeature_->getMap();

	//precompute
	queryLength_=getLength(queryMap_);
	//compute score for query Image for GIFT-TFIDF normalization
	queryScore_=scoreFeatureSet(queryFeature_);

}

//score a featureset according to some scoring function
double BM25Distance::scoreFeatureSet(const SparseHistogramFeature * featureset) {

	MapTypeDouble tmpMap=featureset->getMap();
	double tmpscore=0;
	for (MapTypeDouble::iterator i=tmpMap.begin(); i!=tmpMap.end(); i++) {
		tmpscore+=scoreFeature(featureset, i, i);
	}
	return tmpscore;
}

double BM25Distance::scoreFeature(const SparseHistogramFeature *featureset,
		MapTypeDouble::iterator F, MapTypeDouble::iterator Q) {

	double cf=collectionFrequencies_[F->first];
	double idf=cf;//log((dataBaseSize_-cf+0.5)/(cf+0.5));

	double K=k1_*((1-b_)+b_*featureset->length()/avgDL_);
	double result=idf*(k1_+1)*F->second*(k3_+1)*Q->second/((K+F->second)*(k3_
			+Q->second));

	return result;
}

double BM25Distance::distance(const BaseFeature* queryFeature,
		const BaseFeature* databaseFeature) {
	double result=0.0;
	const SparseHistogramFeature * data=
			dynamic_cast<const SparseHistogramFeature*>(databaseFeature);

	//get map
	MapTypeDouble databaseMap=data->getMap();
	//  queryScore_=0;
	//iterate over query features
	for (MapTypeDouble::iterator i=queryMap_.begin(); i!=queryMap_.end(); i++) {

		//look for that feature in document
		MapTypeDouble::iterator j=databaseMap.find(i->first);
		if (j!=databaseMap.end()) {
			double tmp;
			if (true) {
				double cf=collectionFrequencies_[j->first];
				double idf=cf;
				double K=k1_*((1-b_)+b_*data->length()/avgDL_);
				double wtq=(k3_+1)*i->second/(k3_+i->second);
				double wtd=(idf*(k1_+1)*j->second/(K+j->second));
				tmp=wtd*wtq;
				//      DBG(10)<<K<<" "<<wtd<<" "<<wtq<<" "<<tmp<<endl;
				//   double tmp=scoreFeature(data,j,i);
				result+=tmp;
			} else {
				result+=scoreFeature(data, j, i);
			}

		}

	}
	result/=queryLength_*getLength(databaseMap);
	// DBG(10)<<result<<endl;
	//  result/=queryScore_;
	return 1-result;

}

void BM25Distance::initialize(Database &db, uint distanceIndex) {
	//initialize collection frequencies
	TFIDFDistance::initialize(db, distanceIndex);

	//compute average document length in database
	avgDL_=0;
	for (uint i=0; i<db.size(); i++) {
		//get ImageContainer
		ImageContainer * ic=db[i];

		//get correct Featuremap
		const SparseHistogramFeature * shf=
				dynamic_cast<const SparseHistogramFeature*>((*ic)[distanceIndex]);
		avgDL_+=shf->length();
	}
	avgDL_/=dataBaseSize_;

}
