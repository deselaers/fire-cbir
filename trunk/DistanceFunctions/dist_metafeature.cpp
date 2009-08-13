#include "dist_metafeature.hpp"

using namespace std;

double MetaFeatureDistance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {  

  // Get features from db and query image
  const MetaFeature* db = dynamic_cast<const MetaFeature*>(databaseFeature);
  const MetaFeature* query = dynamic_cast<const MetaFeature*>(queryFeature);

  if(db && query) {

    double dist = 0.0;
    map<string, string>::const_iterator mf_it, curr_key;

    for(mf_it = query->values().begin(); mf_it != query->values().end(); ++mf_it) {
      // The dist is increased if db has no such key or the value of
      // this key is different
      curr_key = db->values().find(mf_it->first);
      if(curr_key == db->values().end()) {
	dist += 1.0;
      } else if(curr_key->second != mf_it->second) {
	dist += 1.0;
      }
    
    }
      
    return dist;
  
  } else {
    ERR << "Features not comparable" << ::std::endl;
    return -1.0;
  }
}
