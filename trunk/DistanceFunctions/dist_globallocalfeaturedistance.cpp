#include "dist_globallocalfeaturedistance.hpp"
#include "localfeatures.hpp"

using namespace std;

double GlobalLocalFeatureDistance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
#ifdef HAVE_KDTREE_LIBRARY
  const LocalFeatures* db=dynamic_cast<const LocalFeatures*>(databaseFeature);
  const LocalFeatures* query=dynamic_cast<const LocalFeatures*>(queryFeature);
  
  if(db && query) {
    double result;
    
    if(query->filename() != calculatedFor) {
      ERR << "System not prepared for this query:" << query->filename() << " beeing prepared for " << calculatedFor  << endl;
      exit(20);
    }
    
    if( hits.find(db->filename()) ==hits.end()) {
      DBG(25) << "No single hit for this database image: " << db->filename() << endl;
      result=hitcounter;
    } else {
      result=hitcounter-hits[db->filename()];
    }
    return result;
  } else {
#endif
    ERR << "Features not comparable" << ::std::endl;
    return -1.0;
#ifdef HAVE_KDTREE_LIBRARY
  }
#endif
}

void GlobalLocalFeatureDistance::start(const BaseFeature *queryFeature) {
#ifdef HAVE_KDTREE_LIBRARY
  const LocalFeatures* query=dynamic_cast<const LocalFeatures*>(queryFeature);
  if(query) {
    calculatedFor=query->filename();
    hits.clear();
    hitcounter=0;
    float *point=new float[query->dim()];
    for(uint i=0;i<query->numberOfFeatures();++i) {
      for(uint j=0;j<query->dim();++j) {
        point[j]=(*query)[i][j];
      }
      knn_kdtree_search_idc(kdt,point,k,epsilon,neighbors);
      for(uint i=0;i<k;++i) {
        hits[string(neighbors[i].labelvec)]++;
        ++hitcounter;
      }
    }
    delete[] point;
    
    DBGI(20,{ for(map<string,uint>::const_iterator i=hits.begin();i!=hits.end();++i) { DBG(20)<< i->first << " " << i->second << endl; }})
    
  } else {
#endif
    ERR << "Not possible to use this feature for this distance!" << endl;
    exit(20);
#ifdef HAVE_KDTREE_LIBRARY 
}
#endif
}

void GlobalLocalFeatureDistance::stop() {
#ifdef HAVE_KDTREE_LIBRARY 
  hits.clear();
  hitcounter=0;
#endif
}


void GlobalLocalFeatureDistance::loadTree(::std::string kdtreefilename) {
#ifdef HAVE_KDTREE_LIBRARY
  DBG(20) << "Loading kdtree from file " << kdtreefilename << endl;
  char *fn=new char[kdtreefilename.size()+1];
  strcpy(fn,kdtreefilename.c_str());
  knn_kdtree_read(fn,0,kdt); // load non binary
  delete fn;
#endif
}
