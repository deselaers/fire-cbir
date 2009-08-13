#include "dist_lfhungarian.hpp"
using namespace std;

void LFHungarianDistance::start(const BaseFeature *) {
}


double LFHungarianDistance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
  double result=0.0;
  
  const LocalFeatures* db=dynamic_cast<const LocalFeatures*>(databaseFeature);
  const LocalFeatures* query=dynamic_cast<const LocalFeatures*>(queryFeature);
  
  if(db && query) {
    double **DistMatrix;
    int **Result;
    
    DBG(20) << "Allocating memory" << endl;
    uint m=db->numberOfFeatures();
    uint n=query->numberOfFeatures();
    DistMatrix=(double**)calloc(sizeof(double*),m);
    Result=(int**)calloc(sizeof(int*),m);
    for(uint i=0;i<m;++i) {
      DistMatrix[i]=(double*)calloc(sizeof(double),n);
      Result[i]=(int*)calloc(sizeof(int),n);
    }
    
    DBG(20) << "Initializing Distance Matrix" << endl;
    for(uint i=0;i<m;++i) {
      for(uint j=0;j<n;++j) {
        DistMatrix[i][j]=eucdist((*db)[i],(*query)[j]);
      }
    }
    
    DBG(20) << "Hungarian: start()" << endl;
    solveMinWeightEdgeCover(DistMatrix,Result,m,n);
    

    DBG(20) << "summing up the distances: ";
    DBG(50) << "Results Array:" << endl;
    for(uint i=0;i<m;++i) {
      for(uint j=0;j<n;++j) {
        BLINK(50) << Result[i][j];
        result+=Result[i][j]*DistMatrix[i][j];
      }
      BLINK(50) << endl;
    }
    BLINK(20) << result << endl;
      
    DBG(20) << "cleaning up" << endl;
    for(uint i=0;i<m;++i) {
      free(DistMatrix[i]);
      free(Result[i]);
    }
    
    free(DistMatrix); 
    free(Result);
    
  } else {
    ERR << "Features not comparable" << endl;
    result=-1.0;
  }
  return result;
}

void LFHungarianDistance::stop(){
}


double LFHungarianDistance::eucdist(const vector<double> &a, const vector<double>& b) const {
  double  result=0.0;
  double tmp;
  for(uint i=0;i<a.size();++i) {
    tmp=a[i]-b[i];
    result+=tmp*tmp;
  }
  return result;
}
