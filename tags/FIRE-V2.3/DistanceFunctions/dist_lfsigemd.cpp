#include "dist_lfsigemd.hpp"
#include "lfsignaturefeature.hpp"
#include "emd.hpp"


using namespace std;


float dist(feature_t *F1, feature_t *F2)
{
  double res=0.0,tmp;

  for(uint d=0;d<F1->D;++d) {
    tmp=F1->vec[d]-F2->vec[d];
    res+=tmp*tmp;
  }
  return sqrt(res); 
}


double LFSignatureEMDistance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
  double result=0.0;
    
  const LFSignatureFeature* db=dynamic_cast<const LFSignatureFeature*>(databaseFeature);
  const LFSignatureFeature* q=dynamic_cast<const LFSignatureFeature*>(queryFeature);
    
  if(db && q) {
    result=0;

    uint dbsize=db->signature().size();
    uint qsize=q->signature().size();

    if(dbsize==0 or qsize==0) return 2.0;

    uint D=db->model()[0].dim;
      
    // make weights for emd;
    float *qweights=new float[qsize];
    float *dbweights=new float[dbsize];
    
    // get size of larger signature
    uint dbsum=0; for(uint i=0;i<db->signature().size();++i) {dbsum+=db->signature()[i];}
    uint qsum=0; for(uint i=0;i<q->signature().size();++i) {qsum+=q->signature()[i];}
      
    for(uint i=0;i<qsize;++i) {
      qweights[i]=float(q->signature()[i])/float(qsum);
    }
    
    for(uint i=0;i<dbsize;++i) {
      dbweights[i]=float(db->signature()[i])/float(dbsum);
    }

    // make feature for emd
    feature_t *dbfeat=new feature_t[dbsize];//
    for(uint i=0;i<dbsize;++i) {
      dbfeat[i].D=D;
      dbfeat[i].vec=new float[D];
      for(uint d=0;d<D;++d) {
        dbfeat[i].vec[d]=db->model()[i].mean[d];
      }
    }
      
    feature_t *qfeat=new feature_t[qsize];
    for(uint i=0;i<qsize;++i) {
      qfeat[i].D=D;
      qfeat[i].vec=new float[D];
      for(uint d=0;d<D;++d) {
        qfeat[i].vec[d]=q->model()[i].mean[d];
      }
    }
      

    signature_t dbsig={db->signature().size(),dbfeat,dbweights};
    signature_t qsig={q->signature().size(),qfeat,qweights};
      
    DBG(20) << "Starting EMD ...";
    result=emd(&dbsig, &qsig, dist,0,0);
    BLINK(20)  << " emd=" <<result << endl;

    //cleaning up
    for(uint i=0;i<dbsize;++i) {delete[] dbfeat[i].vec;} delete[] dbfeat;
    for(uint i=0;i<qsize;++i) {delete[] qfeat[i].vec;} delete[] qfeat;
    delete[] qweights;
    delete[] dbweights;
    
    DBG(100) << VAR(result) << endl;
    return result;
      
  } else {
    ERR << "Features not comparable" << endl;
    return -1.0;
  }
}
