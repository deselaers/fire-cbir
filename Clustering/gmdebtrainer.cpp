#include "gmdebtrainer.hpp"
#include "factory.hpp"

double BaseCriterion::fprime(const FeatureVector& X, const uint k, const UntiedGaussMixtureDensityClassifier& classifier) {
  vector<double> log_p_xc;
  classifier->classify(X,log_p_xc);
  for(uint i=0;i<log_p_xc.size();++i) {
    log_p_xc[i]*=classifier.priors_[i];
  }
  
  double infprime=log_p_xc[k];
  double sum=0.0;
  for(uint c=0;c<log_p_xc.size();++c) {
    sum+=log_p_xc[c];
  }
  infprime/=sum;
  
  return fprime(infprime);
}


UntUntiedGaussMixtureDensityEBTrainer::UntiedGaussMixtureDensityEBTrainer(UntiedGaussMixtureDensityClassifier &classifier, const ::std::string criterion) : classifier_(classifier) {
  Factory<BaseCriterion,BaseCriterion* (*)() ,string> criterionfactory;
  Factor.registerClass("mmi",BaseScoring::create<MMICriterion>);
  
  criterion_=factory.getObject(criterion);

  C_=classifier_.mixtures_.size();
  D_=classifier_.mixtures_[0].dim_;
}

void UntiedGaussMixtureDensityEBTrainer::retrain(const DataSet& traindata) {
  vector<uint>   maxApproxIdx(N_);
  
  for(uint it=0;it<maxIterations_;++it) {
    UntiedGaussMixtureDensityClassifier oldClassifier(classifier_);
    DBG(10) << "Iteration " << it+1 << "/" << maxIterations_ <<"." << endl;
    
    for(uint c=0;c<C_;++c) {
      DBG(10) << "Iteration " << it+1 << "/" << maxIterations_ <<" class " << c << "/" << C_ << "." << std::endl;
      
      uint I=classifier_.mixtures_[c].numberOfClusters();
      GaussMixtureDensity& gmdi=classifier_.mixtures_[c];
      
      vector<double> gamma1(I*D_,0.0);
      vector<double> gammaX(I*D_,0.0);
      vector<double> gammaX2(I*D_,0.0);
      
      DBG(15) << "Getting gammas" << endl;
      for(uint n=0;n<N_;++n) {
        double maxLogPx=-numeric_limits<double>::max();
        
        for(uint i=0;i<I;++i) {
          double logpx=gmdi.logp_x(traindata[n]);
          if(logpx>maxLogPx) {
            maxLogPx=logpx;
            maxApproxIdx[n]=i;
          } 
          
          vector<double> p_cn_xn; classifier.p_c_x(traindata[n],p_cn_xn);
          double gamma=1.0/double(N_)*criterion_->fprime(traindata[n],c,classifier_)*(delta(c,traindata[n].cls())-p_cn_xn[c]);
          uint in=maxApproxIdx[n];
          
          for(uint d=0;d<D_;++d) {
            uint idx=d*I+in;
            double xd=traindata[n][d];
            
            gamma1[idx]+=gamma;
            gammaX[idx]+=gamma*xd;
            gammaX2[idx]+=gamma*xd*xd;
          }
        }



      }
    } // for c
  } // for it
}
    
