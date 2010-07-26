#include <vector>
#include <cmath>
#include <limits>
#include "gmd.hpp"
#include "diag.hpp"
#include "streamparser.hpp"

using namespace std;

/* ----------------------------------------
   Implementation GaussDensity
   ---------------------------------------- */

double GaussDensity::calcClusterNormalizationValue()  {
  clusterNormalizationValue_=0.0;
  for (uint d=0; d<dim_;++d) {
    clusterNormalizationValue_+=log(2.0*M_PI*variance_[d]);
  }
  return clusterNormalizationValue_;
}


double GaussDensity::logp_x(const FeatureVector &x) {
  double tmp;
  double sum=0.0;
  for(uint d=0;d<dim_;++d) {
    tmp=x[d]-mean_[d];
    sum+=(tmp*tmp)/variance_[d];
  }
  return -0.5*(sum+clusterNormalizationValue_);
}

void GaussDensity::estimate(const DataSet& data, const ::std::vector<int>& clusterMemberships, const int clusterId, const double minVar) {
  elements_=0;
  mean_=vector<ftype>(dim_,0.0);
  variance_=vector<ftype>(dim_,0.0);
  
  for(uint n=0;n<data.size();++n) {
    if(clusterMemberships[n]==clusterId) {
      FeatureVector* x=data[n];
      ++elements_;
      for(uint d=0;d<dim_;++d) {
        mean_[d]+=(*x)[d];
      }
    }
  }
  
  for(uint d=0;d<dim_;++d) {
    mean_[d]/=double(elements_);
  }
  
  double tmp;
  for(uint n=0;n<data.size();++n) {
    if(clusterMemberships[n]==clusterId) {
      FeatureVector* x=data[n];
      for(uint d=0;d<dim_;++d) {
        tmp=mean_[d]-(*x)[d];
        variance_[d]+=tmp*tmp;
      }
    }
  }
  
  for(uint d=0;d<dim_;++d) {
    variance_[d]/=double(elements_);
  }
  
  unzeroVariance(minVar);
  calcClusterNormalizationValue();
  DBG(50) << "Cluster center from " << elements_ << " estimated:" << endl
          << *this;
}

void GaussDensity::unzeroVariance(double minVar) {
  for(uint d=0;d<dim_;++d) {
    if(variance_[d]<minVar) variance_[d]=minVar;
  }
}

void GaussDensity::write(::std::ostream &os) const {
  os << "dim " << dim_ << ::std::endl
     << "elements " << elements_ << ::std::endl
     << "clusterNormalizationValue " << clusterNormalizationValue_ << ::std::endl
     << "mean " << mean_ << ::std::endl
     << "variance " << variance_ << ::std::endl;
}

void GaussDensity::read(::std::istream &is) {
  StreamParser sp(is);
  sp.getFromLine("dim",dim_);
  sp.getFromLine("elements",elements_);
  sp.getFromLine("clusterNormalizationValue",clusterNormalizationValue_);
  sp.getVectorFromLine("mean",dim_,mean_);
  sp.getVectorFromLine("variance",dim_,variance_);
}

/* ----------------------------------------
   Implementation GaussMixtureDensity
   ---------------------------------------- */
void GaussMixtureDensity::reestimate(const DataSet& data) {
  DBG(25) << "Finding closest cluster center for each observation." << endl;
  vector<int> clusterMemberships(data.size(),0);
  
  for(uint n=0;n<data.size();++n) {
    uint maxc=clusters_.size();
    double maxlogpx=-numeric_limits<double>::max();
    DBG(35) << "logpx: " ;
    for(uint c=0;c<clusters_.size();++c) {
      double logpx=clusters_[c].logp_x(*(data[n]));
      BLINK(35) << logpx << " ";
      if(logpx>maxlogpx) {
        maxlogpx=logpx;
        maxc=c;
      }
    }
    BLINK(35) << "-> " << VAR(maxc) << endl;
    if(maxc!=clusters_.size()) {
      clusterMemberships[n]=maxc;
    }  else {
      ERR << "Strange: No closest cluster center could be found for observation " << VAR(n) << endl;
    }
  }
  DBG(20) << VAR(clusterMemberships) << endl;
  
  DBG(25) << "Reestimating cluster centers and cluster weights:" ;
  DBG(16) << "Members: ";
  for(uint c=0;c<clusters_.size();++c) {
    clusters_[c].estimate(data,clusterMemberships,c,settings_.minVar);
    clusterWeights_[c]=double(clusters_[c].elements())/double(data.size());
    BLINK(16) << " " << clusters_[c].elements();
  }
  BLINK(16) << endl;
  DBG(25) << "Finished reestimation." << endl;
}

GaussDensity GaussMixtureDensity::splitDensity(uint clusterId) {
  GaussDensity &g1=clusters_[clusterId];
  GaussDensity g2(clusters_[clusterId]);
  double epsilon=settings_.epsilon;
  
  DBG(50) << "Splitting density " << clusterId << endl;
  DBG(50) << "Before splitting: " << endl
          << clusters_[clusterId];
  
  switch(settings_.disturbMode) {
  case DMT_VARIANCE:
    DBG(35) << "Disturb by variance x " << VAR(epsilon) << "."<< endl;
    for(uint d=0;d<dim_;++d) {
      g1.mean()[d]+=epsilon*g1.variance()[d];
      g2.mean()[d]-=epsilon*g1.variance()[d];
    }
    break;
  case DMT_MEAN:
    DBG(35) << "Disturb by +mean x " << VAR(epsilon) << "."<< endl;
    for(uint d=0;d<dim_;++d) {
      g1.mean()[d]+=epsilon*g1.mean()[d];
      g2.mean()[d]-=epsilon*g1.mean()[d];
    }
    break;
  case DMT_MEAN2:
    DBG(35) << "Disturb by ±mean x " << VAR(epsilon) << "."<< endl;
    uint d;
    for(d=0;d<dim_/2;++d) {
      g1.mean()[d]+=epsilon*g1.mean()[d];
      g2.mean()[d]-=epsilon*g1.mean()[d];
    }
    for(;d<dim_;++d) {
      g1.mean()[d]-=epsilon*g1.mean()[d];
      g2.mean()[d]+=epsilon*g1.mean()[d];
    }
    break;
  case DMT_CONST:
    DBG(35) << "Disturb by " << VAR(epsilon) << "."<< endl;
    for(uint d=0;d<dim_;++d) {
      g1.mean()[d]+=epsilon;
      g2.mean()[d]-=epsilon;
    }
    break;
  default:
    ERR << "Unknown " << VAR(settings_.disturbMode) << endl;
  }

  DBG(50) << "After splitting (original):" << endl
          << g1 
          << "After splitting (new):" << endl
          << g2;

  return g2;
}

void GaussMixtureDensity::splitDensities() {
  switch(settings_.splitMode) {
  case SMT_ALL:
    DBG(20) << "Splitting all densities" << endl;
    for(int c=clusters_.size()-1;c>=0;--c) {
      if(clusters_[c].elements()>settings_.minObsForSplit) {
        clusters_.push_back(splitDensity(c));
        clusterWeights_[c]/=2.0;
        clusterWeights_.push_back(clusterWeights_[c]);
      }
    }
    break;
  case SMT_LARGEST: {
    DBG(20) << "Splitting only largest density" << endl;
    uint maxElements=0; uint maxc=clusters_.size();
    for(uint c=0;c<clusters_.size();++c) {
      if(clusters_[c].elements()>maxElements) {
        maxElements=clusters_[c].elements();
        maxc=c;
      }
    }
    if(maxc!=clusters_.size()) {
      if(clusters_[maxc].elements()>settings_.minObsForSplit) {
        clusters_.push_back(splitDensity(maxc));
        clusterWeights_[maxc]/=2.0;
        clusterWeights_.push_back(clusterWeights_[maxc]);
      }
    } else {
      ERR << "No largest cluster found. Not splitting any density." << endl;
    }
    break;
  }
  case SMT_VARIANCE: {
    DBG(20) << "Splitting cluster with highest variance only" << endl;
    double maxVar=-numeric_limits<double>::max(); uint maxc=clusters_.size();
    for(uint c=0;c<clusters_.size();++c) {
      double var=0.0;
      for(uint d=0;d<dim_;++d) {
        var+=clusters_[c].variance()[d]*clusters_[c].variance()[d];
      }
      if(var>maxVar) { 
        maxVar=var;
        maxc=c;
      }
    }
    if(maxc!=clusters_.size()) {
      if(clusters_[maxc].elements()>settings_.minObsForSplit) {
        clusters_.push_back(splitDensity(maxc));
        clusterWeights_[maxc]/=2.0;
        clusterWeights_.push_back(clusterWeights_[maxc]);
      }
    } else {
      ERR << "No cluster with largest variance found. Not splitting any density." << endl;
    }
    break;
  }
  default:
    ERR << "Unknown " << VAR(settings_.splitMode) << endl;
  }
}

void GaussMixtureDensity::poolVariances() {
  switch(settings_.poolMode) {
  case  PMT_NOPOOLING:
    DBG(25) << "NoPooling is easy :-) (sigma_cd)" << endl;
    break;
  case  PMT_CLUSTERPOOLING:
    DBG(25) << "Pooling over clusters: sigma_d" << endl;
    for(uint d=0;d<dim_;++d) {
      double sum=0.0;
      for(uint c=0;c<clusters_.size();++c) {
        sum+=clusters_[c].variance()[d];
      }
      sum/=double(clusters_.size());
      for(uint c=0;c<clusters_.size();++c) {
        clusters_[c].variance()[d]=sum;
      }
    }
    for(uint c=0;c<clusters_.size();++c) {
      clusters_[c].calcClusterNormalizationValue();
    }
    break;
  case  PMT_DIMENSIONPOOLING:
    DBG(25) << "Pooling over dimensions: sigma_c" << endl;
    for(uint c=0;c<clusters_.size();++c) {
      double sum=0.0;
      for(uint d=0;d<dim_;++d) {
        sum+=clusters_[c].variance()[d];
      }
      sum/=double(dim_);
      for(uint d=0;d<dim_;++d) {
        clusters_[c].variance()[d]=sum;
      }
    }
    for(uint c=0;c<clusters_.size();++c) {
      clusters_[c].calcClusterNormalizationValue();
    }
    break;
  case  PMT_ALLPOOLING: {
    DBG(25) << "Pooling over dimensions and cluster: sigma" << endl;
    double sum=0.0;
    for(uint c=0;c<clusters_.size();++c) {
      for(uint d=0;d<dim_;++d) {
        sum+=clusters_[c].variance()[d];
      }
    }
    sum/=double(dim_*clusters_.size());
    for(uint c=0;c<clusters_.size();++c) {
      for(uint d=0;d<dim_;++d) {
        clusters_[c].variance()[d]=sum;
      }
    }
    for(uint c=0;c<clusters_.size();++c) {
      clusters_[c].calcClusterNormalizationValue();
    }
    break;
  }
  default:
    ERR << "Unknown " << VAR(settings_.poolMode) << endl;
  }
}

void GaussMixtureDensity::init(const DataSet &data) {
  vector<int> clusterMemberships(data.size(),0);
  GaussDensity tmpCluster(data[0]->size());
  tmpCluster.estimate(data,clusterMemberships,0,settings_.minVar);
  clusters_.push_back(tmpCluster);
  clusterWeights_.push_back(1.0);
}

void GaussMixtureDensity::deleteToSmallClusters() {
  uint c=0;
  vector<GaussDensity>::iterator aktCluster=clusters_.begin();
  while(aktCluster!=clusters_.end()) {
    if(aktCluster->elements()<settings_.minObs) {
      DBG(15) << "Deleting density with " << aktCluster->elements() << " observations." << endl;
      clusters_.erase(aktCluster);
    } else {
      ++aktCluster;
      ++c;
    }
  }

  DBG(25) << "After deleting clusters: reestimate cluster weights: " ;
  clusterWeights_.resize(clusters_.size());
  uint totalElements=0;
  for(uint c=0;c<clusters_.size();++c) { totalElements+=clusters_[c].elements(); }
  for(uint c=0;c<clusters_.size();++c) { 
    clusterWeights_[c]=double(clusters_[c].elements())/double(totalElements); 
    BLINK(25) << " "<< clusterWeights_[c];
  }
  BLINK(25) << endl;
}

void GaussMixtureDensity::train(const DataSet& data) {
  uint D=data[0]->size();
  uint N=data.size();
  DBG(20) << "Starting to estimate: " << VAR(N) << " " << VAR(D) << endl;
  dim_=D;

  if(clusters_.size()==0) {
    DBG(15) << "Initialization: single Gaussian" << endl;
    init(data);
  } else {
    DBG(15) << "Initialization not necessary. Clustermodel exists already" << endl;
  }
  
  for(uint splitIter=0;splitIter<settings_.maxSplits;++splitIter) {
    DBG(10) << "Splitting " << splitIter+1 << "/" << settings_.maxSplits << ":" << clusters_.size() << " clusters." << endl;
    splitDensities();
        
    for(uint reestimationIter=0;reestimationIter<settings_.iterationsBetweenSplits;++reestimationIter) {
      DBG(15) << "Reestimating " << clusters_.size() << " densities. Iteration " << reestimationIter+1 << "/" << settings_.iterationsBetweenSplits << "." << endl;
      reestimate(data);
      poolVariances();
      deleteToSmallClusters();
      DBG(25) << "After split " << splitIter+1 << " after iteration " << reestimationIter+1 << ":" << endl << *this;
    }
  }
  
  DBG(20) << "Finished estimation" << endl;
  DBG(10) << *this;
}

double GaussMixtureDensity::logp_x(const FeatureVector& x) {
  return logp_x(x,clusterWeights_);
}

double GaussMixtureDensity::logp_x(const FeatureVector& x, const vector<double>& clusterWeights) {
  uint C=clusters_.size();
  vector<double> log_p_x_c_p_c(C);
  
  double logpmax=-numeric_limits<double>::max();
  uint cmax=C;
  
  for(uint c=0;c<C;++c) {
    log_p_x_c_p_c[c]=clusters_[c].logp_x(x)+log(clusterWeights[c]);
    if(log_p_x_c_p_c[c]>logpmax) {
      logpmax=log_p_x_c_p_c[c];
      cmax=c;
    }
  }
  
  double logsum=0.0;
  for(uint c=0;c<C;++c) {
    if(c!=cmax) {
      logsum+=exp(log_p_x_c_p_c[c]-logpmax);
    }
  }
  DBG(35) << logpmax+log1p(logsum) << endl;
  return logpmax+log1p(logsum);
}


bool GaussMixtureDensity::read(::std::istream &is) {
  StreamParser sp(is);
  string tmp;
  uint clusterNo;
  int nOfClusters;
  if(not sp.getFromLine("nOfClusters",nOfClusters)) return false;
  if(not sp.getLine(tmp)) return false; if(tmp!="settings") return false;
  is>>settings_;
  if(not sp.getVectorFromLine("clusterWeights",nOfClusters,clusterWeights_)) return false;
  clusters_.resize(nOfClusters);
  for(uint i=0;i<clusters_.size();++i) {
    if(not sp.getFromLine("cluster",clusterNo)) return false;
    if(clusterNo!=i) return false;
    clusters_[i].read(is);
  }
  return true;
}

bool GaussMixtureDensity::write(::std::ostream &os) const {
  os << "nOfClusters " << clusters_.size() << endl
     << "settings" << endl
     << settings_ 
     << "clusterWeights " << clusterWeights_ << endl;
  for(uint i=0;i<clusters_.size();++i) {
    os << "cluster " << i << endl
       << clusters_[i];
  }
  return true;
}



/* ----------------------------------------
   Implementation GaussMixtureDensityClassifier
   ---------------------------------------- */
int GaussMixtureDensityClassifier::classify(const FeatureVector &x) {
  vector<ftype> p_x_c(C_);
  return classify(x,p_x_c);
}

bool GaussMixtureDensityClassifier::load(const ::std:: string& filename) {
  igzstream is;
  is.open(filename.c_str());
  if(!is.good()) {
    ERR << "Cannot open '" << filename << "' for reading." << endl;
    return false;
  } else {
    if(not this->read(is)) {
      ERR << "Problem when reading '" << filename << "'." << endl;
      return false;
    }
    is.close();
    return true;
  }
}


bool GaussMixtureDensityClassifier::save(const ::std:: string& filename) {
  ogzstream os;
  os.open(filename.c_str());
  if(!os.good()) {
    ERR << "Cannot open '" << filename << "' for writing." << endl;
    return false;
  } else {
    if(not this->write(os)) {
      ERR << "Problem when writing '" << filename << "'." << endl;
      return false;
    }
    os.close();
    return true;
  }
}

void GaussMixtureDensityClassifier::p_c_x(const FeatureVector &x, ::std::vector<double>& p_c_x) {
  vector<double> p_x_c;
  p_c_x.resize(C_);
  classify(x,p_x_c);
  
  double sum=0.0;
  for(uint c=0;c<C_;++c) {
    p_c_x[c]=priors_[c]*p_x_c[c];
    sum+=p_c_x[c];
  }
  for(uint c=0;c<C_;++c) {
    p_c_x[c]/=sum;
  }
}

/* ----------------------------------------
   Implementation UntiedGaussianMixtureDensityClassifier
   ---------------------------------------- */

void UntiedGaussMixtureDensityClassifier::train(const DataSet& traindata) {
  // find number of classes
  C_=0;
  for(DataSet::const_iterator i=traindata.begin();i!=traindata.end();++i) {
    if((*i)->cls()>C_) C_=(*i)->cls();
  }
  ++C_;
  
  uint D=traindata[0]->size();
  
  mixtures_=vector<GaussMixtureDensity>(C_,GaussMixtureDensity(D));
  priors_=vector<double>(C_,0.0);
  //now split traindata classwise and train mixtures

  
  for(uint c=0;c<C_;++c) {
    DataSet trainc;
    for(DataSet::const_iterator i=traindata.begin();i!=traindata.end();++i) {
      if((*i)->cls()==c) {
        trainc.push_back(*i);
      }
    }
    DBG(15) << trainc.size() << " observations for class " << c << "." << endl;
    
    priors_[c]=double(trainc.size())/double(traindata.size());
    
    mixtures_[c].settings()=gmdSettings_;
    mixtures_[c].train(trainc);
    DBG(15) << "Training for class " << c << " finished." << endl;
  }
}

int UntiedGaussMixtureDensityClassifier::classify(const FeatureVector &x, vector<ftype>& log_p_x_c) {
  double max_p_x_c=-numeric_limits<double>::max();
  uint cmax=C_;
  log_p_x_c.resize(C_);
  
  DBG(25) << "log_p_x_c:";
  for(uint c=0;c<C_;++c) {
    log_p_x_c[c]=mixtures_[c].logp_x(x);
    BLINK(25) << " " << log_p_x_c[c] ;
    if(log_p_x_c[c]>max_p_x_c) {
      max_p_x_c=log_p_x_c[c];
      cmax=c;
    }
  }
  BLINK(25) << endl;
  return cmax;
}

bool UntiedGaussMixtureDensityClassifier::read(::std::istream &is) {
  StreamParser sp(is);
  uint tmp;
  
  if(not sp.getFromLine("mixtures",C_)) {return false;}
  if(not sp.getVectorFromLine("priors",C_,priors_)) {return false;}
  mixtures_.resize(C_);
  for(uint i=0;i<C_;++i) {
    DBG(15) << "reading mixture " << i << endl;
    if(not sp.getFromLine("mixture",tmp)) {return false;}
    if(tmp!=i) {ERR << "Strange: Expected mixture " << i << " but got " << tmp << "." << endl; return false;}
    is >> mixtures_[i];
  }
  return true;
}

bool UntiedGaussMixtureDensityClassifier::write(::std::ostream &os) const {
  os << "mixtures " << mixtures_.size() << endl
     << "priors " << priors_ << endl;
  for(uint i=0;i<mixtures_.size();++i) {
    os << "mixture " << i << endl
       << mixtures_[i];
  }
  return true;
}


/* ----------------------------------------
   Implementation TiedGaussMixtureDensityClassifier
   ---------------------------------------- */
void TiedGaussMixtureDensityClassifier::train(const DataSet & traindata) {
  // find number of classes
  C_=0;
  for(DataSet::const_iterator i=traindata.begin();i!=traindata.end();++i) {
    if((*i)->cls()>C_) C_=(*i)->cls();
  }
  ++C_;
  

  DBG(25) << "Estimating a joint mixture." << endl;
  uint dim=traindata[0]->size();
  mixture_=GaussMixtureDensity(dim);
  mixture_.settings()=gmdSettings_;
  mixture_.train(traindata);
  
  uint numberClusters=mixture_.clusters().size();
  
  DBG(25) << "Estimating cluster priors." << endl;
  vector<int> classCounts(C_);
  weights_.resize(C_); for(uint c=0;c<C_;++c) { weights_[c].resize(numberClusters); }
  
  for(uint n=0;n<traindata.size();++n) {
    uint cls=traindata[n]->cls();
    double logpmax=-numeric_limits<double>::max(); uint cmax=numberClusters;
    
    for(uint c=0;c<numberClusters;++c) {
      GaussDensity& gmd=mixture_.clusters()[c];
      double logpx=gmd.logp_x(*(traindata[n]));
      if(logpx>logpmax) {
        logpmax=logpx;
        cmax=c;
      }
    }
    if(cmax!=numberClusters) {
      weights_[cls][cmax]+=1.0;
      ++classCounts[cls];
    } else {
      ERR << "Observation has no closest cluster center :-( That should not happen." << endl;
    }
  }
  
  for(uint c=0;c<C_;++c) {
    for(uint k=0;k<numberClusters;++k) {
      weights_[c][k]/=double(classCounts[c]);
    }
  }
}

int TiedGaussMixtureDensityClassifier::classify(const FeatureVector &x, ::std::vector<ftype>& log_p_x_c) {
  double max_p_x_c=-numeric_limits<double>::max();
  uint cmax=C_;
  log_p_x_c.resize(C_);
  
  DBG(25) << "log_p_x_c:";
  for(uint c=0;c<C_;++c) {
    log_p_x_c[c]=mixture_.logp_x(x,weights_[c]);
    BLINK(25) << " " << log_p_x_c[c] ;
    if(log_p_x_c[c]>max_p_x_c) {
      max_p_x_c=log_p_x_c[c];
      cmax=c;
    }
  }
  BLINK(25) << endl;
  return cmax;
}

bool TiedGaussMixtureDensityClassifier::read(::std::istream &is) {
  StreamParser sp(is);
  uint tmp,nofclusters;
  string line;
  if(not sp.getFromLine("weightvectors",C_)) {return false;}
  if(not sp.getVectorFromLine("priors",C_,priors_)) {return false;}
  if(not sp.getFromLine("nofclusters",nofclusters)) {return false;}
  weights_.resize(C_);
  for(uint c=0;c<C_;++c) {
    if(not sp.getFromLine("weights",tmp)) {return false;}
    if(tmp!=c) {ERR << "Expected weight " << c << " but gut " << tmp << "." << endl; return false;}
    if(not sp.getVectorFromLine("weights",nofclusters,weights_[c])) {return false;}
  }
  if(not sp.getLine(line)) {return false;}
  if(line != "mixture") {ERR << "Expected 'mixture', got '"<< line << "'." << endl; return false;}
  return mixture_.read(is);
}

bool TiedGaussMixtureDensityClassifier::write(::std::ostream &os) const {
  os << "weightvectors " << weights_.size() << endl
     << "priors " << priors_ << endl
     << "nofclusters " << mixture_.clusters().size() << endl;
  for(uint c=0;c<weights_.size();++c) {
    os << "weights " << c << endl
       << "weights " << weights_[c] << endl;
  }
  os << "mixture" << endl
     << mixture_;
  return true;
}
