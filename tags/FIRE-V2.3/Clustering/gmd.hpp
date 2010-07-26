#include <vector>
#include "diag.hpp"
#include "streamparser.hpp"

#ifndef __gmd_hpp__
#define __gdm_hpp__

typedef double ftype;

/** container for observations.
 * this class is a very simple feature vector class derived from a ::std::vector
 * could easily be replaced or extended
 * the only addition is class information
 */
class FeatureVector : public ::std::vector<ftype> {
public:
  FeatureVector() : ::std::vector<ftype>() {}
  FeatureVector(size_type n) : ::std::vector<ftype>(n) {}
  FeatureVector(size_type n, const ftype& t) : ::std::vector<ftype>(n,t) {}
  FeatureVector(const ::std::vector<ftype>& v) : ::std::vector<ftype>(v) {}

  uint& cls() {return cls_;}
  const uint & cls() const {return cls_;}

protected:
  uint cls_;
};

/** a vector of pointers to FeatureVector-s 
 */
typedef ::std::vector<FeatureVector*> DataSet;


/** a Gaussian density 
 *
 */
class GaussDensity {
public:
  
  /** constructor to create a D variate Gaussian density */
  GaussDensity(int D=0) : dim_(D), mean_(D), variance_(D), elements_(0), clusterNormalizationValue_(0.0) {}

  /** return the dimensionality of the Gaussian density. */
  uint &dim() {return dim_;}

  /** how many elements were used to estimate this density */
  uint &elements() {return elements_;}

  /** return the dimensionality of the Gaussian density. */
  const uint &dim() const {return dim_;}

  /** how many elements were used to estimate this density */
  const uint &elements() const {return elements_;}
  
  /** the mean vector of this Gaussian */
  ::std::vector<ftype> &mean() {return mean_;}

  /** the mean vector of this Gaussian */
  const ::std::vector<ftype> &mean() const {return mean_;}

  /** the variance vector of this Gaussian (diagonal covariance) */
  ::std::vector<ftype> &variance() {return variance_;}

  /** the variance vector of this Gaussian (diagonal covariance) */
  const ::std::vector<ftype> &variance() const {return variance_;}
  
  /** calculate the normalization Value of this Gaussian density this
      function has to be called before logp_x if the variance was
      changed */
  double calcClusterNormalizationValue();

  /** calculate log(p(x)) for this Gaussian given the feature vector x */
  double logp_x(const FeatureVector &x);

  /** make sure that the smallest variance value in this Gaussian is
      at least minVar */
  void unzeroVariance(double minVar=1E-6);

  /** estimate the parameters of this Gaussian from the given DataSet
      taking into account all observations data[i] where
      clusterMembership[i]=clusterId
      
      make sure that afterward all variance values are at least minVar
  */
  void estimate(const DataSet& data, const ::std::vector<int>& clusterMemberships, const int clusterId, const double minVar=1E-6);


  /** return the clusterNormalizationValue (\sum_d 2*M_PI*variance_d)
      of this cluster */
  double& clusterNormalizationValue() {return clusterNormalizationValue_;}

  /** return the clusterNormalizationValue (\sum_d 2*M_PI*variance_d)
      of this Gaussian */
  const double& clusterNormalizationValue() const {return clusterNormalizationValue_;}

  /** write this Gaussian to the stream os in a line-based format */
  void write(::std::ostream &os) const;

  /** read a Gaussian from the stream is in the same format that
      write(::std::ostream &os) writes*/
  void read(::std::istream &is);
  
protected:
  /// the dimensionality of the data
  uint dim_;

  /// the mean vector
  ::std::vector<ftype> mean_;

  /// the variance vector
  ::std::vector<ftype> variance_;

  ///the number of observations that were used to estimate this gaussian
  uint elements_;
  
  /// the cluster normalization value (basically the part before exp(...)
  double clusterNormalizationValue_;
};



/* ----------------------------------------
 * Settings for clustering
 * ---------------------------------------- */

/** types for settings of a GaussMixtureDensity */
typedef int PoolModeType;
typedef int DisturbModeType;
typedef int SplitModeType;

/** set the pool mode of a GaussMixtureDensity */
static const PoolModeType PMT_NOPOOLING=1;
static const PoolModeType PMT_CLUSTERPOOLING=2;
static const PoolModeType PMT_DIMENSIONPOOLING=3;
static const PoolModeType PMT_ALLPOOLING=4;

/** set the disturb mode for splitting of GaussMixtureDensity */
static const DisturbModeType DMT_VARIANCE=1;
static const DisturbModeType DMT_MEAN=2;
static const DisturbModeType DMT_MEAN2=3;
static const DisturbModeType DMT_CONST=4;

/** set the split mode (which densities are split) of a GaussMixtureDensity */
static const SplitModeType SMT_ALL=1;
static const SplitModeType SMT_LARGEST=2;
static const SplitModeType SMT_VARIANCE=3;


/** data structure to parameterize a GaussMixtureDensity 
 * the names of the variables speak for themselves
 */
struct GaussMixtureDensitySettings {

  /** how are variances pooled */
  PoolModeType poolMode;
  
  /** how are means disturbed when splitting*/
  DisturbModeType disturbMode;
  
  /** which densities are split when splitting */
  SplitModeType splitMode;

  /** how often are densities split. leads to a maximum of 2^maxSplits
      densities */
  uint maxSplits;

  /** how many observations must a cluster have to be split */
  uint minObsForSplit;

  /** how many observations must a cluster have... if it has less it is deleted */
  uint minObs;
  
  /** how many reestimation iterations are done between splits */
  uint iterationsBetweenSplits;

  /** epsilon that is used for disturbing */
  ftype epsilon;

  /** minimum allowed variance (if the variance is smaller than this
      it is set to this value */
  double minVar;

  /** write a GaussMixtureDensitySettings to a stream */
  bool write(::std::ostream &os) const {
    os << "poolMode " << poolMode << ::std::endl
       << "disturbMode " << disturbMode << ::std::endl
       << "splitMode " << splitMode << ::std::endl
       << "maxSplits " << maxSplits << ::std::endl
       << "minObsForSplit " << minObsForSplit << ::std::endl
       << "minObs " << minObs << ::std::endl
       << "iterationsBetweenSplits " << iterationsBetweenSplits << ::std::endl
       << "epsilon " << epsilon << ::std::endl
       << "minVar " << minVar << ::std::endl;
    return true;
  }

  PoolModeType String2PoolModeType(const ::std::string& s) {
    ::std::string s1=s;
    string_tolower(s1);
    if(s1=="nopooling") {
      return PMT_NOPOOLING;
    } else if (s1=="clusterpooling") {
      return PMT_CLUSTERPOOLING;
    } else if(s1=="dimensionpooling") {
      return PMT_DIMENSIONPOOLING;
    } else if(s1=="allpooling") {
      return PMT_ALLPOOLING;
    } else {
      ERR << "Invarlid poolmode '" << s1 << "'. Defaulting to nopooling." << ::std::endl;
      return PMT_NOPOOLING;
    }
  }

  DisturbModeType String2DisturbModeType(const ::std::string &s) {
    ::std::string s1=s;
    string_tolower(s1);
    if(s1=="variance") {
      return DMT_VARIANCE;
    } else if (s1=="mean") {
      return DMT_MEAN;
    } else if(s1=="mean2") {
      return DMT_MEAN2;
    } else if(s1=="const") {
      return DMT_CONST;
    } else {
      ERR << "Invarlid disturb mode '" << s1 << "'. Defaulting to variance." << ::std::endl;
      return DMT_VARIANCE;
    }
  }

  SplitModeType String2SplitModeType(const ::std::string &s) {
    ::std::string s1=s;
    string_tolower(s1);
    if(s1=="all") {
      return SMT_ALL;
    } else if (s1=="largest") {
      return SMT_LARGEST;
    } else if(s1=="variance") {
      return SMT_VARIANCE;
    } else {
      ERR << "Invarlid split mode '" << s1 << "'. Defaulting to all." << ::std::endl;
      return SMT_ALL;
    }
  }

  ::std::string PoolModeType2String(const PoolModeType p) {
    switch(p) {
    case PMT_NOPOOLING:
      return "nopooling";
    case PMT_CLUSTERPOOLING:
      return "clusterpooling";
    case PMT_DIMENSIONPOOLING:
      return "dimensionpooling";
    case PMT_ALLPOOLING:
      return "allpooling";
    default:
      ERR << "Invalid pool mode: " << p << ". Returning nopooling" << ::std::endl;
      return "nopooling";
    }
  }

  ::std::string DisturbModeType2String(const PoolModeType p) {
    switch(p) {
    case DMT_VARIANCE:
      return "variance";
    case DMT_MEAN:
      return "mean";
    case DMT_MEAN2:
      return "mean2";
    case DMT_CONST:
      return "const";
    default:
      ERR << "Invalid disturb mode: " << p << ". Returning variance" << ::std::endl;
      return "variance";
    }
  }

  ::std::string SplitModeType2String(const PoolModeType p) {
    switch(p) {
    case SMT_ALL:
      return "all";
    case SMT_LARGEST:
      return "largest";
    case SMT_VARIANCE:
      return "variance";
    default:
      ERR << "Invalid split mode: " << p << ". Returning all" << ::std::endl;
      return "all";
    }
  }


  /** read a GaussMixtureDensitySettings from a stream, returns true
      if sucessful, false otherwise */
  bool read(::std::istream &is) {
    StreamParser sp(is);
    if(not sp.getFromLine("poolMode", poolMode)) { return false;}
    if(not sp.getFromLine("disturbMode", disturbMode)) { return false;}
    if(not sp.getFromLine("splitMode", splitMode)) { return false;} 
    if(not sp.getFromLine("maxSplits",maxSplits)) { return false;} 
    if(not sp.getFromLine("minObsForSplit",minObsForSplit)) { return false;} 
    if(not sp.getFromLine("minObs",minObs)) { return false;} 
    if(not sp.getFromLine("iterationsBetweenSplits",iterationsBetweenSplits)) { return false;}
    if(not sp.getFromLine("epsilon",epsilon)) { return false;}
    if(not sp.getFromLine("minVar",minVar)) { return false;}
    return true;
  }
};


/** a GaussMixtureDensity with
 * \item EM training with maximum approximation
 * \item diagonal covariance matrices
 */
class GaussMixtureDensity {
public:  

  /** constructor to for dim-dimensional data */
  GaussMixtureDensity(uint dim=0) : dim_(dim),  clusterWeights_(0), clusters_(0), settings_() {};
  
  /** train a Gaussian mixture density */
  void train(const DataSet& data);
  
  /** init a Gaussian mixture density (single Gaussian) */
  void init(const DataSet &data);

  /** reestimation from given training data */
  void reestimate(const DataSet& data);
  
  /** split all densities according to settings */
  void splitDensities();


  const uint numberOfClusters() const {return clusters_.size();}

  /** split the density with number clusterID according to the
      settings */
  GaussDensity splitDensity(uint clusterId);

  /** pool the variances according to the settings */
  void poolVariances();
  
  /** calculate log( p( x | this GMD ) ) 
   * 
   * this is done numerically stable completely in the logarithmic
   * domain. Use the given clusterWeights for this (this is necessary
   * for tied mixture densities
   */
  double logp_x(const FeatureVector& x, const ::std::vector<double>& clusterWeights);

  /** calculate log( p( x | this GMD ) ) 
   * 
   * this is done numerically stable completely in the logarithmic
   * domain. Use the estimate clusterWeights for this (this is necessary
   * for untied mixture densities
   */
  double logp_x(const FeatureVector& x);

  /** return the vector of densities */
  const ::std::vector<GaussDensity>& clusters() const  {return clusters_;}

  /** return the vector of densities */
  ::std::vector<GaussDensity>& clusters() {return clusters_;}

  /** return the vector of density-priors */
  const ::std::vector<ftype>& clusterWeights() const {return clusterWeights_;}

  /** return the vector of density-priors */
  ::std::vector<ftype>& clusterWeights() {return clusterWeights_;}

  /** return the settings object */
  GaussMixtureDensitySettings & settings() {return settings_;}

  /** return the settings object */
  const GaussMixtureDensitySettings & settings() const {return settings_;}
  
  /** delete mixtures which are too small according to the settings */
  void deleteToSmallClusters();
  
  /** write the complete GaussMixtureDensity to a stream */
  bool write(::std::ostream &os) const;
  
  /** read the complete GaussMixtureDensity from a stream in the same
      format that write writes */
  bool read(::std::istream &is);

  
protected:
  /** dimensionality of the data handled */
  uint dim_;

  /** priors for the clusters */
  ::std::vector<ftype> clusterWeights_;
  
  /** the densities themsevles */
  ::std::vector<GaussDensity> clusters_;
  
  /** the settings object */
  GaussMixtureDensitySettings settings_;
};


/** abstract class for a classifier based on GaussMixtureDensities */
class GaussMixtureDensityClassifier {
public:
  GaussMixtureDensityClassifier(): gmdSettings_(),  C_(0) {}
  
  virtual ~GaussMixtureDensityClassifier() {}
  
  /** train the classifier */
  virtual void train(const DataSet& traindata)=0;
  
  /** classify an observation */
  virtual int classify(const FeatureVector &x);
  
  /** classify an observation and return the log( p(x|c) ) for all
      classes c in the argument vector*/
  virtual int classify(const FeatureVector &x, ::std::vector<ftype>& log_p_x_c)=0;
  
  /** calculate
   *                      p(c) p(x|c)
   *     p(c | x) = --------------------------
   *                   \sum_c p(c) p(x|c)
    */
  virtual void p_c_x(const FeatureVector &x, ::std::vector<double>& p_c_x);
  

  /** load a GaussMixtureDensityClassifier from a file */
  virtual bool load(const ::std:: string&);
  /** save a GaussMixtureDensityClassifier to a file */
  virtual bool save(const ::std:: string&);

  /** read a GaussMixtureDensityClassifier from a stream */
  virtual bool read(::std::istream &is)=0;
  
  /** write a GaussMixtureDensityClassifier to a stream */
  virtual bool write(::std::ostream &os) const =0 ;

  /** return the settings for the embedded GaussMixtureDensities*/
  GaussMixtureDensitySettings & settings() {return gmdSettings_;}
  /** return the settings for the embedded GaussMixtureDensityies */
  const GaussMixtureDensitySettings & settings() const {return gmdSettings_;}
  
  const uint C() const {return C_;}
  uint & C() {return C_;}

protected:
  /** the settings objects */
  GaussMixtureDensitySettings gmdSettings_;

  /** a vector for the prior probabilities of the classes */
  ::std::vector<double> priors_;
  
  /** the number of classes */
  uint C_;
};

/** a classifier for untied Gaussian mixture densites */
class UntiedGaussMixtureDensityClassifier : public GaussMixtureDensityClassifier {
  //friend UntiedGaussMixtureDensityEBTrainer;

public:

  UntiedGaussMixtureDensityClassifier() : GaussMixtureDensityClassifier(), mixtures_() {}

  /** training */
  void train(const DataSet& traindata);

  /** classification and returning of the log ( p(x | c) ) for each class c */
  virtual int classify(const FeatureVector &x, ::std::vector<ftype>& log_p_x_c);
  
  /** read from a stream */
  virtual bool read(::std::istream &is);
  
  /** write to a stream */
  virtual bool write(::std::ostream &os) const ;
  
protected:
  /** a vector containing a GaussMixtureDensity per class */
  ::std::vector<GaussMixtureDensity> mixtures_;
};


/** a (not yet implemented classifier for tied Gaussian mixture densities */
class TiedGaussMixtureDensityClassifier : public GaussMixtureDensityClassifier {
  void train(const DataSet& traindata);
  virtual int classify(const FeatureVector &x, ::std::vector<ftype>& log_p_x_c);
  
  virtual bool read(::std::istream &is);
  virtual bool write(::std::ostream &os) const;
  
protected:
  GaussMixtureDensity mixture_;

  ::std::vector< ::std::vector<ftype> > weights_;
};

/** operators to write the objects of this file to streams */

/** operator to write a vector<ftype> to a stream */
inline ::std::ostream& operator<<(::std::ostream& os, const ::std::vector<ftype> & src) {
  for(::std::vector<ftype>::const_iterator i=src.begin();i!=src.end();++i) os << *i << " " ;
  return os;
}

/** operator to write a vector<int> to a stream */
inline ::std::ostream& operator<<(::std::ostream& os, const ::std::vector<int> & src) {
  for(::std::vector<int>::const_iterator i=src.begin();i!=src.end();++i) os << *i << " " ;
  return os;
}

/** operator to write a GaussDensity to a stream */
inline ::std::ostream& operator<<(::std::ostream& os, const GaussDensity & src) {
  src.write(os);
  return os;
}

/** operator to read a GaussDensity to a stream */
inline ::std::istream& operator>>(::std::istream& is,  GaussDensity & dst) {
  dst.read(is);
  return is;
}


/** operator to write a GaussMixtureDensitySettings object to a
    stream */
inline ::std::ostream& operator<<(::std::ostream& os, const GaussMixtureDensitySettings & src) {
  src.write(os);
  return os;
}

/** operator to read a GaussMixtureDensitySettings from a stream */
inline ::std::istream& operator>>(::std::istream& is,  GaussMixtureDensitySettings & dst) {
  dst.read(is);
  return is;
}

/** operator to write a GaussMixtureDensity to a stream */
inline ::std::ostream& operator<<(::std::ostream& os, const GaussMixtureDensity & src) {
  src.write(os);
  return os;
}

/** operator to read a GaussMixtureDensity from a stream */
inline ::std::istream& operator>>(::std::istream& is,  GaussMixtureDensity & dst) {
  dst.read(is);
  return is;
}

#endif
