#ifndef __gmdebtrainer_hpp__
#define __gmdebtrainer_hpp__

#include "gmd.hpp"
#include "factory.hpp"


public BaseCriterion() {
 public:
  virtual double f(const double x) const=0;
  virtual double fprime(const double x) const=0;
  
  virtual void initAlternativeClasses(const uint c, const uint C)=0;
  virtual uint nextAlternativeClass()=0;
  virtual bool alternativeClassesEnded()=0;
  
  virtual double fprime(const FeatureVector& X, const uint k, const UntiedGaussMixtureDensityClassifier& classifier); 


  // this is necessary for da factory
  template<class T>
    static BaseCriterion* create(const ::std::string&) {
    return new T(configfile, numberOfDistances);
  }

};

class MMICriterion : public BaseCriterion {
public:
  virtual double f(const double x) const { return x;}
  virtual double fprime(const double) const {return 1;}
  
  virtual void initAlternativeClasses(const uint k, const uint C) {
    c_=0;
    k_=k;
    C_=C;
  }
  
  virtual const uint nextAlternativeClass() {
    return c_++;
  }
  
  virtual const bool alternativeClassesEnded() {
    return c_>=C_;
  }
  
private:
  uint C_,c_,k_;
};


class UntiedGaussMixtureDensityEBTrainer {
public:
  UntiedGaussMixtureDensityEBTrainer(UntiedGaussMixtureDensityClassifier &classifier, const ::std::string criterion="mmi"); 
  void retrain(const DataSet& traindata);
  


protected:
  UntiedGaussMixtureDensityClassifier& classifier_;
  
  uint N_, C_, D_;
  

};
#endif
