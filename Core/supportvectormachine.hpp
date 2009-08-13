#ifndef __supportvectormachine__hh__
#define __supportvectormachine__hh__

#include <vector>
#ifdef HAVE_LIBSVM
#include "svm.h"
#else
#warning SVM will not work without libsvm being configured
#endif 
#include "diag.hpp"


/** a svm classifier using libsvm in derived from Classifier
 */
class SupportVectorMachine  {
public:
  
  /**
   * default constructor
   * kernel type:
   *        0 -- linear: u'*v
   *        1 -- polynomial: (gamma*u'*v + coef0)^degree
   *        2 -- radial basis function: exp(-gamma*|u-v|^2)
   *        3 -- sigmoid: tanh(gamma*u'*v + coef0)
   * degree: degree in kernel function
   * double
   */
  
  SupportVectorMachine(int kernel_type=2, int degree=3, double gamma=0.01, double coef0=0, double cost=4);
  
  virtual ~SupportVectorMachine();
  virtual void load(const std::string& filename);
  virtual int classify(const DoubleVector& x, DoubleVector &scores, int index_offset=0);
  virtual void train(const ::std::vector<DoubleVector>& trainVectors, const std::vector<int>& classes);
  virtual bool setParameter(::std::vector<double> &parameterList);
  virtual void setDefaultParam();

  	virtual int C() {
#ifdef HAVE_LIBSVM
  		return svm_get_nr_class(model);
#else
  		return 0;
#endif	
	  }
  
protected:
#ifdef HAVE_LIBSVM
  struct svm_parameter param; // set by the constructor
  struct svm_model *model;    // calculated in train
  struct svm_problem prob;    // set in train
  struct svm_node *x_space;   // set in train
#endif
};

#endif
