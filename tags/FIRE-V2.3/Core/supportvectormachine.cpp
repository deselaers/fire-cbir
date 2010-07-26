#include "supportvectormachine.hpp"
#include <limits>

using namespace std;

SupportVectorMachine::SupportVectorMachine(int kernel_type, 
                                           int degree, 
                                           double gamma, 
                                           double coef0, 
                                           double cost) {
#ifdef HAVE_LIBSVM
	// default values
  param.svm_type = C_SVC;
  param.nu = 0.5;
  param.cache_size = 40;
  param.eps = 1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;

  // set by the user
  param.kernel_type=kernel_type;
  param.degree=degree;
  param.gamma=gamma;
  param.coef0=coef0;
  param.C=cost;
#endif
}

void SupportVectorMachine::load(const std::string& filename) {
#ifdef HAVE_LIBSVM
	model=svm_load_model(filename.c_str());
  if(model==0) {
    ERR << "Cannot load svm from "<< filename << endl;
  }
#endif
}

SupportVectorMachine::~SupportVectorMachine() {
#ifdef HAVE_LIBSVM
  svm_destroy_model(model);
  svm_destroy_param(&param);
  delete[] x_space;
#endif
}

void SupportVectorMachine::train(const ::std::vector<DoubleVector>& trainVectors, const std::vector<int>& classes) {
#ifdef HAVE_LIBSVM
  uint N=trainVectors.size();
  uint D=trainVectors[0].size();
  
  //DBG(10) << VAR(N) << " " <<VAR(D) << endl;
  
  x_space=new svm_node[N*(D+1)]; 
  prob.l=N;
  //prob.y=new double[N];// (double*)malloc(sizeof(double)*N);
  //prob.x=new struct svm_node*[N];//malloc(sizeof(struct svm_node *)*N);

  prob.y=(double*)malloc(sizeof(double)*N);
  prob.x=(struct svm_node**)malloc(sizeof(struct svm_node*)*N);
  //DBG(10) << VAR(N*(D+1)) << endl;
  int xcount=0;
  for (uint n=0; n<N; ++n) {
    prob.y[n]=classes[n];
    DBG(25) << VAR(classes[n]) << endl;
    
    //prob.x[n]=new struct svm_node[D+1];//malloc(sizeof(struct svm_node)*D+1);
    prob.x[n]=&(x_space[xcount]);
    for (uint d=0; d<trainVectors[n].size(); ++d) {
      x_space[xcount].index=d;
      x_space[xcount].value=trainVectors[n][d];
      //prob.x[n][d].index=d;
      //prob.x[n][d].value=trainVectors[n][d];
      ++xcount;
    }
    x_space[xcount++].index=-1;
    //DBG(10)<<VAR(xcount)<< endl;
    //prob.x[n][D].index=-1;
  }
  
  DBG(10) << "converted, going to train" << endl;
  
  model=svm_train(&prob, &param);
#else
  ERR << "SVM cannot be trained. FIRE was compiled without support for it" << std::endl;
#endif
}

int SupportVectorMachine::classify(const DoubleVector& x, DoubleVector& scores, int index_offset) {
#ifdef HAVE_LIBSVM
	struct svm_node* toclassify=(struct svm_node*)malloc(sizeof(struct svm_node)*x.size()+1);
  
  for(uint d=0;d<x.size();++d) {
    toclassify[d].index=d+index_offset;
    toclassify[d].value=x[d];
  }
  toclassify[x.size()].index=-1;

  int C=this->C();
  double *s=new double[C];
  
  svm_predict_values(model,toclassify,s);
  
  scores=DoubleVector(C,0.0);
  double maxScore=-std::numeric_limits<double>::max();
  int argmaxScore=-1;
  for(int c=0;c<C;++c) { 
    scores[c]=s[c]; 
    if(scores[c]>maxScore) {
      maxScore=scores[c];
      argmaxScore=c;
    }
  }
  
  free(toclassify); delete(s);
  return argmaxScore;
#else
  return 0;
#endif
}


// change the user specific settings even after instantiation
bool SupportVectorMachine::setParameter(::std::vector<double> &parameterList){
  bool retrval = false;
#ifdef HAVE_LIBSVM
  if(parameterList.size() == 5){
    // unfortunately two parameters have to be of type int
    // namely the kernel and degree of the kernel function if it is a polynomial
    param.kernel_type=(int)parameterList[0];
    param.degree=(int)parameterList[1];
    param.gamma=parameterList[2];
    param.coef0=parameterList[3];
    param.C=parameterList[4];
    retrval = true;
  }
#endif
  return retrval;
}

void SupportVectorMachine::setDefaultParam(){
#ifdef HAVE_LIBSVM
    param.kernel_type=2;
    param.degree=3;
    param.gamma=0.01;
    param.coef0=0.0;
    param.C=4.0;
#endif
}
