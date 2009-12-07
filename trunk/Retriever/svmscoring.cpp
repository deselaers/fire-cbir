#include "svmscoring.hpp"
#include "diag.hpp"
#include <limits>
#include <vector>
#include <sstream>

using namespace std;

double SvmScoring::getScore(const ::std::vector<double>& dists){
  
  DoubleVector scores(1,0);
  svm_.classify(dists,scores,2);
  return exp(scores[0]);

}

const ::std::string SvmScoring::settings(){

  ostringstream oss;
  oss<<"nothing here"<<std::endl;
  return oss.str();

}
