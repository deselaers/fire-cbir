#ifndef __fuzzy3histo_hpp__
#define __fuzzy3histo_hpp__

#include <vector>
#include <cmath>
#include <iostream>


class Fuzzy3BinHisto : public ::std::vector<double> {
public:
  Fuzzy3BinHisto() : ::std::vector<double>(3,0.0), counter_(0),stepsize_(1.0/3.0) {}
  
  void feed(const double v) {
    ++counter_; 
    int bin=int(v/stepsize_);
    if (bin>2) {
      bin = 2;
    }
    
    /*::std::cout << "v = " << v << ::std::endl;
    ::std::cout << "stepsize_ = " << stepsize_ << " counter_ = " << counter_ << ::std::endl;
    ::std::cout << "bin = " << bin << ::std::endl;*/
    
    double dist1=v-double(bin)*stepsize_;
    double dist2=double(bin+1)*stepsize_-v;
    
    //::std::cout << "dist1 = " << dist1 << " dist2 = " << dist2 << ::std::endl;
    
    double c=stepsize_/2;
    double h=1/c;
    double m=h/c;
    
    if(dist1<dist2) { // closer to the lower border than to the upper border
      double val=(c-dist1);val*=val;val*=m;val/=2;
      if(bin-1>=0) {
        (*this)[bin-1]+=val;
        (*this)[bin]+=(1-val);
      } else {
        (*this)[bin]+=1;
      }
    } else {
      double val=(c-dist2);val*=val;val*=m;val/=2;
      if(bin+1<3) {
        (*this)[bin+1]+=val;
        (*this)[bin]+=(1-val);
      } else {
        (*this)[bin]+=1;
      }
    }
  }
  
  const ::std::vector<double> normalized() const {
    ::std::vector<double> result(3,0.0);
    /*double sum = 0;
    for (int i = 0; i < (*this).size(); i++) {
      sum += (*this).at(i);
    }*/
    if (counter_==0) {ERR << "counter==0" << ::std::endl;}
    for (uint i = 0; i < (*this).size(); i++) {
      result[i] = (*this).at(i) / counter_;
    }
    return (result);
  }


private:
  uint counter_;
  double stepsize_;
};
#endif
