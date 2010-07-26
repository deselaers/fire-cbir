#ifndef __monomial2kernelfunction_hpp__
#define __monomial2kernelfunction_hpp__
#include <math.h>
#include <cmath>
#include "kernelfunction.hpp"
#include <iostream>

class Monomial2KernelFunction : public KernelFunction {
  
public:
  Monomial2KernelFunction(int x1, int y1, int x2, int y2) 
    : x1_(x1),y1_(y1),x2_(x2),y2_(y2) {
    DBG(35) << "Monomial2KernelFunction(" << x1 << ","<< y1 << "," << x2 << "," << y2 << ")" << ::std::endl;
  };

  virtual double calc(const InterpolatingImage &img,
                      const uint t0, const uint t1,
                      const uint r, const uint R, uint z) {

    double phi=2*M_PI*r/R; double sinphi=sin(phi); double cosphi=cos(phi);
    
    float pos1=((double) x1_)*cosphi   + ((double) y1_)*sinphi+t0;
    float pos2=(- ((double) x1_))*sinphi+((double) y1_)*cosphi+t1;
    float pos3=((double) x2_)*cosphi   + ((double) y2_)*sinphi+t0;
    float pos4=(- ((double) x2_))*sinphi+((double) y2_)*cosphi+t1;
    
    double w1=img(pos1,pos2,z); if(w1<0) w1=0;
    double w2=img(pos3,pos4,z); if(w2<0) w2=0;
    
    //DBG(100) << pos1 << " " << pos2  << " " << pos3  << " " << pos4  << " -> " << w1 << " " << w2 << " -> " << sqrt(w1*w2)<<::std::endl;
    
    return sqrt(w1*w2);
  }  
  
private:
  int x1_,y1_,x2_,y2_;
  
};

#endif
