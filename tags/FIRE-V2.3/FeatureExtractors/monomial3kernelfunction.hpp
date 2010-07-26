#ifndef __monomial3kernelfunction_hpp__
#define __monomial3kernelfunction_hpp__
#include <math.h>
#include <cmath>

#include "kernelfunction.hpp"

class Monomial3KernelFunction : public KernelFunction {
  
public:
  Monomial3KernelFunction(int x1, int y1, int x2, int y2, int x3, int y3) 
    : x1_(x1),y1_(y1),x2_(x2),y2_(y2),x3_(x3),y3_(y3) {};

  virtual double calc(const InterpolatingImage &img,
              const uint t0, const uint t1,
              const uint r, const uint R,uint z) {

    double phi=2*M_PI*r/R;    
    double pos1=((double) x1_)*cos(phi)   + ((double) y1_)*sin(phi)+t0;
    double pos2=(- ((double) x1_))*sin(phi)+((double) y1_)*cos(phi)+t1;
    double pos3=((double) x2_)*cos(phi)   + ((double) y2_)*sin(phi)+t0;
    double pos4=(- ((double) x2_))*sin(phi)+ ((double) y2_)*cos(phi)+t1;
    double pos5=((double) x3_)*cos(phi)   + ((double) y3_)*sin(phi)+t0;
    double pos6=(- ((double) x3_))*sin(phi)+ ((double) y3_)*cos(phi)+t1;
    
    double result=exp(log(img(pos1,pos2,z)*img(pos3,pos4,z)*img(pos5,pos6,z))/3);
    return result;
  }  

private:
  int x1_,y1_,x2_,y2_,x3_,y3_;

};

#endif
