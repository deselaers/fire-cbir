#ifndef __monomial1kernelfunction_hpp__
#define __monomial1kernelfunction_hpp__
#include <math.h>
#include <cmath>
#include "kernelfunction.hpp"


class Monomial1KernelFunction : public KernelFunction {
  
public:
  Monomial1KernelFunction(int x1, int y1) 
    : x1_(x1),y1_(y1) {
    DBG(35) << "Monomial1KernelFunction(" << x1 << ","<< y1 << ")" << ::std::endl;
  };

  virtual double calc(const InterpolatingImage &img,
                      const uint t0, const uint t1,
                      const uint r, const uint R,uint z) {
    
    double phi=2*M_PI*r/R;    
    double sinphi=sin(phi);
    double cosphi=cos(phi);
    
    //    DBG(50) << "sinphi=" << sinphi << " cosphi=" << cosphi << ::std::endl;

    float pos1=((double) x1_)*sinphi   +((double) y1_)*sinphi+t0;
    float pos2=(- ((double) x1_))*sinphi+((double) y1_)*cosphi+t1;
    
    return img(pos1,pos2,z);
  }  

private:
  int x1_,y1_,x2_,y2_;
  
};

#endif
