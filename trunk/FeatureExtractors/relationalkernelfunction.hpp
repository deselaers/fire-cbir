#ifndef __relationalkernelfunction_hpp__
#define __relationalkernelfunction_hpp__
#include <math.h>
#include <cmath>

#include <string>
#include <sstream>
#include "kernelfunction.hpp"

// gute werte:  X(0,0)-X(4,0)

class RelationalKernelFunction : public KernelFunction {
  
public:
  RelationalKernelFunction(uint x1, uint y1, uint x2, uint y2) 
    : x1_(x1),y1_(y1),x2_(x2),y2_(y2) {};

  double rel(const double x) const {
    double c=25;
    if(x< (-c)) {
        return 1;
    }
    else if(x>c) {
        return 0;
    }
    else {
        return 1/(2*c)*(c-x);
    }
  }
  
  virtual double calc(const InterpolatingImage &img,
                      const uint t0, const uint t1,
                      const uint r, const uint R, uint z) {
    
    double phi=2*M_PI*r/R; double sinphi=sin(phi); double cosphi=cos(phi);
    
    float pos1=((double) x1_)*cosphi   + ((double) y1_)*sinphi+t1;
    float pos3=((double) x2_)*cosphi   + ((double) y2_)*sinphi+t1;
  
    float pos2=(- ((double) x1_))*sinphi+((double) y1_)*cosphi+t0;
    float pos4=(- ((double) x2_))*sinphi+((double) y2_)*cosphi+t0;
    
    DBG(25) << "positions: " << pos1 << " " << pos2 << " " << pos3 << " " << pos4 << ::std::endl;
    
    double w1=img(pos1,pos2,z); if(w1<0) w1=0;
    double w2=img(pos3,pos4,z); if(w2<0) w2=0;
    
    return rel((w1-w2)  * 255.0);
  }  

private:
  uint x1_,y1_,x2_,y2_;
  
};

#endif
