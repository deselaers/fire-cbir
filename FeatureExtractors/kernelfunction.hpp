#ifndef __kernelfunction_hpp__
#define __kernelfunction_hpp__

#include "interpolatingimage.hpp"
#include <string>
#include <cmath>

class KernelFunction {
public:
  virtual double calc(const InterpolatingImage &img,
                      const uint t0, const uint t1,
                      const uint r, const uint R, const uint z)=0;

  virtual ~KernelFunction() {}

};

#endif
