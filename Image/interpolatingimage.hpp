#ifndef __interpolatingimage_hpp__
#define __interpolatingimage_hpp__


#include "imagefeature.hpp"

class InterpolatingImage {
public:
  InterpolatingImage(const ImageFeature& image, uint splinedegree=3);
  ~InterpolatingImage();
    
  
  double operator()(double x, double y, uint z) const;

  uint zsize() const;


private:
  ImageFeature sourceImage_;
  uint splinedegree_;
#ifdef HAVE_INTERPOL_LIBRARY
  float **coeff_;

#endif
};


#endif
