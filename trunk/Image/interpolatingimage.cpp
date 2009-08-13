#ifndef __interpolationimage__hpp__
#define __interpolationimage__hpp__

#include "interpolatingimage.hpp"

#ifdef HAVE_INTERPOL_LIBRARY
#include "coeff.h"
#include "interpol.h"
#endif

InterpolatingImage::~InterpolatingImage() {
#ifdef HAVE_INTERPOL_LIBRARY
  for(uint i=0;i<sourceImage_.zsize();++i) {
    delete[] coeff_[i];
  }
  delete[] coeff_;
#endif
}

InterpolatingImage::InterpolatingImage(const ImageFeature& image, uint splinedegree) : sourceImage_(image), splinedegree_(splinedegree) {
#ifdef HAVE_INTERPOL_LIBRARY
  coeff_=new float*[image.zsize()];
  for(uint i=0;i<image.zsize();++i) {
    coeff_[i]=new float[image.xsize()*image.ysize()];
    for(uint x=0;x<image.xsize();++x) {
      for(uint y=0;y<image.ysize();++y) {
        coeff_[i][y*image.xsize()+x]=float(image(x,y,i));
      }
    }
    SamplesToCoefficients(coeff_[i],image.xsize(),image.ysize(),splinedegree);
  }
#endif
  

}

double InterpolatingImage::operator()(double x, double y, uint z) const {
  double result;
#ifdef HAVE_INTERPOL_LIBRARY
  result=InterpolatedValue((coeff_[z]),sourceImage_.xsize(),sourceImage_.ysize(),x,y,splinedegree_);
#else
  
  int xx=int(x); int XX=(int(x)+1)%sourceImage_.xsize();
  int yy=int(y); int YY=(int(y)+1)%sourceImage_.ysize();
  
  double t=(x-xx)/(XX-xx);
  double u=(y-yy)/(YY-yy);
  
  result=
     (1-t)*(1-u)*sourceImage_(xx,yy,z)
    +t*(1-u)*sourceImage_(XX,yy,z)
    +t*u*sourceImage_(XX,YY,z)
    +(1-t)*u*sourceImage_(xx,YY,z);

#endif
  return result;
}

#endif


uint InterpolatingImage::zsize() const {
  return sourceImage_.zsize();
}
