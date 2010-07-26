#include <string>
#include <sstream>
#include "relationalfeaturehistogram.hpp"
#include "kernelfunctionmaker.hpp"
#include "fuzzy3binhisto.hpp"
#include "imagelib.hpp"

using namespace std;

HistogramFeature getRelationalFeatureHistogramPart(const ImageFeature &image, 
                                                   const ::std::string &kernelFunction,
                                                   const int startx,
                                                   const int stopx,
                                                   const int starty,
                                                   const int stopy) {
  int R=36;

  ImageFeature grayImage;
  grayImage = makeGray(image);
                                                 
  KernelFunction *f=getNewKernelFunction(kernelFunction);                                                 
  InterpolatingImage iimg(grayImage); 
  
  vector<uint> stepsVec(3,8); 
  HistogramFeature result(stepsVec);  
  for (uint i=0; i < 3;++i) {
    result.min()[i] = 0.0; 
    result.max()[i] = 1.0;
  }
  result.initStepsize();
         
  DBG(30) << "startx: " << startx << " stopx: " << stopx << " starty: " << starty << " stopy: " << stopy << endl;
                       
  // check if the given parameters are within range  
  if ((startx >= stopx) || (starty >= stopy)) {
    DBG(30) << "start value greater than stop value" << endl;
    return (result);
  }
  if ((startx < 0) || (stopx > (int) grayImage.xsize()) || (starty < 0) || (stopy > (int) grayImage.ysize())) {
    DBG(30) << "start or stopvalue out of range" << endl;
    return (result);
    
  }
                     
  Fuzzy3BinHisto f3h;  
  for(uint x = (uint) startx; x < (uint) stopx; ++x) {
    for(uint y = (uint) starty; y < (uint) stopy; ++y) {
      f3h = Fuzzy3BinHisto();
      for(int r=0;r<R;++r) {
        double temp;
        temp = f->calc(iimg,x,y,r,R,0);
        DBG(50) << "r = " << r << " fcalc = " << temp << endl;
        f3h.feed(temp);
      }
      DBG(30) << "f3h = " << f3h[0] << " " << f3h[1] << " " << f3h[2] << " " << endl;
      result.feed(f3h.normalized());
    }
  }
  
  return result;
                                                                     
                                                   
}


HistogramFeature getRelationalFeatureHistogram(const ImageFeature &image, 
                                               const string &kernelFunction) {
												  
  return (getRelationalFeatureHistogramPart(image, kernelFunction, 0, image.xsize(), 0, image.ysize()));

}


