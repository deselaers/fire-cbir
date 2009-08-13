#ifndef __invariantfeaturehistogramm__hpp__
#define __invariantfeaturehistogramm__hpp__

#include <string>

#include "histogramfeature.hpp"
#include "imagefeature.hpp"
#include "kernelfunctionmaker.hpp"


/* extract invariant feature histogram from image
   kernelFunction: description string of a kernelfunction
   samples: number of samples for the monte carlo method, none = -1
   steps: steps of the histogram (per dimension)
   R: number of steps in the rotation
 */
HistogramFeature getInvariantFeatureHistogram(const ImageFeature &image, 
                                              const ::std::string &kernelFunction, 
                                              const int samples,
                                              const int steps,
                                              const int R);

HistogramFeature getInvariantFeatureHistogramPart(const ImageFeature &image, 
                                              const ::std::string &kernelFunction, 
                                              const int samples,
                                              const int steps,
                                              const int R, 
                                              const int startx,
                                              const int stopx,
                                              const int starty,
                                              const int stopy);



#endif
