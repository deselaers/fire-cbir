#ifndef __relationalfeaturehistogramm__hpp__
#define __relationalfeaturehistogramm__hpp__

#include <string>

#include "histogramfeature.hpp"
#include "imagefeature.hpp"
#include "kernelfunction.hpp"
#include "imagelib.hpp"


/* extract invariant feature histogram from image
   kernelFunction: description string of a kernelfunction
   samples: number of samples for the monte carlo method, none = -1
   steps: steps of the histogram (per dimension)
   R: number of steps in the rotation
 */
HistogramFeature getRelationalFeatureHistogram(const ImageFeature &image, 
                                               const ::std::string &kernelFunction);

HistogramFeature getRelationalFeatureHistogramPart(const ImageFeature &image, 
                                                   const ::std::string &kernelFunction,
                                                   const int startx,
                                                   const int stopx,
                                                   const int starty,
                                                   const int stopy);


#endif
