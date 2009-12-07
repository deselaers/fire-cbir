/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

FIRE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include "getpot.hpp"
#include "imagelib.hpp"
#include "colorhsv.hpp"

#ifdef HAVE_FFT_LIBRARY
extern "C" {
  #include FFTW_INCLUDE
}
#endif

struct InterestPoint {
  uint x, y, scale;
  double saliency;
};

class DifferenceOfGaussian {

public:
  
  DifferenceOfGaussian(const ImageFeature& image);

  // use this method if you want to get the a list of 'numPoints' pixel locations
  // of interest in the current image. NOTE: identical locations might be found
  // several times with different patch sizes. The maximal length of the list
  // is bounded by the number of extrema found in the guassian scale space
  // (cf. Lowe: "distinctive image features from scale-invariant keypoints")
  // having at least a magnitude of 'MAGNITUDE_THRESHOLD'.
  std::vector<InterestPoint> getInterestPoints(int numPoints);

  // use this method if you want to determine for each pixel (x, y) the optimal patch size
  // NOTE: the patch this is for each pixel at least 'MIN_PATCH_SIZE', unless the pixels
  // falls close to the border, in which case it might even be 1 (for boderline pixels)
  std::vector<uint> getAllScales();

  static const int MIN_PATCH_SIZE;
  static const int MAX_PATCH_SIZE;
  static const int DEFAULT_PATCH_SIZE;
  static const int MIN_SPACE_BETWEEN_POINTS;// = 1; 

  static const double INIT_SIGMA;// = 0.6; // initial variance for gaussian filter

private:

#ifdef HAVE_FFT_LIBRARY
  void fftwImage(const ImageFeature& img, int& paddedSize, fftw_complex* &hsTransformed, fftw_complex* &vTransformed, bool fixedSize);
  ImageFeature makeFilter(double sigma);
  void applyFilter(fftw_complex* &hsTransformedImage, fftw_complex* &vTransformedImage, fftw_complex* hsTransformedFilter, 
		   fftw_complex* vTransformedFilter, int imgSize);
  ImageFeature fftwBackImage(fftw_complex* &hsTransformedImage, fftw_complex* &vTransformedImage, 
			     fftw_complex* &hsTransformedFilter, fftw_complex* &vTransformedFilter, 
			     int size, int origX, int origY);
  ImageFeature difference(const ImageFeature& img1, const ImageFeature& img2);
  void convolution(ImageFeature& image, double sigma);


#endif
  ImageFeature img;

};


