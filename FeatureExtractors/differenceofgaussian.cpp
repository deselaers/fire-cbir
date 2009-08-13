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

#include "differenceofgaussian.hpp"

const int DifferenceOfGaussian::MIN_PATCH_SIZE = 7;
const int DifferenceOfGaussian::MAX_PATCH_SIZE = 31;
const int DifferenceOfGaussian::DEFAULT_PATCH_SIZE = 11;     
const int DifferenceOfGaussian::MIN_SPACE_BETWEEN_POINTS = 5;
const double DifferenceOfGaussian::INIT_SIGMA = 0.6; 

using namespace std;

struct InterestPointsOrdering : public binary_function<InterestPoint, InterestPoint, bool> {
  bool operator()(InterestPoint ip1, InterestPoint ip2) { return fabs(ip1.saliency) > fabs(ip2.saliency); }
};

DifferenceOfGaussian::DifferenceOfGaussian(const ImageFeature& image) {
  img = image;
}

#ifdef HAVE_FFT_LIBRARY

void DifferenceOfGaussian::fftwImage(const ImageFeature& img, int& paddedSize, fftw_complex* &hsTransformed, fftw_complex* &vTransformed, 
		     bool fixedSize) {
  // Set width = height = nearest power of 2
  int width = (int) img.xsize();
  int height = (int) img.ysize();
  vector<double> complPix(3);
  double rValue, gValue, bValue;

  int maxnn, paddedHeight, paddedWidth, dim;    

  if (!fixedSize) {
    // get larger dimension
    int max = width;
    if (height > max) { 
      max = height; 
    }
    // find next largest power of 2
    maxnn = 1;
    while (maxnn < max) { 
      maxnn *= 2; 
    }
    paddedHeight = maxnn;
    paddedWidth = maxnn;
    paddedSize = maxnn;
  } else {
    maxnn = paddedSize;
    paddedHeight = paddedSize;
    paddedWidth = paddedSize;
  }

  dim = paddedHeight * paddedWidth;
  //hh = height / 2;
  //hb = width / 2;

  hsTransformed = new fftw_complex[dim];
  vTransformed = new fftw_complex[dim];     
  // Allocate memory for temp arrays
  fftw_complex* hsIn = new fftw_complex[dim];
  fftw_complex* vIn = new fftw_complex[dim];

  // Fill vectors with padded image
  //int xoffset = paddedWidth / 2 - width / 2;
  //int yoffset = paddedHeight / 2 - height / 2;
  unsigned int idx;

  ColorHSV(0, 0, 0).complexPixel(complPix);
  for(int x = 0; x < paddedWidth; x++) {
    for(int y = 0; y < paddedHeight; y++) {
      idx = y * paddedWidth + x;
      hsIn[idx].re = complPix[0];
      hsIn[idx].im = complPix[1];
      vIn[idx].re = complPix[2];
      vIn[idx].im = 0.0;
    }
  }

  // parse the data of the image and
  // set the complex arrays for the fourier
  // transformation accordingly
  for(int x = 0; x < width; x++) {
    for(int y = 0; y < height; y++) {
      if (img.zsize() == 1) {
	rValue = img(x, y, 0);
	gValue = bValue = rValue;
      } else  {
	rValue = img(x, y, 0);
	gValue = img(x, y, 1);
	bValue = img(x, y, 2);
      }
      ColorHSV(rValue, gValue, bValue).complexPixel(complPix);
      //idx = (y + yoffset) * paddedWidth + (x + xoffset);
      idx = y * paddedWidth + x;
      hsIn[idx].re = complPix[0];
      hsIn[idx].im = complPix[1];
      vIn[idx].re = complPix[2];
      vIn[idx].im = 0.0;
    }
  }
  

  // Fourier-Transformation
  fftwnd_plan plan = fftw2d_create_plan(paddedWidth, paddedHeight, FFTW_FORWARD, FFTW_ESTIMATE); 
  fftwnd_one(plan, hsIn, hsTransformed);
  fftwnd_one(plan, vIn, vTransformed);
  fftwnd_destroy_plan(plan);
  
  // we do not need the complex input data anymore
  delete[] hsIn;
  delete[] vIn;
}


ImageFeature DifferenceOfGaussian::makeFilter(double sigma) {
  DBG(20) << "making filter for sigma= " << sigma << endl;

  int dim = 1 + 2 * ((int)  (3.0 * sigma));
  if (dim % 2 == 0) {
    dim++; // we want odd dimension (to have a unique patch center pixel)
  }
  ImageFeature mask(dim, 1, 1);

  double sigma2sq = 2 * sigma * sigma;
  double normalize = 1.0 / (sqrt(2.0 * M_PI) * sigma);

  for (int i = 0 ; i < dim; i++) {
    int rel = i - (dim - 1) / 2;
    mask(i, 0, 0) = exp(-double(rel * rel) / sigma2sq) * normalize;
  }

  return mask;
}

void DifferenceOfGaussian::convolution(ImageFeature& image, double sigma) {
  ImageFeature hMask = makeFilter(sigma);
  // horizontal convolution
  convolve(image, hMask);

  // vertical convolution
  ImageFeature vMask(1, hMask.xsize(), 1);
  for (int y = 0; y < (int) vMask.ysize(); y++) {
    vMask(0, y, 0) = hMask(y, 0, 0);
  }
  convolve(image, vMask);
}

void DifferenceOfGaussian::applyFilter(fftw_complex* &hsTransformedImage, fftw_complex* &vTransformedImage,
		       fftw_complex* hsTransformedFilter, fftw_complex* vTransformedFilter, int imgSize) {
  // Apply filter
  for(int y = 0; y < imgSize; y++) {
    for(int x = 0; x < imgSize; x++) {
      fftw_complex HSttImg = hsTransformedImage[y * imgSize + x];
      fftw_complex VttImg = vTransformedImage[y * imgSize + x];
      fftw_complex HSttFilter = hsTransformedFilter[y * imgSize + x];
      fftw_complex VttFilter = vTransformedFilter[y * imgSize + x];
      // complex multiplication
      hsTransformedImage[y * imgSize + x].re = (HSttImg.re * HSttFilter.re) - (HSttImg.im * HSttFilter.im);
      hsTransformedImage[y * imgSize + x].im = (HSttImg.re * HSttFilter.im) - (HSttImg.im * HSttFilter.re);
      // complex multiplication
      vTransformedImage[y * imgSize + x].re = (VttImg.re * VttFilter.re) - (VttImg.im * VttFilter.im);
      vTransformedImage[y * imgSize + x].im = (VttImg.re * VttFilter.im) - (VttImg.im * VttFilter.re);
    }
  }
}

ImageFeature DifferenceOfGaussian::fftwBackImage(fftw_complex* &hsTransformedImage, fftw_complex* &vTransformedImage, 
				 fftw_complex* &hsTransformedFilter, fftw_complex* &vTransformedFilter, 
				 int size, int origX, int origY) {
  fftwnd_plan plan;
  fftw_complex* hsResult = new fftw_complex[size * size];
  fftw_complex* vResult = new fftw_complex[size * size];
  plan = fftw2d_create_plan(size, size, FFTW_BACKWARD, FFTW_ESTIMATE); 
  fftwnd_one(plan, hsTransformedImage, hsResult);
  fftwnd_one(plan, vTransformedImage, vResult);
  //fftwnd_one(plan, hsTransformedFilter, hsResult);
  //fftwnd_one(plan, vTransformedFilter, vResult);

  /*
  int HSno, Vno;
  if(GABOR_USE_HS) {
    Vno = d*2;
    HSno = Vno + 1;
  } else {
    Vno = d;
    HSno = 0;
  }
  */

  ImageFeature newImage(origX, origY, 1);
  int xoffset = 0;// size / 2 - origX / 2;
  int yoffset = 0;//size / 2 - origY / 2;
  int idx;
  for(int y = 0; y < origY; y++) {
    for(int x = 0; x < origX; x++) {
      idx = (y + yoffset) * size + (x + xoffset);
      newImage(x, y, 0) = sqrt(vResult[idx].re * vResult[idx].re + vResult[idx].im * vResult[idx].im);
      /*
      if (GABOR_USE_HS) {
        (*this)(x - horizontalMargin, y - verticalMargin, HSno) = sqrt(hsResult[idx].re * hsResult[idx].re + hsResult[idx].im * hsResult[idx].im);
      }
      */
    }
  }

  fftwnd_destroy_plan(plan);
  delete[] hsResult;
  delete[] vResult;
  delete[] hsTransformedImage;
  delete[] vTransformedImage;
  delete[] hsTransformedFilter;
  delete[] vTransformedFilter;
  normalize(newImage);
  return newImage;
}

ImageFeature DifferenceOfGaussian::difference(const ImageFeature& img1, const ImageFeature& img2) {
  ImageFeature newImage = ImageFeature(img1.xsize(), img1.ysize(), 1);
  for (int x = 0; x < (int) newImage.xsize(); x++) {
    for (int y = 0; y < (int) newImage.ysize(); y++) {
      newImage(x, y, 0) = img1(x, y, 0) - img2(x, y, 0);
    }
  }
  return newImage;
}
#endif

vector<uint> DifferenceOfGaussian::getAllScales() {
#ifdef HAVE_FFT_LIBRARY
  vector<uint> scales;

  double sigma = INIT_SIGMA;
  
  ImageFeature prevImage = img;
  ImageFeature newImage = img;

  // calculate number of difference images
  double tmpSigma = INIT_SIGMA;
  int num = 0;
  for (int i = 0; ; i++) {
    tmpSigma *= 1.2;
    int dim = 1 + ((int) (6.0 * tmpSigma));
    if (dim % 2 == 0) {
      dim++; 
    }
    num++;
    if (dim > MAX_PATCH_SIZE) {
      break;
    }
  }

  ImageFeature* diffImages = new ImageFeature[num];
  int* sizes = new int[num];
  
  for (int i = 0; i < num; i++) {
    /*
      int paddedImageSize;
      fftw_complex *hsTransformedImage = NULL, *vTransformedImage = NULL;
      fftwImage(img, paddedImageSize, hsTransformedImage, vTransformedImage, false);
      ImageFeature filterImage = makeFilter(i);
      fftw_complex *hsTransformedFilter = NULL, *vTransformedFilter = NULL;
      fftwImage(filterImage, paddedImageSize, hsTransformedFilter, vTransformedFilter, true);
      applyFilter(hsTransformedImage, vTransformedImage, hsTransformedFilter, vTransformedFilter, paddedImageSize);
      ImageFeature newImage = fftwBackImage(hsTransformedImage, vTransformedImage, hsTransformedFilter, vTransformedFilter, 
      paddedImageSize, (int) img.xsize(), (int) img.ysize());
      normalize(newImage);
    */
    
    int dim = 1 + ((int) (6.0 * sigma));
    if (dim % 2 == 0) {
      dim++; // we want odd dimension (to have a unique patch center pixel)
    }

    sigma *= 1.2;
    convolution(newImage, sigma);
    diffImages[i] = difference(prevImage, newImage);
    sizes[i] = dim;
    prevImage = newImage;
  }

  for (int x = 0; x < (int) img.xsize(); x++) {
    for (int y = 0; y < (int) img.ysize(); y++) {
      int extrIndex = -1;
      double extrValue = 0.0;
      for (int i = 1; i < num - 1; i++) {
	if (((diffImages[i](x, y, 0) > diffImages[i - 1](x, y, 0)) &&
	     (diffImages[i](x, y, 0) > diffImages[i + 1](x, y, 0))) ||
	    ((diffImages[i](x, y, 0) < diffImages[i - 1](x, y, 0)) &&
	     (diffImages[i](x, y, 0) < diffImages[i + 1](x, y, 0)))) {
	  if (extrIndex == -1) {
	    if ((x - sizes[i] / 2 >= 0) && (int(img.xsize()) - x > sizes[i] / 2) && (y - sizes[i] / 2 >= 0) && (int(img.ysize()) - y > sizes[i] / 2)) {
	      extrValue = fabs(diffImages[i](x, y, 0));
	      extrIndex = i;
	    }
	  } else {
	    if (fabs(diffImages[i](x, y, 0)) > extrValue) {
	      if ((x - sizes[i] / 2 >= 0) && (int(img.xsize()) - x > sizes[i] / 2) && (y - sizes[i] / 2 >= 0) && (int(img.ysize()) - y > sizes[i] / 2)) {
	      	extrValue = fabs(diffImages[i](x, y, 0));
	      	extrIndex = i;
	      }
	    }
	  }
	}
      }
      int optSize = 1;
      if (extrIndex != -1) {
	optSize = sizes[extrIndex];
      } else {
	int distToBorder = img.xsize();;
	if (x < distToBorder) {
	  distToBorder = x;
	}
	if (int(img.xsize()) - x - 1 < distToBorder) {
	  distToBorder = int(img.xsize()) - x - 1;
	}
	if (y < distToBorder) {
	  distToBorder = y;
	}
	if (int(img.ysize()) - y - 1< distToBorder) {
	  distToBorder = int(img.ysize()) - y - 1;
	}
	if (distToBorder * 2 + 1 >= DEFAULT_PATCH_SIZE) {
	  optSize = DEFAULT_PATCH_SIZE;
	} else {
	  optSize = distToBorder * 2 + 1;
	}
      }
      scales.push_back(optSize);
    }
  }

  delete[] diffImages;
  delete[] sizes;  

  return scales;
#endif

}

#ifdef HAVE_FFT_LIBRARY
bool isExtremum(const ImageFeature& imgPrev, const ImageFeature& imgCur, const ImageFeature& imgNext, int x, int y) {
  // test for local maximum
  bool extr = true;
  for (int xDiff = x - 1; xDiff <= x + 1; xDiff++) {
    for (int yDiff = y - 1; yDiff <= y + 1; yDiff++) {
      extr = (imgCur(x, y, 0) > imgPrev(xDiff, yDiff, 0)) && (imgCur(x, y, 0) > imgNext(xDiff, yDiff, 0));
      if (!extr) break;
    }
    if (!extr) break;
  }
  if (extr) {
    for (int xDiff = x - 1; xDiff <= x + 1; xDiff++) {
      for (int yDiff = y - 1; yDiff <= y + 1; yDiff++) {
	if ((xDiff != x) && (yDiff != y)) {
	  extr = (imgCur(x, y, 0) > imgCur(xDiff, yDiff, 0));
	}
	if (!extr) break;
      }
      if (!extr) break;
    }
  }
  if (extr) return true;
  // test for local minimum
  extr = true;
  for (int xDiff = x - 1; xDiff <= x + 1; xDiff++) {
    for (int yDiff = y - 1; yDiff <= y + 1; yDiff++) {
      extr = (imgCur(x, y, 0) < imgPrev(xDiff, yDiff, 0)) && (imgCur(x, y, 0) < imgNext(xDiff, yDiff, 0));
      if (!extr) break;
    }
    if (!extr) break;
  }
  if (extr) {
    for (int xDiff = x - 1; xDiff <= x + 1; xDiff++) {
      for (int yDiff = y - 1; yDiff <= y + 1; yDiff++) {
	if ((xDiff != x) && (yDiff != y)) {
	  extr = (imgCur(x, y, 0) < imgCur(xDiff, yDiff, 0));
	}
	if (!extr) break;
      }
      if (!extr) break;
    }
  }
  return extr;
}

#endif
vector<InterestPoint> DifferenceOfGaussian::getInterestPoints(int numPoints) {
#ifdef HAVE_FFT_LIBRARY

  vector<InterestPoint> interestPoints;
  double sigma = INIT_SIGMA;
  ImageFeature prevImage = img;
  ImageFeature newImage = img;

  // calculate number of difference images
  double tmpSigma = INIT_SIGMA;
  int num = 0;
  for (int i = 0; ; i++) {
    tmpSigma *= 1.2;
    int dim = 1 + ((int) round(8.0 * tmpSigma));
    if (dim % 2 == 0) {
      dim++; 
    }
    num++;
    if (dim > MAX_PATCH_SIZE) {
      break;
    }
  }

  ImageFeature* diffImages = new ImageFeature[num];
  int* sizes = new int[num];
  
  for (int i = 0; i < num; i++) {
    
    int dim = 1 + ((int) round(8.0 * sigma));
    if (dim % 2 == 0) {
      dim++; // we want odd dimension (to have a unique patch center pixel)
    }

    sigma *= 1.2;
    convolution(newImage, sigma);
    diffImages[i] = difference(prevImage, newImage);
    sizes[i] = dim;
    prevImage = newImage;
  }

  for (int x = 1; x < (int) img.xsize() - 1; x++) {
    for (int y = 1; y < (int) img.ysize() - 1; y++) {
      for (int i = 1; i < num - 1; i++) {
	if (isExtremum(diffImages[i - 1], diffImages[i], diffImages[i + 1], x, y)) {
	  if ((sizes[i] >= MIN_PATCH_SIZE) && (sizes[i] <= MAX_PATCH_SIZE)) {
	    if ((x - sizes[i] / 2 >= 0) && (int(img.xsize()) - x > sizes[i] / 2) && (y - sizes[i] / 2 >= 0) && (int(img.ysize()) - y > sizes[i] / 2)) {
	      double extrValue = fabs(diffImages[i](x, y, 0));
	      InterestPoint ip;
	      ip.x = x;
	      ip.y = y;
	      ip.scale = sizes[i];
	      ip.saliency = extrValue;
	      interestPoints.push_back(ip);
	    } else {
	      break;
	    }
	  }
	}
      }
    }
  }

  DBG(20) << "Found " << interestPoints.size() << " Difference-of-Gaussian points" << endl;
  if ((int) interestPoints.size() > numPoints) {
    // sort interest points
    sort(interestPoints.begin(), interestPoints.end(), InterestPointsOrdering());
  }

  // take only points which do not lie too close to each other
  map<int, InterestPoint> selectedPointsMap;
  map<int, InterestPoint> selectedInterestPointsMap;

  int width = (int) img.xsize();
  int i = 0;
  int selected = 0;

  DBG(20) << "Choosing the best points..." << endl;
  do {
    // get the next best interest point
    InterestPoint nextInterestPoint = interestPoints[i];
    
    DBG(30) << "Current considered point at (" << nextInterestPoint.x << "," << nextInterestPoint.y << "), saliency " << nextInterestPoint.saliency << endl;
    // check if there is another point too close to this one
    bool otherPointFound = false;
    for (int xDiff = -MIN_SPACE_BETWEEN_POINTS; xDiff <= +MIN_SPACE_BETWEEN_POINTS; xDiff++) {
      for (int yDiff = -MIN_SPACE_BETWEEN_POINTS; yDiff <= +MIN_SPACE_BETWEEN_POINTS; yDiff++) {
	int point = (nextInterestPoint.y + yDiff) * width + nextInterestPoint.x + xDiff;
	if (selectedInterestPointsMap.find(point) != selectedInterestPointsMap.end()) {
    int curSize = (selectedInterestPointsMap.find(point))->second.scale;
    if (abs(curSize - (int) nextInterestPoint.scale) < 10) {
		  double curInterest = (selectedInterestPointsMap.find(point))->second.saliency;
      otherPointFound = true;
      DBG(40) << "  point (" << (nextInterestPoint.x + xDiff) << "," << (nextInterestPoint.y + yDiff) << ") with saliency " << curInterest << " is too close" << endl;
      break;
    }
	}
      }
      if (otherPointFound) {
	break;
      }
    }
    if (!otherPointFound) {
      // there is no other interest point in the neighborhood of this interest point
      // therefore, we can safely add this point
      int point = nextInterestPoint.y * width + nextInterestPoint.x;
      selectedInterestPointsMap[point] = nextInterestPoint;
      selected++;
      DBG(30) << " taking current point, now having " << selected << endl;
    } else {
      // in the neighborhood of this interest point, there is another one with even higher saliency
      // therefore, do not add this point
      DBG(30) << " skipping current point" << endl;
    }

    i++;
  } while ((selected < numPoints) && (i < (int) interestPoints.size()));
  DBG(30) << "choice finished " << endl;

  delete[] diffImages;
  delete[] sizes;  
  interestPoints.clear();

  // collect all interest points from the map and return them
  for (map<int, InterestPoint>::const_iterator iterator = selectedInterestPointsMap.begin(); iterator != selectedInterestPointsMap.end(); iterator++) {
    interestPoints.push_back(iterator->second);
    DBG(30) << "using interest point (" << iterator->second.x << "," << iterator->second.y << ") with saliency " << iterator->second.saliency << endl;
  }
  DBG(20) << "returning " << interestPoints.size() << " interest points" << endl;

  return interestPoints;
#endif
}


