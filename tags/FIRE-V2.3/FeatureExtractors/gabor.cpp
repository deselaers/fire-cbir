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

#include "gabor.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
//#include <values.h>
#include "colorhsv.hpp"
#include "imagelib.hpp"

using namespace std;

Gabor::Gabor(const ImageFeature &img) {
  originalImage = img;
  if ((img.zsize() != 1) && (img.zsize() != 3)) {
    cout << "Gabor (constructor): original ImageFeature must either be a RGB or a gray-value image !" << endl;
  }
}

Gabor::~Gabor() {
}

// public method called from outside to start the gabor feature extraction
// controls the full process of extracting numPha_ * numFreq_ gabor features.
// returns no result. results are available via the getImage(...) function
// If the horizontal resp. vertical margin is set to a value > 0,
// the original image is enlarged by a black border of the given size,
// to decrease nonlinear effects at the image's edges. Anyway, this
// only makes sense if the image itself has a block background color.
void Gabor::calculate(int numPha_, int numFreq_, uint horizonalMargin_, uint verticalMargin_) {
  numPha = numPha_;
  numFreq = numFreq_;
  horizontalMargin = horizonalMargin_;
  verticalMargin = verticalMargin_;
  
  // adjust underlying imagefeature so that it can hold the 
  // gabor features to be calculated later
  xsize_ = originalImage.xsize();
  ysize_ = originalImage.ysize();
  zsize_ = numPha * numFreq * (GABOR_USE_HS ? 2 : 1);
  data_.resize(zsize_);
  for (uint c = 0; c < zsize_; ++c) {
    data_[c].resize(xsize_ * ysize_);
  }
  
#ifdef HAVE_FFT_LIBRARY
  
  // fourier transformed image data, set by the init(...) method
  fftw_complex *vTransformed = NULL, *hsTransformed = NULL;
  
  init(hsTransformed, vTransformed);
  extractGabor(hsTransformed, vTransformed);
  
  // final cleanup of allocated memory
  delete[] hsTransformed;
  delete[] vTransformed;
  
#else
  
  DBG(10) << "This class is compiled without the fftw library. Gabor feature extraction will not work." << endl;

#endif
}


// returns the 'layer'-th gabor feature. for a frequency f and phase p
// the feature number is calculated as "f * # phases + p"
// the result is undefined if the call to this method is not
// preceded by a call to "calculate()"
ImageFeature Gabor::getImage(int layer) {
  ImageFeature img(this->xsize(), this->ysize(), (GABOR_USE_HS? 2: 1));
  for (uint x = 0; x < this->xsize(); x++)
    for (uint y = 0; y < this->ysize(); y++) {
      if (GABOR_USE_HS) {
        img(x, y, 0) = (*this)(x, y, 2 * layer);
        img(x, y, 1) = (*this)(x, y, 2 * layer + 1);
      } else {
        img(x, y, 0) = (*this)(x, y, 2 * layer);
      }
    }
  if (GABOR_USE_HS) {
    return makeGray(img);
  } else {
    return img;
  }
}


#ifdef HAVE_FFT_LIBRARY

// initializes the gabor feature extraction
// sets up the data structure and performs the fourier transformation
void Gabor::init(fftw_complex* &hsTransformed, fftw_complex* &vTransformed) {
  // Set width = height = nearest power of 2
  uint width = originalImage.xsize() + 2 * horizontalMargin;
  uint height = originalImage.ysize() + 2 * verticalMargin;
  vector<double> complPix(3);
  double rValue, gValue, bValue;
    
  // get larger dimension
  uint max = width;
  if (height > max) { max = height; }
  
  // find next largest power of 2
  maxnn = 1;
  while (maxnn < max) { maxnn *= 2; }
  
  paddedHeight = maxnn;
  paddedWidth = maxnn;
  dim = paddedHeight * paddedWidth;

  hh = height / 2;
  hb = width / 2;

  hsTransformed = new fftw_complex[dim];
  vTransformed = new fftw_complex[dim];     
  // Allocate memory for temp arrays
  fftw_complex* hsIn = new fftw_complex[dim];
  fftw_complex* vIn = new fftw_complex[dim];

  // Fill vectors with padded image
  int xoffset=paddedWidth/2-width/2;
  int yoffset=paddedHeight/2-height/2;
  unsigned int idx;

  ColorHSV(0, 0, 0).complexPixel(complPix);
  for(uint x = 0; x < paddedWidth; x++) {
    for(uint y = 0; y < paddedHeight; y++) {
      idx = y * paddedWidth + x;
      hsIn[idx].re = complPix[0];
      hsIn[idx].im = complPix[1];
      vIn[idx].re = complPix[2];
      vIn[idx].im = 0.0;
    }
  }

  // parse the data of the underlying image and
  // set the complex arrays for the fourier
  // transformation accordingly
  for(uint x = 0; x < width; x++) {
    for(uint y = 0; y < height; y++) {
      // consider the margins 
      if ((x < horizontalMargin) || (x >= width - horizontalMargin) ||
          (y < verticalMargin) || (y >= height - verticalMargin)) {
        rValue = bValue = gValue = 0;
      } else {
        if (originalImage.zsize() == 1) {
          rValue = originalImage(x-horizontalMargin, y-verticalMargin, 0);
          gValue = bValue = rValue;
        } else  {
          rValue = originalImage(x-horizontalMargin, y-verticalMargin, 0);
          gValue = originalImage(x-horizontalMargin, y-verticalMargin, 1);
          bValue = originalImage(x-horizontalMargin, y-verticalMargin, 2);
        }
      }
      ColorHSV(rValue, gValue, bValue).complexPixel(complPix);
      idx = (y + yoffset) * paddedWidth + (x + xoffset);
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

// calculates a filter for the given phase and frequency
void Gabor::calcFilter(complex<double>* &curFilter, unsigned int curPha, int curFreq) {
  // Erase old filter
  for (uint i = 0; i < dim; i++) {
    curFilter[i] = complex<double>(0.0, 0.0);
  }
  
  double alpha = sqrt(log(2.0)/2.0);
  double phasendiff = M_PI / ((double)GABOR_NUM_PHASES);
  double sigmafactor = alpha / (tan(phasendiff/2.0)*M_PI);
  double ratio = (3.0 * alpha) / (M_PI * sigmafactor);
  curFreq += GABOR_FREQ_START;

  double freq = sqrt(2.0)*pow(2.0,double(curFreq));

  double pha = (double(curPha)) * (M_PI / (double) numPha);

  double u0 = cos(pha)*freq;
  double v0 = sin(pha)*freq;
  double U = u0 * cos(pha) + v0 * sin(pha);
  double V = - u0 * sin(pha) + v0 * cos(pha);
  double deviation = sigmafactor / freq;

  double s = -2*M_PI*M_PI*deviation*deviation;

  double u,v,delta1,delta2;
  
  for (int y = -hh ; y< hh; y++) {
    for (int x = -hb; x<hb; x++) {
      u = x * cos(pha) + y * sin(pha);
      v = -x* sin(pha) + y * cos(pha);
      delta1= ratio*(u - U);
      delta2= v - V;
      
      int idxx, idxy;
      idxx=(x+paddedWidth)%paddedWidth;
      idxy=(y+paddedHeight)%paddedHeight;
      curFilter[(idxy*paddedWidth)+idxx] = exp(std::complex<double>(s*(delta1*delta1+delta2*delta2),0.0));
    }
  }

  curFilter[hh * paddedWidth + hb] = 0;
}


// extract the gabor features
// for each phase and each frequency, a filter is computed,
// which is applied to the fourier-transformed data.
// afterwards the inverse fourier transformation is performed,
// and the result data is copied back into the image
void Gabor::extractGabor(fftw_complex* hsTransformed, fftw_complex* vTransformed) {
  
  fftw_complex* hsResult = new fftw_complex[dim];
  fftw_complex* vResult = new fftw_complex[dim];
  fftw_complex* hsSave = new fftw_complex[dim];
  fftw_complex* vSave = new fftw_complex[dim];
  complex<double>* curFilter = new complex<double>[dim];

  fftwnd_plan plan;
  plan = fftw2d_create_plan(paddedWidth, paddedHeight, FFTW_BACKWARD, FFTW_ESTIMATE); 
  
  for(uint aktFreq=0; aktFreq<numFreq; aktFreq++) {
    for(unsigned int aktPha=0; aktPha<numPha; aktPha++)  {

      // Calculate filter
      calcFilter(curFilter, aktPha, aktFreq);

      // copy transformed data
      memcpy(&(hsSave[0]), &(hsTransformed[0]), sizeof(fftw_complex) * dim);
      memcpy(&(vSave[0]), &(vTransformed[0]), sizeof(fftw_complex) * dim);

      // Apply filter
      for(unsigned int y=0;y<paddedHeight;++y) {
        for(unsigned int x=0;x<paddedWidth;++x) {
	  
          fftw_complex HStt = hsTransformed[y * maxnn + x];
          fftw_complex Vtt = vTransformed[y * maxnn + x];
          fftw_complex filter;
          filter.re = curFilter[y * paddedWidth + x].real();
          filter.im = curFilter[y * paddedWidth + x].imag();
	  
          // complex multiplication of 'HStt' with 'filter'
          hsTransformed[y * maxnn + x].re = (HStt.re * filter.re) - (HStt.im * filter.im);
          hsTransformed[y * maxnn + x].im = (HStt.re * filter.im) + (HStt.im * filter.re); // here was a minus instead of the plus
          // complex multiplication of 'Vtt' with 'filter'
          vTransformed[y * maxnn + x].re = (Vtt.re * filter.re) - (Vtt.im * filter.im);
          vTransformed[y * maxnn + x].im = (Vtt.re * filter.im) + (Vtt.im * filter.re); // here was a minus instead of the plus
        }
      }

      // Fourier-Transformation backwards
      fftwnd_one(plan, hsTransformed, hsResult);
      fftwnd_one(plan, vTransformed, vResult);
      
      // Save calculated feature
      copyToResult(hsResult, vResult, aktFreq * numPha + aktPha);

      // restore transformed data, needed for next extraction step
      memcpy(&(hsTransformed[0]), &(hsSave[0]), sizeof(fftw_complex) * dim);
      memcpy(&(vTransformed[0]), &(vSave[0]), sizeof(fftw_complex) * dim);

    }
  }

  fftwnd_destroy_plan(plan);

  // cleanup allocated memory
  delete[] curFilter;
  delete[] hsSave;
  delete[] vSave;
  delete[] hsResult;
  delete[] vResult;
}

// copies the output store in the hsResult and vResult arrays into
// the underlying imagefeature object for the given dimension
void Gabor::copyToResult(fftw_complex* hsResult, fftw_complex* vResult, int d)
{
  int HSno, Vno;
  if(GABOR_USE_HS) {
    Vno = d*2;
    HSno = Vno + 1;
  } else {
    Vno = d;
    HSno = 0;
  }
  
  int width = originalImage.xsize() + 2 * horizontalMargin;
  int height = originalImage.ysize() + 2 * verticalMargin;
  int xoffset=maxnn/2-width/2;
  int yoffset=maxnn/2-height/2;
  int idx;
  
  for(uint y = verticalMargin; y < height - verticalMargin; ++y) {
    for(uint x = horizontalMargin; x < width - horizontalMargin; ++x) {
      // consider the margins
      idx = (y + yoffset) * maxnn + (x + xoffset);
      (*this)(x - horizontalMargin, y - verticalMargin, Vno) = sqrt(vResult[idx].re * vResult[idx].re + vResult[idx].im * vResult[idx].im);
      if (GABOR_USE_HS) {
        (*this)(x - horizontalMargin, y - verticalMargin, HSno) = sqrt(hsResult[idx].re * hsResult[idx].re + hsResult[idx].im * hsResult[idx].im);

      }
    }
  }
}

#endif
