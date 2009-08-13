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
#ifndef __imagelib_hpp__
#define __imagelib_hpp__


#include "imagefeature.hpp"
#include "affinetransformation.hpp"
#include <string>


void sobelv(ImageFeature &img);
void sobelh(ImageFeature &img);
void sobeldiagonal1(ImageFeature &img);
void sobeldiagonal2(ImageFeature &img);
void laplace(ImageFeature &img);
void convolve(ImageFeature &img, const ImageFeature &filter);
void fftconvolve(ImageFeature &img, const ImageFeature &filter);
void normalize(ImageFeature &img);
void gammaCorrection(ImageFeature &img, double gamma);

/// multiply to images, pixelwise
void multiply(ImageFeature &img1, const ImageFeature& img2);

/// multiply an image with a scalar, pixelwise
void multiply(ImageFeature &img, const double s);

void square(ImageFeature &img1);
void power(ImageFeature &img, int p);
void signedPower(ImageFeature &img, int p);
void absolute(ImageFeature &img);

// add plane 0 from b to plane z from a
void add(ImageFeature &a, const ImageFeature &b, const uint z);

//divide each value from a by d
void divide(ImageFeature &a, const double &d);


// get the min/max of the layer c of the image
double minimum(const ImageFeature &img, uint c=0);
double maximum(const ImageFeature &img, uint c=0);


::std::pair<uint,uint> argmin(const ImageFeature &img,uint c=0);
::std::pair<uint,uint> argmax(const ImageFeature &img,uint c=0);


// make an image grayscale
// type=0 -> maximum
// type=1 -> mean
ImageFeature makeGray(const ImageFeature &img, const uint type=0);

//rotate an image 90 degree clockwise
ImageFeature rotate90(const ImageFeature &img);

//flip an image (horizontally?)
ImageFeature flip(const ImageFeature &img); 


//this is linear gray value modification can be used for normalization.
void contrast(ImageFeature &img, const double min, const double max);

//save one layer of the image (normalized) to a gray value image.
void saveLayer(uint idx, const ImageFeature &img, const ::std::string & filename);

//local intensity normalization... originally by roberto paredes
// was used in the ikv-fibers project
//parameters: wndsize= windowsize, 
//            bnd= graylevel bound for equalization
//            c= image layer to be normalized
void equalize(ImageFeature &image, uint wndsize, double bnd, uint c, double fact_eq=0.5, uint offset=5);


ImageFeature zeropad(const ImageFeature& image, uint xpad, uint ypad);
ImageFeature affineTransformation(const ImageFeature &image, const AffineTransformation& aff);

ImageFeature localvariance(const ImageFeature &image, const uint winsize);
double localvariance(const ImageFeature &image, const uint winsize, const uint x, const uint y);

ImageFeature getPatch(const ImageFeature &image, const uint x, const uint y, const uint winsize);
ImageFeature getPatch(const ImageFeature &image, const uint left, const uint top, const uint right, const uint bottom);
void setPatch(ImageFeature& image, const uint xpos, const uint ypos, const ImageFeature &patch);

// convert the image to a vector
::std::vector<double> toVector(const ImageFeature &image);

double localEntropy(const ImageFeature &image, const uint xpos, const uint ypos, const uint zpos, const uint winsize);
ImageFeature localEntropy(const ImageFeature &image, const uint winsize);

//make everything below threshold black, everything above white
ImageFeature threshold(const ImageFeature &image, const double threshold);

// scale using simple bresenham
ImageFeature scale(const ImageFeature& image, const uint newWidth, const uint newHeight);

// scale using fft
ImageFeature fftscale(const ImageFeature& image, const uint newWidth, const uint newHeight);


void fftshift(ImageFeature &img);


ImageFeature fft(const ImageFeature &img);

// return true if the specified point is within the allowed range for the image, false otherwise
bool inImage(const ImageFeature &img, const int x, const int y, const int z=0);

// mark the specifed point with a "+" in the specified color with line length len
void cross(ImageFeature &img, const uint x, const uint y,const ::std::vector<double>& color,const uint len=2);

// mark the specifed point with a "x" in the specified color with line length len
void diagcross(ImageFeature &img, const uint x, const uint y,const ::std::vector<double>& color,const uint len=2);

//paint a box around the point (centerx,centery) with winsize and color given
void box(ImageFeature &img, const uint centerx, const uint centery, const ::std::vector<double> &color, const uint winsize);

//paint a rectangle into the image
void rect(ImageFeature &img, const uint x, const uint y, const uint width, const uint height, const ::std::vector<double> &color, const bool dash = false);

// get mean and variance for the specified layer of the image
void meanandvariance(const ImageFeature &img, double &mean, double &variance, const uint layer=0);

// mean and variance normalization for all layers of the input image
// default parameters will results in standard normal distribution
void meanAndVarianceNormalization(ImageFeature &img, const double mu=0.0, const double sigma=1.0);

// cut values outside the given range to the min and max values given
void cutoff(ImageFeature &img, const double min=0.0, const double max=1.0);

// add offset to all values in the image
void shift(ImageFeature &img, const double offset=0.0);

// apply histogram normalization to image
void histogramNormalization(ImageFeature &img);

//append vertically
// images have to have the same width
void appendVertically(ImageFeature &imga, const ImageFeature &imgb);

//append vertically
// images have to have the same height
void appendHorizontally(ImageFeature &imga, const ImageFeature &imgb);

// remove shaver many pixels at each edge of the image, leading to a smaller image
void shave(ImageFeature &imga, const uint shaver);

void RGBtoHSV(ImageFeature &img);
void HSVtoRGB(ImageFeature &img);

void DCT(ImageFeature &src, ImageFeature &dest);
void iDCT(ImageFeature &src, ImageFeature &dest);


#endif
