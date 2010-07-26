/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

FIRE is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software Foundation, Inc.,
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/
#ifndef __imagefeature_hpp__
#define __imagefeature_hpp__

#include <vector>
#include <iostream>
#include <iomanip>
#include "vectorfeature.hpp"
#include "diag.hpp"

#ifdef HAVE_IMAGE_MAGICK
#include "Magick++.h"
#else
#warning "Functionality of imagefeature is reduced without ImageMagick: no loading, saving, displaying." << endl;
#endif

/**
   Structure for representing a pixel in the HSV color space
**/
typedef struct __HSVPixel {  double h,s,v; } HSVPixel;

/**
 * Class for processing and storing of images in the FIRE framework.
 * Many helper functions are available in the imagelib.hpp/cpp files
 *
 */
class ImageFeature : public VectorFeature {
protected:
  /// the dimensions of the image, xsize=width, ysize=height, zsize=depth (color channels etc.)
  uint xsize_, ysize_, zsize_;

  /// the data itself, here the outer vector represents the layers,
  /// the inner vector is the vector representation of one layer each.
  ::std::vector< ::std::vector<double> > data_;


#ifdef HAVE_IMAGE_MAGICK
  /// convert the image to the data structure used by image
  /// magick. This is needed for saving and displaying of images.
  Magick::Image makeMagickImage(const uint &idx1=0, const uint &idx2=1, const uint& idx3=2) const;
#else
#warning "no makeMagickImage without ImageMagick"
#endif  

public:

  /*------------------------------------------------------------
    Constructor/Destructor
    ------------------------------------------------------------*/

  /// do nothing constructor
  ImageFeature();

  /// construct empty image of given size
  ImageFeature(uint xsize, uint ysize,uint zsize);


  /// make a new image from a vector
  ImageFeature(const ::std::vector<double> &vec, uint x, uint y, uint z=1) {
    xsize_=x;
    ysize_=y;
    zsize_=z;
    data_.resize(z);
    for(uint c=0;c<z;++c) {
      data_[c].resize(x*y);
      for(uint x=0;x<xsize_;++x) {
        for(uint y=0;y<ysize_;++y) {
          data_[c][y*xsize_+x]=vec[y*xsize_*zsize_+x*zsize_+c];
        }
      }
    }
  }

  /// deconstruct image
  virtual ~ImageFeature();

  virtual ImageFeature *clone() const {return new ImageFeature(*this);}
  /*------------------------------------------------------------
    Loading/Saving
    ------------------------------------------------------------*/

  ///inherited from BaseFeature, but overwriten (here we load images,
  ///not text files), not using the read method
  virtual bool load(const ::std::string &filename);
  virtual bool load(const ::std::string& filename, bool forceGray);

  /// here an image is written in an image file format, not using the
  /// write method. This is inherited from base feature but
  /// overwritten.
  virtual void save(const ::std::string& filename) {save(filename,0,1,2);}
  
  /// here an image is written in an image file format, not using the
  /// write method
  virtual void save(const ::std::string& filename, const uint& idx1, const uint& idx2, const uint& idx3);

  /// read a plain text version of the image. This is the inverse to
  /// write. This is only used in large feature files et. al. The load
  /// method here uses other means of reading images
  virtual bool read(::std::istream & is);

  /// read a binary mode version of the image from the given stream.
  /// this is the inverse to writeBinary and is only used in large
  /// binary feature files.
  virtual bool readBinary(::std::istream & is);

  /// write a plain text version of the image. This is the inverse to
  /// read. This is only used in large feature files et. al. The save
  /// method here uses other means of reading images
  virtual void write(::std::ostream & os); 

  /// write a binary version of the image to the given stream. this is the
  /// inverse to readBinary and is only used in large binary feature files
  /// et al.
  virtual void writeBinary(::std::ostream & os);

  virtual void createFromJF(DoubleVector* imageData);
  virtual DoubleVector* toJF(int channel = 0);
  virtual void createFromPixelset(int width, int height, int* pixels);
  virtual void createFromPixelset(int width, int height, const ::std::vector<int>& pixels);

  /// display the image on the X-screen in a image magick window
  virtual void display(const uint& idx1=0, const uint& idx2=1, const uint& idx3=2) const;


  /*------------------------------------------------------------
    Access to the data
    ------------------------------------------------------------*/
  
  /// how many entries in total, derived from vectorfeature 
  virtual const uint size() const {  return xsize_*ysize_*zsize_;}
  
  /// calculate the binary size of an image 
  virtual const unsigned long int calcBinarySize() const;

  ///width
  virtual const uint xsize() const {  return xsize_;}

  ///height
  virtual const uint ysize() const {  return ysize_;}
  
  ///depth (number of (color) layers)
  virtual const uint zsize() const {  return zsize_;}

  /// get the idx-th pixel (const)
  virtual double operator[](uint idx) const {  return data_[idx%zsize_][idx/zsize_];}
  /// get the idx-th pixel
  virtual double& operator[](const uint idx) {  return data_[idx%zsize_][idx/zsize_];}

  /// return the values at position x,y for all layers in a vector
  virtual const ::std::vector<double> operator()(uint x, uint y) const;

  /// access pixel at position (x,y), layer c
  virtual double& operator()(uint x, uint y, uint c) {  return data_[c][y*xsize_+x];}

  /// access pixel at position (x,y), layer c (const)
  virtual const double& operator()(uint x, uint y, uint c) const {  return data_[c][y*xsize_+x];}

  /// append the layers from the given ImageFeature to the current one
  virtual void append(const ImageFeature& img);

  /// return the corresponding layer as imagefeature
  virtual const ImageFeature layer(const uint i) const;


  // return the value of pixel (x, y) as a HSV pixel
  HSVPixel hsvPixel(int x, int y);


  ::std::string& filename() {return filename_;}
  const ::std::string& filename() const {return filename_;}
  

  virtual ::std::vector<double> layerVector(const uint i) const {return data_[i];}

  /// change the size of the image. contents are destroyed by this
  virtual void resize(const uint width, const uint height, const uint depth) {
    xsize_=width; ysize_=height; zsize_=depth;
    
    data_.resize(zsize_);
    for(uint c=0;c<zsize_;++c) {
      data_[c].resize(xsize_*ysize_);
    }
  }

private:
  ::std::string filename_;


};

/// print an image feature to a stream
inline ::std::ostream& operator<<(::std::ostream &os, const ImageFeature& src) {
  os << "ImageFeature: X=" << src.xsize() << " Y=" << src.ysize() << " Z=" << src.zsize() << " => SIZE=" << src.size() << ::std::endl;

  for(uint y=0;y<src.ysize();++y) {
    for(uint x=0;x<src.xsize();++x) {
      for(uint c=0;c<src.zsize();++c) {
        os << ::std::setw(5) << ::std::setprecision(4) << double(src(x,y,c)) << " ";
      }
      os << "\t";
    }
    os << ::std::endl;
  }
  os << ::std::endl;
  return os;
}

#endif
