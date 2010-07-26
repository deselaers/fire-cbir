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
#include <algorithm>
#include "imagecontainer.hpp"
#include "vectorfeature.hpp"
#include "basefeature.hpp"
#include "imagefeature.hpp"


using namespace std;

ImageContainer::ImageContainer() : features_() {
}

ImageContainer::ImageContainer(const string& filename, const uint nof) : basename_(filename), features_(nof) {
}
 

ImageContainer::~ImageContainer() {
  for(uint i=0;i<features_.size();++i) {
    delete features_[i];
  }
}


BaseFeature*& ImageContainer::operator[](uint idx) {
  return features_[idx];
}

const BaseFeature* ImageContainer::operator[](uint idx) const {
  return features_[idx];
}


const ::std::string& ImageContainer::basename() const {
  return basename_;
}


const uint& ImageContainer::clas() const {
  return class_;
}

uint& ImageContainer::clas() {
  return class_;
}


const DescriptionSet& ImageContainer::description() const {
  return description_;
}

DescriptionSet& ImageContainer::description() {
  return description_;
}

::std::vector<double> ImageContainer::asVector() {
  ::std::vector<double> vec;
  for (uint i = 0; i < features_.size(); i++) {
    ImageFeature *imf=dynamic_cast<ImageFeature *>(features_[i]);
    if(imf) {
      DBG(25) << "Image feature!" << endl;
      for  (uint j = 0; j < imf->size(); j++) {
        vec.push_back((*imf)[j]);
      }
    } else {
      VectorFeature *vf=dynamic_cast<VectorFeature *>(features_[i]);
      if (vf) {
        const ::std::vector<double> &featurevec = vf->data();
        for  (uint j = 0; j < featurevec.size(); j++) {
          vec.push_back(featurevec[j]);
        }
      } else {
        ERR << "Features can not be vectorized" << endl;
        exit(20);
      }
    }
  }
  return vec;
}

ImageContainer ImageContainer::operator-(const ImageContainer & img) const {
  
  ImageContainer result=ImageContainer(*this);
  for (uint i = 0; i < features_.size(); i++) {
    VectorFeature *vf=dynamic_cast<VectorFeature *>(features_[i]);

    if (vf) {
      (*(result.features_[i])).operator-=(*vf);
    }
    else{
      ERR<<"Feature cant be vectorized!"<<endl;
      exit(20);
    }
  }
  return result;
}
  
struct CopyFromPointer {
  template<class T> T* operator()(const T *r) const {
    return r->clone();
  }
};                                                                                                        


ImageContainer::ImageContainer(const ImageContainer& src) :
  basename_(src.basename_), features_(src.features_.size()), class_(src.class_), description_(src.description_) {
  std::transform(src.features_.begin(), src.features_.end(), features_.begin(), CopyFromPointer());
}

/// clears the imagecontainer
void ImageContainer::clear(void) {
  for (uint j=0; j<features_.size(); ++j)
    delete features_[j];
}
