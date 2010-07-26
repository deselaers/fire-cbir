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

#include "sparsehistogramfeature.hpp"
#include "createsparsehisto.hpp"
#include "imagefeature.hpp"
#include "getpot.hpp"
#include "jflib.hpp"
#include <vector>
#include <sstream>
#include <iostream>

using namespace std;

CreateSparseHisto::CreateSparseHisto(uint stepsPerDimension) {
  stepsPerDimension_ = stepsPerDimension;
  width_ = -1;
  height_ = -1;
}

void CreateSparseHisto::addLayers(const ImageFeature& img) {
  if (((width_ != -1) && ((int) img.xsize() != width_)) ||
      ((height_ != -1) && ((int) img.ysize() != height_))) {
    ERR << "size of image to be added does not match size of previously added image!" << endl;
  } else if ((width_ == -1) && (height_ == -1)) {
    width_ = img.xsize();
    height_ = img.ysize();
    layers_ = ImageFeature(width_, height_, 0);
    layers_.append(img);
  } else {
    layers_.append(img);
  }
}

void CreateSparseHisto::write(string filename) {
  ifstream is(filename.c_str());
  SparseHistogramFeature shf;
  if (is == NULL) {
    DBG(10) << "Creating histogram '" << filename << "'" << endl;
    vector<uint> size(layers_.zsize());
    vector<double> max(layers_.zsize());
    vector<double> min(layers_.zsize());
    for (uint i = 0; i < layers_.zsize(); i++) {
      size[i] = stepsPerDimension_;
      max[i] = 1.0;
      min[i] = 0.0;
    }
    shf = SparseHistogramFeature(size);
    shf.max() = max;
    shf.min() = min;
    shf.initStepsize();
  } else {
    DBG(10) << "Histogram '" << filename << "' already present, remove it first." << endl;
    return;
  }
  for (int x = 0; x < width_; x++) {
    for (int y = 0; y < height_; y++) {
      SHPoint point(layers_.zsize());
      for (uint d = 0; d < layers_.zsize(); d++) {
        point[d] = layers_(x, y, d);
      }
      shf.feed(point);
    }
  }
  shf.calcRelativeFrequency();
  shf.save(filename);
}
