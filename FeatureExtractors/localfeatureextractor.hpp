#ifndef __localfeatureextractor_hpp__
#define __localfeatureextractor_hpp__
/*
  This file is part of the FIRE -- Flexible Image Retrieval System

  FIRE is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  FIRE is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU General Public License
  along with FIRE; if not, write to the Free Software Foundation,
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include <string>

#include "pca.hpp"
#include "diag.hpp"
#include "localfeatures.hpp"
#include "database.hpp"
#include "imagefeature.hpp"

struct LocalFeatureExtractorSettings {
  // extractionpoints

  uint waveletPoints;
  uint dogPoints;
  uint randomPoints;
  uint gridX;
  uint gridY;
  bool alignGrid;
  bool alignedGrid;

  // preprocessing
  bool forceGray;
  uint padding;

  //patches
  bool patches;
  ::std::vector<uint> extractionSizes;
  uint winsize;

  //sift
  bool sift;

  //color histograms
  bool histo;
  int histosteps;

  // postprocessing
  bool pca;
  bool loadpca;
  uint pcadim;
  ::std::string pcafile;
  bool pcaonly;
  ::std::string tmppath;
  ::std::string addonsuffix;

  ::std::string suffix() const;
  int dim() const;
  bool consistent() const;

  bool read(::std::istream& is) const;
  bool write(::std::ostream &os) const;
};

class LocalFeatureExtractor {
public:

  LocalFeatures& extractFromImage(const ImageFeature &img, LocalFeatures& lf);

  void extractFromDatabase(const Database &db);

  void preprocess(ImageFeature &img);
  void getExtractionPositions(const ImageFeature &img, ::std::vector<FeatureExtractionPosition>& extractionPositions);
  void extractPatches(ImageFeature &img, const ::std::vector<FeatureExtractionPosition> &positions, LocalFeatures &lf);
  void extractSIFT(ImageFeature &img, const ::std::vector<FeatureExtractionPosition> &positions, LocalFeatures &lf);
  void extractHisto(ImageFeature &img, const ::std::vector<FeatureExtractionPosition> &positions, LocalFeatures &lf);

  void pcaTransform(const PCA &pca, LocalFeatures &lf);


  LocalFeatureExtractorSettings& settings() {return settings_;}
  const LocalFeatureExtractorSettings& settings() const {return settings_;}

protected:
  LocalFeatureExtractorSettings settings_;

};

/** operator to write a GaussMixtureDensity to a stream */
inline ::std::ostream& operator<<(::std::ostream& os, const LocalFeatureExtractorSettings & src) {
  src.write(os);
  return os;
}



#endif
