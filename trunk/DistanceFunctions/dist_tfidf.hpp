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
#ifndef __dist_tfidf_hpp__
#define __dist_tfidf_hpp__

#include "vectorfeature.hpp"
#include "sparsehistogramfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include <iostream>



class TFIDFDistance : public BaseDistance {
public:
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature);

  virtual ::std::string name() {return "tfidf";}
  
  //this initializes the term frequencies of the features contained in the query image
  virtual void start(const BaseFeature * queryFeature);

  //clears term frequencies
  virtual void stop();

  //this initializes the collection frequencies of all features by iterating over the database
  virtual void initialize(Database &db,uint distanceIndex);

protected:
  
  MapTypeDouble queryMap_;
  double queryScore_;
  MapTypeDouble queryTermFrequencies_;
  double queryLength_;
  MapTypeDouble collectionFrequencies_;
  const SparseHistogramFeature* queryFeature_;
  uint  dataBaseSize_;
  

  //compute termfrequencies for some featuremap
  MapTypeDouble getTermFrequencies(MapTypeDouble inQuery);

  //compute euclidean length of features
  double getLength(const MapTypeDouble  inMap);

  //calculate some score given a featureset
  double scoreFeatureSet(const SparseHistogramFeature * featureset);
  double scoreFeature(const SparseHistogramFeature * featureset, MapTypeDouble::iterator F);
  
};

#endif
