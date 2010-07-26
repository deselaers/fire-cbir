/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

FIRE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the ilied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
#ifndef __dist_bm25_hpp__
#define __dist_bm25_hpp__

#include "vectorfeature.hpp"
#include "sparsehistogramfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include "dist_tfidf.hpp"
#include <iostream>

/* This class implements BM25 weighting(scoring) as presented by Beaulie et al, "Okapi at TREC-5" */

class BM25Distance : public TFIDFDistance {
public:

  BM25Distance(double k1,double k3, double b):k1_(k1),k3_(k3),b_(b){
  }
 
 //distance computation
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature);
  
  //--
  virtual ::std::string name() {return "BM25";}
  
  //compute collection frequencies and average document length
  virtual  void initialize(Database &db,uint distanceIndex);

  virtual void start(const BaseFeature * queryFeature);

private:

  //bm25 parameters
  double k1_,k3_,b_;

  //average document length
  double avgDL_;

  //scoring function
  double scoreFeature(const SparseHistogramFeature * featureset, MapTypeDouble::iterator F, MapTypeDouble::iterator Q);
  double scoreFeatureSet(const SparseHistogramFeature * featureset);

};
#endif
