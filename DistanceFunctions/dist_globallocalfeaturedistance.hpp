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
#ifndef __dist_globallocalfeaturedistance_hpp__
#define __dist_globallocalfeaturedistance_hpp__

#include "vectorfeature.hpp"
#include "vectorvectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include <iostream>
#include <map>
#include <string>
#include <stdlib.h>

#ifdef HAVE_KDTREE_LIBRARY
#include "knn_api.h"
#else
#warning "KDTREE LIBRARY NOT AVAILABLE. GlobalLocalFeatureDistance will not work."
#endif

class GlobalLocalFeatureDistance : public BaseDistance {
public:
  
  GlobalLocalFeatureDistance(::std::string fn="tree.kdt", uint ka=10, double eps=0.1): k(ka), epsilon(eps), filename(fn) {
#ifdef HAVE_KDTREE_LIBRARY
    kdt=(t_knn_kdtree *)malloc(1*sizeof(t_knn_kdtree));

    DBG(10) << "tree=" << fn << " k=" << k << " epsilon=" << epsilon << ::std::endl;

    DBG(10) << "Loading Tree from " << filename;
    loadTree(filename);
    BLINK(10) << "... done" << ::std::endl;

    neighbors=new t_knn[k];
    uint dim=kdt->dim;
    for (uint i=0; i<k; ++i) {
      neighbors[i].featvec=(float *)malloc(dim*sizeof(float));
      neighbors[i].labelvec=(char *)calloc(LIBKNN_LABEL_SIZE,sizeof(char));
    }
#endif
  }
  
  virtual ~GlobalLocalFeatureDistance() {
#ifdef HAVE_KDTREE_LIBRARY
    knn_kdtree_free(kdt);
    delete[] neighbors;
    free(kdt);
#endif
  }
  
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature);  
  virtual ::std::string name() {return "globallocalfeaturedistance";}
  virtual void start(const BaseFeature *);
  virtual void stop();
  
private:
  //load a tree into this distance object
  virtual void loadTree(::std::string filename);
  
  // distances are currently calculated for the query image with the name in this string:
  ::std::string calculatedFor;
  // here the hits for each image are stored (the image is not contained if no hit was obtained)
  ::std::map< ::std::string, uint > hits;
  // here we count the total number of hits.
  uint hitcounter;

  // the k for k-Nearest neighbor
  uint k;
  
  // the epsilon for the approximate search
  double epsilon;
  
  // the filename of the search tree
  ::std::string filename;
  
  // here is the tree
#ifdef HAVE_KDTREE_LIBRARY
  t_knn_kdtree *kdt;
  t_knn *neighbors;
#endif
  
};
#endif
