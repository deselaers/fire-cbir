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
#ifndef __distancemaker_hpp__
#define __distancemaker_hpp__

#include <map>
#include <string>
#include "diag.hpp"
#include "basedistance.hpp"
#include "dist_euclidean.hpp"
#include "dist_crosscorrelation.hpp"
#include "dist_l1.hpp"
#include "dist_jsd.hpp"
#include "dist_kld.hpp"
#include "dist_chisquare.hpp"
#include "dist_histogramintersection.hpp"
#include "dist_reldev.hpp"
#include "dist_relbindev.hpp"
#include "dist_oneminusfidelity.hpp"
#include "dist_sqrtoneminusfidelity.hpp"
#include "dist_logtwominusfidelity.hpp"
#include "dist_arccosfidelity.hpp"
#include "dist_oneminusfidelitysquare.hpp"
#include "dist_idm.hpp"
#include "dist_binaryfeature.hpp"
#include "dist_globallocalfeaturedistance.hpp"
#include "dist_metafeature.hpp"
#include "dist_textfeature.hpp"
#include "dist_mpeg7.hpp"
#include "dist_histogrampair.hpp"
#include "dist_distfile.hpp"
#include "dist_weightedeuclidean.hpp"
#include "dist_lfhungarian.hpp"
#include "dist_facefeature.hpp"
#include "dist_lfsigemd.hpp"
#include "dist_rast.hpp"
#include "dist_tfidf.hpp"
#include "dist_bm25.hpp"
#include "dist_smart2.hpp"
#include "basefeature.hpp"
#include "dist_weightedl1.hpp"

class DistanceMaker {
private:

  /// a map to make a DistanceType from a DistanceName. This is
  /// initialized once in the constructor and used in
  /// makeDistance(const string..)
  ::std::map<const ::std::string, DistanceType> distanceNames_;

  /// a map to know which is the default distance for a given feature
  /// type. This is intialized once in the constructor and used in
  /// getDefaultDistance(...)
  ::std::map<const FeatureType, DistanceType> defaultDistances_;
public:
  DistanceMaker();

  /// given a distancename, split it into name and parameteres (by
  /// ':') and create the distance object with the given parameters
  BaseDistance* makeDistance(const ::std::string& distancename);
  
  /// given a FeatureType find out which is the default distance, and
  /// return this.
  BaseDistance* getDefaultDistance(const FeatureType featType);
  
  /// given a DistanceType and a (possibly empty) parameter set,
  /// create the appropriate distance object and return it
  BaseDistance* makeDistance(const DistanceType distType, const ::std::string &par="");
  
  /// return a list of all available distances. this is generated from
  /// the distanceNames_ map
  ::std::string availableDistances() const;
};

#endif
