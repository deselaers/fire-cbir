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
#ifndef __basedistance_hpp__
#define __basedistance_hpp__

#include <string>
#include <cmath>
#include "basefeature.hpp"
#include "database.hpp"

typedef uint DistanceType;
static const DistanceType DT_BASE=1000;
static const DistanceType DT_EUCLIDEAN=1001;
static const DistanceType DT_L1=1002;
static const DistanceType DT_CHISQUARE=1003;
static const DistanceType DT_JSD=1004;
static const DistanceType DT_KLD=1005;
static const DistanceType DT_HISTOGRAMINTERSECTION=1006;
static const DistanceType DT_RELDEV=1007;
static const DistanceType DT_RELBINDEV=1008;
static const DistanceType DT_ONEMINUSFIDELITY=1009;
static const DistanceType DT_SQRTONEMINUSFIDELITY=1010;
static const DistanceType DT_LOGTWOMINUSFIDELITY=1011;
static const DistanceType DT_ARCCOSFIDELITY=1012;
static const DistanceType DT_SINFIDELITY=1013;
static const DistanceType DT_QUADRATICFORM=1014;
static const DistanceType DT_EMD=1015;
static const DistanceType DT_TIMEWARP=1016;
static const DistanceType DT_TANGENT=1017;
static const DistanceType DT_IDM=1018;
static const DistanceType DT_LFDIREKT=1019;
static const DistanceType DT_LFIDM=1020;
static const DistanceType DT_IRM=1021;
static const DistanceType DT_HUNGARIANREGIONMATCHING=1022;
static const DistanceType DT_CROSSCORRELATION=1023;
static const DistanceType DT_BINARYFEATURE=1024;
static const DistanceType DT_GLFD=1025;
static const DistanceType DT_METAFEATURE=1026;
static const DistanceType DT_MPEG7=1027;
static const DistanceType DT_HISTOPAIR=1028;
static const DistanceType DT_TEXTFEATURE=1029;
static const DistanceType DT_DISTFILE=1030;
static const DistanceType DT_LFHUNGARIAN=1031;
static const DistanceType DT_FACEFEAT=1032;
static const DistanceType DT_WEIGHTEDEUCLIDEAN=1033;
static const DistanceType DT_LFSIGEMD=1034;
static const DistanceType DT_RAST=1035;
static const DistanceType DT_TFIDF=1036;
static const DistanceType DT_BM25=1037;
static const DistanceType DT_SMART2=1038;
static const DistanceType DT_WEIGHTEDL1=1039;

class BaseDistance {
protected:
  uint distanceIndex_;
public:
  virtual ::std::string name() {return "base";}
  virtual ~BaseDistance() {}
  virtual double distance(const BaseFeature*, const BaseFeature*) {return 0.0;}
  virtual void initialize(Database &, uint) {};
  virtual void start(const BaseFeature*) {}
  virtual void stop(){}
    /** tune the parameters of the distance function given a set of positive and negative queries, e.g. after relevance feedback. */
  virtual void tune(const std::vector<const BaseFeature*>&, const std::vector<const BaseFeature*>&) {}
};

#endif
