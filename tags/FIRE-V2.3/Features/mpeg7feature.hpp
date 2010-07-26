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
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#ifndef __mpeg7feature_hpp__
#define __mpeg7feature_hpp__

#include <iostream>
#include <string>
#include <map>
#include "basefeature.hpp"


/** here we define some magic numbers for mpeg 7 features. These magic
    numbers are used to distinguish between the different types of
    MPEG 7 features which are all handeled the same in fire. But the
    information is necessary to make the MPEG7Distance able to know
    how to run the MPEG7 software and which MPEG7 features have been
    extracted for a certain database.
*/
typedef uint MPEG7Type;
static const MPEG7Type ScalableColor     =10001;
static const MPEG7Type ColorStructure    =10002;
static const MPEG7Type DominantColor     =10003; 
static const MPEG7Type ColorLayout       =10004;
static const MPEG7Type TextureBrowsing   =10005;
static const MPEG7Type HomogeneousTexture=10006;
static const MPEG7Type LocalEdgeHistogram=10007;
static const MPEG7Type ShapeHistogram    =10008;
static const MPEG7Type RegionBasedShape  =10009;
static const MPEG7Type ContourBasedShape =10010;

/**
 * Feature for MPEG7 features.
 * 
 * Idea: We do not want to handle all MPEG7 features on our own, but
 * use the MPEG7 reference software for everything related to MPEG 7
 * features. That is:
 * 
 * a) the MPEG 7 reference software is used for feature extraction
 * b) the MPEG 7 reference software is used for feature comparison
 * 
 * Due to this, here we have only dummy features which are configured
 * by the filelist to make fire know which MPEG 7 features are
 * available.  The comparison of these features is done by the
 * MPEG7Distance (which is implemented in dist_mpeg7.*) and this
 * distance calls the appropriate programs from the MPEG 7 reference
 * software for query. Beforehand the MPEG7 software has to have been
 * used for feature extraction. A set of scripts which does this in a
 * fire compatible manner will be provided in the tools directory.
 *
 * More information on this topic can be found in the dist_mpeg7
 * implementation.
 * 
 * This implementation can also be seen as a reference if you want to
 * combine your information retrieval system with FIRE.
 */
class MPEG7Feature : public BaseFeature {
private:
  MPEG7Type mpegtype_;
  ::std::string basename_;


public:
  ///this is initialized in mpeg7feature.cpp (static)
  static ::std::map< ::std::string , MPEG7Type > subtypeMap_;
  

  /// return a reference to the basename of the image of this feature
  ::std::string & basename() {return basename_;}

  /// return a const reference to the basename of the image of this feature
  const ::std::string & basename() const {return basename_;}

  
  virtual ~MPEG7Feature() {subtypeMap_.clear();}

  virtual MPEG7Feature* clone() const {return new MPEG7Feature(*this);}
  /// constructor, here the subtype can be given, if none is given: set
  /// to zero, which is not a valid type
  MPEG7Feature(MPEG7Type t=0);
  
  /// the inverse of suffix2MPEG7Type
  ::std::string MPEG7Type2Name(const MPEG7Type type) const;
    
  /// get the mpeg7 type string of this feature
  const ::std::string  getMPEG7Type() const { return MPEG7Type2Name(mpegtype_); }
  
  /// given a suffix, determine the mpeg7 type
  MPEG7Type suffix2MPEG7Type(const ::std::string &suffix);
  
  /// set the mpeg7 type of this feature
  void setMPEG7Type(const ::std::string &suffix) { mpegtype_=suffix2MPEG7Type(suffix); }
  
  ///do nothing, this is a virtual feature, the feature is handeled by
  ///the MPEG7 reference software
  bool read(::std::istream &) {return true;}
  
  ///do nothing, this is a virtual feature, the feature is handeled by
  ///the MPEG7 reference software
  void write(::std::ostream &) {}
};

#endif
