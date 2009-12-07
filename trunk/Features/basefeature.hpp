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
#ifndef __basefeature_hpp__
#define __basefeature_hpp__
#include "diag.hpp"
#include <iostream>
#include <string>

typedef uint FeatureType;


/** Variables to have unique identifiers for the different feature
    types. These are used e.g. in FeatureLoader. Each newly defined
    feature has to be given an identifier here! */
static const FeatureType FT_BASE=100;
static const FeatureType FT_IMG=101;
static const FeatureType FT_VEC=102;
static const FeatureType FT_HISTO=103;
static const FeatureType FT_LF=104;
static const FeatureType FT_GABOR=105;  //currently not used
static const FeatureType FT_BLOBS=106;  //currently not used
static const FeatureType FT_REGIONS=107; //currently not used
static const FeatureType FT_BINARY=108; 
static const FeatureType FT_VECTORVECTOR=109; //currently not used
static const FeatureType FT_META=110;

/** this feature type is a "imaginary feature". The suffix just tells
    which MPEG-7 features are extracted and the MPEG7Distance calls
    the MPEG-7 reference software to obtain the distances. This is
    made to simplify the process of using MPEG-7 features by just
    using the reference software as an external program. */ 
static const FeatureType FT_MPEG7=111;
static const FeatureType FT_HISTOPAIR=112;
static const FeatureType FT_SPARSEHISTO=113;
static const FeatureType FT_TEXT=114;
static const FeatureType FT_TEXT_EN=116;
static const FeatureType FT_TEXT_FR=117;
static const FeatureType FT_TEXT_GE=118;

/** this feature type is basically necessary to use precalculated
 * distances.  This feature loads safed distance files and the
 * comparison is then based on these files instead of the real
 * features. This feature has to be compared using
 * DistanceFileDistance. 
 */
static const FeatureType FT_DISTFILE=115;

/** FIRE FaceFeature .face files this feature is a special inheritance
 * from vectorvectorfeature supporting multiple detected faces per
 * image.
 */
static const FeatureType FT_FACEFEAT=116;


/** featuretype for PASCAL annotation converted using
    pascalannotation.py */
static const FeatureType FT_PASCALANNOTATION=117;


static const FeatureType FT_LFSIGNATURE=118;

/** featcls files written by lfclustering
    will be used for the dist_RAST */
static const FeatureType FT_LFPOSCLSIDFEAT=119;

/** this feature type is deprecated, you can never write a file in
    this format, but the old fire has files in this format, and for
    compatibility reasons, we can read it */
static const FeatureType FT_OLDHISTO=151;



/** Base class, mainly abstract, for features. Every feature to be
 * used in this framework should be derived from this class. Here only
 * the absolute basic functionality is given. For most features it is
 * probably more usefull to derive from VectorFeature.
*/
class BaseFeature {
protected:
  FeatureType type_;
  
public:
  /** constructor setting the feature type, anyway this is not
      necessary as we cannot create an object of this type as we have
      read and write purely abstract.
   */
  BaseFeature() :type_(FT_BASE) {}

  virtual ~BaseFeature() {}

  virtual BaseFeature *clone() const=0;

  /** get the information of what type this feature is. here, the
      static const FeatureTypes from above are allowed */
  const FeatureType& type() const { return type_;}

  /** get and set the information of what type this feature is. here,
      the static const FeatureTypes from above are allowed */
  FeatureType& type() { return type_;}

  /** load a feature from a (possibly gzipped) file. This function is
      implemented here, but does nothing apart from opening a file and
      call read. This allows for only implementing the read function
      in derived classes. */
  virtual bool load(const ::std::string &filename);

  /** save a feature to a (probably gzipped) file. This function is
      implemented here, but does nothing apart from opening a file and
      call write. This allows for only implementing the write function
      in derived classes. */
  virtual void save(const ::std::string &filename);

  /** read a feature from a given istream. This should be well
      defined, as it is used in LargeFeatureFile as well, as in load,
      and possibly in other places, too. */ 
  virtual bool read(::std:: istream & is)=0;
  
  /** read a feature from a given istream in binary mode. This should be well
      defined, as it is used in LargeBinaryFeatureFile as well, as in loadBinary,
      and possibly in other places, too. */ 
  virtual bool readBinary(::std:: istream &){ERR << "Not supported for this featuretype." <<::std::endl; return true;}  // CHANGE
  
  /** save a feature to a given ostream. This should be well defined,
      as it is used in LargeFeatureFile as well, as in save and
      possibly in other places, too.*/ 
  virtual void write(::std::ostream & os)=0; 
  
  /** save a feature to a given ostream in binary mode. This should be well defined,
      as it is used in LargeFeatureFile as well, as in save and
      possibly in other places, too.*/ 
  virtual void writeBinary(::std::ostream &) {ERR << "Not supported for this featuretype." <<::std::endl;}

  /** calculate the binary size of one feature of the given type
   */
  virtual const unsigned long int calcBinarySize() const { return 0;}

  BaseFeature & operator-=(const BaseFeature &){ 
    ERR << "not implemented" << std::endl;
    return (*this);
  }
  
  // this is necessary for da factory
  template<class T>
  static BaseFeature* create() {
    return new T();
  }
  
  


};

#endif
