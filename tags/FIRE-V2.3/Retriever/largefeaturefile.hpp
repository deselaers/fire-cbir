/* This file is part of the FIRE -- Flexible Image Retrieval System
  
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
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#ifndef __largefeaturefile__hpp__
#define __largefeaturefile__hpp__

#include <string>
#include "diag.hpp"
#include "imagecontainer.hpp"
#include "vectorfeature.hpp"
#include "gzstream.hpp"
#include "featureloader.hpp"

/**
 * this is a class to allow for large files containing one features of
 * one type for a whole database to reduce filesystem cluttering.
 *
 * files are organized as follows:
 *
 * FIRE_largefeaturefile 
 * suffix <suffix>    here the type of the features has to be determined
 * filename <filename>
 * feature information
 * filename <filename>
 * feature information
 *
 * the feature information is read/written by the read/write methods
 * from the features, thus make sure that this work sufficiently
 * stable.
 */

class LargeFeatureFile {
private:
  ::std::string filename_,suffix_,lastSuffix_;
  igzstream ifs_;
  ogzstream ofs_;
  FeatureType type_;
  FeatureLoader fl_;
  
  bool reading_, writing_;

public:

  /*---------------------------------------------------------
    reading
    ---------------------------------------------------------*/
  // initialize a file for reading, that is read the header
  LargeFeatureFile(::std::string filename, ::std::string relevantSuffix);
  
  // read one feature into the j-th feature of the given
  // ImageContainer
  void readNext(ImageContainer *img, uint j);

  // close the read file
  void closeReading();

  /*---------------------------------------------------------
    writing
    ---------------------------------------------------------*/

  // initalize a file for writing. That is, write the header
  LargeFeatureFile(::std::string filename, ::std::string suffix, ::std::string comment);
  
  // write the next feature into the LargeFeature file, that is, write
  // the j-th feature from the given ImageContainer.
  void writeNext(ImageContainer *img, uint j);
  
  // close the write file
  void closeWriting();
  

};
  


#endif
