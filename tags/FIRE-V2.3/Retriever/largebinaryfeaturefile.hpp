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
  
#ifndef _largebinaryfeaturefile_hpp_
#define _largebinaryfeaturefile_hpp_

#include<string>
#include<fstream>
#include"diag.hpp"
#include"basefeature.hpp"
#include"imagecontainer.hpp"


/**
 * this is a class to allow for large files containing one features of
 * one type for a whole database to reduce filesystem cluttering. 
 * this file is binary and contains no encoded newlines or other spaces
 * besides as a binary file no comments are included
 *
 * files are organized as follows:
 *
 * FIRE_largebinaryfeaturefile [Type Char[28]]
 * <suffixtype> [Type unsingend int] here the type of the features has to be determined
 * <maximum length of a image file name> [Type uint]
 * <number of saved features> [Type long unsingend int]
 * <featuresize of single feature in bytes> [Type unsingend long int]
 * <different sizes> [Type bool] indicate if the features differ in binary size 
 * <filename> [Type Char[B_FILENAMESIZE]] // the default value is 50
 * feature information
 * <filename>
 * feature information
 *
 * the feature information is read/written by the readBinary/writeBinary methods
 * from the features, thus make sure that this work sufficiently
 * stable.
 */


static const unsigned long int B_HEADERSIZE = sizeof(char[28])+sizeof(uint)+sizeof(FeatureType)+2*sizeof(unsigned long int)+sizeof(bool);
const uint B_FILENAMESIZE = 50; 

class LargeBinaryFeatureFile{

private:

  FeatureType suffixtype_;
  unsigned long int numsaved_;
  unsigned long int featuresize_;
  uint filenamesize_;
  bool differ_;
 
  std::ifstream ifs_;
  std::ofstream ofs_;
 
  bool reading_, writing_;
  // if loaded is true if the header and at least the feature data of one image was read
  bool loaded_;
 
public:

  ~LargeBinaryFeatureFile() {}

  /*-------------------------------------------------------- 
    reading
    ---------------------------------------------------------*/
  // initalize a file for reading, i.e. read the header including
  // the number of saved features, the featuresize and the filename length
  LargeBinaryFeatureFile(::std::string filename);

  // read one feature into the j-th feature of the given
  // ImageContainer
  bool readNext(ImageContainer *img, uint j); 

  // close the read file
  void closeReading();
  
  /*---------------------------------------------------------
    writing
    ---------------------------------------------------------*/
  
  // initalize a file for writing. That is, write the header
  LargeBinaryFeatureFile(::std::string filename, uint suffixtype,unsigned long int numsaved, unsigned long int featuresize,bool differ=false,uint filenamelength=B_FILENAMESIZE);
  
  // write the next feature into the LargeFeature file, that is, write
  // the j-th feature from the given ImageContainer.
  void writeNext(ImageContainer *img, uint j);
  
  // close the write file
  void closeWriting();
  
  /*---------------------------------------------------------
   * helper functions
   *--------------------------------------------------------*/
   
  void seekreading(unsigned long int& seekpos){
    ifs_.clear();
    ifs_.seekg(B_HEADERSIZE+seekpos); 
  }
  
  void seekwriting(unsigned long int& seekpos){
    ofs_.clear(); 
    ofs_.seekp(B_HEADERSIZE+seekpos);
  }
  
  void setFilenamesize(uint size){ filenamesize_=size; }
  
  uint getFilenamesize() { return filenamesize_;}
  
  const FeatureType getFeatureType(){ return suffixtype_;}
  
  const bool loaded(){ return loaded_;}

  const unsigned long int getFeaturesize() { return featuresize_;}

  const unsigned long int getNumSaved() { return numsaved_; }

};
#endif
