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
  
#include <string>
#include <iostream>
#include <fstream>
#include "largebinaryfeaturefile.hpp"
#include "imagecontainer.hpp"
#include "basefeature.hpp"
#include "vectorfeature.hpp"
#include "binaryfeature.hpp"
#include "imagefeature.hpp"
#include "histogramfeature.hpp"
#include "sparsehistogramfeature.hpp"
#include "histogrampairfeature.hpp"

using namespace std;

LargeBinaryFeatureFile::LargeBinaryFeatureFile(string filename){
  ifs_.open(filename.c_str(),ios::in | ios::binary);
  if(!ifs_.good() || !ifs_){
    ERR << "Cannot open LargeBinaryFeatureFile '" <<filename  << "'. Aborting." << endl;
    exit(20);
  } else {
    reading_ = true;
    writing_ = false;
    loaded_ = false;
    char magic[28];
    
    // read the MagicNumber
    ifs_.read(magic,sizeof(char[28]));
    if(strcmp(magic,"FIRE_largebinaryfeaturefile")!= 0){
      ERR << filename << " is not a FIRE_largefeaturefile. Aborting" << endl;
      exit(20);
    }
    // read Feature Type
    ifs_.read((char*)&suffixtype_,sizeof(FeatureType));
    DBG(10) << VAR(suffixtype_) << endl;
    // read the filename length
    ifs_.read((char*)&filenamesize_,sizeof(uint));
    DBG(10) << VAR(filenamesize_) << endl;
    // read how many features are saved in this file
    ifs_.read((char*)&numsaved_,sizeof(unsigned long int));
    DBG(10) << VAR(numsaved_) << endl;
    // read how much space in the file one features requires
    ifs_.read((char*)&featuresize_,sizeof(unsigned long int));
    DBG(10) << VAR(featuresize_) << endl;
    // read the flag indicating whether the features differ in size
    ifs_.read((char*)&differ_,sizeof(bool));
    DBG(10) << VAR(differ_) << endl;
  } // end else
}

LargeBinaryFeatureFile::LargeBinaryFeatureFile(string filename, uint suffixtype, unsigned long int numsaved, unsigned long int featuresize,bool differ,uint filenamelength){
  
  ofs_.open(filename.c_str(),ios::out|ios::binary);
  if(!ofs_.good()){
    ERR << "Cannot open LargeBinaryFeatureFile '" << filename << "' for writing. Aborting!" << endl;
    exit(20);
  } else {
    reading_ = false;
    // write MagicNumber
    char magic[28] = "FIRE_largebinaryfeaturefile";
    ofs_.write(magic,sizeof(char[28]));
    // write featuretype
    ofs_.write((char*)&suffixtype,sizeof(FeatureType));
    suffixtype_=suffixtype;
    // write the length of the filenames that is beeing used
    ofs_.write((char*)&filenamelength,sizeof(uint));
    filenamesize_=filenamelength;
    // write how many features will be saved in this file
    ofs_.write((char*)&numsaved,sizeof(unsigned long int));
    numsaved_=numsaved;
    // write how much space one feature needs
    // this must be adapted by the size of the filenames
    featuresize+= filenamelength;
    ofs_.write((char*)&featuresize,sizeof(unsigned long int)); 
    featuresize_=featuresize;
    // write if the features differ in size
    ofs_.write((char*)&differ,sizeof(bool));
    differ_=differ;
    writing_ = true;
    loaded_= false;
  }  
} 

bool LargeBinaryFeatureFile::readNext(ImageContainer *img, uint j){
  if(reading_){
    // read the file name
    char file[filenamesize_];
        
    ifs_.read(file,sizeof(char[filenamesize_]));
    
    //DBG(10) << VAR(file) << endl;
    //BG(10) << VAR(ifs_.good()) << " " << VAR(ifs_.bad()) << " " << VAR(ifs_.eof()) << " " << VAR(ifs_.fail()) << endl;
  	
    // remove the possible inserted ';' for length consistency
    // note that we assume that the symbol ';' occurs only at the end and
    // is not a normal part of the filename
    string filename(file);
    
    // hier gibt es einen fehler, wenn:
    // - der dateiname ein Seminkolon enthält
    // - oder: der dateiname genau 50 zeichen lang ist.
    
    uint p;
    for(p=filename.size()-1; p>=0 and filename[p]==';' ;--p);
    filename.erase(p+1);

    // now read the feature
    if(filename!=img->basename()){ 
      ERR << "Expected feature for file '" << img->basename() << "', got '" << filename << "'." << endl;
      exit(20);
    } else {
      BaseFeature* feat = NULL;
      switch(suffixtype_){
      case FT_VEC: feat = new VectorFeature(); break;
      case FT_BINARY: feat = new BinaryFeature(); break; 
        //  case FT_DISTFILE: feat = new DistanceFileFeature(); break;
        //  case FT_FACEFEAT: feat = new FaceFeature(); break;
      case FT_IMG: feat = new ImageFeature(); break; 
        //  case FT_OLDHISTO: // depreciated
      case FT_HISTO: feat = new HistogramFeature(); break;
      case FT_SPARSEHISTO: feat = new SparseHistogramFeature(); break; 
        //    case FT_HISTOPAIR: feat = new HistogramPairFeature(); break;
        //	  case FT_LF: feat = new LocalFeatures(); break;
        //	  case FT_LFPOSCLSIDFEAT: feat = new LFPositionClusterIdFeature(); break;
        //	  case FT_LFSIGNATURE: feat = new LFSignatureFeature(); break;
        //	  case FT_MPEG7: feat = new MPEG7Feature(); break;
        //	  case FT_META: feat = new MetaFeature(); break;
        //	  case FT_TEXT_EN:
        //	  case FT_TEXT_FR:
        //	  case FT_TEXT_GE:
        //	  case FT_TEXT: feat = new TextFeature(); break;
        //	  case FT_PASCALANNOTATION: feat = new PascalAnnotationFeature(); break;
        // not yet implemented
      case FT_GABOR:
      case FT_BLOBS:
      case FT_REGIONS:
        ERR << "This feature type is not yet implemented: " << suffixtype_ << endl;
        return false;
      default:
        ERR << "This feature type is not yet known: " << suffixtype_ << endl;
        return false;
      } 
          
      bool readBool = feat->readBinary(ifs_); 
      if(!readBool){
        return false;
      }
      img->operator[](j)=feat;
      // if the features differ in size the padded zeros have to be skipped
      if(differ_){
        long unsigned int local = img->operator[](j)->calcBinarySize();
        long unsigned int currPos = ifs_.tellg();
        ifs_.seekg(currPos+(featuresize_-local-B_FILENAMESIZE));
      }
    }
  } else {
    ERR << "file not in reading mode" << endl;
    return false;
  }
  loaded_=true;
  return true;
}

void LargeBinaryFeatureFile::writeNext(ImageContainer *img, uint j){
  if(writing_){
    // write the filename
    string filen = img->basename();
    // stretch the filenname to a length of filenamesize_
    uint fsize = (filenamesize_)-filen.size();
    filen.append(fsize-1,';'); 
    ofs_.write(filen.c_str(),sizeof(char[filenamesize_]));
    // write the feature
    DBG(105) << "write Data for file: " << img->basename() << endl;
    img->operator[](j)->writeBinary(ofs_);
    
    // padd the feature by zeros if the feauters of this type differ in size
    long unsigned int localSize = 0;
    bool fill = false;
    if(differ_){
      localSize = img->operator[](j)->calcBinarySize();
      for(long unsigned int i=0;i<(featuresize_-localSize-B_FILENAMESIZE);i++){
	ofs_.write((char*)&fill,sizeof(bool));
      }
    }
     
  } else {
    ERR << "File not in writing mode" << endl;
  }
}

void LargeBinaryFeatureFile::closeReading(){
  ifs_.close();
}

void LargeBinaryFeatureFile::closeWriting(){
  ofs_.close();
}

  
