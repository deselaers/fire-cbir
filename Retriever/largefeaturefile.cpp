#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "largefeaturefile.hpp"
#include "basefeature.hpp"
#include "histogramfeature.hpp"
#include "localfeatures.hpp"
#include "imagefeature.hpp"
#include "featureloader.hpp"
#include "vectorfeature.hpp"
#include "textfeature.hpp"

using namespace std;

LargeFeatureFile::LargeFeatureFile(::std::string filename, string suffix) :  lastSuffix_(suffix),fl_() {
  istringstream iss; 
  string keyword,line;
  reading_=true; writing_=false;

  ifs_.open(filename.c_str());
  if(!ifs_.good() || !ifs_) {
    ERR << "Cannot open LargeFeatureFile '" <<filename  << "'. Aborting." << endl;
    exit(20);
  } else {
    getline(ifs_,line);
    if(line!="FIRE_largefeaturefile") {
      cout << "Not a FIRE_largefeaturefile. Aborting" << endl;
      exit(20);
    }

    getline(ifs_,line); iss.clear(); iss.str(line); iss >> keyword;
    if(keyword=="suffix") {
      iss >> suffix_;
    } else {
      ERR << "Not a valid FIRE LargeFeatureFile, expected 'suffix', got '" << keyword << "'." << endl;
      exit(20);
    }
    
    // eat comments
    while(ifs_.peek() =='#') {
      getline(ifs_,line);
    }
  }
}
  
void LargeFeatureFile::readNext(ImageContainer *img, uint j) {
  if(reading_) {
    string line,keyword,filename;
    istringstream iss;
    getline(ifs_,line);
    iss.clear(); iss.str(line);
    iss >> keyword;
    if(keyword!="filename") {
      ERR << "Expected 'filename', got '" << keyword << "'." << endl;
      exit(20);
    } else {
      iss >> filename;
      if(filename!=img->basename()) {
        ERR << "Expected feature for file '" << img->basename() << "', got '" << filename << "'." << endl;
      }
      BaseFeature *feat=fl_.makeNewFeature(lastSuffix_);
      if(!feat) {
        ERR << "Unknown feature type: " << VAR(lastSuffix_) << " " << VAR(suffix_) << endl;
      }
      
      if(not feat->read(ifs_)) {
        ERR << "Reading feature for file '" << img->basename() << "' failed." << endl;
        exit(20);
      }
      img->operator[](j)=feat;
    }
  } else {
    ERR << "File not in reading mode" << endl;
    exit(20);
  }
}


void LargeFeatureFile::closeReading() {
  ifs_.close();
}

void LargeFeatureFile::closeWriting() {
  ofs_.close();
}

LargeFeatureFile::LargeFeatureFile(::std::string filename, ::std::string suffix, ::std::string comment) {
  istringstream iss; 
  string keyword,line;
  reading_=false; writing_=true;

  ofs_.open(filename.c_str());
  if(!ofs_.good() || !ofs_) {
    ERR << "Cannot open LargeFeatureFile '" << filename << "' for writing. Aborting!" << endl;
    exit(20);
  } else {
    ofs_ << "FIRE_largefeaturefile" << endl
        << "suffix " << suffix << endl
        << "# Comment: " << comment << endl;
  }
}

void LargeFeatureFile::writeNext(ImageContainer *img, uint j) {
  if(writing_) {
    ofs_ << "filename " << img->basename() << endl;
    img->operator[](j)->write(ofs_);
  } else {
    ERR << "File not in writing mode" << endl;
  }
}
