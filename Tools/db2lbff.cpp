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
  59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

/** 
 * A program to convert any image retrieval database in FIRE format to
 * LBFF (large binary feature file). But only vectorfeatures (e.g. vectors,
 * images, histograms are considered)
 *
 */

#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include"basefeature.hpp"
#include"vectorfeature.hpp"
#include"histogramfeature.hpp"
#include"binaryfeature.hpp"
#include"imagefeature.hpp"
#include"sparsehistogramfeature.hpp"
#include"largebinaryfeaturefile.hpp"
#include"database.hpp"
#include"getpot.hpp"

using namespace std;

void usage(){
  cout << "Usage():" << endl
       << "-h, --help                     give help" << endl
       << "-f, --filelist <file>          specify FIRE filelist to be converted" << endl
       << "                               to large binary feature file format. this option must be set" << endl
       << "-t, --targetdirectory <path>   specify where the large binary feature files should be saved" << endl
       << "                               this is optional. when not used the large binary feature files" << endl
       << "                               will be created in the directory as specified by the path included" << endl
       << "                               in the given FIRE filelist." << endl
       << endl;   
  exit(20);
}

bool calcMaxFeatSize(const Database& db,unsigned long int & ftsize,uint suffIdx ){
  bool differ = false;
  unsigned long int help = 0;
  for(uint imgIdx=1;imgIdx<db.size();imgIdx++){
    help = ((*(db[imgIdx]))[suffIdx])->calcBinarySize();
    if(ftsize!=help){
      differ = true;
      ftsize = max(ftsize,help);
    }
  }
  return differ;
}

int main(int argc, char** argv){
  GetPot cl(argc,argv);
  string path;
  string filelist;
  bool pathset = false;
	
  //parse commandline via getpot
  vector<string> ufos = cl.unidentified_options(6,"-h","--help","-f","--filelist","-t","--targetdirectory"); //6
	
  if(ufos.size()!=0) {
    for(vector<string>::const_iterator i=ufos.begin();i!=ufos.end();++i) {
      cout << "Unknown option detected: " << *i << endl;
    }
    usage();
  }
  	
  if(cl.search(2,"-h","--help")){
    usage();
  }
  	
  if(cl.search(2,"-f","--filelist")){
    filelist = cl.follow("filelist",2,"-f","--filelist");
  } else {
    ERR << "No fielist specified for converting" << endl;
    usage();
  }
  	
  if(cl.search(2,"-t","--targetdirectory")){
    path = cl.follow("~",2,"-t","--targetdirectory");
    pathset = true;
  }

  // create database and loadfilelist
  Database db;
  DBG(10) << "filelist = " << filelist << endl; 
  if(db.loadFileList(filelist) != 0){
    DBG(10) << "loaded filelist" << endl;
    db.loadFeatures();
    DBG(10) << "loaded features" << endl;
    // write lbff files
    for(uint i = 0; i< db.numberOfSuffices();++i){
      string filename;
      if(pathset){
        filename=path+"/"+db.suffix(i)+".lbff";
      } else {
        filename=db.path()+"/"+db.suffix(i)+".lbff";
      }
      // get feature type
      FeatureType ftype = db.featureType(i);
      // calculate how large in binary one feature is
      unsigned long int ftsize = 0;
      // determine if features of equal type differ in size
      bool differ = false;
      switch(ftype){
      case FT_HISTO: 
        ftsize = ((*(db[0]))[i])->calcBinarySize();
        differ = calcMaxFeatSize(db,ftsize,i);
        DBG(10) <<  "type = FT_HISTO" << endl;
        if(differ){
          DBG(10) << "features differ in size " << endl;
        }
        break;
      case FT_IMG: 
        ftsize = ((*(db[0]))[i])->calcBinarySize();
        differ = calcMaxFeatSize(db,ftsize,i);
        DBG(10) << "type = FT_IMG" << endl;
        if(differ){
          DBG(10) << "features differ in size " << endl;
        }
        break;
      case FT_VEC:
        ftsize = ((*(db[0]))[i])->calcBinarySize();
        differ = calcMaxFeatSize(db,ftsize,i);
        DBG(10) << "type = FT_VEC" << endl;
        if(differ){
          DBG(10) << "features differ in size " << endl;
        }
        break;
      case FT_SPARSEHISTO:
        ftsize = ((*(db[0]))[i])->calcBinarySize();
        differ = calcMaxFeatSize(db,ftsize,i);
        DBG(10) << "type = FT_SPARSEHISTO" << endl;
        if(differ){
          DBG(10) << "features differ in size " << endl;
        }
        break;
      case FT_BINARY:
        ftsize = ((*(db[0]))[i])->calcBinarySize();
        DBG(10) << "type = FT_BINARY" << endl;
        break;
      default:
        ERR << "unknown feature type "<<db.suffix(i)<<" in FIRE filelist present" << endl;
        exit(20);
      }
      DBG(10) << "got size " << ftsize << endl;
      /*//remove possible .gz from the filename to be created
      uint gzpos = filename.rfind(".gz");
      if (gzpos != string::npos){
        filename.erase(gzpos,3);
      }*/
      // note that the length of the filename is added to the feature size in the constructor of the largebinaryfeaturefiles
      LargeBinaryFeatureFile lbff(filename,ftype,(unsigned long int)db.size(),ftsize,differ);
      DBG(10) << "fileheader written" << endl;
      for(uint j = 0; j< db.size();++j){
        lbff.writeNext(db[j],i);
      } 
      lbff.closeWriting(); 
      DBG(10) << "features written" << endl;
      DBG(10) << "information written to file " << filename << endl;
    } 
  } else {
    ERR << "Error loading FIRE filelist; exiting" << endl;
    exit(20); 
  }
  exit(0);
  	
}

