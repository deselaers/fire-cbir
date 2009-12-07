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
using namespace std;

#include "imagecomparator.hpp"
#include "diag.hpp"
#include <iostream>
#include <csignal>
#include <sstream>

using namespace std;

ImageComparator::ImageComparator() : distances_(0), cacheActive_(false)
{}

ImageComparator::ImageComparator(uint size)  :distances_(size), cacheActive_(false) {};

void ImageComparator::initialize(Database &db) {
  for(uint i=0;i<distances_.size();++i) {
    distances_[i]->initialize(db,i);
  }
}


ImageComparator::~ImageComparator() {
  for(uint i=0;i<distances_.size();++i) {
    delete distances_[i];
  }
  if(cacheActive_) closeCache();
}

void ImageComparator::start(const ImageContainer *query) {
  for(uint i=0;i<distances_.size();++i) {
    distances_[i]->start((*query)[i]);
  }
}

void ImageComparator::start(const ImageContainer *query, uint distanceID) {
  if (distanceID < distances_.size()) {
    distances_[distanceID]->start((*query)[distanceID]);
  }
}

vector<double> ImageComparator::compare(const ImageContainer *queryImage, const ImageContainer* databaseImage) {
  vector<double> result(distances_.size());
  DBG(35) << "Comparing " << queryImage->basename() << " with " << databaseImage->basename() << endl;
  for(uint i=0;i<distances_.size();++i) {
    result[i]=this->compare(queryImage, databaseImage,i);
    //result[i]=distances_[i]->distance((*queryImage)[i],(*databaseImage)[i]);
  }
  return result;
}

double ImageComparator::compare(const ImageContainer *queryImage, const ImageContainer *databaseImage, const uint distanceID) {
  if (distanceID >= distances_.size()) { ERR << "distanceID larger than distances_.size " << endl; }
  
  DBG(35) << "Comparing " << queryImage->basename() << " with " << databaseImage->basename() << " with regard to distance "<< distances_[distanceID]->name() <<endl;
  
  double d,result;
#ifdef HAVE_SQLITE3
  if(cacheActive_ and getFromCache(databaseImage->basename(), queryImage->basename(), distances_[distanceID]->name(),d)) {
    result=d;
  } else {
#endif
    result = distances_[distanceID]->distance((*queryImage)[distanceID],(*databaseImage)[distanceID]);
#ifdef HAVE_SQLITE3
    if(cacheActive_) { setInCache(databaseImage->basename(), queryImage->basename(), distances_[distanceID]->name(),result); }
  }
#endif
  return result;
}

void ImageComparator::stop() {
  for(uint i=0;i<distances_.size();++i) {
    distances_[i]->stop();
  }
}

void ImageComparator::stop(uint distanceID) {
  if(distanceID < distances_.size())
  {
    distances_[distanceID]->stop();
  }
}

void ImageComparator::distance(const uint& idx, BaseDistance* d) {
  if(distances_.size() < idx+1) {
    distances_.resize(idx+1);
  }
  if(d!=distances_[idx]) {delete distances_[idx];}

  distances_[idx]=d;
}

BaseDistance *ImageComparator::distance(const uint& idx) const {
  return distances_[idx];
}


uint ImageComparator::size() const {
  return distances_.size();
}

void ImageComparator::setCache(const std::string& cachefile) {
  cacheFileName_=cachefile;
}


void ImageComparator::openCache() {
#ifdef HAVE_SQLITE3
  int rc=sqlite3_open(cacheFileName_.c_str(), &sqliteDB_);
  if(rc) {
    ERR << "Cannot open cache database in file '" << cacheFileName_ << "'." << endl;
    ERR << "Ignoring and continuing without..." << endl;
    sqlite3_close(sqliteDB_);
  } else {
    DBG(10) << "Database successfully opened." << endl;
    cacheActive_=true;
  }
#else
  cacheActive_=false;
#warning "SQLite distance caching will not be avilable"
#endif
  }

void ImageComparator::closeCache() {
#ifdef HAVE_SQLITE3
  sqlite3_close(sqliteDB_);
#endif
}

void ImageComparator::createTable() {
  string create="create table cache(dbimg varchar, qimg varchar, distfct varchar, dist double);";
  string createindex="create index dbimgidx on cache (dbimg,qimg,dist);";
}

#ifdef HAVE_SQLITE3
static int callbackInsert(void *, int argc, char **argv, char **azColName){
  int i;
  for(i=0; i<argc; i++){ printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL"); }
  printf("\n");
  return 0;
}
#endif

#ifdef HAVE_SQLITE3
static int callbackDistance(void *dist, int argc, char **argv, char **azColName) {
  int i;

  double v;
  if(argc!=1) { DBG(10) << "Found " << argc << " distances for this combination of (dbimg,qimg,dist)" << endl; }
  for(i=0; i<argc; i++){
    if(string(azColName[i])==string("dist")) {
      std::istringstream iss(argv[i]);
      iss >> v;
    }
    *((double*)dist)=v;
    DBG(50) << VAR(v) << endl;
  }
  return 0;
}
#endif

bool ImageComparator::getFromCache(const std::string& dbimg, const std::string &qimg, const std::string& distname, double &dist) {
#ifdef HAVE_SQLITE3
  //string query="select dist from cache where dbimg=\"asdf\" and qimg=\"sdfg\" and distfct=\"idm\";";
 
  ostringstream oss;
  oss << "select dist from cache where dbimg=\"" << dbimg << "\" and qimg=\"" << qimg << "\" and distfct=\"" << distname << "\";";
 
  string query=oss.str();
  char *zErrMsg=0;
  double distFromDB=-1000;
  int rc=sqlite3_exec(sqliteDB_,query.c_str(),callbackDistance,&distFromDB,&zErrMsg);
  if( rc!=SQLITE_OK ){
    ERR << "SQL error: " << zErrMsg << endl;
    sqlite3_free(zErrMsg);
  } 
  dist=distFromDB;
  DBG(50) << "cache ->" << VAR(distFromDB) << endl;

  return dist>-1;
#else
#warning "Caching SQLite distance caching will not be avilable"
  return false;
#endif
}

bool ImageComparator::setInCache(const std::string& dbimg, const std::string &qimg, const std::string& distname, const double &dist) {
#ifdef HAVE_SQLITE3
  DBG(50) << "-> cache" << endl;
  if(cacheActive_) {  
    // string query="insert into cache (\"asdf\",\"sdfg\",\"idm\",7.8);";
    ostringstream oss;
    oss << "insert into cache values(\"" << dbimg << "\",\"" << qimg << "\",\"" << distname << "\","<<dist<<");";
    string query=oss.str(); 
    double distFromDB;
    char *zErrMsg=0;


    int rc=sqlite3_exec(sqliteDB_,query.c_str(),callbackInsert,&distFromDB,&zErrMsg);
    if( rc!=SQLITE_OK ){
      ERR << "SQL error: " << zErrMsg << endl;
      sqlite3_free(zErrMsg);
      return false;
    } 
    return true;
  
  } else {
    ERR << "Caching not active... not saved " << endl;
    return false;
  }
#else
#warning "Caching SQLite distance caching will not be avilable"
  return false;
#endif

}

void ImageComparator::tuneDistances(const vector<ImageContainer*>& posQueries, const vector<ImageContainer*>& negQueries) {
  DBG(10) << "Tuning distances" << endl;  
  for(uint i=0;i<distances_.size();++i) {
    vector<const BaseFeature*> posFeat(posQueries.size()), negFeat(negQueries.size());
    for(uint n=0;n<posQueries.size();++n) {posFeat[n]=(*posQueries[n])[i];}
    for(uint n=0;n<negQueries.size();++n) {negFeat[n]=(*negQueries[n])[i];}    
    distances_[i]->tune(posFeat,negFeat);
  }
}
