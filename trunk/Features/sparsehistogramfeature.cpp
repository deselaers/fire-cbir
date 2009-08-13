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

#include "sparsehistogramfeature.hpp"
#include <map>
#include <vector>
#include <sstream>
#include <math.h>

using namespace std;


SparseHistogramFeature::SparseHistogramFeature() {
  type_ = FT_SPARSEHISTO;
  counter_  = 0;
  dimensions_ = 0;
}

SparseHistogramFeature::SparseHistogramFeature(vector<uint> steps) {
  type_ = FT_SPARSEHISTO;
  steps_ = steps;
  counter_ = 0;
  dimensions_ = steps.size();
}

SparseHistogramFeature::SparseHistogramFeature(std::vector<uint> steps,const ::std::vector<double>& min, const ::std::vector<double>& max) {
  type_ = FT_SPARSEHISTO;
  steps_ = steps;
  counter_ = 0;
  dimensions_ = steps.size();
  min_=min;
  max_=max;
  this->initStepsize();
}


SparseHistogramFeature::SparseHistogramFeature(uint dimensions, uint steps) {
  type_ = FT_SPARSEHISTO;
  steps_ = vector<uint>(dimensions, steps);
  counter_ = 0;
  dimensions_ = dimensions;
}

SparseHistogramFeature::~SparseHistogramFeature() {
  data_.clear();
};


const uint SparseHistogramFeature::size() const {
  return data_.size();
}

const uint SparseHistogramFeature::length() const{
  return length_;
}

const unsigned long int SparseHistogramFeature::calcBinarySize() const{
  unsigned long int bsize = 0;
  // the data sizes
  // dim,  bins size and counter are of type uint
  bsize+= 3*sizeof(uint);
  // the min and max vectors
  bsize+= 2*(unsigned long int)min_.size() *(unsigned long int)sizeof(double);
  // the bins hashmap
  bsize+= (unsigned long int)bins_.size() * ((unsigned long int)sizeof(Position)+sizeof(uint));
  return bsize;
}

const uint SparseHistogramFeature::counter() const {
  return counter_;
}

uint SparseHistogramFeature::binCount(const Position& pos) {
  if (bins_.count(pos) == 0) {
    return 0;
  } else {
    MapTypeInt::iterator i = bins_.find(pos);
    return i->second;
  }
}

double SparseHistogramFeature::operator()(const Position& pos) {
  if (data_.count(pos) == 0) {
    return 0.0;
  } else {
    MapTypeDouble::iterator i = data_.find(pos);
    return i->second;
  }
}

MapTypeDouble SparseHistogramFeature::getMap() const {
  return data_;
}

MapTypeInt SparseHistogramFeature::getBinMap() const{
  return bins_;
}

void SparseHistogramFeature::calcRelativeFrequency() {
  for (MapTypeInt::iterator it = bins_.begin(); it != bins_.end(); it++) {
    data_[it->first] = double(bins_[it->first]) / double(counter_);
  }
}

void SparseHistogramFeature::feedbin(const Position& pos) {
  bins_[pos]++;
  counter_++;
  data_[pos] = double(bins_[pos]) / double(counter_);
}

Position SparseHistogramFeature::feed(const SHPoint& point) {
  Position pos=pointToPos(point);
  feedbin(pos);
  return pos;
}

vector<double>& SparseHistogramFeature::max() {
  return max_;
}
const vector<double>& SparseHistogramFeature::max() const {
  return max_;
}

vector<double>& SparseHistogramFeature::min() {
  return min_;
}

const vector<double>& SparseHistogramFeature::min() const {
  return min_;
}

void SparseHistogramFeature::initStepsize() {
  stepSize_ = vector<double>(dimensions_);
  for (uint i = 0; i < dimensions_; i++) {
    stepSize_[i] = (max_[i] - min_[i]) / double(steps_[i]);
  }
}

Position SparseHistogramFeature::pointToPos(const SHPoint& point) const {
  Position posString;
  for(uint i = 0; i < dimensions_; i++) {
    if (((point[i] - min_[i]) / stepSize_[i]) >= steps_[i]) {
      posString += char(steps_[i] - 1 + 48);
    } else if (((point[i] - min_[i]) / stepSize_[i]) < 0) {
      posString += char(48);
    } else if (((point[i] - min_[i]) / stepSize_[i]) < steps_[i]) {
      posString += char((point[i] - min_[i]) / stepSize_[i] + 48);
      //  } else {
      //    ERR << "Invalid element to feed into the histogram" << endl;
      //    ERR << "For dimension " << i << " with " << steps_[i] << " steps: " << point[i] << endl;
      //    exit(1);
    }
  }
  return posString;
}

void SparseHistogramFeature::expand(MapTypeDouble& expandedMap, const double smoothFactor) const {
  
  DBG(30) << "expanding with smoothFactor=" << smoothFactor << endl;

  // resize expanded map appropriately
  expandedMap.resize(data_.bucket_count() * 2 * dimensions_);
  
  Position adjustPosition;
  for (MapTypeDouble::const_iterator iBegin = data_.begin(); iBegin != data_.end(); iBegin++) {
    double value = iBegin->second;
    Position pos = iBegin->first;
    
    // first copy the value at the current position into the expanded map
    // but only (1.0 - smootheFactor) of it
    MapTypeDouble::iterator found = expandedMap.find(pos);
    if (found != expandedMap.end()) {
      found->second += value * (1.0 - smoothFactor);
    } else {
      Position newPosition = pos;
      expandedMap[newPosition] = value * (1.0 - smoothFactor);
    }
    
    // then count neighbors over all dimensions
    int numNeighbors = 0;
    for (uint adjust = 0; adjust < dimensions_; adjust++) {
      if (pos[adjust] > 1) {
	numNeighbors++;
      }
      if ((uint) pos[adjust] < (uint) steps_[adjust]) {
	numNeighbors++;
      }
    }
    
    // finally adjust neighbor positions in the expanded map
    // iterate over all dimensions
    double smoothPortion = (value * smoothFactor) / numNeighbors;
    
    for (uint adjust = 0; adjust < dimensions_; adjust++) {
      
      adjustPosition = pos;
      uint stepValue = adjustPosition.c_str()[adjust];
      // decrease position in dimension adjust by 1
      // first neighbor
      if (stepValue > 48) {
	string lessStr;
	lessStr += char(stepValue - 1);
	adjustPosition.replace(adjust, 1, lessStr);
	MapTypeDouble::iterator found = expandedMap.find(adjustPosition);
	if (found != expandedMap.end()) {
	  found->second+=smoothPortion;
	} else {
	  Position newPosition = adjustPosition;;
	  expandedMap[newPosition] = smoothPortion;
	}
      }
      
      // decrease position in dimension adjust by 1
      // other neighbor
      if (stepValue < 48 + (uint) steps_[adjust] - 1) {
	string moreStr;
	moreStr += char(stepValue + 1);
	adjustPosition.replace(adjust, 1, moreStr);
	MapTypeDouble::iterator found = expandedMap.find(adjustPosition);
	if (found != expandedMap.end()) {
	  found->second += smoothPortion;
	} else {
	  Position newPosition = adjustPosition;
	  expandedMap[newPosition] = smoothPortion;
	}
      }
    }
  }
  DBG(50) << "Expansion done, expanded map contains " << expandedMap.size() << " bins." << endl;
}



uint SparseHistogramFeature::prunePoorBins(uint threshold) {
  int pruned = 0;
  MapTypeInt::iterator it = bins_.begin();
  while (it != bins_.end()) {
    if (it->second < threshold) {
      pruned++;
      counter_ -= it->second;
      it++;
      bins_.erase(it->first);
      data_.erase(it->first);
    } else {
      it++;
    }
  }
  if (threshold > 0) {
    for (MapTypeInt::iterator it = bins_.begin(); it != bins_.end(); it++) {
      data_[it->first] = double(bins_[it->first]) / double(counter_);
    }
    DBG(10) << "pruned " << pruned << " bins below threshold " << threshold << endl;
  }
  return pruned;
}

uint SparseHistogramFeature::prunePoorBins(double threshold) {
  int pruned = 0;
  MapTypeDouble::iterator it = data_.begin();
  while (it != data_.end()) {
    if (it->second < threshold) {
      pruned++;
      counter_ -= bins_[it->first];
      it++;
      bins_.erase(it->first);
      data_.erase(it->first);
    } else {
      it++;
    }
  }
  for (MapTypeInt::iterator it = bins_.begin(); it != bins_.end(); it++) {
    data_[it->first] = double(bins_[it->first]) / double(counter_);
  }
  if (threshold > 0) {
    DBG(10) << "pruned " << pruned << " bins below threshold " << threshold << endl;
  }
  return pruned;
}


MapTypeDouble SparseHistogramFeature::getDifferenceHistogram(const MapTypeDouble& otherMap) {
  MapTypeDouble diffHisto;
  for (MapTypeDouble::iterator it = data_.begin(); it != data_.end(); it++) {
    diffHisto[it->first] = it->second;
  }
  for (MapTypeDouble::const_iterator it = otherMap.begin(); it != otherMap.end(); it++) {
    if (diffHisto.find(it->first) != diffHisto.end()) {
      diffHisto[it->first] = fabs(diffHisto[it->first] - it->second);
    } else {
      diffHisto[it->first] = it->second;
    }
  }
  return diffHisto;
}

void SparseHistogramFeature::getFilledNeighborBins(const Position& , vector<Position>& ) {
  /*
  for (uint i = 0; i < dimensions_; i++) {
    Position neighborPos = pos;
    uint stepValue = pos.c_str()[i];
    // decrease position in dimension i by 1
    // first neighbor
    if (stepValue > 48) {
      string lessStr;
      lessStr += char(stepValue - 1);
      neighborPosition.replace(i, 1, lessStr);
      MapTypeInt::iterator found = data_.find(neighborPosition);
      if (found != data_.end()) {
	if (found->second > 0) {
	  neighbors.push_back(neighborPosition);
	}
      }
    }
    // decrease position in dimension adjust by 1
    // other neighbor
    if (stepValue < 48 + (uint) steps_[i] - 1) {
      string moreStr;
      moreStr += char(stepValue + 1);
      neighborPosition.replace(i, 1, moreStr);
      MapTypeInt::iterator found = data_.find(neighborPosition);
      if (found != data_.end()) {
	if (found->second > 0) {
	  neighbors.push_back(neighborPosition);
	}
      }
    }
  }
  */
}



bool SparseHistogramFeature::read(istream &is) {
  istringstream iss;
  string line,token;

  // MagickNumber
  getline(is,line);
  if(line!="FIRE_sparse_histogram") {
    ERR << "Not reading a valid histogram, expected 'FIRE_sparse_histogram', got '"<< line <<"'." << endl;
    return false;
  }

  // Comments
  while('#'==is.peek()) { // comment lines
    getline(is,line); 
  }
  
  // dim
  token=""; getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="dim") {
    iss >> dimensions_;
  } else {
    ERR << "Expected 'dim', got '" << line <<"'." << endl;
    return false;
  }

  // counter
  token=""; getline(is,line); iss.clear(); iss.str(line); iss >> token;
  uint tmpCounter;
  if(token=="counter"){
    iss >> tmpCounter;
  } else {
    ERR << "Expected 'counter', got '" << token <<"'." << endl;
    return false;
  }

  // steps
  getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="steps") {
    steps_=vector<uint>(0);
    while(!iss.eof()) {
      uint tmp;
      iss >> tmp;
      steps_.push_back(tmp);
    }
    if (dimensions_ != steps_.size()) {
      ERR << "steps " << steps_.size() << " and dim " << dimensions_ << " inconsistent." << endl;
    }
  } else {
    ERR << "Expected 'steps', got '" << line <<"'." << endl;
    return false;
  }

  // min
  getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if (token=="min") {
    min_ = vector<double>(0);
    while (!iss.eof()) {
      double tmp;
      iss >> tmp;
      min_.push_back(tmp);
    }
    if (dimensions_ != min_.size()) {
      ERR << "min "<<min_.size() << "  and dim " << dimensions_ << " inconsistent." << endl;
      return false;
    }
  } else {
    ERR << "Expected 'min', got '" << line <<"'." << endl;
    return false;
  }

  //max
  getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="max") {
    max_=vector<double>(0);
    while(!iss.eof()) {
      double tmp;
      iss >> tmp;
      max_.push_back(tmp);
    }
    if (dimensions_ != max_.size()) {
      ERR << "max " << max_.size() << " and dim " << dimensions_ << " inconsistent." << endl;
      return false;
    }
  } else {
    ERR << "Expected 'max', got '" << line <<"'." << endl;
    return false;
  }

  //bins
  uint bins;
  getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="bins") {
    iss >> bins;
  } else {
    ERR << "Expected 'bins', got '" << line <<"'." << endl;
    return false;
  }

  // data
  bins_ = MapTypeInt((int) (bins * 1.5));
  data_ = MapTypeDouble((int) (bins * 1.5));
  for (uint i = 0; i < bins; i++) {
    getline(is,line); iss.clear(); iss.str(line); iss >> token;
    if (strcmp(token.c_str(), "data") == 0) {
      Position pos;
      iss >> pos;
      
      uint value;
      iss >> value;

      DBG(100) << VAR(pos) << " " << VAR(value) << endl;
      bins_[pos] = value;
      data_[pos] = double(value) / double(tmpCounter);
      counter_ += value;
      
    } else {
      ERR << "Expected 'data', got '" << line <<"'." << endl;
      return false;
    }
  }
  if (tmpCounter != counter_) {
    ERR << "Read wrong value for counter: " << tmpCounter << "(actual value is " << counter_ << ", " VAR(bins) <<")" << endl;
  }
  DBG(50) << "Sparse histogram contains " << bins << " bins." << endl;

  this->initStepsize();
  iss.clear();

  //calculate total document length
  length_=0;
  for (MapTypeInt::iterator i=bins_.begin();i!=bins_.end();i++){
    length_+=i->second;
  }
  
  return true;
}

bool SparseHistogramFeature::readBinary(istream &is){
  if(!is.good() || is.eof()){
  	return false;
  }
  // read dimension information
  is.read((char*)&(dimensions_),sizeof(uint));
  if(is.fail() || is.eof()){
  	return false;
  }

  // counter
  is.read((char*)&(counter_),sizeof(uint));
  if(is.fail() || is.eof()){
  	return false;
  }

  // steps
  // due to binary format steps must be of dimension dimensions_
  steps_.resize(dimensions_);
  for(uint i=0;i<dimensions_;++i){
    is.read((char*)&(steps_[i]),sizeof(uint));
    if(is.fail() || (i<dimensions_-1 && is.eof())){
    	return false;
    }
  }
   
  // min
  min_.resize(dimensions_);
  for(uint i=0;i<dimensions_;++i){
    is.read((char*)&(min_[i]),sizeof(double));
    if(is.fail() || (i<dimensions_-1 && is.eof())){
    	return false;
    }
  }
  
  // max
  max_.resize(dimensions_);
  for(uint i=0;i<dimensions_;++i){
    is.read((char*)&(max_[i]),sizeof(double));
    if(is.fail() || (i<dimensions_-1 && is.eof())){
    	return false;
    }
  }
 
  // bins
  uint binscount;
  is.read((char*)&(binscount),sizeof(uint));
  if(is.fail() || is.eof()){
  	return false;
  }
  
  // data
  bins_= MapTypeInt((int) (binscount * 1.5));
  data_= MapTypeDouble((int) (binscount * 1.5));
  for(uint i=0;i<binscount;++i){
    Position pos;
    is.read((char*)&pos,sizeof(Position));
    if(is.fail() || is.eof()){
    	return false;
    }
    uint value;
    is.read((char*)&value,sizeof(uint));
    if(is.fail() || (i<binscount-1 && is.eof())){
    	return false;
    }
    bins_[pos] = value;
    data_[pos] = (double)value/(double)counter_;
  }

  DBG(50) << "Sparse histogram contains " << binscount << " bins." << endl;
  this->initStepsize();

  length_=0;
  for (MapTypeInt::iterator i=bins_.begin();i!=bins_.end();i++){
    length_+=i->second;
  }
  
  return true;
}

void SparseHistogramFeature::write(ostream &os) {
  os << "FIRE_sparse_histogram" << endl
     << "# Sparse histogram file saved for Fire V2" << endl
     << "dim " << dimensions_ << endl
     << "counter " << counter_ << endl
     << "steps";
  for (uint i = 0; i < steps_.size(); ++i) {
    os << " "<< steps_[i];
  } 
  os << endl;
  os << "min";
  for (uint i = 0; i < min_.size(); ++i) {
    os << " "<< min_[i];
  } 
  os << endl;
  os << "max";
  for (uint i = 0;i < max_.size();++i) {
    os << " "<< max_[i]; 
  } 
  os << endl;
  os << "bins " << data_.size() << endl;

  MapTypeInt::const_iterator iBegin = bins_.begin();
  MapTypeInt::const_iterator iEnd = bins_.end();

  while (iBegin != iEnd) {
    uint value = iBegin->second;
    const Position pos = iBegin->first;
    os << "data ";
    os << pos << " "; 
    os << value << endl;
    iBegin++;
  }
}

void SparseHistogramFeature::writeBinary(ostream &os){
  os.write((char*)&(dimensions_),sizeof(uint));
  os.write((char*)&(counter_),sizeof(uint));
  // now write the vectors
  // first the steps
  for(uint i=0;i<steps_.size();++i){
    os.write((char*)&(steps_[i]),sizeof(uint));
  }
  // now the vectors for min and max
  for(uint i=0;i<min_.size();++i){
    os.write((char*)&(min_[i]),sizeof(double));
  }
  for(uint i=0;i<max_.size();++i){
    os.write((char*)&(max_[i]),sizeof(double));
  }
  // write bins
  os.write("bins",sizeof(char[5]));
  uint dsize = data_.size();
  os.write((char*)&dsize,sizeof(uint));
  // write the data
  MapTypeInt::const_iterator it=bins_.begin();
  MapTypeInt::const_iterator end = bins_.end();
  while(it != end){
    uint value = it->second;
    const Position pos = it->first;
    // beware! we skip here the magic 'data' word !
    os.write((char*)&pos,sizeof(Position));
    os.write((char*)&value,sizeof(uint));
    it++;
  }
}

void SparseHistogramFeature::collectBinValues(const MapTypeDouble& theMap, vector<double>& binValues) const {
  int numBins = 1;
  Position pos = "";
  for (int i = 0; i < (int) steps_.size(); i++) {
    numBins *= steps_[i];
    pos.append("0");
  }
  binValues.resize(numBins);
  if (theMap.find(pos) != theMap.end()) {
    binValues[0] = theMap.find(pos)->second;
  } else {
    binValues[0] = 0.0;
  }
  //DBG(10) << VAR(pos) << endl;

  //cout << "adding bin value for position " << pos << endl;
  for (int i = 1; i < numBins; i++) {
    int index = 0;
    bool done = false;
    while (!done) {
      char character = pos.c_str()[index] - 48;
      character++;
      if (character >= (int) steps_[index]) {
	string charStr;
	charStr += char(48);
	pos.replace(index, 1, charStr);
	index++;
      } else {
	string charStr;
	charStr += char(character + 48);
	pos.replace(index, 1, charStr);
	done = true;
      }
    }
    
//    DBG(10) << VAR(pos) << endl;

    if (theMap.find(pos) != theMap.end()) {
      binValues[i] = theMap.find(pos)->second;
    } else {
      binValues[i] = 0.0;
    }
    //cout << "adding bin value for position " << pos << endl;
  }
}


void SparseHistogramFeature::getAllBinValues(vector<double>& allBinValues) const {
  collectBinValues(data_, allBinValues);
}

void SparseHistogramFeature::getAllSmoothedBinValues(vector<double>& allSmoothedBinValues, double smootheFactor) const {
  MapTypeDouble smoothedMap;
  expand(smoothedMap, smootheFactor);
  collectBinValues(smoothedMap, allSmoothedBinValues);
}


