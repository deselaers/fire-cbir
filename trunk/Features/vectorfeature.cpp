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
#include <sstream>
#include "vectorfeature.hpp"
#include "diag.hpp"

using namespace std;

bool VectorFeature::read(istream & is) {
  string line;
  getline(is,line);
  istringstream iss;
  if(!is.good()) {
    ERR << "Error reading" << endl;
    return false;
  }
  
  if(line!="FIRE_vectorfeature") {
    DBG(30) << "Magic number not found, ignoring, because this could be an old file" << endl;
  }
  
  getline(is,line); iss.str(line); string keyword; iss >> keyword;
  if("dim"==keyword) {
    uint size;
    iss >> size;
    data_.resize(size);
  } else {
    ERR << "Expected 'dim', got " << line << endl;
    return false;
  }

  getline(is,line); iss.clear(); iss.str(line); iss >> keyword;
  if ("data"==keyword) {
    for(uint i=0;i<data_.size();++i) {
      iss >> data_[i];
    }
  } else {
    ERR << "Expected 'data', got " << line << endl;
    return false;
  }
  return true;
}

bool VectorFeature::readBinary(istream &is){
	if(!is.good()){
      ERR << "Error reading" << endl;
      return false;      	
	}

	uint size = 0;
    is.read((char*)&size,sizeof(uint));
	if(is.fail() || is.eof()){
		return false;
	}
	data_.resize(size);
	
    // read the data
    for(uint i = 0; i < data_.size();++i){
      is.read((char*)&(data_[i]),sizeof(double));
      if(is.fail() || (i<size-1 && is.eof())){
      	return false;
      }
    }
    return true;
}

void VectorFeature::write(ostream &os) {
  os <<"FIRE_vectorfeature" << endl;
  os << "dim " << data_.size() << endl
     << "data ";
  for(uint i=0;i<data_.size();++i) {
    os << data_[i] << " ";
  }
  os << endl;
}

void VectorFeature::writeBinary(ostream &os){
	uint size = data_.size();
	os.write((char*)&size,sizeof(uint));
	for(uint i = 0; i< data_.size();++i) { 
      os.write((char*)&(data_[i]),sizeof(double)); } 
}

const uint VectorFeature::size() const {
  return data_.size();
}

const unsigned long int VectorFeature::calcBinarySize() const{
  unsigned long int bsize = 0;
  // add the size of the dimension number and the data
  bsize += sizeof(uint);
  bsize += (unsigned long int)data_.size() * (unsigned long int)sizeof(double);
  return bsize;
}

VectorFeature VectorFeature::operator-(const VectorFeature & v ) const
  
  {
    if (v.size() != size()){
      DBG(10)<<"THIS SHOULD NEVER HAPPEN!"<<std::endl;
    }
    VectorFeature result = VectorFeature((*this));
    for (size_t i=0; i < size(); ++i)
      result[i] = (*this)[i] - v[i];
    return result;
  }

VectorFeature & VectorFeature::operator-=(const VectorFeature & v ) 
  {
    if (v.size() != size()){
      DBG(10)<<"THIS SHOULD NEVER HAPPEN!"<<std::endl;
    }
    for (size_t i=0; i < size(); ++i)
     (*this)[i] -= v[i];
    return (*this);
  }
