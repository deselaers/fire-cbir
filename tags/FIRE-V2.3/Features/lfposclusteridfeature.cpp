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

#include <iostream>
#include <sstream>
#include "diag.hpp"
#include "lfposclusteridfeature.hpp"

using namespace std;

bool LFPositionClusterIdFeature::read(istream &is) {
  uint x, y;
  string line;
  int clsid;
  while(not is.eof()) {
    vector<int> tmp;
  
    getline(is,line);
    if(!is.eof()) {
      istringstream iss(line);
      iss >> x >> y;
      while( not iss.eof()) {
        iss >> clsid;
        tmp.push_back(clsid);
      }
      positions_.push_back(PositionPair(x,y));
      data_.push_back(tmp);
    }
  }
  dim_=data_[0].size();
  DBG(15) << "Read " << data_.size() << " feature vectors of dim " << dim_ << endl;
  return true;
}

void LFPositionClusterIdFeature::write(ostream &os) {
  for(uint i=0;i<data_.size();++i) {
    os << position(i).first << " " << position(i).second;
    for(uint d=0;d<dim_;++d) {
      os << " " << data_[i][d];
    }
    os << endl;
  }
}

void LFPositionClusterIdFeature::push_back(const PositionPair& pos, const ::std::vector<int>& feat) {
  data_.push_back(feat);
  positions_.push_back(pos);
}
