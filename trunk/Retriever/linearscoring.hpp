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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef __linearscoring_hpp__
#define __linearscoring_hpp__
#include "diag.hpp"
#include "imagecomparator.hpp"
#include "basescoring.hpp"
#include <sstream>
#include "gzstream.hpp"

class LinearScoring : public BaseScoring {
private:
  ::std::vector<double> weights_;

public:

  LinearScoring(uint no=0) :
    weights_(no,1.0) {
    type_="linear";
  }
  LinearScoring(const ::std::string& filename) {
    type_="linear";
    load(filename);
  }
  LinearScoring(const ::std::string& filename, const uint& numberOfDistances) {
    type_="linear";
    if (filename=="") {
      //if theres is no config file present, assume there is one distance and set that weight to 1
      if (numberOfDistances==0) {
        weights_.push_back(1.0);
      } else {
        weights_.resize(numberOfDistances,1.0);
      }
    } else {
      load(filename);
    }
  }

  virtual void load(const ::std::string& filename) {
    igzstream is;
    is.open(filename.c_str());
    if (!is.good()) {
      ERR << "Cannot open '" << filename << "' for reading." << ::std::endl;
    } else {
      uint N;
      ::std::istringstream iss;
      ::std::string line, token;

      while ('#'==is.peek()) {
        getline(is, line);
      }

      getline(is, line);
      iss.clear();
      iss.str(line);
      iss >> token;
      if (token!="nofweights") {
        ERR << "not a valid weights file" << ::std::endl << " got '" << token << "' instead of 'nofweights'." << ::std::endl;
        exit(10);
      } else {
        iss >> N;
      }

      weights_=::std::vector<double>(N);

      getline(is,line);iss.clear(); iss.str(line); iss >> token;
      if(token!="weights") {
        ERR << "not a valid weights file" << ::std::endl
        << " got '" << token << "' instead of 'weights'." << ::std::endl;
        exit(10);
      } else {
        for(uint i=0;i<N;++i) {
          iss >> weights_[i];
        }
      }
      is.close();
    }

    DBG(10) << "Read " << weights_.size() << " weights from '" << filename <<"': ";
    for(uint i=0;i<weights_.size();++i) BLINK(10) << " " <<weights_[i]; BLINK(10) << ::std::endl;
  }

  virtual double getScore(const ::std::vector<double>& dists) {
    double result=0.0;
    uint M=dists.size();
    if (M>weights_.size()) {weights_.resize(M,0.0);}
    //DBG(10)<< "SCORE :";
    for(uint i=0;i<M;++i) {
      result+=weights_[i]*dists[i];
      //BLINK(10)<<" "<<dists[i];
    }
    //BLINK(10)<<std::endl;
    return exp(-result);
  }

  virtual double& weight(const uint idx) {
    return weights_[idx];
  }

  virtual const double& weight(const uint idx) const {
    return weights_[idx];
  }

  virtual const ::std::string settings() {
    ::std::ostringstream oss;
    oss << "size " << weights_.size();
    for(uint i=0;i<weights_.size();++i) {
      oss << " weight " << i << " " << weights_[i];
    }
    return oss.str();
  }

};

#endif
