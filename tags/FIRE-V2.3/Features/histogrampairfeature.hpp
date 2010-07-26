/*
  This file is part of the FIRE -- Flexible Image Retrieval System

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
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#ifndef __histogrampairfeature_hpp__
#define __histogrampairfeature_hpp__

#include <string>
#include <iostream>
#include <vector>
#include "histogramfeature.hpp"
#include <sstream>

class HistogramPairFeature : public BaseFeature {
private:
  HistogramFeature completeHisto_;
  ::std::vector<uint> counters_;
  ::std::vector< ::std::vector<uint> > data_;

public:

  virtual ~HistogramPairFeature() {}
  HistogramPairFeature(){ type_=FT_HISTOPAIR;}
  HistogramPairFeature(const HistogramFeature & h) :completeHisto_(h) {type_=FT_HISTOPAIR;}

  uint nOfHistograms() const {return data_.size();}

  virtual HistogramPairFeature* clone() const {return new HistogramPairFeature(*this);}
  
  const HistogramFeature& completeHisto() const {return completeHisto_;}
  HistogramFeature& completeHisto() {return completeHisto_;}
  
  const HistogramFeature bghisto(uint j) const {
    uint counter=completeHisto_.counter()-counters_[j];
    ::std::vector<uint> bins=completeHisto_.bins();
    for(uint i=0;i<completeHisto_.size();++i) {
      bins[i]-=data_[j][i];
    }
    return HistogramFeature(counter,bins);
  }
  
  const HistogramFeature histo(uint j) const {
    return HistogramFeature(counters_[j], data_[j]);
  }
  

  virtual void push_back(const HistogramFeature& h) {
    data_.push_back(h.bins());
    counters_.push_back(h.counter());
  }

  /** read a feature from a given istream. This should be well
      defined, as it is used in LargeFeatureFile as well, as in load,
      and possibly in other places, too. */ 
  virtual bool read(::std:: istream & is) {
    ::std::string tmp, line;
  
  
    getline(is,tmp);
    if(tmp!="FIRE_HistogramFeaturePairs") {
      ERR << "expecting FIRE_HistogramFeaturePairs, got: " << tmp << ::std::endl;
      return false;
    }
    
    getline(is,tmp);
    if(tmp!="completehistogram") {
      ERR << "expecting completehistogram, got: " << tmp << ::std::endl;
      return false;
    }
    
    completeHisto_=HistogramFeature();
    DBG(25) << "Reading the complete histo ...";
    if(not completeHisto_.read(is)) return false;
    BLINK(25) << "finished" << ::std::endl;
    
    uint nofhistos;
    getline(is,tmp);
    ::std::istringstream iss(tmp);
    iss >> tmp >> nofhistos;
    if(tmp!="numberOfVectors") {
      ERR << "expecting numberOfVectors, got: " << tmp << ::std::endl;
      return false;
    }
    
    data_=::std::vector< ::std::vector<uint> >(nofhistos,::std::vector<uint>(completeHisto_.size()));
    counters_=::std::vector<uint>(nofhistos);
    uint aktHisto=0;
    for(uint i=0;i<nofhistos;++i) {
      is >> tmp >> aktHisto >> counters_[i];
      if (aktHisto != i) {ERR << "aktHisto(" << aktHisto<< ") != i(" << i << ")" << ::std::endl;}
      
      for(uint j=0;j<completeHisto_.size();++j) {
        is >> data_[i][j];
      }
    }
    return true;
  }
  
  /** save a feature to a given ostream. This should be well defined,
      as it is used in LargeFeatureFile as well, as in save and
      possibly in other places, too.*/ 
  virtual void write(::std::ostream & os) {
    os << "FIRE_HistogramFeaturePairs" << ::std::endl
       << "completehistogram" << ::std::endl;
    completeHisto_.write(os);
    os << "numberOfVectors " << data_.size() << ::std::endl;
    for(uint i=0;i<data_.size();++i) {
      os << "histo " << i << " " << counters_[i];
      for(uint j=0;j<data_[i].size();++j) {
        os << " " << data_[i][j];
      }
      os << ::std::endl;
    }
  }
  
};

#endif
