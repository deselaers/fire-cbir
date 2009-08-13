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
#ifndef __lfsignaturefeature_hpp__
#define __lfsignaturefeature_hpp__

#include <string>
#include <iostream>
#include <sstream>
#include "basefeature.hpp"
#include "gaussiandensity.hpp"
                   

class LFSignatureFeature : public BaseFeature {
private:
  uint nOfBins_;
  ::std::vector< GaussianDensity > model_;
  ::std::vector< uint > signature_;
  
public:
 
  /// constructor, if no parameter (or 0) given -> defaults to false,
  /// else to true
  LFSignatureFeature() : nOfBins_(0), model_(), signature_()  {
    type_=FT_LFSIGNATURE;
  }
  
  LFSignatureFeature * clone() const {return new LFSignatureFeature(*this);}

  /// derived from base feature, read from stream
  ::std::vector<GaussianDensity>& model() {return model_;}
  const ::std::vector<GaussianDensity>& model() const {return model_;}
  ::std::vector<uint>& signature() {return signature_;}
  const ::std::vector<uint>& signature() const {return signature_;}

  bool read(::std::istream &is);
  void write(::std::ostream &os);

};

#endif
