#include "lfsignaturefeature.hpp"

bool LFSignatureFeature::read(::std::istream &is) {
    ::std::string line;
    getline(is,line);
    if (line != "FIRE_LFSignature") {
      ERR << "Not a FIRE_LFSignature." << ::std::endl;
      return false;
    }
    while(!is.eof()) {
      getline(is,line);
      if(!is.eof()) {
        ::std::istringstream iss(line);
        ::std::string keyword;
        iss >> keyword;
        if(keyword == "gaussians") {
          iss >> nOfBins_;
          model_=::std::vector<GaussianDensity>(nOfBins_);
          signature_=::std::vector<uint>(nOfBins_);
        } else if(keyword=="signature") {
          for(uint i=0;i<nOfBins_;++i) {
            iss >> signature_[i];
          }
        } else if(keyword=="density") {
          uint no;
          iss >> no;
          iss >> keyword;
          if(keyword=="elements") {
            iss >> model_[no].elements;
          } else if(keyword=="dim") {
            iss >> model_[no].dim;
          } else if(keyword=="mean") {
            GaussianDensity& c=model_[no];
            c.mean=::std::vector<double>(c.dim);
            for(uint i=0;i<c.dim;++i) {
              iss >> c.mean[i];
            }
          } else if(keyword=="variance") {
            GaussianDensity& c=model_[no];
            c.sigma=::std::vector<double>(c.dim);
            for(uint i=0;i<c.dim;++i) {
              iss >> c.sigma[i];
            }
          } else {
            ERR << "Reading density received unknown keyword in position 3: '" << keyword << "'." << ::std::endl;
            return false;
          }
        } else {
          ERR << "Unknown keyword '" << keyword << "'." << ::std::endl;
          return false;
        }
      }
    }
    return true;
  }

  /// derived from base feature, write to stream
  void LFSignatureFeature::write(::std::ostream &os) {
    os << "FIRE_LFSignature" << ::std::endl
       << "gaussians " << model_.size() << ::std::endl
       << "signature" ;
    for(uint i=0;i<signature_.size();++i) os << " " << signature_[i];
    os << ::std::endl;
    for(uint i=0;i<model_.size();++i) {
      GaussianDensity& c=model_[i];
      os << "density " << i << " elements " << c.elements << ::std::endl;
      os << "density " << i << " dim " << c.dim << ::std::endl;
      os << "density " << i << " mean";
      for(uint k=0;k<c.dim;++k) { os << " " << c.mean[k]; } os << ::std::endl;
      
      os << "density " << i << " variance";
      for(uint k=0;k<c.dim;++k) { os << " " << c.sigma[k]; } os << ::std::endl;
    }
  }

