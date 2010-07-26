#include "mpeg7feature.hpp"
#include <map>
#include <string>

using namespace std;

/// this is necessary to initialize the static member
map<string, MPEG7Type> MPEG7Feature::subtypeMap_;

/// constructor, here the subtype can be given, if none is given: set
/// to zero, which is not a valid type
MPEG7Feature::MPEG7Feature(MPEG7Type t) : mpegtype_(t) {
  type_=FT_MPEG7;
  if(subtypeMap_.size()==0) {
    subtypeMap_["ScalableColor"]=ScalableColor;
    subtypeMap_["ColorStructure"]=ColorStructure;
    subtypeMap_["DominantColor"]=DominantColor;
    subtypeMap_["ColorLayout"]=ColorLayout;
    subtypeMap_["TextureBrowsing"]=TextureBrowsing;
    subtypeMap_["HomogeneousTexture"]=HomogeneousTexture;
    subtypeMap_["EdgeHistogram"]=LocalEdgeHistogram;
    subtypeMap_["3DShapeSpectrum"]=ShapeHistogram;
    subtypeMap_["RegionShape"]=RegionBasedShape;
    subtypeMap_["ContourShape"]=ContourBasedShape;
  }
}


string MPEG7Feature::MPEG7Type2Name(const MPEG7Type type) const {
  string result="UnknownType";
  for(map< string, MPEG7Type >::const_iterator i=subtypeMap_.begin();i!=subtypeMap_.end();++i) {
    if(i->second==type) {
      result=i->first;
      break;
    }
  }
  return result;
}

/// given a suffix, determine the mpeg7 type
MPEG7Type MPEG7Feature::suffix2MPEG7Type(const string &suffix) {
  string lastSuffix=suffix.substr(suffix.rfind(".")+1,suffix.size());
  MPEG7Type result;
  if(lastSuffix=="mp7") { 
    uint posPoint=suffix.rfind(".");
    uint posPoint2=suffix.rfind(".",posPoint-1);
    lastSuffix=suffix.substr(posPoint2+1,posPoint-posPoint2-1);
  }
  if(subtypeMap_.find(lastSuffix)==subtypeMap_.end()) {
    ERR << "unknown MPEG7 subtype: " << lastSuffix << endl;
    result=0;
  } else {
    result=subtypeMap_[lastSuffix];
  }
  return result;
}
