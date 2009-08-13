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
#include "featureloader.hpp"
#include "vectorfeature.hpp"
#include "imagefeature.hpp"
#include "histogramfeature.hpp"
#include "sparsehistogramfeature.hpp"
#include "binaryfeature.hpp"
#include "localfeatures.hpp"
#include "metafeature.hpp"
#include "textfeature.hpp"
#include "histogrampairfeature.hpp"
#include "mpeg7feature.hpp"
#include "distancefilefeature.hpp"
#include "facefeature.hpp"
#include "pasannotation.hpp"
#include "lfsignaturefeature.hpp"                   
#include "lfposclusteridfeature.hpp"

using namespace std;

FeatureLoader::FeatureLoader() {
  map_["png"]=FT_IMG;
  map_["pgm"]=FT_IMG;
  map_["jpg"]=FT_IMG;
  map_["jpg"]=FT_IMG;
  map_["vec"]=FT_VEC;
  map_["histo"]=FT_HISTO;
  map_["lf"]=FT_LF;
  map_["sift"]=FT_LF;
  map_["bin"]=FT_BINARY;
  map_["txt"]=FT_META;
  map_["mp7"]=FT_MPEG7;
  map_["histopair"]=FT_HISTOPAIR;
  map_["histopairs"]=FT_HISTOPAIR;
  map_["sparsehisto"]=FT_SPARSEHISTO;
  map_["sparsehistogram"]=FT_SPARSEHISTO;
  map_["textID"]=FT_TEXT;
  map_["dists"]=FT_DISTFILE;
  map_["facefeat"]=FT_FACEFEAT;
  map_["faces"]=FT_FACEFEAT;
  map_["annotation"]=FT_PASCALANNOTATION;
  map_["lfsig"]=FT_LFSIGNATURE;
  map_["featcls"]=FT_LFPOSCLSIDFEAT;

  // not yet in factory
  map_["gab"]=FT_GABOR;
  map_["blb"]=FT_BLOBS;
  map_["reg"]=FT_REGIONS;
  map_["oldhisto"]=FT_OLDHISTO;
  map_["oldvec"]=FT_VEC;

  featureFactory_.registerClass("png",BaseFeature::create<ImageFeature>);
  featureFactory_.registerClass("jpg",BaseFeature::create<ImageFeature>);
  featureFactory_.registerClass("pgm",BaseFeature::create<ImageFeature>);
  featureFactory_.registerClass("histo",BaseFeature::create<HistogramFeature>);
  featureFactory_.registerClass("oldhisto",BaseFeature::create<HistogramFeature>);
  featureFactory_.registerClass("vec",BaseFeature::create<VectorFeature>);
  featureFactory_.registerClass("oldvec",BaseFeature::create<VectorFeature>);
  featureFactory_.registerClass("lf",BaseFeature::create<LocalFeatures>);
  featureFactory_.registerClass("sift",BaseFeature::create<LocalFeatures>);
  featureFactory_.registerClass("bin",BaseFeature::create<BinaryFeature>);
  featureFactory_.registerClass("txt",BaseFeature::create<MetaFeature>);
  featureFactory_.registerClass("meta",BaseFeature::create<MetaFeature>);
  featureFactory_.registerClass("mp7",BaseFeature::create<MPEG7Feature>);
  featureFactory_.registerClass("histopair",BaseFeature::create<HistogramPairFeature>);
  featureFactory_.registerClass("histopairs",BaseFeature::create<HistogramPairFeature>);
  featureFactory_.registerClass("sparsehisto",BaseFeature::create<SparseHistogramFeature>);
  featureFactory_.registerClass("sparsehistogram",BaseFeature::create<SparseHistogramFeature>);
  featureFactory_.registerClass("textID",BaseFeature::create<TextFeature>);
  featureFactory_.registerClass("dists",BaseFeature::create<DistanceFileFeature>);
  featureFactory_.registerClass("facefeat",BaseFeature::create<FaceFeature>);
  featureFactory_.registerClass("faces",BaseFeature::create<FaceFeature>);
  featureFactory_.registerClass("annotation",BaseFeature::create<PascalAnnotationFeature>);
  featureFactory_.registerClass("lfsig",BaseFeature::create<LFSignatureFeature>);
  featureFactory_.registerClass("featcls",BaseFeature::create<LFPositionClusterIdFeature>);
}
  
FeatureType FeatureLoader::suffix2Type(const ::std::string& suffix) const {
  FeatureType result=0;
  if(map_.find(suffix) == map_.end()) {
    ERR << "No filetype assiciated with this suffix: " << suffix << endl;
  } else {
    result=map_.find(suffix)->second;    
  }
  DBG(50) << VAR(result) << endl;
  return result;
}


BaseFeature *FeatureLoader::makeNewFeature(const ::std::string& suffix) const {
  BaseFeature *result=NULL;
  result=featureFactory_.getObject(suffix);
  if(!result) {
    ERR << "No feature type associated with this suffix: " << suffix << endl;
  } 
  return result;
}

BaseFeature* FeatureLoader::load(const ::std::string& basename, const ::std::string& suffix, const ::std::string& lastSuffix, const ::std::string& path) {
    
  BaseFeature *result=NULL;
  string filename=path+"/"+basename+"."+suffix;
  
  FeatureType type=suffix2Type(lastSuffix);
  DBG(50) << VAR(type) << endl;
  
  result=makeNewFeature(lastSuffix);

  // now: special cases. Default is at the end...
  if(suffix=="oldhisto") {
    DBG(35) << "Loading OldHistogram from '" << filename << "'." << endl;
    dynamic_cast<HistogramFeature*>(result)->loadOld(filename);
    DBG(35) << "Loaded OldHistogram from '" << filename << "'." << endl;
  } else if (suffix=="mp7") {
    DBG(35) << "Set MPEG7 properties" << endl;
    dynamic_cast<MPEG7Feature*>(result)->setMPEG7Type(suffix);
    dynamic_cast<MPEG7Feature*>(result)->basename()=basename;
  } else { // default
    DBG(35) << "Loading from '" << filename << "'." << endl;
    result->load(filename);
    DBG(35) << "Loaded from '" << filename << "'." << endl;
  }
  return result;

}
