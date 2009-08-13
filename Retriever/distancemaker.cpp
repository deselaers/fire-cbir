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
#include <string>
#include <sstream>
#include "distancemaker.hpp"
#include "stringparser.hpp"

using namespace std;

DistanceMaker::DistanceMaker() {
  distanceNames_["base"]=DT_BASE;
  distanceNames_["basedist"]=DT_BASE;
  distanceNames_["euclidean"]=DT_EUCLIDEAN;
  distanceNames_["l1"]=DT_L1;
  distanceNames_["cityblock"]=DT_L1;
  distanceNames_["jsd"]=DT_JSD;
  distanceNames_["jd"]=DT_JSD;
  distanceNames_["kld"]=DT_KLD;
  distanceNames_["chisquare"]=DT_CHISQUARE;
  distanceNames_["histogramintersection"]=DT_HISTOGRAMINTERSECTION;
  distanceNames_["his"]=DT_HISTOGRAMINTERSECTION;
  distanceNames_["reldev"]=DT_RELDEV;
  distanceNames_["relativedeviation"]=DT_RELDEV;
  distanceNames_["relbindev"]=DT_RELBINDEV;
  distanceNames_["relativebindeviation"]=DT_RELBINDEV;
  distanceNames_["crosscorrelation"]=DT_CROSSCORRELATION;
  distanceNames_["fidelity"]=DT_ONEMINUSFIDELITY;
  distanceNames_["omf"]=DT_ONEMINUSFIDELITY;
  distanceNames_["oneminusfidelity"]=DT_ONEMINUSFIDELITY;
  distanceNames_["sqrtoneminusfidelity"]=DT_SQRTONEMINUSFIDELITY;
  distanceNames_["logtwominusfidelity"]=DT_LOGTWOMINUSFIDELITY;
  distanceNames_["arccosfidelity"]=DT_ARCCOSFIDELITY;
  distanceNames_["sinfidelity"]=DT_SINFIDELITY;
  distanceNames_["oneminusfidelitysquare"]=DT_SINFIDELITY;
  distanceNames_["idm"]=DT_IDM;
  distanceNames_["imagedistortionmodel"]=DT_IDM;
  distanceNames_["binary"]=DT_BINARYFEATURE;
  distanceNames_["binaryfeature"]=DT_BINARYFEATURE;
  distanceNames_["globallocalfeaturedistance"]=DT_GLFD;
  distanceNames_["metafeature"]=DT_METAFEATURE;
  distanceNames_["textfeature"]=DT_TEXTFEATURE;
  distanceNames_["glfd"]=DT_GLFD;
  distanceNames_["mpeg7"]=DT_MPEG7;
  distanceNames_["histpair"]=DT_HISTOPAIR;
  distanceNames_["lfhungarian"]=DT_LFHUNGARIAN;
  distanceNames_["distfile"]=DT_DISTFILE;
  distanceNames_["weightedeuclidean"]=DT_WEIGHTEDEUCLIDEAN;
  distanceNames_["facefeature"]=DT_FACEFEAT;
  distanceNames_["facefeat"]=DT_FACEFEAT;
  distanceNames_["lfsigemd"]=DT_LFSIGEMD;
  distanceNames_["rast"]=DT_RAST;
  distanceNames_["tfidf"]=DT_TFIDF;
  distanceNames_["bm25"]=DT_BM25;
  distanceNames_["smart2"]=DT_SMART2;
  distanceNames_["weightedl1"]=DT_WEIGHTEDL1;
  //tbc
  
  defaultDistances_[FT_BASE]=DT_BASE;
  defaultDistances_[FT_IMG]=DT_EUCLIDEAN;
  defaultDistances_[FT_VEC]=DT_EUCLIDEAN;
  defaultDistances_[FT_HISTO]=DT_JSD;
  defaultDistances_[FT_LF]=DT_GLFD;
  defaultDistances_[FT_BINARY]=DT_BINARYFEATURE;
  defaultDistances_[FT_META]=DT_METAFEATURE;
  defaultDistances_[FT_MPEG7]=DT_MPEG7;
  defaultDistances_[FT_HISTOPAIR]=DT_HISTOPAIR;
  defaultDistances_[FT_SPARSEHISTO]=DT_JSD;
  defaultDistances_[FT_OLDHISTO]=DT_JSD;
  defaultDistances_[FT_META]=DT_METAFEATURE;
  defaultDistances_[FT_TEXT]=DT_TEXTFEATURE;
  defaultDistances_[FT_DISTFILE]=DT_DISTFILE;
  defaultDistances_[FT_FACEFEAT]=DT_FACEFEAT;
  defaultDistances_[FT_LFSIGNATURE]=DT_LFSIGEMD;
  defaultDistances_[FT_LFPOSCLSIDFEAT]=DT_RAST;
  //tbc
}


BaseDistance* DistanceMaker::makeDistance(const ::std::string& distancename) {
  BaseDistance* result; 
  string dn, par;

  uint pos=distancename.find(":");
  if(pos<distancename.size()) {
    dn=distancename.substr(0,pos);
    par=distancename.substr(pos+1,distancename.size()-1);
  } else {
    dn=distancename;
    par="";
  }
  
  //  cout << pos << " " << dn <<" "<< par<<  endl;
  
  if(distanceNames_.find(dn)==distanceNames_.end()) {
    ERR << "No such distance defined: " << distancename << ::std::endl;
    result=new BaseDistance();
  } else {
    result=makeDistance(distanceNames_[dn],par);
  }
  return result;
}


BaseDistance* DistanceMaker::getDefaultDistance(const FeatureType featType) {
  BaseDistance* result; 
  if(defaultDistances_.find(featType)==defaultDistances_.end()) {
    ERR << "No default distance for FeatureType " << featType << " defined." << ::std::endl;
    result=new BaseDistance();
  } else {
    result=makeDistance(defaultDistances_[featType]);
  }
  return result;
}


BaseDistance* DistanceMaker::makeDistance(const DistanceType distType, const string &par) {
  BaseDistance* result; 
  switch(distType) {
  case DT_BASE:
    result=new BaseDistance();
    break;
  case DT_EUCLIDEAN:
    result=new EuclideanDistance();
    break;
  case DT_CROSSCORRELATION: {
    int d=getIntAfter(par,"D=",4);
    result=new CrosscorrelationDistance(d);
    break;
  }
  case DT_L1:
    result=new L1Distance();
    break;
  case DT_JSD:{
    double smoothFactor=getDoubleAfter(par,"SF=",0.5);
    result=new JSDDistance(smoothFactor);
    break;}
  case DT_KLD:
    result=new KLDDistance();
    break;
  case DT_CHISQUARE:
    result=new ChisquareDistance();
    break;
  case DT_HISTOGRAMINTERSECTION:
    result=new HistogramintersectionDistance();
    break;
  case DT_RELDEV:
    result=new RelativeDeviationDistance();
    break;
  case DT_RELBINDEV:
    result=new RelativeBinDeviationDistance();
    break;
  case DT_ONEMINUSFIDELITY:
    result=new OneMinusFidelityDistance();
    break;
  case DT_SQRTONEMINUSFIDELITY:{
    result=new SqrtOneMinusFidelityDistance();
    break;}
  case DT_LOGTWOMINUSFIDELITY:{
    result=new LogTwoMinusFidelityDistance();
    break;}
  case DT_ARCCOSFIDELITY:{
    result=new ArccosFidelityDistance();
    break;}
  case DT_SINFIDELITY:{
    result=new SqrtOneMinusFidelityDistance();
    break;}
  case DT_BINARYFEATURE:{
    result=new BinaryFeatureDistance();
    break;}
  case DT_METAFEATURE:{
    result=new MetaFeatureDistance();
    break;}
  case DT_TEXTFEATURE:{
    string serverstr=getStringAfter(par,"SERVER=","localhost");
    int portno=getIntAfter(par,"PORT=",4242);
    string language=getStringAfter(par,"LANG=","None");
    DBG(10) << "textfeature " << VAR(serverstr) << " "<< VAR(portno) << " " << VAR(language) << endl;
    result=new TextFeatureDistance(serverstr, portno, language);
    break;}
  case DT_IDM: {
    int wr1=getIntAfter(par,"WR1=",3);
    int wr2=getIntAfter(par,"WR2=",1);
    double threshold=getDoubleAfter(par,"TH=",0.05);
    bool sobel=not getBooleanString(par,"NOSOBEL");
    
    result=new ImageDistortionModelDistance(wr1,wr2,threshold,sobel);
    break;
  }

  case DT_GLFD: {
    string treename=getStringAfter(par,"TREE=","tree.kdt");
    uint k=getIntAfter(par,"K=",10);
    double epsilon=getDoubleAfter(par,"EPS=",0.1);
    
    result=new GlobalLocalFeatureDistance(treename,k,epsilon);
    break;}

  case DT_MPEG7:{
    string xmmainpath=getStringAfter(par,"XMMAIN=","./XMMain.exe");
    string mpegdatapath=getStringAfter(par,"MPD=","mpegdata");
    
    result=new MPEG7Distance(xmmainpath,mpegdatapath);
    break;}
  case DT_HISTOPAIR:{
    string grounddistname=getStringAfter(par,"GDIST=","jsd");
    double centerweight=getDoubleAfter(par,"CWEIGHT=",0.5);
    
    result=new HistogramPairDistance(makeDistance(grounddistname), centerweight);
    break;}
  case DT_DISTFILE:{
    string scoringname=getStringAfter(par,"SCORING=","linear",' ');
    bool clearing= not getBooleanString(par,"NOCLEARING");
    result=new DistanceFileDistance(scoringname,clearing);
    break;
  }
  case DT_WEIGHTEDEUCLIDEAN: {
    double alpha = getDoubleAfter(par, "weight=", 0.5);
    result = new WeightedEuclideanDistance(alpha);
    break;
  }
 case DT_LFHUNGARIAN:{
    result=new LFHungarianDistance();
    break;
  }
  case DT_FACEFEAT: {
    string allOrOne=getStringAfter(par,"USE=","one");
    result=new FaceFeatureDistance(allOrOne);
    break;
  }
  case DT_LFSIGEMD: {
    result=new LFSignatureEMDistance();
    break;
  }
  case DT_RAST: {
    int eps=getIntAfter(par,"EPS=",4);
    double tolerance=getDoubleAfter(par,"TOLERANCE=",0.05);
    int mindx=getIntAfter(par,"MINDX=",-200);
    int maxdx=getIntAfter(par,"MAXDX=",200);
    int mindy=getIntAfter(par,"MINDY=",-100);
    int maxdy=getIntAfter(par,"MAXDY=",100);
    int minq=getIntAfter(par,"MINQ=",10);
    double amin=getDoubleAfter(par,"AMIN=",-0.1);
    double amax=getDoubleAfter(par,"AMAX=",0.1);
    double minscale=getDoubleAfter(par,"MINSCALE=",0.8);
    double maxscale=getDoubleAfter(par,"MAXSCALE=",1.2);
    result=new RASTDistance(eps,tolerance,mindx,maxdx,mindy,maxdy,minq,amin,amax,minscale,maxscale);
    break;
  }
  case DT_TFIDF: {

    DBG(10)<< "Creating TFIDFDistance"<<endl;
    result=new TFIDFDistance();
    break;
  }
  case DT_BM25: {

    DBG(10)<< "Creating BM25 Distance"<<endl;
    double k1=1.5;
    double k3=8;
    double b=0.68;
    result=new BM25Distance(k1,k3,b);
    break;
  }

  case DT_SMART2:{
    DBG(10)<< "Creating SMART2Distance"<<endl;
    result=new SMART2Distance();
    break;
  }
  
  case DT_WEIGHTEDL1:{
    DBG(10) << "Creating Weighted L1 Distance" << endl;
    int maxiter=getIntAfter(par,"MAXITER=",100);
    double stepwidth=getDoubleAfter(par,"STEPWIDTH=",0.1);
    double regularisationWeight=getDoubleAfter(par,"REGWEIGHT=",0.0);
    result=new WeightedL1Distance(maxiter,stepwidth,regularisationWeight);
    break;
  }
    
  default:
    ERR << "Distancetype not known: " << distType << endl;
    result=new BaseDistance();
    break;
  }
  return result;
}

string DistanceMaker::availableDistances() const {
  ostringstream oss("");
    
  for(map<const string,DistanceType>::const_iterator i=distanceNames_.begin();i!=distanceNames_.end();++i) {
    oss << i->first << " ";
  }
  return oss.str();
}
