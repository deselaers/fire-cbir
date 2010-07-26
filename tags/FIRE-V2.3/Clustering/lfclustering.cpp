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
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include "diag.hpp"
#include <libgen.h>
#include "baseclusterer.hpp"
#include "em.hpp"
#include "getpot.hpp"
#include "localfeatures.hpp"
#include "clusterlocalfeatures.hpp"
#include "distancemaker.hpp"
#include "histogramfeature.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include "gzstream.hpp"
#include "database.hpp"
using namespace std;

typedef struct marker_ {
  marker_() : color(3,0) {}
  char type;
  vector<double> color;
} Marker;

map<uint,Marker> readHowToMarkFile(string filename) {
  //  file is structured as follows:
  //  no (+|x) r g b
  map<uint, Marker> result;
  igzstream ifs; ifs.open(filename.c_str());
  if(!ifs.good() || !ifs) {
    ERR << "Cannot open markfile " << filename << ". Aborting." << endl;
    exit(20);
  }
  
  string line;
  Marker tmp;
  while(!ifs.eof() && filename!="") {
    getline(ifs,line);      
    if(line!="") {
      istringstream iss(line);
      uint no;
      tmp=Marker();
      iss >> no >> tmp.type >> tmp.color[0] >> tmp.color[1] >> tmp.color[2];
      result[no]=tmp;
    }
  }
  ifs.close();
  return result;
}

void USAGE() {
  cout << "USAGE: lfclustering - cluster local features" << endl
       << "  lfclustering <options> --[fire]filelist <list of localfeature files>" << endl
       << "  Available options are: " << endl
       << "   CLUSTEROPTION" << endl
       << "    -position            to add position data when clustering" << endl
       << "    -splitMode <splitmodename> in {allSplit,largestSplit,varianceSplit}" <<endl
       << "    -poolMode <poolmodename> in {*noPooling*, clusterPooling, dimensionPooling, bothPooling}" << endl
       << "    -disturbMode <distmodename> in {constDisturb, meanDisturb, meanDisturb2, varianceDisturb}" << endl
       << "    -maxSplit <maximal number of splits> " << endl
       << "    -stopWithNClusters <maximal number of clusters>" << endl
       << "    -dontSplitBelow <minimal number of observations in a cluster to be split>" << endl
       << "    -iterationsBetweenSplits <number of reestimate iterations after each split>" << endl
       << "    -minObservationsPerCluster <minimal number of observation in a cluster not to be purged>" << endl
       << "    -epsilon <epsilon to disturb means>" << endl
       << "    -dist <distance name, if not specified mahalanobis>" << endl
       << "    --saveModel <filename> to save the model of the clustering to the specified file after clustering" << endl
       << "    --saveBeforeSplits to save the model before a split" << endl
       << "    --startWith <filename> to load a model which is used as start point" << endl
       << endl
       << "   HISTOGRAMIZATION (probably only suitable if --filelist was used" << endl
       << "    --noClustering <filename> don't cluster the local features, but load the model from the given file" << endl
       << "    --histogramization   to create local feature histograms for each of the clustered files" << endl
       << "    --histoOutputPath  where to save the histograms" << endl
       << "    --suffix             to change the suffix, usually 'histo.gz'" << endl
       << "    --fuzzy              to use fuzzy assignments for histogramization" << endl
       << "       fuzzy is only supported if used with --noClustering." << endl
       << "    --fuzzyalpha <0.5>   kernel size for fuzzy histogram" <<endl
       << "    --fuzzysuffix <fhisto.vec.gz> suffix for fuzzy histos (safed as vector features)" << endl
       << endl
       << "   SPECIALFUNCTIONS (DANGEROUS)" << endl
       << "     take care, these functions are HACKED! " << endl
       << "     these functions can only be applied together with --noclustering" << endl
       << "    --savePositionImage to save an image, where at position where a local feature was extracted"  << endl
       << "                     a mark denoting the cluster is created. Which mark to use for which cluster " << endl
       << "                     is specified with the --howtomark file. " << endl
       << "    --markSuffix <suffix> (default .positionmark.png) this suffix is appended to the name of the file" << endl
       << "                     containing the local features to be histogramized" << endl
       << "    --howtomark <filename>   in this file it is specified how to mark the positions where local features " << endl
       << "                       have been extracted" << endl
       << "    --cropSuffix <N>  how many characters have to be removed from the tail of the lf-file to obtain the image to be loaded and marked" << endl
       << "    --savePatches <clnr> to save all patches which are assigned to the given clusternumber" << endl
       << "               combination with marking is strange!!!" << endl
       << "    --discardDimension <d> to set dimension d to zero in each lf" << endl
       << "    --discardRange <start> <end> to delete all dimension in the range [start, end]" << endl
       << "    --writeFeatClsFiles  to write featcls files" << endl
       << "      option to this: --accIdPos accumulate identical position" << endl
       << "    --featclssuffix <suffix>  default: featcls" << endl
       << endl;
}

int main(int argc, char **argv) {
  GetPot cl(argc, argv);
  
  if(!(cl.search("--filelist") || cl.search("--firefilelist")) || cl.search("-h")) { USAGE(); exit(20);}
  
  
  bool savePatches=cl.search("--savePatches");

  bool marking=cl.search("--savePositionImage");
  uint cropCharacters=cl.follow(10,"--cropSuffix");

  uint discardDimension = cl.follow(0,"--discardDimension");
  bool discardDim = cl.search("--discardDimension");
  uint discardRangeStart = 0, discardRangeEnd = 0;
  bool discardRange = cl.search("--discardRange");
  if (discardRange) {
    discardRangeStart = cl.next(0);
    discardRangeEnd = cl.next(0);
  }

  map<uint, Marker> markers;
  if(marking|| savePatches) {
    markers = readHowToMarkFile(cl.follow("markerfile","--howtomark"));
  }
  
  bool fuzzy=cl.search("--fuzzy");
  bool withPosition = cl.search("-position");
  string fuzzysuffix=cl.follow("fhisto.vec.gz","--fuzzysuffix");
  double fuzzyAlpha=cl.follow(0.5,"--fuzzyalpha");
  
  bool saveBeforeSplits=cl.search("--saveBeforeSplits");
  string histopath=cl.follow("","--histoOutputPath");

  string splitMode = cl.follow("allSplit","-splitMode");
  string disturbMode = cl.follow("meanDisturb","-disturbMode");
  string poolMode = cl.follow("noPooling","-poolMode");
  uint maxSplit = cl.follow(10,"-maxSplit");
  uint stopWithNClusters = cl.follow(999999999,"-stopWithNClusters");
  uint dontSplitBelow = cl.follow(8,"-dontSplitBelow");
  uint iterationsBetweenSplits = cl.follow(10,"-iterationsBetweenSplits");
  uint minObservationsPerCluster = cl.follow(5,"-minObservationsPerCluster");
  double epsilon = cl.follow(0.1,"-epsilon");
  string distanceMaker = cl.follow("","-dist");
  string modelfile=cl.follow("","--saveModel");
  vector<string> filenames;
  bool writeFeatClsFiles=cl.search("--writeFeatClsFiles");
  bool accIdPos=cl.search("--accIdPos");
  string featclssuffix=cl.follow("featcls","--featclssuffix");
    
  ResultVector clusterinformation;
    
  //------------------------------------------------
  // read filelist  & local features

  if(cl.search("--filelist")) {
    cl.search("--filelist");
    string filename=cl.next("-");
    while(filename != "-") {
      igzstream ifs; ifs.open(filename.c_str());
      if(!ifs.good() || !ifs) {
        ERR << "Cannot open filelist " << filename << ". Aborting." << endl;
        exit(20);
      } else {
        DBG(10) << "Reading filelist: " << filename << endl;
        
        while(!ifs.eof() && filename!="") {
          getline(ifs,filename);      
          if(filename!="") {
            filenames.push_back(filename);
          }
        }
        ifs.close();
      }
      filename=cl.next("-");
    }
    DBG(10) << "Filelist reading done... going to process " << filenames.size() << " files in the following." << endl;
  } else if (cl.search("--firefilelist")) {
    string filename=cl.next("");
    Database db;
    DBG(10) << "Reading filelist from FIRE filelist" << endl;
    db.loadFileList(filename);
    

    for(uint n=0;n<db.size();++n) {
      filenames.push_back(db.path()+"/"+db.filename(n)+"."+db.suffix(0));
    }
    DBG(10) << "Filelist reading done... going to process " << filenames.size() << " files in the following." << endl;

    
  }
  // filelist reading done
  //------------------------------------------------
  

  if(!cl.search("--noClustering")) { // here clustering is applied
    //------------------------------------------------
    DBG(10) << "Loading local features from filelist" << endl;
    vector<LocalFeatures*> locfeat;
    LocalFeatures* lf;

    uint dim=0;
    uint dimchanges=0;
    uint totalNumber=0;

    
    uint modulo=filenames.size()/10;
    
    for(uint i=0;i<filenames.size();++i) {
      if (modulo != 0) {
        if (i%modulo==0) {DBG(10) << "Read " << i << "/" << filenames.size() << "." << endl;}
      }
      string filename=filenames[i];
      lf=new LocalFeatures();
      lf->load(filename);
      if(discardDim) {
        lf->discardDimension(discardDimension);
      }
      if (discardRange) {
        lf->discardRange(discardRangeStart, discardRangeEnd);
      }
      if(lf->dim()!=dim) {
        dim=lf->dim();
        DBG(10) << "Setting dim to " << dim << endl;
        ++dimchanges;
      }

      totalNumber+=lf->numberOfFeatures();
      locfeat.push_back(lf);
      totalNumber++;
    }
    DBG(10) << totalNumber << " local features loaded." << endl;
    // local features are loaded
    //------------------------------------------------
    // now clustering
    DBG(10) << "Starting to cluster" << endl;
    clusterinformation = cluster(locfeat, withPosition, splitMode, maxSplit, stopWithNClusters, disturbMode, poolMode, dontSplitBelow, 
                                 iterationsBetweenSplits, minObservationsPerCluster, epsilon, distanceMaker,modelfile,saveBeforeSplits,cl.follow("","--startWith"));
    DBG(10) << "Clustering done" << endl;
    //------------------------------------------------
    
    if(cl.search("--histogramization")) {


      if(fuzzy) {
        ERR << "Fuzzy not supported for this case. Only supported if run with --noClustering." << endl;
        exit(10);
      }
      
      DBG(10) << "Starting histogramization" << endl;
      string suffix=cl.follow("histo.gz","--suffix");
      
      uint currentPosition=0;
      uint nOfClusters=(*max_element(clusterinformation.begin(),clusterinformation.end()))+1;
      
      for(uint i=0;i<filenames.size();++i) {
        DBG(15) << "Histogramization of " << filenames[i] ;
        HistogramFeature histo(nOfClusters); histo.min()=vector<double>(1,0.0); histo.max()=vector<double>(1,nOfClusters);
        for(uint j=0;j<locfeat[i]->numberOfFeatures();++j) {
          histo.feedbin(clusterinformation[currentPosition]);
          ++currentPosition;
        }
        
        if(histopath!="") {
          string bname=basename(const_cast<char*>(filenames[i].c_str()));
          histo.save(histopath+"/"+bname+"."+suffix);
        } else {
          histo.save(filenames[i]+"."+suffix);
        }
        BLINK(15) << " done" << endl; 
      }
      
      DBG(15) << "currentPosition=" << currentPosition << ", clusterinformation.size()=" << clusterinformation.size() << endl;
      DBG(10) << "Histogramization finished" << endl;
    }
  } else { // option --noClustering was given, thus loading model
    //here I have to load the model
    string filename=cl.follow("model.clustermodel","--noClustering");
    LocalFeatures lf;
    
    DBG(10) << "Loading clustermodel from '" << filename << "'." ;
    BaseClusterer* clusterer=new EM();
    clusterer->loadModel(filename);
    BLINK(10) << "..done" << endl;
    
    ImageFeature img; 
    
    if(cl.search("--histogramization")) {
      DBG(10) << "Starting histogramization" << endl;
      string suffix=cl.follow("histo.gz","--suffix");
      
      uint nOfClusters=clusterer->numberOfClusters();
      
      for(uint i=0;i<filenames.size();++i) {
        
        string fn;
        if(marking|| savePatches) {
          fn=string(filenames[i],0,filenames[i].size()-cropCharacters);
          DBG(10) << "Loading image " << fn << " for marking or patching." << endl;
          img.load(fn);
        }
        

        DoubleVector dists; // these are needed only if fuzzy,
                            // otherwise this is ignored
        VectorFeature fuzzyHisto(nOfClusters);
        
        
        DBG(10) << "Histogramization of " << filenames[i] ;
        HistogramFeature histo(nOfClusters); histo.min()=vector<double>(1,0.0); histo.max()=vector<double>(1,nOfClusters);
        lf=LocalFeatures();
        lf.load(filenames[i]);
        if(discardDim) {lf.discardDimension(discardDimension);}

        
        ogzstream rasts;
        if(writeFeatClsFiles) {

          string rastfn=string(filenames[i])+string(".")+string(cl.follow("featcls","--featclssuffix"));
          DBG(10) << VAR(rastfn) << endl;
          rasts.open(rastfn.c_str());
        }
        
        FeatureExtractionPosition oldpos; oldpos.x=-1; oldpos.y=-1;
        
        
        for(uint j=0;j<lf.numberOfFeatures();++j) {
          uint cls=clusterer->classify(lf[j],dists);
         
          if(writeFeatClsFiles) {
            if(!accIdPos) {
              rasts << lf.position(j).x << " " << lf.position(j).y << " " << cls << endl;
            } else {
              if(oldpos.x==int(lf.position(j).x) and oldpos.y==int(lf.position(j).y)) {
                rasts << " " <<cls ; 
              } else {
                if(oldpos.x!=-1) {
                  rasts << endl << lf.position(j).x << " " << lf.position(j).y << " " << cls;
                } else {
                  rasts  << lf.position(j).x << " " << lf.position(j).y << " " << cls;
                }
                oldpos.x=lf.position(j).x; oldpos.y=lf.position(j).y;
              }
            }
          }

          if(marking||savePatches) {
            if(markers.count(cls) >0) {
              ostringstream oss;
              Marker m=markers[cls];
              if(!savePatches) {
                if(m.type=='x') { diagcross(img,lf.position(j).x, lf.position(j).y, m.color);
                } else if(m.type=='+') { cross(img,lf.position(j).x, lf.position(j).y, m.color);
                } else { ERR << "Unknown marking type" << endl; }
                box(img,lf.position(j).x, lf.position(j).y,m.color,lf.extractionSize(j));
              } else {
                if(savePatches && (markers[cls].type=='+' || markers[cls].type=='x')) {
                  oss.str();
                  oss <<filenames[i] << "-" << cls << "-" << j << ".png";
                  ImageFeature patch=getPatch(img,lf.position(j).x, lf.position(j).y,lf.winsize());
                  patch.save(oss.str());
                  DBG(10) << "Patch saved to " << oss.str() << endl;
                }
              }
            }
          }
          histo.feedbin(cls);
          if(fuzzy) {
            double minDist=*min_element(dists.begin(),dists.end());
            for(DoubleVector::iterator i=dists.begin();i!=dists.end();++i) {*i-=minDist;}
            for(DoubleVector::iterator i=dists.begin();i!=dists.end();++i) {*i=exp(-(*i)/fuzzyAlpha);}
            double sum=0.0;
            for(DoubleVector::iterator i=dists.begin();i!=dists.end();++i) {sum+=*i;}
            for(DoubleVector::iterator i=dists.begin();i!=dists.end();++i) {*i/=sum;}
            
            for(uint i=0;i<nOfClusters;++i) fuzzyHisto[i]+=dists[i];
          }
        }
        
        if(writeFeatClsFiles) {
          if(accIdPos) {rasts << endl;}
          rasts.close();
        }

        if(fuzzy) {
          double sum=0.0;
          for(uint j=0;j<nOfClusters;++j) {sum+=fuzzyHisto[j];}
          for(uint j=0;j<nOfClusters;++j) {fuzzyHisto[j]/=sum;}
          
          fuzzyHisto.save(filenames[i]+"."+fuzzysuffix);
        } else {
          if(histopath!="") {
            string bname=basename(const_cast<char*>(filenames[i].c_str()));
            histo.save(histopath+"/"+bname+"."+suffix);
          } else {
            histo.save(filenames[i]+"."+suffix);
          }
        }
        


        BLINK(10) << " done" << endl;
        if(marking) {
          string fn=filenames[i]+string(".")+string(cl.follow("positionmark.png","--markSuffix"));
          DBG(10) << "Markimage is saved to " << fn << endl;
          img.save(fn);
        }
      }
      DBG(10) << "Histogramization finished" << endl;
    }
  }
  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}

