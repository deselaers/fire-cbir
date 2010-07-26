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
#include <vector>
#include "getpot.hpp"
#include "gzstream.hpp"
#include "diag.hpp"
#include "facefeature.hpp"
#include "pca.hpp"
#include "imagelib.hpp"
#include "imagefeature.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "processfacefeatures [options] (--features filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "    -h, --help   show this help" << endl
       << "    -s, --suffix suffix of output files" << endl
       << endl
       << "    -estimatepca <pcafile> load featurefiles and estimate pca matrix" << endl
       << endl
       << "    -transformpca <pcafile> load pcafile and transform features" << endl
       << "        needs: -dim <dim> for dimension to which we want to transform" << endl
       << endl
       << "    -normalize        normalize all images" << endl 
       << "    -histnormalize    histogram normalize all features" << endl
       << "    -visualize        save face features as images" << endl
       << "    -crop <cropdim>   crop <cropdim> pixels in each direction" << endl
       << "       need: -width <w> -height <h>" << endl
       << endl
       << "    -mark to mark where faces have been extracted" << endl
       << "       needs -width <w> -height <h> -feature <featurefile> -image <imagefile>" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  //command line parsing
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();
  
 
  //get list of files to be processed
  vector<string> infiles;

  string suffix; // set default value here !
  if (cl.search(2, "--suffix", "-s")) {
    suffix = cl.next(" ");
  }

  if(cl.search("-mark")) {
    string ffile=cl.follow("feature.facefeat.gz","-feature");
    string ifile=cl.follow("image.png","-image");
    
    ImageFeature img; img.load(ifile);
    FaceFeature feat; feat.load(ffile);
    vector<double> color(3); color[0]=color[1]=color[2]=1.0;
    
    for(uint i=0;i<feat.nOfFaces();++i) {
      Face &f =feat.face(i);
      
      rect(img,f.posx(),f.posy(),f.width(),f.height(),color);
    }
    img.save(ifile+"."+suffix);
    exit(0);
  }

  if(cl.search("--images")) {
    string filename=cl.next(" ");;
    while(filename!=" ") {
      infiles.push_back(filename);
      filename=cl.next(" ");
    }
  } else if (cl.search("--filelist")) {
    string filename="test";
    igzstream ifs; ifs.open(cl.follow("list","--filelist"));
    if(!ifs.good() || !ifs) {
      ERR << "Cannot open filelist " <<cl.follow("list","--filelist")  << ". Aborting." << endl;
      exit(20);
    }
    while(!ifs.eof() && filename!="") {
      getline(ifs,filename);
      if(filename!="") {
        infiles.push_back(filename);
      }
    }
    ifs.close();
  } else {
    USAGE();
    exit(20);
  }
  
  if(cl.search("-estimatepca")) {

    FaceFeature f; f.load(infiles[0]);
    PCA pca(f.nOfCoefficients());

    // processing the files
    for(uint i=0;i<infiles.size();++i) {
      string filename=infiles[i];
      DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;
      FaceFeature tmp; tmp.load(filename);
      for(uint i=0;i<tmp.nOfFaces();++i) {
        pca.putData(tmp[i]);
      }
      DBG(20) << "Finished with '" << filename << "'." << endl;
    }
    pca.dataEnd();
    pca.calcPCA();
    pca.save(cl.follow("pca.pca","-estimatepca"));
    DBG(10) << "PCA saved to " << cl.follow("pca.pca","-estimatepca") << endl;
  } else if(cl.search("-transformpca")) {
    PCA pca;
    pca.load(cl.follow("pca.pca","-transformpca"));
    DBG(10) << "PCA loaded from " << cl.follow("pca.pca","-transformpca") << endl;


    uint dim=cl.follow(20,"-dim");

    // processing the files
    for(uint i=0;i<infiles.size();++i) {
      string filename=infiles[i];
      DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;
      FaceFeature tmp; tmp.load(filename);
      for(uint i=0;i<tmp.nOfFaces();++i) {
        tmp[i]=pca.transform(tmp[i],dim);
      }
      tmp.nOfCoefficients()=dim;
      tmp.save(filename+"."+suffix);
      DBG(20) << "Finished with '" << filename << "'." << endl;
    }
    
  } else if(cl.search("-normalize")) {
    uint width=cl.follow(16,"-width");
    uint height=cl.follow(16,"-height");
    for(uint i=0;i<infiles.size();++i) {
      string filename=infiles[i];
      DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;
      FaceFeature tmp; tmp.load(filename);
      for(uint i=0;i<tmp.nOfFaces();++i) {
        ImageFeature img(tmp[i],width,height);
        normalize(img);
        tmp[i]=img.layerVector(0);
      }
      tmp.save(filename+"."+suffix);
      DBG(20) << "Finished with '" << filename << "'." << endl;
    }
    
  } else if(cl.search("-histnormalize")) {
    uint width=cl.follow(16,"-width");
    uint height=cl.follow(16,"-height");
    for(uint i=0;i<infiles.size();++i) {
      string filename=infiles[i];
      DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;
      FaceFeature tmp; tmp.load(filename);
      for(uint i=0;i<tmp.nOfFaces();++i) {
        ImageFeature img(tmp[i],width,height);
        histogramNormalization(img);
        tmp[i]=img.layerVector(0);
      }
      tmp.save(filename+"."+suffix);
      DBG(20) << "Finished with '" << filename << "'." << endl;
    }
  } else if(cl.search("-visualize")) {
    uint width=cl.follow(16,"-width");
    uint height=cl.follow(16,"-height");
    for(uint i=0;i<infiles.size();++i) {
      string filename=infiles[i];
      DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;
      FaceFeature tmp; tmp.load(filename);
      ImageFeature img(0,height,1);
      for(uint i=0;i<tmp.nOfFaces();++i) {
        ImageFeature tmpimg(tmp[i],width,height);
        appendHorizontally(img,tmpimg);
      }
      img.save(filename+".png");
      DBG(20) << "Finished with '" << filename << "'." << endl;
    }

  } else if(cl.search("-crop")) {
    uint width=cl.follow(16,"-width");
    uint height=cl.follow(16,"-height");
    uint cropdim=cl.follow(3,"-crop");
    for(uint i=0;i<infiles.size();++i) {
      string filename=infiles[i];
      DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;
      FaceFeature tmp; tmp.load(filename);
      ImageFeature img(0,height,1);
      tmp.nOfCoefficients()=(width-2*cropdim)*(height-2*cropdim);
      for(uint i=0;i<tmp.nOfFaces();++i) {
        ImageFeature tmpimg(tmp[i],width,height);
        shave(tmpimg,cropdim);
        tmp[i]=tmpimg.layerVector(0);

      }
      tmp.save(filename+"."+suffix);
      DBG(20) << "Finished with '" << filename << "'." << endl;
    }
  } else {
    USAGE(); 
    exit(20);
  }



  

  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
