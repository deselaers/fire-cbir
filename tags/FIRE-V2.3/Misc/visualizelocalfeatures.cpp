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
#include "diag.hpp"
#include "getpot.hpp"
#include "localfeatures.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include "pca.hpp"
using namespace std;


void USAGE() {
  cout << "USAGE: visualizelocalfeatures <options>" << endl
       << "   -lf <localfeaturesfile>" << endl
       << "   -layerwise" << endl
       << "   -nonorm           (performs no normalization)" << endl
       << "   -scale <n>        (scales the local feature by the integer factor n)" << endl
       << "   -save <filename>  (saves the images to <filename>_$i.png)" << endl
       << "   -unpca <filename> (backtransforms image with the given pca file)" << endl
       << "   -dontshow         (suppresses display)" << endl
       << "   -h,--help    show this help" << endl
       << endl;
}

int main(int argc, char**argv) {
  GetPot cl(argc, argv);

  if(cl.search(2,"-h","--help") or !cl.search("-lf")) {USAGE(); exit(0);}

  bool nonorm = cl.search(1, "-nonorm");
  LocalFeatures lf; lf.load(cl.follow("test.lf.gz","-lf"));
  uint sc = cl.follow(1, "-scale");
  uint w=uint(lf.winsize()*2+1);
  ImageFeature img(w,w,lf.zsize());

  string savefilename = cl.follow("", "-save");

  bool backtransform = false;
  PCA pca;
  if(cl.search("-unpca"))
  {
	  if( !pca.load(cl.follow("", "-unpca")) )
	  {
		  ERR << "Error loading PCA file!" << endl;
		  abort();
	  }
	  backtransform = true;
  }

  bool dontshow = cl.search("-dontshow");

  for(uint i=0;i<lf.size();++i) {
    DBG(10) << "lf " << i << endl;
    if(!backtransform) {
    	for(uint j=0;j<lf.dim();++j) {
    		img[j] = lf.getData()[i][j];
    	}
    } else {
    	vector<double> backtransformed = pca.backTransform(lf.getData()[i]);
    	for(uint j=0;j<(uint)pca.dim();++j) {
    		img[j] =  backtransformed[j];
    	}
    }
    if (sc != 1) {
      img = scale(img, w * sc, w * sc);
    }

    if (!nonorm) {
      normalize(img);
    }
    if(!dontshow) {
    	if(cl.search("-layerwise")) {

    		for(uint c=0;c<img.zsize();++c) {
    			BLINK(10) << " " << c;
    			img.display(c,c,c);
    		}
    		BLINK(10) << endl;
    	} else {
    		BLINK(10) << "all layers" << endl;;
    		img.display();
    	}
    }

    if (savefilename != "") {
    	std::ostringstream filenamewithnr;
    	filenamewithnr << savefilename << "_" << setw(3) << setfill('0') << i << ".png";
    	img.save(filenamewithnr.str());
    }
  }
}
