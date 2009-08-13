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
#include "diag.hpp"
#include "getpot.hpp"

#include "pca.hpp"
#include "vectorfeature.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"

using namespace std;


void USAGE() {
  cout << "USAGE: pcavectortoimage <options>" << endl
       << "   -h,--help            show this help" << endl
       << "   -p,--pca <filename>  load the pca transformation from the specified file" << endl
       << "   -b,--base <n>        visualize the n-th pca-vector" << endl
       << "   -M,--mean            visualize the pca mean vector" << endl
       << "   -m,--minMean         subtract mean after backtransformation" << endl
       << "                        if visualizing pca vectors this should be used" << endl
       << "   -v,--vec <filename>  load the vector to be transformed from the specified file" << endl
       << "      -b, -v, and -M are mutually exclusive" << endl
       << "   -i,--img <filename>  save the backtransformed vector to the given image " << endl
       << "                                 (squared if no dimensions specified)" << endl
       << "   -n,--norm            normalize the output image" << endl
       << "   -w,--width <width>   the output image has the specified width" << endl
       << "   -h,--height <height> the output image has the specified height" << endl
       << "   -c,--color           the output image is a color image (depth=3) otherwise gray (depth=1)" << endl
       << "   --nopca              the vector just has to be converted to an image, and is NOT pca transformed" << endl
       << "        (take care, this has to be consistent with the size of the pca, i.e. w*h=dim" << endl
       << "           otherwise strange things may happen.)" << endl
       << endl;
}

int main(int argc, char**argv) {
  GetPot cl(argc, argv);
  
  if(cl.search(2,"-h","--help")) {USAGE(); exit(0);}
  string pcafile=cl.follow("pca.pca",2,"-p","--pca");
  string imagefile=cl.follow("image",2,"-i","--img");

  VectorFeature vec;
  PCA pca; pca.load(pcafile);
  
  uint width, height, depth;
  if(cl.search(2,"-c","--color")) {
    depth=3;
  } else {
    depth=1;
  }
    
  if(!cl.search(2,"-w","--width") && !cl.search(2,"-h","--height")) {
    width=uint(sqrt(double(pca.dim())));
    height=width;
    if(int(height*width)!=pca.dim()) {
      ERR << "pca not for squared images, specify width or height" << endl
          << "height=" << height << "* width=" << width << "!= size=" << pca.dim() << endl;
      exit(20);
    } 
  } else {
    if(cl.search(2,"-w","--width") && !cl.search(2,"-h","--height")) {
      width=cl.follow(10,2,"-w","--width");
      height=pca.dim()/width;
      if(int(height*width)!=pca.dim()) {
        ERR << "pca images of this width, specify valid values" << endl
            << "height=" << height << "* width=" << width << "!= size=" << pca.dim() << endl;
        exit(20);
      }
    } else if(!cl.search(2,"-w","--width") && cl.search(2,"-h","--height")) {
      height=cl.follow(10,2,"-j","--height");
      width=pca.dim()/height;
      if(int(height*width)!=pca.dim()) {
        ERR << "pca images of this height, specify valid values" << endl
            << "height=" << height << "* width=" << width << "!= size=" << pca.dim() << endl;
        exit(20);
      }
    } else {
      height=cl.follow(10,2,"-j","--height");
      width=cl.follow(10,2,"-w","--width");
      if(int(height*width)!=pca.dim()) {
        ERR << "pca images of this height and width, specify valid values" << endl
            << "height=" << height << "* width=" << width << "!= size=" << pca.dim() << endl;
        exit(20);
      }
    }
  }



  vector<double> backtransformed;
  
  if(cl.search(2,"-v","--vec") && !cl.search(2,"-b","--base") && !cl.search(2,"-M","--mean")) {
    string vecfile=cl.follow("vec.vec",2,"-v","--vec");
    DBG(10) << "Loading Vectorfile " << vecfile << endl;
    vec.load(vecfile);    
    DBG(10) << "Vector to be backtransformed" ;
    for(uint i=0;i<vec.size();++i) {
      BLINK(10) << " "<< vec[i];
    }
    BLINK(10) << endl;

    if(cl.search("-n1st"))  vec[0]=0;

    if(cl.search("--nopca")) {
      backtransformed=vec.data();
    } else {
      backtransformed=pca.backTransform(vec.data());
    }
    
    DBG(10) << "Backtransformed Vector" ;
    for(uint i=0;i<backtransformed.size();++i) {
      BLINK(10) << " " <<backtransformed[i];
    }
    BLINK(10) << endl;

  } else if(cl.search(2,"-b","--base") && !cl.search(2,"-v","--vec") && !cl.search(2,"-M","--mean")) {
    uint base=cl.follow(0,2,"-b","--base");
    backtransformed=pca.eigenvector(base);
  } else if(!cl.search(2,"-b","--base") && !cl.search(2,"-v","--vec") && cl.search(2,"-M","--mean")) {
    backtransformed=pca.mean();
  } else {
    USAGE();
    exit(20);
  }
  
  ImageFeature image(width,height,depth);
  
  if(cl.search(2,"-m","--minMean")) {
    DBG(10) << "Subtracting mean" << endl;
    for(uint i=0;i<backtransformed.size();++i) {
      backtransformed[i]-pca.mean()[i];
    }
  }
  
  for(uint i=0;i<backtransformed.size();++i) {
    image[i]=backtransformed[i];
  }

  
  if(cl.search(2,"-n","--norm")) {
    DBG(10) << "normalization" << endl;
    normalize(image);
  }
  
  
//  //make values positive
//  shift(image,-minimum(image,0));
//
  cutoff(image);
//

  DBG(10) << "Going to save:" ;
  for(uint i=0;i<image.size();++i) {
    BLINK(10) << " " <<image[i];
  }
  BLINK(10) << endl;



  image.save(imagefile);
  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
