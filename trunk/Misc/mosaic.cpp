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
#include "diag.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include <vector>
#include "stdlib.h"
#include "getpot.hpp"
using namespace std;

void USAGE() {
  cout << "USAGE: mosaic [options] --images image1 [image2 [image3 ... ] ]" << endl
       << "  available options: " << endl
       << "     --winx   how many windows in x direction (10)" << endl
       << "     --winy   how many windows in y direction (10)" << endl
       << "     --borderwidth  how much border between windows in output image (0)" << endl
       << "     --positionDeviation  how strong can patches in original image deviate from calculated position (20)" << endl
       << "     --affineDeviation  how much can the affine transformation deviate from identity (0.1) " << endl
       << "     --output  where is the output saved (output.png)" << endl
       << endl;
}


double randomDouble(const double max) {
  return (max*rand())/double(RAND_MAX+1);
}

int randomInt(const int max) {
  uint r=rand();
  double rmax=double(RAND_MAX)+1;
  r=int(double(r)*double(max)/double(rmax));
  return r;
}

int main(int argc, char**argv) {
  GetPot cl(argc,argv);
  srand((unsigned)time(0));
 
  vector<ImageFeature> inputimages;
  if (!cl.search("--images")) {
    USAGE();
    exit(20);
  }

  cl.search("--images");
  string filename=cl.next(" ");
  while(filename!=" ") {
    DBG(10) << "Loading " << filename << endl;
    ImageFeature tmp;
    tmp.load(filename);
    inputimages.push_back(tmp);
    filename=cl.next(" ");
  }
  
  uint dimx=inputimages[0].xsize();
  uint dimy=inputimages[0].ysize();
  uint dimz=inputimages[0].zsize();
  
 
  uint WINX=cl.follow(10,"--winx");
  uint WINY=cl.follow(10,"--winy");
  uint BORDERWIDTH=cl.follow(0,"--borderwidth");
  double AFFDEV=cl.follow(0.1,"--affineDeviation");
  uint winWidth=dimx/WINX;
  uint winHeight=dimy/WINY;
  uint POSDEV=cl.follow(20,"--positionDeviation");
  uint POSDEV2=POSDEV/2;
  if(POSDEV > winWidth/2 || POSDEV>winHeight/2) {
    POSDEV=min(winWidth/2,winHeight/2);
    DBG(10) << "Setting posdev=" << POSDEV << endl;
  }
  
  ImageFeature outputimage(dimx+(WINX+1)*BORDERWIDTH, dimy+(WINY+1)*BORDERWIDTH, dimz);

  DBG(10) "Patching the new image" << endl;
  for(uint winx=0;winx<WINX;++winx) {
    for(uint winy=0;winy<WINY;++winy) {
      
      uint sourceImage=randomInt(inputimages.size());
      BLINK(30) << " " << sourceImage << flush;
      int top=max(int(randomInt(POSDEV)+int(winy*winHeight)-int(POSDEV2)),0);
      int bottom=max(randomInt(POSDEV)+int((winy+1)*winHeight)-int(POSDEV2),0);
      int left=max(randomInt(POSDEV)+int(winx*winWidth)-int(POSDEV2),0);
      int right=max(randomInt(POSDEV)+int((winx+1)*winWidth)-int(POSDEV2),0);
      
      ImageFeature patch=getPatch(inputimages[sourceImage],left,top,right,bottom);
      
      AffineTransformation aff(1+randomDouble(AFFDEV),randomDouble(AFFDEV), randomDouble(AFFDEV), 1+randomDouble(AFFDEV),0,0);
      patch=affineTransformation(patch,aff);
      
      setPatch(outputimage,left+(winx+1)*BORDERWIDTH,top+(winy+1)*BORDERWIDTH,patch);
    }
    BLINK(30) << endl;
  }

  outputimage.save(cl.follow("output.png","--output"));
}
