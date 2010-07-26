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
#include "imagefeature.hpp"
#include "getpot.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "  colororgray [options] filename"
       << "    -- to check whether the given file is a color or a gray image." << endl
       << "  Available options are:"<< endl
       << "    -q,--quick <nr (default=20)> check only nr pixels" << endl
       << endl;
}


bool checkPixel(const ImageFeature& im,uint x, uint y) {
  bool result=true;
  for(uint z=1;z<im.zsize()&&result;++z) {
    if(im(x,y,z)!=im(x,y,0)) result=false;
  }
  return result;
}

int main(int argc, char**argv) {
  GetPot cl(argc, argv);
  
  ImageFeature im;
  im.load(cl[cl.size()-1]);
  
  bool result=true;

  if(cl.search(2,"--quick","-q")) {
    srand((unsigned)time(0));

    uint nr=cl.follow(20,2,"--quick","-q");
    for(uint i=0;i<nr&&result;++i) {
      uint x=rand()%im.xsize();
      uint y=rand()%im.ysize();
      
      if(!checkPixel(im,x,y)) result=false;
    }
  } else {
    for(uint x=0;x<im.xsize()&&result;++x) {
      for(uint y=0;y<im.ysize()&&result;++y) {
        if(!checkPixel(im,x,y)) result=false;
      }
    }
  }
  
  if(result) {
    cout << "GRAY" << endl;
  } else {
    cout << "COLOR" << endl;
  }
}
