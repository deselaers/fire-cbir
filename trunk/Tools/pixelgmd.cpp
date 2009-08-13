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
#include "jflib.hpp"
#include <cmath>
#include "gmd.hpp"
#include "getpot.hpp"

using namespace std;



int main(int argc , char** argv) {
  GetPot cl(argc,argv);

  ImageFeature img1;
  img1.load(argv[1],true);


  GaussMixtureDensitySettings settings;
  settings.poolMode=settings.String2PoolModeType(cl.follow("nopooling","--poolmode"));
  settings.disturbMode=settings.String2DisturbModeType(cl.follow("variance","--disturbmode"));
  settings.splitMode=settings.String2SplitModeType(cl.follow("all","--splitmode"));
  settings.maxSplits=cl.follow(2,"--maxsplits");
  settings.minObsForSplit=cl.follow(8,"--minobsforsplit");
  settings.minObs=cl.follow(4,"--minobs");
  settings.iterationsBetweenSplits=cl.follow(10,"--iterationsbetweensplits");
  settings.minVar=cl.follow(0.1,"--minvariance");
  settings.epsilon=cl.follow(0.01,"--epsilon");

  GaussMixtureDensity gmd;
  gmd.settings()=settings;
  
  DataSet data;
  
  uint dim=3;
  for(uint x=0;x<img1.xsize();++x) {
    for(uint y=0;y<img1.ysize();++y) {
      FeatureVector *f=new FeatureVector(dim);
      (*f)[0]=x;
      (*f)[1]=y;
      (*f)[2]=img1(x,y,0);
      data.push_back(f);
    }
  }
  
  gmd.train(data);

  
}
