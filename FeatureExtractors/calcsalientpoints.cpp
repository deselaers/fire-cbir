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

#include "salientpoints.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include "point.hpp"
#include <stdlib.h>
#include "stdio.h"

using namespace std;

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("Usage: salienpoints <image file> <# of salient points>\n");
    exit(1);
  }
  ImageFeature imgFeature;
  imgFeature.load(argv[1]);
  SalientPoints sp(imgFeature);
  
  vector<Point> points = sp.getSalientPoints(atoi(argv[2]));


  vector<double> col(3,1.0);
  

  DBG(10) << points.size() << " salient points extracted, wanted " << atoi(argv[2]) << endl;
  for (uint i = 0; i < points.size(); i++) {
    Point p = points[i];
    DBG(10) <<"(" << p.x << ", " << p.y << ")" << endl;
    cross(imgFeature,p.x,p.y,col,2);
  } 
  
  imgFeature.display();

}
