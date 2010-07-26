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

#include <limits>
#include "dist_idm.hpp"
#include "imagelib.hpp"
#include "dist_euclidean.hpp"
using namespace std;

double ImageDistortionModelDistance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
  double dist=-10.0;
    
  const ImageFeature* db=dynamic_cast<const ImageFeature*>(databaseFeature);
  const ImageFeature* query=dynamic_cast<const ImageFeature*>(queryFeature);
    
  if(db && query) {
    dist=0.0;
    const ImageFeature DB=*db;
    const ImageFeature Q=*query;
    double bestDist;
    
    ImageFeature db=DB;
    ImageFeature query=Q;
    ImageFeature tmpImg;
    
    if(sobel_) {
      DBG(35) << "Sobelfiltering...";
      
      tmpImg=db;
      sobelv(db);
      sobelh(tmpImg);
      db.append(tmpImg);        
      
      tmpImg=query;
      sobelv(query);
      sobelh(tmpImg);
      query.append(tmpImg);
      
      for(uint x=0;x<db.xsize();++x) {
        for(uint y=0;y<db.ysize();++y) {
          db(x,y,0)+=4.0; db(x,y,0)/=8.0;
          db(x,y,1)+=4.0; db(x,y,1)/=8.0;}}
      
      for(uint x=0;x<query.xsize();++x) {
        for(uint y=0;y<query.ysize();++y) {
          query(x,y,0)+=4.0; query(x,y,0)/=8.0;
          query(x,y,1)+=4.0; query(x,y,1)/=8.0;}}
      

      BLINK(35) << "done" << endl;
    }

    // deformation of reference image, that is search for optimum
    // position in reference image

    double tmp;

    if(db.xsize()==query.xsize() && db.ysize()==query.ysize()) {
      DBG(35) << "comparing images of same size" << endl;

      for(int y=0;y<int(query.ysize());++y) {
        for(int x=0;x<int(query.xsize());++x) {
          
          bestDist=numeric_limits<double>::max();
          
          for(int xx=x-int(wr1_);xx<=x+int(wr1_);++xx) {
            if(xx>=0 && xx<int(db.xsize())) {
              for(int yy=y-int(wr1_);yy<=y+int(wr1_);++yy) {
                if(yy >= 0 && yy < int(db.ysize())) {
                  tmp=pixelDist(query,db,x,y,xx,yy);
                  if(tmp<bestDist || (tmp==bestDist && xx==x && yy==y)) {
                    bestDist=tmp;
                  }
                }
              }
            }
          }
          dist+=bestDist;
        }
      }
    } else {
      DBG(35) << "comparing images of different sizes" << endl;
      for(int y=0;y<int(query.ysize());++y) {
        int db_y=int( double((db.ysize()-1)*y)/double(query.ysize()-1)+0.5);
        for(int x=0;x<int(query.xsize());++x) {
          int db_x=int( double((db.xsize()-1)*x)/double(query.xsize()-1)+0.5);
          bestDist=numeric_limits<double>::max();

          for(int xx=db_x-int(wr1_);xx<=db_x+int(wr1_);++xx) {
            if(xx>=0 && xx<int(db.xsize())) {
              for(int yy=db_y-int(wr1_);yy<=db_y+int(wr1_);++yy) {
                if(yy >= 0 && yy < int(db.ysize())) {
                  tmp=pixelDist(query,db,x,y,xx,yy);
                  if(tmp<bestDist || (tmp==bestDist && xx==x && yy==y)) {
                    bestDist=tmp;
                  }
                }
              }
            }
          }
          dist+=bestDist;
        }
      }
    }
  }
  return dist;
}

double ImageDistortionModelDistance::pixelDist(const ImageFeature& query,
                                               const ImageFeature& db, 
                                               const uint xquery, const uint yquery, const uint xdb, const uint ydb) const {
  
  double dist=0.0, tmpDist, aktDist;
  int queryx, queryy, dbx, dby;
  int cnt;
  
  for(uint c=0;c<query.zsize();++c) {
    tmpDist=0.0;
    cnt=0;

    for(int i=-int(wr2_);i<=int(wr2_);++i) {
      dby=ydb+i;  
      queryy=yquery+i;
      if(dby>=0 && dby <int(db.ysize()) && queryy>=0 && queryy<int(query.ysize())) {
        for(int j=-int(wr2_);j<=int(wr2_);++j) {
          dbx=xdb+j;
          queryx=xquery+j;
          if(dbx>=0 && dbx <int(db.xsize()) &&  queryx>=0 && queryx <int(query.xsize())) {
            aktDist=((query)(queryx,queryy,c)-(db)(dbx,dby,c));
            aktDist*=aktDist;
            tmpDist+=aktDist;
            cnt+=1;
          }
        }
      }
    }
    tmpDist/=cnt;
    dist+=tmpDist;
  }

  if(threshold_>0) {
    return min(threshold_, dist);
  } else {
    return dist;
  }
}
