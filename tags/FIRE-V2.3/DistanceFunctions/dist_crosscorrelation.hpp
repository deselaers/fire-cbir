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
#ifndef __dist_crosscorrelation_hpp__
#define __dist_crosscorrelation_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include "imagefeature.hpp"
#include <iostream>
#include <math.h>

class CrosscorrelationDistance : public BaseDistance {
private: 
  int d_;
public:
  CrosscorrelationDistance(int d=4) : d_(d){
  }

  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
    double result=0.0;
    double maxCorel=0.0;
    
    const ImageFeature* db=dynamic_cast<const ImageFeature*>(databaseFeature);
    const ImageFeature* query=dynamic_cast<const ImageFeature*>(queryFeature);
    
    const ImageFeature DB=*db;
    const ImageFeature Q=*query;

    if(db && query) {
      if(db->xsize() ==query->xsize() && db->ysize() == query->ysize()) {
        if(db->xsize() == db->ysize()) {
          //compare images
          int dim=db->xsize();
          
          int aX,aY;
          double avgQ=0.0, avgDB=0.0;
    
          //get image averages
          for(int y=0;y<dim;++y) {
            for(int x=0;x<dim;++x) {
              avgQ+=Q(x,y,0);
              avgDB+=DB(x,y,0);
            }
          }
          
          avgQ/=double(dim*dim);
          avgDB/=double(dim*dim);
          
          double sumQ, sumDB;
          double tmp;
          
          for(int n=-d_;n<=d_;++n) {  // maximising over m and n
            for(int m=-d_;m<=d_;++m) {
              
              double corel=0.0; sumDB=0.0; sumQ=0.0;
              
              for(int y=0;y<dim;++y) {
                aY=y+n;
                for(int x=0;x<dim;++x) {
                  aX=x+m;
                  
                  tmp=Q(y,x,0)-avgQ;
                  tmp*=tmp;
                  sumQ+=tmp;
                  
                  if(aX>=0&&aX<dim && aY>=0&&aY<dim) {
                    corel+=(DB(aY,aX,0)-avgDB)*(Q(y,x,0)-avgQ);
                    
                    tmp=DB(aY,aX,0)-avgDB;
                    tmp*=tmp;
                    sumDB+=tmp;
                  }
                }
              }
              
              corel/=sqrt(sumDB*sumQ);
              if(maxCorel<corel) {
                maxCorel=corel;
              }
            }
          }
          result=1-maxCorel;
          
        } else {
          ERR << "only square images." << ::std::endl;
          result=-1.0;
        }
      } else {
        ERR << "only images of same size." << ::std::endl;
        result=-1.0;
      }
    } else {
      ERR << "only image features" << ::std::endl;
      result=-1.0;
    }
    return result;
  }

  virtual ::std::string name() {return "crosscorrelation";}
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
};

#endif
