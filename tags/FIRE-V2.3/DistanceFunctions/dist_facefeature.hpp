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
#ifndef __dist_facefeature_hpp__
#define __dist_facefeature_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include "histogrampairfeature.hpp"
#include <iostream>
#include <limits>
#include "facefeature.hpp"

class FaceFeatureDistance : public BaseDistance {

private:
  bool useAll;

public:
  FaceFeatureDistance(::std::string allOrOne="one") {
    if(allOrOne=="one") {
      useAll=false;
      DBG(10) << "USE=one" << ::std::endl;
    } else if(allOrOne=="all") {
      useAll=true;
      DBG(10) << "USE=all" << ::std::endl;
    } else {
      ERR << "Unknown parameter: " << allOrOne << " assuming one." << ::std::endl;
      useAll=false;
      DBG(10) << "USE=one" << ::std::endl;
    }
  }
  
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
    double result;
    
    const FaceFeature* db=dynamic_cast<const FaceFeature*>(databaseFeature);
    const FaceFeature* query=dynamic_cast<const FaceFeature*>(queryFeature);
    
    if(db && query) {
      DBG(55) << VAR(query->numberOfVectors()) << " " << VAR(db->numberOfVectors()) << ::std::endl;
      if(query->numberOfVectors() ==0) {
        // no query face
        result=0.0;
      } else {
        
        //using all faces in query image

        if(useAll) {
          if(db->numberOfVectors()==0) {
            //current database image has no face
            result=0.0;
            for(uint i=0;i<query->numberOfVectors();++i) {
              double tmp=0.0;
              const ::std::vector<double> &f=query->operator[](i);
              for(uint i=0;i<f.size();++i) {
                tmp+=f[i]*f[i];
              }
              result+=tmp;
            }
          } else {
            //database and query image have faces

            result=0.0;
            for(uint i=0;i<query->numberOfVectors();++i) {
              const ::std::vector<double> &f=query->operator[](i);

              double tmpres=::std::numeric_limits<double>::max();
              
              for(uint j=0;j<db->numberOfVectors();++j) {
                double tmp=0.0;              
                
                const ::std::vector<double> &fd=db->operator[](j);
                for(uint d=0;d<f.size();++d) {
                  double t=f[d]-fd[d];
                  tmp+=t*t;
                }
                if(tmp<tmpres) tmpres=tmp;
              }
              result+=tmpres;
            }
          }          
        } else {
          // using only best match from query image
          if(db->numberOfVectors()==0) {
            //current database image has no face
            result=0.0;
            for(uint i=0;i<query->numberOfVectors();++i) {
              double tmp=0.0;
              const ::std::vector<double> &f=query->operator[](i);
              for(uint i=0;i<f.size();++i) {
                tmp+=f[i]*f[i];
              }
              if(tmp>result) result=tmp;
            }
          } else {
            //database and query image have faces
            result=::std::numeric_limits<double>::max();
            for(uint i=0;i<query->numberOfVectors();++i) {
              const ::std::vector<double> &f=query->operator[](i);
              
              for(uint j=0;j<db->numberOfVectors();++j) {
                double tmp=0.0;              
                
                const ::std::vector<double> &fd=db->operator[](j);
                for(uint d=0;d<f.size();++d) {
                  double t=f[d]-fd[d];
                  tmp+=t*t;
                }
                if(tmp<result) result=tmp;
              }
            }
          }
        }
      }
      return result;
    } else {
      ERR << "Features not comparable" << ::std::endl;
      return -1 ;
    }
  }
  
  virtual ::std::string name() {return "facefeature";}
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
};
  
#endif
