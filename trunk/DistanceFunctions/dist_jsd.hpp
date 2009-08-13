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
#ifndef __dist_jsd_hpp__
#define __dist_jsd_hpp__

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include "sparsehistogramfeature.hpp"
#include <iostream>
#include <cmath>
using namespace std;

class JSDDistance : public BaseDistance {
private:
  double smoothFactor_;

public:
  
  JSDDistance(double sF=0.5) : smoothFactor_(sF) {}
  
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
    
    DBG(30) << "distance function of jsd distance called..." << std::endl;

    double result=0.0;
    
    const VectorFeature* db=dynamic_cast<const VectorFeature*>(databaseFeature);
    const VectorFeature* query=dynamic_cast<const VectorFeature*>(queryFeature);

    if(db && query) {


#ifdef OLD_SLOW_JSD
      double tmp;
      double n1,n2;
      if (query->size() != db->size()) ERR << "qsize != dbsize" << ::std::endl;
      uint s = query->size();

      for(uint i=0;i<s;++i) {
        n1=(db->operator[](i));
        n2=(query->operator[](i));

        if(n1!=0) {
          tmp=2*n1;
          tmp/=(n1+n2);
          tmp=log(tmp);
          tmp*=n1;
          result+=tmp;
        }
        
        if(n2!=0) {
          tmp=2*n2;
          tmp/=(n1+n2);
          tmp=log(tmp);
          tmp*=n2;
          result+=tmp;
        }
      }
#else
      float tmp1,tmp2;
      float n1,n2,by_n;
      if (query->size() != db->size()) ERR << "qsize != dbsize (lengths of feature vectors are not the same)" << ::std::endl;
      uint s = query->size();                                                                                                
                             
      for(uint i=0;i<s;++i) {
        n1=(float)(db->operator[](i));
        n2=(float)(query->operator[](i));
        if(n1==0 && n2==0) continue;
        by_n=2.0f/(n1+n2);

        if(n1==0){
            tmp2=by_n*n2;
            tmp2=log(tmp2);
            tmp2*=n2;
            result+=tmp2;
        }
        else
            if(n2==0){
                tmp1=by_n*n1;
                tmp1=log(tmp1);
                tmp1*=n1;
                result+=tmp1;
            }
            else{ // interleave independent paths
                tmp2=by_n*n2;
                tmp1=by_n*n1;
                tmp2=log(tmp2);
                tmp1=log(tmp1);
                tmp2*=n2;
                tmp1*=n1;
                result+=tmp1;
                result+=tmp2;
            }
        }
#endif
      return result;
      
    } else {

      const SparseHistogramFeature* db=dynamic_cast<const SparseHistogramFeature*>(databaseFeature);
      const SparseHistogramFeature* query=dynamic_cast<const SparseHistogramFeature*>(queryFeature);
      
      if (db && query) {

	if (smoothFactor_ < 0) {

	  // not implemented yet
	  ERR << "smoothFactor < 0 not implemented yet" << std::endl;
	  exit(1);

	} else if (smoothFactor_ == 0.0) {
	  
	  double tmp,n1,n2; 
	  Position pos; 
	  MapTypeDouble dbMap = db->getMap();
	  MapTypeDouble queryMap = query->getMap();

	  /// iterate over the map from the database histogram,
	  /// get the corresponding values from the query histogram, add
	  /// the according term to the result and delete the entry from the 
	  /// query histogram to be sure to not treat it a second time.
	  for(MapTypeDouble::iterator i = dbMap.begin(); i != dbMap.end(); ++i) {
	    n1 = i->second;
	    pos = i->first;
	    MapTypeDouble::iterator j = queryMap.find(pos);
	    if (j != queryMap.end()) {
	      n2 = j->second;
	      queryMap.erase(j->first);
	    } else {
	      n2 = 0.0;
	    }
	    
	    if (n1 != 0) {
	      tmp = 2 * n1;
	      tmp /= (n1 + n2);
	      tmp = log(tmp);
	      tmp *= n1;
	      result += tmp;
	    }
	    
	    if (n2 != 0) {
	      tmp = 2 * n2;
	      tmp /= (n1 + n2);
	      tmp = log(tmp);
	      tmp *= n2;
	      result += tmp;
	    }
	  }
	  
	  n1 = 0.0;
	  for(MapTypeDouble::iterator j = queryMap.begin(); j != queryMap.end(); ++j) {
	    n2 = j->second;
	    tmp = 2 * n2;
	    tmp /= (n1 + n2);
	    tmp = log(tmp);
	    tmp *= n2;
	    result += tmp;
	  }
	  return result;

	} else {

	  // expand the sparse histogram features and compare them
	  MapTypeDouble expandedDb, expandedQuery;
	  db->expand(expandedDb,smoothFactor_);
	  query->expand(expandedQuery,smoothFactor_);
	  
	  double tmp, n1, n2;
	  Position pos; 
	  MapTypeDouble::iterator j;
	  
	  /// iterate over the expanded map from the database histogram,
	  /// get the corresponding values from the query histogram, add
	  /// the according term to the result and delete the entry from 
	  /// the query histogram to be sure to not treat it a second time
	  for(MapTypeDouble::iterator i = expandedDb.begin(); i != expandedDb.end(); ++i) {
	    n1 = i->second;
	    pos = i->first;
	    j = expandedQuery.find(pos);
	    if (j != expandedQuery.end()) {
	      n2 = j->second;
	      expandedQuery.erase(j->first);
	    } else {
	      n2 = 0.0;
	    }
	    
	    if(n1 != 0) {
	      tmp = 2 * n1;
	      tmp /= (n1 + n2);
	      tmp =log(tmp);
	      tmp *= n1;
	      result += tmp;
	    }
	    
	    if(n2 != 0) {
	      tmp = 2 * n2;
	      tmp /= (n1 + n2);
	      tmp = log(tmp);
	      tmp *= n2;
	      result += tmp;
	    }
	  }
	  
	  n1 = 0.0;
	  for(MapTypeDouble::iterator j = expandedQuery.begin(); j != expandedQuery.end(); ++j) {
	    n2 = j->second;
	    tmp = 2 * n2;
	    tmp /= (n1 + n2);
	    tmp = log(tmp);
	    tmp *= n2;
	    result += tmp;
	  }
	  return result;
	}

      } else {
	ERR << "Features not comparable" << ::std::endl;
	return -1.0;
      }
	  
    }
  }

  virtual ::std::string name() {return "jsd";}
  virtual void start(const BaseFeature *) {}
  virtual void stop(){}
};

#endif
