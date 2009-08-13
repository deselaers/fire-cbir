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

#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
#include <iostream>
#include "lfposclusteridfeature.hpp"
#include "dist_rast.hpp"
#ifdef HAVE_RAST_LIBRARY
#include "misc.h"
#include "struct.h"
#include "geo.h"
#include "vecmat.h"
#include "rast.h"
#else
#warning RAST matching will not work. library missing.
#endif


double RASTDistance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
#ifdef HAVE_RAST_LIBRARY    
  double result=0.0;
  const LFPositionClusterIdFeature* db=dynamic_cast<const LFPositionClusterIdFeature*>(databaseFeature);
  const LFPositionClusterIdFeature* query=dynamic_cast<const LFPositionClusterIdFeature*>(queryFeature);
    
  if(db and query and db->dim() == query->dim()) {
      
    // setting up the rast matching
    autodel<RastF2D> rast(makeRastF2D());
    rast->set_maxresults(1);
    rast->set_verbose(0);
    rast->set_tolerance(tolerance_);
    rast->set_score_threshold(0.1);
    rast->set_min_q(minq_);
    rast->set_xrange(mindx_,maxdx_);
    rast->set_yrange(mindy_,maxdy_);
    rast->set_arange(amin_,amax_);
    rast->set_srange(minscale_,maxscale_);
    rast->set_lsq(0);
    // if verbosity level is >= 20 make rast verbose, too
    DBGI(20,rast->set_verbose(1));
      
    // provide the query features as model
    for(uint i=0;i<query->size();++i) {
      Features features;
      for(uint d=0;d<query->dim();++d) {
        features.push((*query)[i][d]);
      }
      rast->add_msource(query->position(i).first, query->position(i).second,eps_,features);
    }
      
    // provide the database features as instance
    for(uint i=0;i<db->size();++i) {
      Features features;
      for(uint d=0;d<db->dim();++d) {
        features.push((*db)[i][d]);
     } 
      rast->add_ipoint(db->position(i).first, db->position(i).second,features);
    }
      
    // start the rast matching
    rast->match();
      
    // if verbosity is high enough print the result
    DBGI(15,{
        for(int i=0;i<rast->nresults();i++) {
          DBG(15) << i << " " 
                  << rast->ubound(i) << " " 
                  << rast->lbound(i) << " " 
                  << rast->translation(i,0)<< " " 
                  << rast->translation(i,1)<< " " 
                  << rast->angle(i) << " " 
                  << rast->scale(i) << ::std::endl; 
        }
      });
      
      
    if(rast->nresults()>0) {
      result=query->size()-rast->ubound(0);
    } else {
      result=query->size();
    }
    // destructor
    return result;
  } else {
    ERR << "Features not comparable" << ::std::endl;
    return -1.0;
  }

#else  // no rast available
  return -1.0;
#endif


}

