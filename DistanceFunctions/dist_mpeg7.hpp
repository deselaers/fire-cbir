/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

FIRE is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software Foundation, Inc.,
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#ifndef __dist_mpeg7_hpp__
#define __dist_mpeg7_hpp__
                   
#include <map>
#include <string>
#include <iostream>
#include "vectorfeature.hpp"
#include "diag.hpp"
#include "basedistance.hpp"
                   
/** MPEG7Distance
 *
 * This is a special distance function. Here not real features from
 * the fire framework are compared and processed but the XM MPEG7
 * reference software ( at the time of writing this software available
 * from
 * http://www.lis.e-technik.tu-muenchen.de/research/bv/topics/mmdb/e_mpeg7.html,
 * 16. November, 2004 is used to compare these features. 
 * 
 * In the tools directory you should find some scripts which do
 * appropriate feature extraction for these features and create the
 * necessary data, that is: the par and bitstream files.
 */

class MPEG7Distance : public BaseDistance {
private:
  /// path for XMMain.exe
  ::std::string xmmain_;
  
  /** path for the data to be read by XMMain.exe (par, list, and
   * bitstreams), these paths are the filenames without suffices, that
   * if this e.g.
   *
   * /work/deselaers/database/mpeg7/collayout
   *
   * it is assumed that we have the three files
   *
   * /work/deselaers/database/mpeg7/collayout.par
   * /work/deselaers/database/mpeg7/collayout.mp7
   * /work/deselaers/database/mpeg7/collayout.lst
   *
   * otherwise the result might be strange.
   */
  ::std::string mpeg7data_;
  
  /** here we save the distances which we get from parsing the output
   * from XMMain in the start routine. The data here are just "handed
   * out" by the distance function.
   */
  ::std::map< ::std::string, double> distances_;
  
public:
  
  /// constructor: set the path of the XMMain.exe binary and of the
  /// mpeg data to be read this is the bitstream and the parfile
  MPEG7Distance(::std::string xmmain, ::std::string mpeg7data) : xmmain_(xmmain), mpeg7data_(mpeg7data)  {}

  virtual ~MPEG7Distance() {}
  
  /** 
   * return the distances between the two passed mpeg7 features. These
   * are not calculated here, but precalculated in start
   */
  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature);
  
  virtual ::std::string name() {return "mpeg7";}
  
  
  /// initialize the query for the given BaseFeature. This runs the
  /// XMMain.exe from the MPEG7 reference software and parses the
  /// output
  virtual void start(const BaseFeature *queryFeature);
  
  /// forget the distances for the last query which was prepared.
  virtual void stop();
};

#endif
