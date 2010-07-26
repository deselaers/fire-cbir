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
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#ifndef __svmscoring_hpp__
#define __svmscoring_hpp__
#include "diag.hpp"
#include "basescoring.hpp"
#include "supportvectormachine.hpp"

class SvmScoring: public BaseScoring {
 
private:

  SupportVectorMachine svm_;

public:
  SvmScoring(const ::std::string &filename){type_="svm"; load(filename);}     
  SvmScoring(const ::std::string &filename, const uint&){
    type_="svm"; if(filename=="") {
      DBG(10) << "need configfile" << ::std::endl;} else{ load(filename);
    }
  }
  
  void load(const ::std::string &filename){
    svm_.load(filename);
  }
  
  virtual double getScore(const ::std::vector<double>& dists);

  virtual const ::std::string settings();
};

#endif
