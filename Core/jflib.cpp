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
#include "jflib.hpp"
#include "gzstream.hpp"
#include <fstream>
#include <iostream>
#include <ios>
using namespace std;

pair<DoubleVectorVector,IntVector> readJF(string filename) {
  igzstream is; is.open(filename.c_str());

  if(!is.good()) {
    ERR << "Reading from " << filename << endl;
    exit(20);
  }
  
  DoubleVectorVector resultData;
  IntVector resultClasses;
  int maxCls;
  uint dim;
  int cls=0;
  DoubleVector *tmpDblVec;
  if(is.good() && is) {
    is >> maxCls >> dim;
    while(cls != -1) {
      double tmpDbl;
      tmpDblVec=new DoubleVector(dim);
      is >> cls;
      if(cls !=-1) {
        for(uint i=0;i<dim;++i) {
          is >> tmpDbl;
          (*tmpDblVec)[i]=tmpDbl;
        }
        resultClasses.push_back(cls);
        resultData.push_back(tmpDblVec);
      }
    }
    is.close();
  }  
  return pair<DoubleVectorVector, IntVector>(resultData, resultClasses);
}

void writeJF(int maxCls, DoubleVectorVector pixelData, IntVector classData, string filename) {

  ogzstream os; os.open(filename.c_str(), ios_base::out);
  if(!os.good()) {
    ERR << "Writing to " << filename << endl;
    exit(20);
  }

  uint vsize = pixelData[0]->size();
  os << maxCls << " " << vsize << endl;

  for (uint j = 0; j < pixelData.size(); j++) {
    os << classData[j];
    DoubleVector* dv = pixelData[j];
    for (uint i = 0; i < dv->size(); i++) {
      os << " "<< (*dv)[i];
    }
    os << endl;
  }

  os << -1 << endl;
  os.close();

}
