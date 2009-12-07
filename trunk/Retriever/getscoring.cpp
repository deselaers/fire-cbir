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

#include "basescoring.hpp"
#include "maxentscoring.hpp"
#include "linearscoring.hpp"
#include "svmscoring.cpp"
#include "getscoring.hpp"
#include "maxentscoringsecondorder.hpp"
#include "maxentscoringfirstandsecondorder.hpp"
#include <string>
#include "factory.hpp"

using namespace std;

BaseScoring * getScoring(const ::std::string& scoringname,const uint numberOfDistances) {
  BaseScoring *result=NULL;
  ::std::string configfile="";

  uint pos=scoringname.find("CONFIG=");
  if(pos<scoringname.size()) {
    for(uint i=pos+7; (scoringname[i]>='0' and scoringname[i]<='9') or (scoringname[i]>='a' and scoringname[i]<='z') or (scoringname[i] >= 'A' and scoringname[i]<='Z') or (scoringname[i]=='.') or (scoringname[i]=='/') or (scoringname[i]=='-') and i<scoringname.size();++i) {
      configfile+=scoringname[i];
    }
  }
  
  Factory<BaseScoring,BaseScoring* (*)(const ::std::string&, const uint&) ,string> factory;
  factory.registerClass("linear",BaseScoring::create<LinearScoring>);
  factory.registerClass("maxent",BaseScoring::create<MaxEntScoring>);
  factory.registerClass("maxent2nd",BaseScoring::create<MaxEntSecondOrderScoring>);
  factory.registerClass("maxent1st2nd",BaseScoring::create<MaxEntFirstAndSecondOrderScoring>);
  factory.registerClass("svm",BaseScoring::create<SvmScoring>);
  
  string realscoringname=scoringname.substr(0,scoringname.find(":"));
  DBG(10) << VAR(realscoringname) << endl;
  result=factory.getObject(realscoringname,configfile,numberOfDistances);
  
  if(!result) {
    ERR << "Scoring '" << realscoringname << "' cannot be resolved by da factory. Returning LinearScoring." << ::std::endl;
    result=new LinearScoring(numberOfDistances);
  }

  DBG(10) << "Scoring: " << result->type() << endl;
  return result;
}

const ::std::string listScorings() {
  return "linear maxent maxent2nd maxent1st2nd svm";
}
