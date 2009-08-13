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
#ifndef __featureloader_hpp__
#define __featureloader_hpp__

#include <map>
#include <string>
#include "diag.hpp"
#include "basefeature.hpp"
#include "factory.hpp"

class FeatureLoader {
private:
  ::std::map<const ::std::string, FeatureType> map_;
  Factory<BaseFeature,BaseFeature* (*)(),::std::string> featureFactory_;

public:
  FeatureLoader();  
  BaseFeature* load(const ::std::string& basename, const ::std::string& suffix, const ::std::string &lastSuffix, const ::std::string& path);
  BaseFeature *makeNewFeature(const ::std::string& suffix) const ;
  ::std::string relevantSuffix(const ::std::string & suffix) const;
  
  FeatureType suffix2Type(const ::std::string& suffix) const;


};


#endif
