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
#ifndef __dist_textfeature_hpp__
#define __dist_textfeature_hpp__

#include "basedistance.hpp"
#include <map>

class TextFeatureDistance : public BaseDistance {
public:

  virtual ~TextFeatureDistance(void) {}

  TextFeatureDistance(::std::string server = "127.0.0.1",
                      unsigned port = 4242,
                      ::std::string language = "None") {
    server_ = server;
    port_ = port;
    language_ = language;
  }

  virtual double distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature);
  virtual void query_wmir(const ::std::string& _query);

 
  virtual ::std::string name() {return "textfeature";}
  virtual ::std::string language() {return language_;}
  virtual void start(const BaseFeature *);
  virtual void stop(){}

  virtual void getServerSettings(::std::string &server, unsigned &port, ::std::string &language);
  
private:
  ::std::map< ::std::string, double> rsv_table_;
  double max_rsv_;
  ::std::string queryfile_;
  ::std::string server_;
  unsigned port_;
  ::std::string language_;
};

#endif
