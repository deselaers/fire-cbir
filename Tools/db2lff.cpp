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

/**
 * A program to convert any image retrieval database in FIRE format to
 * JF (joerg file). But only vectorfeatures (e.g. vectors, images,
 * histograms are considered)
 *
 * This program is necessary if e.g. experiments with maximum entropy
 * have to be made with features from FIRE.
 */


#include <iostream>
#include <string>

#include "vectorfeature.hpp"
#include "getpot.hpp"
#include "database.hpp"
#include "largefeaturefile.hpp"
#include "jflib.hpp"

using namespace std;

void USAGE()
{
  cout << "USAGE for db2lff -- convert any FIRE database to a FIRE database in LargeFeatureFile-Format" << endl
  << "    db2lff <options> --filelist <filelist>" << endl
  << "  Options: " << endl
  << "    -h   show this help and exit" << endl
  << endl;
}

int main(int argc, char **argv)
{
  GetPot cl(argc,argv);

  if(cl.search("-h") || !cl.search("--filelist")) {USAGE(); exit(20);}

  string filelistfilename=cl.follow("filelist","--filelist");


  Database db;
  db.loadFileList(filelistfilename);
  db.loadFeatures();

  for(uint j=0;j<db.numberOfSuffices();++j)
  {
    string filename=db.path()+"/"+db.suffix(j)+".lff";
    DBG(10) << "Writing '" << filename << "'." << endl;
    LargeFeatureFile lff(filename,db.suffix(j),"kein Kommentar");
    for(uint i=0;i<db.size();++i)
    {
      lff.writeNext(db[i],j);
    }
    lff.closeWriting();
  }
}
