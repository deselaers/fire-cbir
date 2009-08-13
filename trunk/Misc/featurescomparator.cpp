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
#include "distancemaker.hpp"

using namespace std;

/* This is a small program that loads a database of different image
   features and compares them to each other. The purpose of this program
   is for example to check different features of different versions of
   a software when reemplementing it.   */

void USAGE() {
  cout << "featurescomparator [options] --filelist <filelist>" << endl
       << "  Options: " << endl
       << "    -d   [jsd|...] distance used to compare features" << endl
       << "         default: jsd " << endl
       << "    -h   show this help and exit" << endl
       << endl;
}

int main(int argc, char **argv) {
	
  GetPot cl(argc,argv);
  
  if(cl.search("-h") || !cl.search("--filelist")) {USAGE(); exit(20);}
  
  // reading distance's name
  string distanceName = "jsd";
  if (cl.search("--d")) {
    distanceName = cl.next(" ");
  }
  
  // reading files
  string filelistfilename=cl.follow("filelist","--filelist");
  Database db;
  db.loadFileList(filelistfilename);
  db.loadFeatures();
  
  // creating the distance
  double dist;
  BaseDistance* distance;
  DistanceMaker distanceMaker;
  distance = distanceMaker.makeDistance(distanceName);
  
  // iterating over all pictures and displaying the distances
  for(uint i=0;i<db.size();++i) {
    cout << "picture no " << i << endl;
    for(uint j=0;j<db.numberOfSuffices();++j) {
      for(uint k=j+1;k<db.numberOfSuffices();++k) {
        dist=distance->distance(db[i]->operator[](j),db[i]->operator[](k));
        cout << "comparing feature " << j << " to feature " << k << " with result: " << dist << endl;        
      }
    }
  }

}
