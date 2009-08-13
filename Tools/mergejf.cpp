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
#include "getpot.hpp"
#include "diag.hpp"
#include "getpot.hpp"
#include "gzstream.hpp"
#include <iostream>

using namespace std;

DoubleVector* collect(const vector<DoubleVector*>& vecs) {
  DoubleVector* dv = new DoubleVector();
  for (uint i = 0; i < vecs.size(); i++) {
    DoubleVector* vec = vecs[i];
    for (uint j = 0; j < vec->size(); j++)
      dv->push_back((*vec)[j]);
  }
  return dv;
}

bool classDataEquals(const IntVector& i1, const IntVector& i2) {
  if (i1.size() != i2.size()) {
    ERR << "Not same number of classes" << endl;
    return false;
  }
  for (uint i = 0; i < i1.size(); i++) {
    if (i1[i] != i2[i]) {
      ERR << "classes different in line " << VAR(i) << ": " << VAR(i1[i]) << " " << VAR(i2[i]) << endl;
      return false;
    }
  }
  return true;
}

void USAGE() {
  cout << "USAGE:" << endl
       << "mergejf [options] (-i inputfile1 [inputfile2 ... ] | -f <filelist>) -o outputfile" << endl
       << "   Options:" << endl
       << "    -h, --help          show this help" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv) {

  GetPot cl(argc, argv);
  string outputJoergfile;

  if (cl.search(2, "-h", "--help")) {
    USAGE();
  }
  vector<string> infiles;

  if(cl.search("-i")) {
    string filename=cl.next(" ");;
    while ((filename!=" ") && (filename != "-o")) {
      infiles.push_back(filename);
      filename=cl.next(" ");
    }
  } else if (cl.search("-f")) {
    string filename="test";
    igzstream ifs; ifs.open(cl.follow("list","--filelist"));
    if(!ifs.good() || !ifs) {
      ERR << "Cannot open filelist " <<cl.follow("list","--filelist")  << ". Aborting." << endl;
      exit(20);
    }
    while(!ifs.eof() && filename!="") {
      getline(ifs,filename);
      if(filename!="") {
        infiles.push_back(filename);
      }
    }
    ifs.close();
  } else {
    USAGE();
  }

  if (cl.search("-o")) {
    outputJoergfile = cl.next(" ");
  } else {
    USAGE();
  }

  DBG(10) << "Merging content of " << infiles.size() << " joergfiles into " << outputJoergfile << "." << endl;


  vector<DoubleVectorVector> imageData;
  IntVector classData;

  for (uint i = 0; i < infiles.size(); i++) {
    DBG(10) << "load joergfile " << infiles[i] << endl;
    pair<DoubleVectorVector, IntVector> image = readJF(infiles[i]);
    DBG(20) << "with " << (image.first).size() << " images and " << (image.first)[0]->size() << " features each." << endl;
    if (imageData.empty()) {
      imageData.push_back(image.first);
      classData = image.second;
    } else {
      if (imageData[0].size() != image.first.size()) {
	ERR << "joergfiles to be merged must have the same number of images !" << endl;
	exit(1);
      }
      if (!classDataEquals(classData, image.second)) {
	ERR << "joergfiles to be merged must have the same class data " << endl;
	exit(2);
      }
      imageData.push_back(image.first);
    }
  }

  DoubleVectorVector newImageData;

  // collect the content
  for (uint i = 0; i < imageData[0].size(); i++) {
    vector<DoubleVector*> content;
    for (uint j = 0; j < imageData.size(); j++) {
      content.push_back(imageData[j][i]);
    }
    newImageData.push_back(collect(content));
  }

  DBG(10) << "write jörgfile " << outputJoergfile << endl;
  DBG(20) << "with " << ((int) newImageData.size()) << " images and " << ((int) newImageData[0]->size()) << " features each." << endl;
  writeJF(10, newImageData, classData, outputJoergfile);

}
