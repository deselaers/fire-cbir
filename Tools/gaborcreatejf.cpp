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

#include "gabor.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include "jflib.hpp"
#include "diag.hpp"
#include "getpot.hpp"
#include <iostream>

using namespace std;

void collectVectors(const vector<DoubleVector*>& vecs, DoubleVector* dv) {
  dv->clear();
  for (uint i = 0; i < vecs.size(); i++) {
    dv->insert(dv->end(), (vecs[i])->begin(), (vecs[i])->end());
  }
}

void USAGE() {
  cout << "USAGE:" << endl
       << "gaborcreatejf [options] -i <input joergfile> -o <output joergfile>" << endl
       << "   Options:" << endl
       << "    -h, --help          show this help" << endl
       << "    -p, --phases        number of phases for the gabor filter (default: 2)" << endl
       << "    -f, --frequencies   number of frequencies for the gabor filter (default: 2)" << endl
       << "    -s, --selection     comma-separated list of gabor layers to consider." << endl
       << "                        each one must be >= 0 and < phases * frequencies" << endl
       << endl;
  exit(20);
}

int main(int argc, char** argv)
{

  const int NUM_PHASES_DEFAULT = 2;
  const int NUM_FREQUENCIES_DEFAULT = 2;
  int numPhases = NUM_PHASES_DEFAULT;
  int numFrequencies = NUM_FREQUENCIES_DEFAULT;

  GetPot cl(argc,argv);
  string inputJoergfile;
  string outputJoergfile;
  bool hasSelection = false;
  vector<int> selection;

  if (cl.search(2, "-h", "--help")) {
    USAGE();
  }

  if (cl.search("-i")) {
    inputJoergfile = cl.next(" ");
  } else {
    USAGE();
  }
  if (cl.search("-o")) {
    outputJoergfile = cl.next(" ");
  } else {
    USAGE();
  }
  if (cl.search(2, "-p", "--phases")) {
    numPhases = cl.next(NUM_PHASES_DEFAULT);
  }
  if (cl.search(2, "-f", "--frequencies")) {
    numFrequencies = cl.next(NUM_FREQUENCIES_DEFAULT);
  }
  if (cl.search(2, "-s", "--selection")) {
    hasSelection = true;
    string param = cl.next(" ");
    char cParam[100];
    strcpy(cParam, param.c_str());
    char* token = strtok(cParam, ",");
    while (token != NULL) {
      selection.push_back(atoi(token));
      token = strtok(NULL, ",");
    }
  }

  DBG(20) << "Loading " << inputJoergfile << "..." <<  endl;
  
  pair<DoubleVectorVector, IntVector> jfs = readJF(inputJoergfile);
  DBG(10) << "Loaded " << jfs.first.size() << " images." << endl;
  DoubleVectorVector gabors(0);
  DoubleVector singleJf;
  ImageFeature img;
  ImageFeature gImg;
  Gabor *gabor;

  DBG(10) << "Calculating gabor features..." << endl;
  for (uint i = 0; i < jfs.first.size(); i++) {

    vector<DoubleVector*> fullData;

    img.createFromJF(jfs.first[i]);
    // display for debugging purposes
    //img.display(); 
    normalize(img);
    //img.display(); 
    gabor = new Gabor(img);
    gabor->calculate(numPhases, numFrequencies);
    for (int l = 0; l < numPhases * numFrequencies; l++) {
      if (!hasSelection) {
	gImg = gabor->getImage(l);
	power(gImg, 2);
	normalize(gImg);
	//gImg.display();
	fullData.push_back(gImg.toJF());
      } else {
	for (vector<int>::iterator it = selection.begin(); it != selection.end(); it++) {
	  if ((*it) == l) {
	    gImg = gabor->getImage(l);
	    power(gImg, 2);
	    normalize(gImg);
	    //gImg.display();
	    fullData.push_back(gImg.toJF());
	  }
	}
      }
    }
    delete gabor;
    DoubleVector *dv = new DoubleVector();
    collectVectors(fullData, dv);
    fullData.clear();
    gabors.push_back(dv);


  }

  DBG(10) << "Writing " << outputJoergfile << endl;
  writeJF(10, gabors, jfs.second, outputJoergfile);

}
