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
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include "getpot.hpp"
#include <iostream>

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "sobelcreatejf [options] -i <input joergfile> -o <output joergfile>" << endl
       << "   Options:" << endl
       << "    -h, --help          show this help" << endl
       << "    -s, --selection     comma-separated list of sobel filters to apply." << endl
       << "                        possible values are:" << endl
       << "                        horizontal, vertical, diagonal, laplace" << endl
       << "                        (default is horizonal,vertical)" << endl
       << endl;
  exit(20);
}

void collectVectors(const vector<DoubleVector*>& vecs, DoubleVector* dv) {
  dv->clear();
  for (uint i = 0; i < vecs.size(); i++) {
    dv->insert(dv->end(), (vecs[i])->begin(), (vecs[i])->end());
  }
}

// main function. reads a JF image, applies a filter and saves the data
// back into another JF.

int main(int argc, char**argv) {

  GetPot cl(argc, argv);
  string inputJoergfile;
  string outputJoergfile;
  bool withHorizontal = true;
  bool withVertical = true;
  bool withDiagonal = false;
  bool withLaplace = false;

  if (cl.search("-h")) {
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

  if (cl.search(2, "-s", "--selection")) {
    withHorizontal = false;
    withVertical = false;
    withDiagonal = false;
    withLaplace = false;
    string param = cl.next(" ");
    char cParam[100];
    strcpy(cParam, param.c_str());
    char* token = strtok(cParam, ",");
    while (token != NULL) {
      if (strcmp(token, "vertical") == 0) {
	withVertical = true;
	DBG(20) << "will apply vertical sobel filter" << endl;
      } else if (strcmp(token, "horizontal") == 0) {
	withHorizontal = true;
	DBG(20) << "will apply horizontal sobel filter" << endl;
      } else if (strcmp(token, "diagonal") == 0) {
	withDiagonal = true;
	DBG(20) << "will apply both diagonal sobel filters" << endl;
      } else if (strcmp(token, "laplace") == 0) {
	withLaplace = true;
	DBG(20) << "will apply laplace filter" << endl;
      } else {
	USAGE();
      }
      token = strtok(NULL, ",");
    }
  }

  DoubleVectorVector imageData;
  IntVector classData;
  DBG(10) << "load jörgfile" << inputJoergfile << "\n";

  pair<DoubleVectorVector, IntVector> image = readJF(inputJoergfile);
  imageData = image.first;
  classData = image.second;

  DBG(20) << "jörgfile has " << ((int) imageData.size()) << " images with " << ((int) imageData[0]->size()) << " features each." << endl;

  DoubleVectorVector newImageData;

  // collect the content
  for (uint i = 0; i < imageData.size(); i++) {
    vector<DoubleVector*> content;
    ImageFeature img;
    if (withHorizontal) {
      img.createFromJF(imageData[i]);
      sobelh(img);
      square(img);
      normalize(img);
      content.push_back(img.toJF());
    }
    if (withVertical) {
      img.createFromJF(imageData[i]);
      sobelv(img);
      square(img);
      normalize(img);
      content.push_back(img.toJF());
    }
    if (withDiagonal) {
      img.createFromJF(imageData[i]);
      sobeldiagonal1(img);
      square(img);
      normalize(img);
      content.push_back(img.toJF());
      img.createFromJF(imageData[i]);
      sobeldiagonal2(img);
      square(img);
      normalize(img);
      content.push_back(img.toJF());
    }
    if (withLaplace) {
      img.createFromJF(imageData[i]);
      laplace(img);
      square(img);
      normalize(img);
      content.push_back(img.toJF());
    }

    DoubleVector* dv = new DoubleVector();
    collectVectors(content, dv);
    content.clear();
    newImageData.push_back(dv);

  }

  DBG(10) << "write jörgfile " << outputJoergfile << endl;
  DBG(20) << "with " << ((int) newImageData.size()) << " images and " << ((int) newImageData[0]->size()) << " features each." << endl;
  writeJF(10, newImageData, classData, outputJoergfile);

}

