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

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include "gzstream.hpp"
#include "localfeatures.hpp"
#include "diag.hpp"
#include "database.hpp"
#include "gaussiandensity.hpp"
#include "getpot.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE: lfcreatejf [--normalize] --train <filelist> [ --clusters <filelist> ] --test <filelist> --jftrain <joergfile name> -jftest <joergfile name>" << endl;
  exit(1);
}

// returns the number of classes found in the given database
// the method does not look for distinct class labels, but
// returns the highest class label + 1
// therefore, classes are assumed to be labeled starting from
// 0 to n-1 
int getNumClasses(const Database& db) {
  int maxClassNr = -1;
  for (int n = 0; n < (int) db.size(); n++) {
    int clazz = db[n]->clas();
    if (clazz < 0) {
      ERR << "found image with class " << clazz << " in database!" << endl;
      exit(1);
    } else {
      if (clazz > maxClassNr) {
      	maxClassNr = clazz;
      }
    }
  }
  return maxClassNr + 1;
}

// loads model for each class k
void loadModel(const string filename, vector<GaussianDensity>& densities) {
  igzstream is; is.open(filename.c_str());
  if (!is || !is.good()) {
    ERR << "Unable to read file '" << filename << "'." << endl;
  } else {
    string line;
    getline(is, line);
    if (!is.good()) {
      ERR << "Error reading file '" << filename << "'." << endl;
    } else {
      if (line != "# EM/LBG clustering model file") {
        ERR << "This is probably not an EM-Model file" << endl
            << "Continuing anyway" << endl;
      }
      while (!is.eof()) {
        getline(is,line);
        if (!is.eof()) {
          istringstream iss(line);
          string keyword;
          iss >> keyword;
          if (keyword == "gaussians" ) {
            uint size;
            iss >> size;
            densities.resize(size);
          } else if (keyword == "density") {
            uint no;
            iss >> no;
            iss >> keyword;
            if (keyword == "elements") {
              iss >> densities[no].elements;
            } else if (keyword == "dim") {
              iss >> densities[no].dim;
            } else if (keyword == "mean") {
              GaussianDensity& c = densities[no];
              c.mean = vector<double>(c.dim);
              for (uint i = 0; i < c.dim; ++i) {
                iss >> c.mean[i];
              }
            } else if (keyword == "variance") {
              GaussianDensity& c = densities[no];
              c.sigma = vector<double>(c.dim);
              for (uint i = 0; i < c.dim; ++i) {
                iss >> c.sigma[i];
              }
            } else {
              ERR << "Reading density received unknown keyword in position 3: '" << keyword << "'." << endl;
            }
	  } else if ((keyword == "splitmode") || (keyword == "maxsplits") || (keyword == "stopWithNClusters") || (keyword == "disturbMode") ||
		     (keyword == "poolMode") || (keyword == "dontSplitBelow") || (keyword == "iterationsBetweenSplits") || (keyword == "epsilon") ||
		     (keyword == "minObservationsPerCluster") || (keyword == "distance")) {
	    // ignore these keywords
          } else {
            ERR << "Unknown keyword '" << keyword << "'." << endl;
          }
        }
      }
    }
    is.close();
  }
  DBG(15) << "read " << densities.size() << " clusters" << endl;
}

void calcVariances(const Database& db, vector<double>& variances) {
  vector<double> tmpvariances;
  DBG(10) << "calculating variances" << endl;
  int numPatches = 0;
  for (int n = 0; n < (int) db.size(); n++) {
    for (int s = 0; s < (int) db.numberOfSuffices(); s++) {
      string filename = db.filename(n) + "." + db.suffix(s);
      LocalFeatures lf;
      lf.load(filename);
      if ((n == 0) && (s == 0)) {
	variances.resize(lf.dim(), 0.0);
	tmpvariances.resize(lf.dim(), 0.0);
      }
      for (int l = 0; l < (int) lf.numberOfFeatures(); l++) {
	vector<double> localFeature = lf[l];
	for (int d = 0; d < (int) localFeature.size(); d++) {
	  variances[d] += localFeature[d] * localFeature[d];
	  tmpvariances[d] += localFeature[d];
	}
	numPatches++;
      }
    }
  }
  for (int d = 0; d < (int) variances.size(); d++) {
    variances[d] /= (double) numPatches;
    tmpvariances[d] /= (double) numPatches;
    variances[d] -= tmpvariances[d] * tmpvariances[d];
    DBG(20) << "variance for dimension " << d << ": " << variances[d] << endl;
  }
}


int main(int argc, char** argv) {

  GetPot cl(argc, argv);

  // process command line arguments
  if (cl.search(2, "-h", "--help")) {
    USAGE();
  }

  string testFilelist = "", trainFilelist = "", joergfileTrain = "", joergfileTest = "",
    clusterFilelist = "";
  if (!cl.search(1, "--test")) {
    USAGE();
  } else {
    testFilelist = cl.next("");
  }
  if (!cl.search(1, "--train")) {
    if (!cl.search(1, "--clusters") || cl.search("--normalize")) {
      USAGE();
    }
  } else {
    trainFilelist = cl.next("");
  }
  bool clusters;
  if (!cl.search(1, "--clusters")) {
    clusters = false;
  } else {
    clusterFilelist = cl.next("");
    clusters = true;
  }
  if (!cl.search(1, "--jftrain")) {
    USAGE();
  } else {
    joergfileTrain = cl.next("");
  }
  if (!cl.search(1, "--jftest")) {
    USAGE();
  } else {
    joergfileTest = cl.next("");
  }

  DBG(10) << "loading testing filelist " << testFilelist << endl;
  Database testDb;
  testDb.loadFileList(testFilelist);
  DBG(10) << "loading training filelist " << trainFilelist << endl;
  Database trainDb;
  trainDb.loadFileList(trainFilelist);
  Database clusterDb;
  if (clusters) {
    DBG(10) << "loading cluster filelist " << trainFilelist << endl;
    clusterDb.loadFileList(clusterFilelist);
  }
  DBG(10) << "writing to training joergfile '" << joergfileTrain << "' and testing joergfile '" << joergfileTest << "'" << endl;
  
  uint numDifferentClasses = getNumClasses(trainDb);
  ogzstream ofsTrain; ofsTrain.open(joergfileTrain.c_str());
  ogzstream ofsTest; ofsTest.open(joergfileTest.c_str());

  vector<double> variances;
  bool normalize = cl.search("--normalize");
  if (normalize) {
    calcVariances(trainDb, variances);
  }

  if (clusters) {
    int lfDimension = -1;
    int numClusters = 0;
    for (int i = 0; i < (int) clusterDb.size(); i++) {
      for (int s = 0; s < (int) clusterDb.numberOfSuffices(); s++) {
	string filename = clusterDb.filename(i) + "." + clusterDb.suffix(s);
	DBG(15) << "loadind clusters from file " << filename << endl;
	vector<GaussianDensity> densities;
	loadModel(filename, densities);
	if ((i == 0) && (s == 0)) {
	  lfDimension = densities[0].dim;
	  DBG(15) << "clusters have " << lfDimension << " dimensions" << endl;
	  ofsTrain << numDifferentClasses << " " << lfDimension << endl;
	} else {
	  if (lfDimension != (int) densities[0].dim) {
	    ERR << "found cluster with " << densities[0].dim << " dimensions, do not match previous ones with " << lfDimension << " dimensions" << endl;
	    ofsTrain.close();
	    ofsTest.close();
	    exit(1);
	  }
	}
	for (int c = 0; c < (int) densities.size(); c++) {
	  ofsTrain << clusterDb[i]->clas();
	  for (int d = 0; d < (int) densities[c].mean.size(); d++) {
	    if (normalize) {
	      ofsTrain << " " << (densities[c].mean[d] / variances[d]);
	    } else {
	      ofsTrain << " " << densities[c].mean[d];
	    }
	  }
	  ofsTrain << endl;
	  numClusters++;
	}
      }
    }
    ofsTrain << "-1" << endl;
    ofsTrain.close();
    DBG(15) << "finished writing to " << joergfileTrain << endl;
    DBG(15) << "wrote " << numClusters << " clusters in total" << endl;

  } else {

    DBG(15) << "loaded " << trainDb.size() << " local feature files with " << trainDb.numberOfSuffices() << " suffices each" << endl;
    DBG(15) << "having features from " << numDifferentClasses << " different classes" << endl;
    int lfDimension = -1;
    int numPatches = 0;
    for (int n = 0; n < (int) trainDb.size(); n++) {
      for (int s = 0; s < (int) trainDb.numberOfSuffices(); s++) {
	string filename = trainDb.filename(n) + "." + trainDb.suffix(s);
	DBG(15) << "processing " << filename << endl;
	LocalFeatures lf;
	lf.load(filename);
	DBG(15) << "loaded " << lf.numberOfFeatures() << " local features" << endl;
	if ((n == 0) && (s == 0)) {
	  lfDimension = lf.dim();
	  DBG(15) << "local features have " << lfDimension << " dimensions" << endl;
	  ofsTrain << numDifferentClasses << " " << lfDimension << endl;
	} else {
	  if (lfDimension != (int) lf.dim()) {
	    ERR << "found local features with " << lf.dim() << " dimensions, do not match previous ones with " << lfDimension << " dimensions" << endl;
	    ofsTest.close();
	    ofsTrain.close();
	    exit(1);
	  }
	}
	for (int l = 0; l < (int) lf.numberOfFeatures(); l++) {
	  ofsTrain << trainDb[n]->clas();
	  vector<double> localFeature = lf[l];
	  for (int d = 0; d < (int) localFeature.size(); d++) {
	    if (normalize) {
	      ofsTrain << " " << (localFeature[d] / variances[d]);
	    } else {
	      ofsTrain << " " << localFeature[d];
	    }
	  }
	  ofsTrain << endl;
	  numPatches++;
	}
      }
    }
    ofsTrain << "-1" << endl;
    ofsTrain.close();
    DBG(15) << "finished writing to " << joergfileTrain << endl;
    DBG(15) << "wrote " << numPatches << " patches in total" << endl;
  }

  // write testing joergfile
  int lfDimension = -1;
  int numPatches = 0;
  for (int n = 0; n < (int) testDb.size(); n++) {
    for (int s = 0; s < (int) testDb.numberOfSuffices(); s++) {
      string filename = testDb.filename(n) + "." + testDb.suffix(s);
      DBG(15) << "processing " << filename << endl;
      LocalFeatures lf;
      lf.load(filename);
      DBG(15) << "loaded " << lf.numberOfFeatures() << " local features" << endl;
      if ((n == 0) && (s == 0)) {
	lfDimension = lf.dim();
	DBG(15) << "local features have " << lfDimension << " dimensions" << endl;
	ofsTest << numDifferentClasses << " " << lfDimension << endl;
      } else {
	if (lfDimension != (int) lf.dim()) {
	  ERR << "found local features with " << lf.dim() << " dimensions, do not match previous ones with " << lfDimension << " dimensions" << endl;
	  ofsTest.close();
	  exit(1);
	}
      }
      for (int l = 0; l < (int) lf.numberOfFeatures(); l++) {
	ofsTest << testDb[n]->clas();
	vector<double> localFeature = lf[l];
	for (int d = 0; d < (int) localFeature.size(); d++) {
	  if (normalize) {
	    ofsTest << " " << (localFeature[d] / variances[d]);
	  } else {
	    ofsTest << " " << localFeature[d];
	  }
	}
	ofsTest << endl;
	numPatches++;
      }
    }
  }
  ofsTest << "-1" << endl;
  ofsTest.close();
  DBG(15) << "finished writing to " << joergfileTest << endl;
  DBG(15) << "wrote " << numPatches << " patches in total" << endl;

  DBG(10) << "cmdline was: "; printCmdline(argc,argv);

}
