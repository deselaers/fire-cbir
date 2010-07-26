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
 * This program is necessary if e.g. experiments with  maximum entropy
 * have to be made with features from FIRE.
 */

#include <iostream>
#include <string>
#include "vectorfeature.hpp"
#include "getpot.hpp"
#include "database.hpp"
#include "jflib.hpp"
#include "sparsehistogramfeature.hpp"
#include "localfeatures.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE for db2jf -- convert any FIRE database to jf" << endl << "    db2jf <options> --filelist <filelist> --outfile <filename>" << endl << "    options:" << endl << "       -h                : show this help" << endl << "       -smoothe <factor> : when converting sparse histograms, smoothe them using" << endl << "                           the given factor [0.0 - 1.0]  before conversion" << endl << "                          this option is ignored, if no sparse histograms are converted" << endl << "        -lf to switch to local feature mode: this only allows one suffix! it might work, but will do weird things with more suffices" << endl << endl;
}

int main(int argc, char **argv) {
  GetPot cl(argc, argv);

  if (cl.search("-h") || !cl.search("--filelist") || !cl.search("--outfile")) {
    USAGE();
    exit(20);
  }

  string filelistfilename=cl.follow("filelist", "--filelist");
  string outfilefilename=cl.follow("outfile.jf", "--outfile");

  Database db;
  db.loadFileList(filelistfilename);
  db.loadFeatures();

  DoubleVectorVector data;

  IntVector classData;
  DoubleVector *tmpFeat;

  double tmp;
  double smootheFactor = cl.follow(0.0, "-smoothe");
  if ((smootheFactor < 0.0) || (smootheFactor > 1.0)) {
    ERR << "invalid smoothe factor given: " << smootheFactor << endl;
    ERR << "0.0 (no smoothing) <= smootheFactor <= 1.0 (total smoothing) required !" << endl;
  }

  if (!cl.search("-lf")) {
    for (uint i=0; i<db.size(); ++i) {
      tmpFeat=new DoubleVector();
      for (uint j=0; j<db.numberOfSuffices(); ++j) {
        const VectorFeature* dbfeat=dynamic_cast<const VectorFeature *>(db[i]->operator[](j));
        if (dbfeat) {
          for (uint k=0; k<dbfeat->size(); ++k) {
            tmp= dbfeat->operator[](k);
            if (isnan(tmp)) {
              ERR << "nan in " << db.filename(i) << "." << endl;
              exit(20);
            }
            tmpFeat->push_back(tmp);
          }
        } else {
          const SparseHistogramFeature* dbfeat=dynamic_cast<const SparseHistogramFeature*>(db[i]->operator[](j));
          if (dbfeat) {
            DoubleVector t;
            if (smootheFactor != 0.0) {
              dbfeat->getAllSmoothedBinValues(t, smootheFactor);
            } else {
              dbfeat->getAllBinValues(t);
            }

            for (uint i=0; i<t.size(); ++i) {
              tmpFeat->push_back(t[i]);
            }

          } else {
            ERR << "Cannot convert this feature type. Maybe try -lf" << endl;
            exit(EXIT_FAILURE);
          }
        }
      }
      data.push_back(tmpFeat);
      classData.push_back(db[i]->clas());
    }
  } else {
    DBG(10) << "in lf mode" << endl;
    for (uint i=0; i<db.size(); ++i) {
      tmpFeat=new DoubleVector();
      for (uint j=0; j<db.numberOfSuffices(); ++j) {
        const LocalFeatures* dbfeat=dynamic_cast<const LocalFeatures*>(db[i]->operator[](j));
        if (dbfeat) {
          for (uint f=0; f<dbfeat->size(); ++f) {
            DoubleVector * feat=new DoubleVector(dbfeat->getData()[f]);
            data.push_back(feat);
            classData.push_back(db[i]->clas());
          }
        }
      }
    }
  }

  DBG(10)<< VAR(data.size()) << endl;
  int maxCls=*(max_element(classData.begin(), classData.end()))+1;
  writeJF(maxCls, data, classData, outfilefilename);
}

