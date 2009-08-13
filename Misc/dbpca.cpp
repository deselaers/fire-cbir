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
 * A program to apply PCA transformation to any image retrieval
 * database in FIRE format. Only vectorfeatures
 * (e.g. vectors, images, histograms) and local features can be used
 */


#include <iostream>
#include <string>
#include "vectorfeature.hpp"
#include "getpot.hpp"
#include "database.hpp"
#include "jflib.hpp"
#include "sparsehistogramfeature.hpp"
#include "localfeatures.hpp"
#include "pca.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE for dbpca -- apply pca on any FIRE database" << endl
       << "    db2jf <options> --filelist <filelist>" << endl
       << "    options:" << endl
       << "       -p <dim>     : specify output dimensionality" << endl
       << "       -h           : show this help" << endl
       << "       -s           : save pcas to PCA-<filelistfilename>-<suffix>.pca.gz" << endl
       << "       -l <pcafile> : load pcas from pcafile" << endl
       << endl;
}

int main(int argc, char **argv) {
  GetPot cl(argc,argv);
  
  if(cl.search("-h") || !cl.search("--filelist")) {USAGE(); exit(20);}
  
  string filelistfilename=cl.follow("filelist","--filelist");
  int PCAdim=cl.follow(40,"-p");

  if(cl.search("-s") && cl.search("-l"))
  {
    cout << "ERROR: Loading the PCA AND saving it again does not make sense!" << endl;
    exit(0);
  }

  bool loadpca = cl.search("-l");
  bool savepca = cl.search("-s");

  string pcafile;
  if(loadpca)
  {
    pcafile = cl.follow("pcafile", "-l");
  }

  Database db;
  FeatureLoader fl;
  ImageContainer img;
  db.loadFileList(filelistfilename);
  
  DoubleVectorVector data;

  IntVector classData;
  
  int dim=0;
  
  for(uint j=0;j<db.numberOfSuffices();++j) {
    
    db[0]->operator[](j) = fl.load(db[0]->basename(),db.suffix(j),db.relevantSuffix(j),db.path());
    
    const VectorFeature* dbfeat=dynamic_cast<const VectorFeature *>(db[0]->operator[](j));
    if(dbfeat) {
      dim=dbfeat->size();
    } else {
      const LocalFeatures* dbfeat=dynamic_cast<const LocalFeatures*>(db[0]->operator[](j));
      if(dbfeat) {
        dim=dbfeat->getData()[0].size();
      }
    }
    
    PCA pca(dim);

    if(loadpca)
    {
      DBG(10) << "Loading PCA from " << pcafile << endl;
      pca.load(pcafile);
    }
    else
    {
      DBG(10) << "Collecting data for  PCA estimation" << endl;
      for(uint i=0;i<db.size();++i) {
        
        // Load Feature
        img = ImageContainer(db.filename(i), db.numberOfSuffices());
        db.loadQuery(db.filename(i), &img);
        
        // Calc PCA
        DBG(10) << "Accumulating for suffix " << j+1 <<"/" << db.numberOfSuffices()  << " image " <<i+1 <<"/" << db.size() << endl;
        const VectorFeature* dbfeat=dynamic_cast<const VectorFeature *>(img[j]);
        if(dbfeat) {
          pca.putData(dbfeat->data());
        } else {
          const LocalFeatures* dbfeat=dynamic_cast<const LocalFeatures*>(img[j]);
          if(dbfeat) {
            for(uint f=0;f<dbfeat->size();++f) {
              pca.putData(dbfeat->getData()[f]);
            }
          } else {
            ERR << "Feature type for suffix " << j << " not supported: " << db.suffix(j)<< endl;
          }
        }
        
        // Unload features
        img.clear();
      }
      
      pca.dataEnd();
      pca.calcPCA();
    }

   
    
    if(savepca)
    {
      cout << "Saving PCA to " << "PCA-"+filelistfilename+"-"+db.suffix(j)+".pca.gz" << endl;
      pca.save("PCA-"+filelistfilename+"-"+db.suffix(j)+".pca.gz");
    }


    for(uint i=0;i<db.size();++i) {
      
      // Load feature
      img = ImageContainer(db.filename(i),db.numberOfSuffices());
      db.loadQuery(db.filename(i), &img);
      
      // PCA Transform
      DBG(10) << "Transforming for suffix " << j+1 <<"/" << db.numberOfSuffices()  << " image " <<i+1 <<"/" << db.size() << endl;
      const VectorFeature* dbfeat=dynamic_cast<const VectorFeature *>(img[j]);
      if(dbfeat) {
        DoubleVector tr=pca.transform(dbfeat->data(),PCAdim);
        VectorFeature toSave(tr);
        toSave.save(db.path()+"/"+db.filename(i)+"."+db.suffix(j)+".pca.vec.gz");
      } else {
        const LocalFeatures* dbfeat=dynamic_cast<const LocalFeatures*>(img[j]);
        if(dbfeat) {
          LocalFeatures toSave=*dbfeat;
          toSave.getData().clear();
          toSave.dim()=PCAdim;
          for(uint f=0;f<dbfeat->size();++f) {
            DoubleVector tr=pca.transform(dbfeat->getData()[f],PCAdim);
            toSave.getData().push_back(tr);
          }
          toSave.save(db.path()+"/"+db.filename(i)+".pca."+db.suffix(j));
        } else {
          ERR << "Feature type for suffix " << j << " not supported: " << db.suffix(j)<< endl;
        }
      }
      
      // Unload features
      img.clear();
    }

    
  }     // for j
}

