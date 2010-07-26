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
#include <string>
#include <vector>
#include "localfeatures.hpp"
#include "getpot.hpp"
#include "gzstream.hpp"
#include "database.hpp"
#ifdef HAVE_KDTREE_LIBRARY
#include "knn_api.h"
#else
#warning "KDTREE LIBRARY NOT AVAILABLE. This program will not make any sense."
#endif


using namespace std;

void USAGE() {
  cout << "USAGE: " << endl
       << "   createlfkdtree (--filelist <filelist>|--firefilelist <firefilelist>) --kdtree <treefile>" << endl
       << "   creates a kdtree from the local feature files specified in filelist" << endl
       << "   Available options:" << endl
       << "     -h, --help       show this help and exit" << endl
       << endl;
}

int main(int argc, char **argv) {
  GetPot cl(argc, argv);
  char *ctmp;  
  string kdtreefilename=cl.follow("tree.kdt","--kdtree");
  string filelistname;
  string filename;
  vector<string> filenames;
  
  if(cl.search("--filelist")) {
    filelistname=cl.next("-");
    while(filelistname != "-") {
      igzstream ifs; ifs.open(filelistname.c_str());
      if(!ifs.good() || !ifs) {
        ERR << "Cannot open filelist " << filelistname << ". Aborting." << endl;
        exit(20);
      } else {
        DBG(10) << "Reading filelist: " << filelistname << endl;
        filename="asdf";
        while(!ifs.eof() && filename!="") {
          getline(ifs,filename);      
          if(filename!="") {
            filenames.push_back(filename);
          }
        }
        ifs.close();
      }
      filelistname=cl.next("-");
    }
    DBG(10) << "Filelist reading done... going to process " << filenames.size() << " files in the following." << endl;
  } else if (cl.search("--firefilelist")){
    filelistname=cl.follow("files.txt","--firefilelist");
    Database db;
    DBG(10) << "Reading filelist from FIRE filelist" << endl;
    db.loadFileList(filelistname);
    
    
    for(uint n=0;n<db.size();++n) {
      filenames.push_back(db.path()+"/"+db.filename(n)+"."+db.suffix(0));
    }
    DBG(10) << "Filelist reading done... going to process " << filenames.size() << " files in the following." << endl;
  }

  if(cl.search(2,"--help","-h")) {USAGE(); exit(20);}
  
  if(!cl.search("--kdtree") || (!cl.search("--filelist") and !cl.search("--firefilelist"))) {USAGE(); exit(20);}
#ifdef HAVE_KDTREE_LIBRARY  
  LocalFeatures lf;
  t_fvec* vectors=new t_fvec;
  uint number=0;
  uint dimchanges=0;

  for(int i=0;i<filenames.size();++i) {
    filename=filenames[i];
    DBG(10) << "Consistency check for features from '" << filename << "'." << endl;
    lf=LocalFeatures();
    lf.load(filename);
    if(vectors->dim!=int(lf.dim())) {
      vectors->dim=lf.dim();
      DBG(10) << "Setting dim to " << vectors->dim << endl;
      ++dimchanges;
    }
    number+=lf.size();
  }
  if(dimchanges==1) {
    DBG(10) << "Loading " << number<< " local features of dimension " << vectors->dim << endl;
  } else {
    DBG(10) << "Probably not consistent: There have been local features of different dimensionalities: Exiting" << endl;
    exit(10);
  }

  
  vectors->featvec=new float*[number];
  vectors->labelvec=new char*[number];
  vectors->nvectors=number;
  int aktvec=0;
  DBG(10) << "Consistency check finished. Now loading features." << endl;

  for(int j=0;j<filenames.size();++j) {
    filename=filenames[j];
    DBG(10) << "Loading " << filename << " " << VAR(aktvec)<< endl;
    lf=LocalFeatures();
    lf.load(filename);
    for(uint i=0;i<lf.size();++i) {
      
      vectors->featvec[aktvec]=new float[vectors->dim];
      for(int j=0;j<vectors->dim;++j) {
        vectors->featvec[aktvec][j]=lf[i][j];
      }
      ctmp=new char[lf.filename().size()+1];
      strcpy(ctmp,lf.filename().c_str());
      
      vectors->labelvec[aktvec]=ctmp;
      DBG(50) << "Copied filename: " <<vectors->labelvec[aktvec]<< endl;
      ++aktvec;
    }
    DBG(30) << "aktvec=" << aktvec << endl;
  }
  
  t_knn_kdtree *kdt=(t_knn_kdtree *)malloc(1*sizeof(t_knn_kdtree)); 
  
  DBG(10) << "Creating kdtree from these data" << endl;
  knn_kdtree_create(vectors,10,kdt);
  DBG(10) << "Saving kdtree to " << kdtreefilename << endl;
  char *fn=new char[kdtreefilename.size()];
  strcpy(fn,kdtreefilename.c_str());
  knn_kdtree_save(kdt,10,0,fn);
  DBG(10) << "cmdline was: " ; printCmdline(argc,argv);
  #endif
}
