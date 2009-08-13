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
#include "localfeatureextractor.hpp"
#include <argtable2.h>

using namespace std;

/** print usage information */
void usage(void **argtable, const char *progname) {
    FILE *fp = stdout;

    fprintf(fp, "Usage: %s ", progname);
    arg_print_syntaxv(fp, argtable, "\n");
    arg_print_glossary(fp, argtable, " %-50s %s\n");
    exit(0);
}


static string ID="$Id: $";

int main(int argc, char **argv) {

  const char *progname="extractlf";
  
  struct arg_int *arg_wavelet     =arg_int0("w","waveletPoints","<wavelet points>","number of wavelet salient points to be extracted");
  struct arg_int *arg_dog         =arg_int0("d","dogPoints","<dog points>","number of difference-of-Gaussian salient points to be extracted");
  struct arg_int *arg_random      =arg_int0("r","randomPoints","<random points>","number of random points to be extracted");
  struct arg_int *arg_gridx       =arg_int0("x","gridx","<grid resolution x>","resolution of a regular grid in x direction");
  struct arg_int *arg_gridy       =arg_int0("y","gridy","<grid resolution y>","resolution of a regular grid in y direction");
  struct arg_lit *arg_gridalign   =arg_lit0("A","alignGrid","calculate extraction size from x,y parameters of x and scale to winsize");
  struct arg_lit *arg_aligned     =arg_lit0("a","alignedGrid","use aligned grid with the defined extraction sizes");
  
  struct arg_int *arg_exSize      =arg_intn("s","extractionSizes","<extraction size>",0,100,"specify an extraction size, can occur up to ten times (applies not to dog)");
  
  struct arg_lit *arg_gray        =arg_lit0("g","gray","make image gray before processing");
  struct arg_lit *arg_padding     =arg_lit0(NULL,"padding","pad images after point determination and before patch extraction");
  
  struct arg_int *arg_winsize      =arg_int0("W","winsize","<winsize>","scale extracted patches to this size");
  
  struct arg_lit *arg_sift        =arg_lit0("S","sift","use sift features instead of patches");
  
  struct arg_lit *arg_pca         =arg_lit0("p","pca","use pca");
  struct arg_int *arg_pcadim      =arg_int0("P","pcadim","<dim>","dimensionality after pca");
  struct arg_file *arg_pcafile    =arg_file0(NULL,"pcafile","<pcafile>","file to read PCA from or to write PCA to");
  struct arg_lit *arg_pcaload     =arg_lit0(NULL,"loadpca","load pca, if not specified, pca is estimated and written to the specified file");
                                            
  struct arg_file *arg_filelist   =arg_file1(NULL,NULL,"<filelist>","database file of images to be processed");
  struct arg_lit *help            =arg_lit0("h", "help", "show help and exit");
  struct arg_end *end             =arg_end(20);
  
  void* argtable[]={arg_wavelet, arg_dog, arg_random, arg_gridx, arg_gridy, arg_gridalign, arg_aligned, 
                    arg_exSize, arg_gray, arg_padding, arg_winsize, arg_sift, 
                    arg_pca, arg_pcadim, arg_pcafile, arg_pcaload, arg_filelist,help,end};
  

  if(arg_nullcheck(argtable)!=0) { ERR << "Insufficient memory for cmdline parsing!" << endl; exit(20);}
  int nerrors=arg_parse(argc,argv,argtable);
  
  if(help->count > 0) {
    usage(argtable, progname);
  }
  if(nerrors > 0) {
    arg_print_errors(stderr, end, progname);
    fprintf(stderr, "Try '%s --help' for more information.\n",progname);
    exit(EXIT_FAILURE);
  }
  
  LocalFeatureExtractorSettings settings;
  
  settings.waveletPoints=arg_wavelet->ival[0];
  settings.dogPoints=arg_dog->ival[0];
  settings.randomPoints=arg_random->ival[0];
  settings.gridX=arg_gridx->ival[0];
  settings.gridY=arg_gridy->ival[0];
  settings.alignedGrid=(arg_aligned->count!=0);
  settings.alignGrid=(arg_gridalign->count!=0);
  
  settings.forceGray=(arg_gray->count!=0);
  settings.padding=(arg_padding->count!=0);
  settings.patches=(arg_sift->count==0);
  settings.sift=(arg_sift->count!=0);

  settings.extractionSizes.resize(arg_exSize->count);
  for(int i=0;i<arg_exSize->count;++i) { settings.extractionSizes[i]=arg_exSize->ival[i]; }
  settings.winsize=arg_winsize->ival[0];
  
  settings.pca=(arg_pca->count!=0);
  settings.pcadim=arg_pcadim->ival[0];
  settings.loadpca=(arg_pcaload->count!=0);
  settings.pcafile=string(arg_pcafile->filename[0]);
  
  string filelist=arg_filelist->filename[0];
  
  LocalFeatureExtractor lfe;
  lfe.settings()=settings;
  
  Database db;
  DBG(10) << "Reading database with files to process." << endl;
  db.loadFileList(filelist);
  
  DBG(10) << "Starting to extract features..." << endl;
  lfe.extractFromDatabase(db);
}
