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
#include "diag.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include <vector>
#include "stdlib.h"
#include "getpot.hpp"

#include <dirent.h>
#include <float.h>

using namespace std;

void solveAssignmentProblemDoubleRect(double **Array, int **Result, int m, int n) {
  // IMPORTANT! The values of Array[size1][size2] are changed in this routine
  //
  // adopted from KNUTH 1994 (The Stanford GraphBase, pp. 74ff)
  // <numbers> refer to the paragraphs in Knuth' text 

  // local variables <14>
  int i;
  int k; // current row of interest
  int l; // current column of interest
  int j; // another interesting column
  int*col_mate;  // column matching a row or -1
  int*row_mate;  // row matching a column or -1
  int*parent_row;  // parent in the forest / ancestor of column's mate or -1
  int*unchosen_row;  // node in the forest
  int*slack_row;  // row where minimum slack[] occurs
  int t; // total number of nodes in the forest
  int q; // total number of explored nodes in the forest
  int unmatched; 
  double s;  // current matrix element of interest
  double*row_dec;  // subtracted from row (\sigma_k)
  double*col_inc;  // added to column (\tau_l)
  double*slack;  // minimum uncovered entry in column
  double del;
  double*Array_k;
  
  if(m>n){
    fprintf(stderr,"cannot have m>n in %s, line %d",__FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  // allocate intermediate data structures <15>
  col_mate=(int*)calloc(sizeof(int),m);
  row_mate=(int*)calloc(sizeof(int),n);
  parent_row=(int*)calloc(sizeof(int),n);
  unchosen_row=(int*)calloc(sizeof(int),m);
  slack_row=(int*)calloc(sizeof(int),n);

  row_dec=(double*)calloc(sizeof(double),m);
  col_inc=(double*)calloc(sizeof(double),n);
  slack=(double*)calloc(sizeof(double),n);
      
  // initialize result matrix to zero                                                                                      
  for(i=0;i<m;++i) {
    for(j=0;j<n;++j) {
      Result[i][j]=0; 
    } 
  }
  
  if(m==n){ // otherwise heuristic does not work
    // subtract column minimum from each column <12>
    for(l=0;l<n;++l){
      s=Array[0][l];
      for(k=1;k<n;++k) {
	if(Array[k][l]<s) {s=Array[k][l];}
      }
      for(k=0;k<n;++k) {
	Array[k][l]-=s;
      }
    }
  }
  // end <12>

  // initialize t, row_mate, parent_row, col_inc, slack <16>
  t=0; // forest starts empty
  for(l=0;l<n;++l){
    row_mate[l]=-1;
    parent_row[l]=-1;
    col_inc[l]=0.0;
    slack[l]=DBL_MAX;
  }

  // choose starting assignment <16>
  for(k=0;k<m;++k){
    s=Array[k][0];
    Array_k=Array[k];
    for(l=1;l<n;++l)
      if(Array_k[l]<s) s=Array_k[l];
    row_dec[k]=s;  // row_dec[k] = row minimum
    for(l=0;l<n;++l) // at minimum set row_mate and col_mate if column has no mate yet
      if((s==Array_k[l])&&(row_mate[l]<0)){
        col_mate[k]=l;  
        row_mate[l]=k;
        goto row_done;
      }
    col_mate[k]=-1; // if column already has a mate, row is unchosen
    unchosen_row[t++]=k;
  row_done:;
  }

  // Hungarian algorithm <18>
  // at most m stages with O(mn) operations -> O(m^2 n) runtime
  if(t==0) goto done;
  unmatched=t;
  while(1){
    q=0;
    while(1){
      while(q<t) {
	// explore node q of the forest <19>
	k=unchosen_row[q]; // k iterates over unchosen rows, indexed by q
	s=row_dec[k];      
	for(l=0;l<n;++l)
	  if(slack[l]>0.0){
	    del=Array[k][l]-s+col_inc[l]; // this is the "real" array value (-dec+inc)
	    // ??? if (del<0.0) del=0.0;
	    if(del<slack[l]){
	      if(del<=0.0){ // we found a new zero at [k][l]  //del==0 -> DBL_EPS??
		// changed == to <= since it can only be smaller than zero due to numerical "problems"
		if(row_mate[l]<0) goto breakthru; // matching can be increased
		slack[l]=0.0; // this column will now be chosen
		parent_row[l]=k;
		unchosen_row[t++]=row_mate[l];
	      }
	      else{
		slack[l]=del;
		slack_row[l]=k;
	      }
	    }
	  } // end <19>
        q++;
      }
      // introduce new zero into the matrix by modifying row_dec and col_inc <21>
      // we have explored the entire forest; none of the unchosen rows
      // has led to a breakthrough
      // an unchosen column with smallest slack will allow further progress
      s=DBL_MAX;
      for(l=0;l<n;++l)
        if(slack[l]>0.0 && slack[l]<s)
          s=slack[l];  // find minimum non-zero slack
      for(q=0;q<t;++q)
        row_dec[unchosen_row[q]]+=s; // and decrease all unchosen rows
      for(l=0;l<n;++l) {
        if(slack[l]>0.0){ // column l is not chosen   
          slack[l]-=s;
          if(slack[l]<=0.0) // slack[l]==0 -> DBL_EPS??  
	    // changed == to <= since it can only be smaller than zero due to numerical "problems"
	    // look at new zero <22>
            {
              k=slack_row[l];
              if(row_mate[l]<0){
                for(j=l+1;j<n;++j)
                  if(slack[j]<=0.0) col_inc[j]+=s; // slack[l]==0 -> DBL_EPS??
                goto breakthru; // matching can be increased 
              }
	      else { // not a breakthrough, but the forest continues to grow
                parent_row[l]=k;
                unchosen_row[t++]=row_mate[l];
              }
            }
        }
	else {
	  col_inc[l]+=s;
	}
      }
      // end <21>
    }
  breakthru:  
    // update matching by pairing row k with column l <20>
    while(1){
      j=col_mate[k];
      col_mate[k]=l;
      row_mate[l]=k;
      if(j<0)break;
      k=parent_row[j];
      l=j;
    }
    // end <20>
    if(--unmatched==0) goto done;
    // prepare for new stage by reinitializing the forest <17>
    t=0;
    for(l=0;l<n;l++){
      parent_row[l]=-1;
      slack[l]=DBL_MAX;
    }
    for(k=0;k<m;++k)
      if(col_mate[k]<0){
        unchosen_row[t++]=k;
      }
    // end <17>
  }
 done:
  // end <18>

  for (i=0;i<m;++i){
    Result[i][col_mate[i]]=1;
  }

  free(col_mate);
  free(row_mate);
  free(parent_row);
  free(unchosen_row);
  free(row_dec);
  free(col_inc);
  free(slack);
  free(slack_row);
}


void USAGE() {
  cout << "USAGE: collage [options] --input image-file --dir images-directory" << endl
       << "  available options: " << endl
       << "     --winx   how many windows in x direction (2)" << endl
       << "     --winy   how many windows in y direction (2)" << endl
       << "     --res    resolution of windows in x direction (128)" << endl
    //    << "     --sc-inp  where is the scaled input image is saved (input-scaled.png)" << endl
       << "     --gammas  number of gamma values used (3)" << endl
       << "     --dupl    number of times each image can be used (1)" << endl
       << "     --output  where is the output saved (output.png)" << endl
       << endl;
}

double eucDistance(ImageFeature &A, ImageFeature &B){
  uint dimX=A.xsize();
  uint dimY=A.ysize();
  uint dimZ=A.zsize();

  if(dimX!=B.xsize() ||
     dimY!=B.ysize() ||
     dimZ!=B.zsize()){
    cerr << "error in eucDistance: image sizes do not match " << endl;
    exit(20);
  }
  double sum=0.0,tmp;
  for (uint z=0; z<dimZ; ++z)
    for (uint y=0; y<dimY; ++y)
      for (uint x=0; x<dimX; ++x){
	tmp=A(x,y,z)-B(x,y,z);
	sum+=tmp*tmp;
      }
  return sum;
}

ImageFeature extractcenterpart(ImageFeature &im, uint resX, uint resY){
  uint top,left,bottom,right;
//  double scale;
  double aspect=((double) resY)/((double) resX);
  uint currentxsize=im.xsize();
  uint currentysize=im.ysize();
  double currentAspect=((double)currentysize) / ((double)currentxsize);
  if (currentAspect>aspect){
    // take middle row
    left=0;
    right=currentxsize;
    uint height=(uint)(aspect*((double)currentxsize));
    top=(currentysize-height)/2;
    bottom=top+height;
//    scale=((double)resY)/((double)height);
  }
  else {
    // take middle column
    top=0;
    bottom=currentysize;
    uint width=(uint)(((double)currentysize)/aspect);
    left=(currentxsize-width)/2;
    right=left+width;
//    scale=((double)resX)/((double)width);
  }
  ImageFeature patch=getPatch(im,left,top,right,bottom);
//  AffineTransformation aff(scale,0,0,scale,0,0);
//  patch=affineTransformation(patch,aff);
//  patch=getPatch(patch,0,0,resX,resY);
  patch=scale(patch,resX+1,resY+1);
  return patch;
}


int main(int argc, char**argv) {
  double gammas[]={1.3,0.8,1.6,0.6,2.0,0.5,2.5,0.4,3.0,0.3,5.0,0.2};
  const uint maxGammas=12;

  GetPot cl(argc,argv);

  // check nonoptional parameters
  if (!(cl.search("--input")&&cl.search("--dir"))) {
    USAGE();
    exit(20);
  }

  // read some optional parameters
  uint winX=cl.follow(2,"--winx");
  uint winY=cl.follow(2,"--winy");
  uint totalsize=(winX*winY);
  uint numgammas=cl.follow(3,"--gammas")-1; // one is the image itself with gamma = 1.0
  if(numgammas>maxGammas){
    numgammas=maxGammas;
    cout << " maximum number of gammas is " << maxGammas+1 << ", setting parameter to that value" << endl;
  }
  uint dupl=cl.follow(1,"--dupl");

  // read input image
  ImageFeature origimage;
  cl.search("--input");
  string filename=cl.next("default");
  DBG(10) << "loading " << filename << endl;
  origimage.load(filename);
  //should now catch: not an image
  uint dimX=origimage.xsize();
  uint dimY=origimage.ysize();
  uint dimZ=origimage.zsize();

  // set resolution of patches
  uint resX=cl.follow(128,"--res");
  uint resY=(resX*dimY*winX)/(dimX*winY);
  DBG(20) << "resX = " << resX << "; resY = " << resY << endl;

  // read all images from directory + extract center part of correct aspect
  cl.search("--dir");
  string dirname=cl.next("default");
  DBG(10) << "processing directory " << dirname << endl;
  vector<ImageFeature> patches;
  DIR* dirp;
  struct dirent* dirt;
  if(!(dirp=opendir(dirname.c_str()))) {
    cerr << "failed to open directory " << dirname << endl;
    exit(20);
  }
  while((dirt=readdir(dirp)) != NULL) {
    if(dirt->d_name[0]!='.'){
      string filename = dirname;
      filename += "/" + (string) (dirt->d_name);
      DBG(10) << "loading + processing " << filename << endl;
      ImageFeature tmp;
      tmp.load(filename);
      // catch: not an image
      tmp=extractcenterpart(tmp,resX,resY);
      patches.push_back(tmp);
      if(numgammas>0){
	ImageFeature tmp2=tmp;
	for(uint actgamma=0;actgamma<numgammas;++actgamma){
	  gammaCorrection(tmp2,gammas[actgamma]);
	  patches.push_back(tmp2);
	  tmp2=tmp;
	}
      }
    }
  }
  closedir(dirp);
  uint numImages=patches.size();

  // only for Hungarian
  if(totalsize>dupl*numImages){
    cerr << "there are not enough images in directory " << dirname << endl;
    //cout << totalsize << " " << winX << " " << winY << " " << numImages << "  " << dupl << " " << numgammas << endl; 
    exit(20);    
  }

  // allocate matrices
  double **distanceMatrix=(double**)calloc(sizeof(double*), totalsize);
  int **assignment=(int**)calloc(sizeof(int*),totalsize);
  for(uint i=0;i<totalsize;++i) {
    distanceMatrix[i]=(double*)calloc(sizeof(double),numImages);
    assignment[i]=(int*)calloc(sizeof(int),numImages);
  }

  // calculate distances between image patches
  // allow brightness adjustment?
  DBG(10) "calculating patch distances" << endl;
  uint newX=resX*winX, newY=resY*winY;
  double scale=((double)(newX))/((double)dimX);
  DBG(20) "scale is " << scale << endl;
  ImageFeature largerOrig(newX,newY,dimZ);
  AffineTransformation aff(scale,0,0,scale,0,0);
  if(scale>=1.0){ // here "largerimage" is correct; the input is enlarged
    setPatch(largerOrig,0,0,origimage);
    // don't care about half pixels for now, superresolution can be added later
    largerOrig=affineTransformation(largerOrig,aff);
    //largerOrig.save(cl.follow("input-scaled.png","--sc-inp"));
  }
  else{ // "largerimage" actually is a smaller image
    largerOrig=affineTransformation(origimage,aff);
  }
  ImageFeature origPatch;
  for(uint winx=0;winx<winX;++winx) {
    for(uint winy=0;winy<winY;++winy) {
      BLINK(15) << ".";
      origPatch=getPatch(largerOrig,winx*resX,winy*resY,winx*resX+resX,winy*resY+resY);
      // assume color images (this should be determined automatically)
      // use Euclidean distance (this should better be L*a*b distance or similar)
      for(uint i=0;i<numImages;++i){
        distanceMatrix[winy*winX+winx][i]=eucDistance(origPatch,patches[i]);
      }
    }
    BLINK(15) << endl;
  }
 
  // determine the best assignment  
  DBG(10) "determining best assignment" << endl;

//  // identity
//  for(uint i=0; i<totalsize; i++){
//    assignment[i][i]=1;
//  }

//  // local minimum
//  for(uint i=0; i<totalsize; i++){
//    uint minind=0;
//    double min=distanceMatrix[i][0];
//    for(uint j=1; j<numImages; ++j){
//      if(distanceMatrix[i][j]<min) {
//	minind=j;
//	min=distanceMatrix[i][j];
//      }
//    }
//    assignment[i][minind]=1;
//  }

  // Hungarian
  if(dupl>1){
    for(uint i=0;i<totalsize;++i) {
      distanceMatrix[i]=(double*)realloc(distanceMatrix[i],sizeof(double)*dupl*numImages);
      assignment[i]=(int*)realloc(assignment[i],sizeof(int)*dupl*numImages);
      for(uint j=1;j<dupl;++j){
	memcpy(&(distanceMatrix[i][j*numImages]),&(distanceMatrix[i][0]),sizeof(double)*numImages);
	memcpy(&(assignment[i][j*numImages]),&(assignment[i][0]),sizeof(int)*numImages); // only zeros here
      }
    }
  }
  solveAssignmentProblemDoubleRect(distanceMatrix, assignment, totalsize, dupl*numImages);

  // determine and write output image
  
  ImageFeature outputimage(newX, newY, dimZ);

  DBG(10) "assembling the new image" << endl;

  uint sourceImage=0;
  for(uint winx=0;winx<winX;++winx) {
    for(uint winy=0;winy<winY;++winy) {
      
      int *currentRow=assignment[winy*winX+winx];
      for(uint i=0;i<dupl*numImages; i++){
	if (currentRow[i]==1) {
	  sourceImage=i%numImages;
	  break;
	}
      }
      // instead of using the pre-extracted patches, we could use the patch 
      // of highest possible resolution here (extract again from image)
      setPatch(outputimage, winx*resX, winy*resY, patches[sourceImage]);
    }
  }

  outputimage.save(cl.follow("output.png","--output"));

  // now we should free some memory...
}



