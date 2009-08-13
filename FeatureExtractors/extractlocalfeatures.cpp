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
#include "imagefeature.hpp"
#include "getpot.hpp"
#include "imagelib.hpp"
#include "gzstream.hpp"
#include "pca.hpp"
#include "differenceofgaussian.hpp"
#include "salientpoints.hpp"
#include "pasannotation.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE:" << endl
       << "extractlocalfeatures [options] (--color|--gray) (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "    -h, --help   show this help" << endl
       << "    --suffix <suffix>     to set the suffix of the output file, default: lf.gz" << endl
       << "    --winsize <winsize>   to set the size of the extracted local features, default: 5" << endl
       << "    --extractionSize <sizes> to specify the sizes of extracted patches to be scaled to winsize later" << endl
       << "    --padding             take lf ranging out of the image, default: no" << endl
       << "    --normalize           to mean/variance normalize extracted features, default: no" << endl
       << "    --histonorm           to apply histogram normalization before local feature extraction" << endl
       << "    --minmaxnorm          to apply minimum/maximum normalization before local feature extraction" << endl
       << "    --tmpforlf <tmp>      to save localfeatures which are not necessary after pca transformation" << endl
       << "    --sobel n             filter values with n x n sobel matrices, n can be 3, 5 or 7" << endl
       << "    --annotation <anno>   use annotation which are stored in the given directory" << endl
       << "    --spatial             add spatial derivatives to patches" << endl
       << "  OPTIONS WHERE FEATURES ARE EXTRACTED" << endl
       << "   if more than one of these methods is specified, the extraction points are added" << endl
       << "    --varpoints <no>      to take the first no features instead of all above a certain variance threshold" << endl
       << "    --varthreshold <t>    to set the variance threshold" << endl
       << "    --salientpoints <no>  to use salient point feature extraction and get <no> features" << endl
       << "    --dogpoints <no>      to use <no> Difference-of-Gaussian interest points" << endl
       << "    --uniformx <x>        to specify a regular grid of feature extraction points" << endl
       << "    --uniformy <y>        to specify a regular grid of feature extraction points" << endl
       << "    --randompoins <no>    to select <no> points randomly distributed over the image" << endl
       << "    --alignedgrid         to select extraction points based on an aligned grid" << endl
       << "    --overlapping     to select extraction points based on an overlapping grid" << endl
       << "  OPTIONS REGARDING PCA" << endl
       << "    --pca <dim>           to apply pca dimensionality reduction right after extraction" << endl
       << "    --loadPCA <filename>  don't calculate pca on the local features created, but load the given one" << endl
       << "    --savePCA <filename>  save the calculated PCA to the given file" << endl
       << "  OPTIONS REGARDING VISUALIZATION" << endl
       << "    --markPoints          with this option, the complete image will be saved with extraction points marked" << endl
       << endl;
  exit(20);
}



bool operator<(const pair<double,Point> &l, const pair<double,Point> &r) { return l.first<r.first;}

void extractPatch(const ImageFeature& img, uint x, uint y, uint w, bool normalize, uint winsize, uint savesize,
                  vector<double>& patchVector) {
  vector<double> tmpVector;
  ImageFeature patch=getPatch(img,x,y,w);
  if(normalize) { 
    meanAndVarianceNormalization(patch,0.5,0.5); 
    cutoff(patch);
  }
  if (w != winsize) { 
    patch=scale(patch,savesize,savesize);
  }
  tmpVector = toVector(patch);
  DBGI(100,{ patch = scale(patch, 100, 100); patch.display(); });
  patchVector.insert(patchVector.end(), tmpVector.begin(), tmpVector.end());
}

LocalFeatures extractPatches(const ImageFeature &img, // the image where we want to extract the patches
                             uint winsize,            // the winsize of the saved patches (2*winsize+1)²
                             bool padding,            // use padding?
                             uint salientpoints,      // how many salientpoints should be used
                             uint dogpoints,          // how many Difference-of-Gaussian interest points should be used
                             uint randompoints,       // how many random points should be used
                             uint varpoints,          // how many points based on local variance should be used
                             uint uniformx,           // how many grid points in x-direction
                             uint uniformy,           // how many grid points in y-direction
                             bool alignedgrid,        // determine extraction points based on an aligned grid
                             bool overlapping,        //     aligned grid is overlapping such that each patch's quarter is covered twice
                             double varthreshold,     // what is the variance threshold to be used for point finding
                             bool markPoints,         // save image with extraction positions marked?
                             string filename,         // the filename of the input image
                             vector<uint> extractionSizes, // different winsizes to be used for extraction
                             // if empty: use winsize only
                             bool mvnormalize,          // normalize mean and variance
                             bool forceGray,          // are we working on gray images?
                             bool useSobel,           // append sobel filtered image data
                             bool spatial,            // if spatial derivative grid points are used
                             const ImageFeature& sobelFilterHorizontal,
                             const ImageFeature& sobelFilterVertical,
                             bool hasAnnotation,      // is there an annotation of this image?
                             const PascalAnnotationFeature& annotation // optinal annotation of the image
                             ) {
  /*------------------------------------------------------------
    init
    ------------------------------------------------------------*/
  
  bool useScaling;
  uint paddingOffset=0; 
  ImageFeature padded;
  ImageFeature varimg;
  Point p;

  LocalFeatures lf;
  lf.winsize()=winsize;
  lf.padding()=padding;
  lf.filename()=filename;
  lf.varthreshold()=varthreshold;
  lf.zsize()=img.zsize();
  lf.imageSizeX()=img.xsize();
  lf.imageSizeY()=img.ysize();


  if(extractionSizes.size()==0) {extractionSizes.push_back(winsize); useScaling=false;}
  else {sort(extractionSizes.begin(), extractionSizes.end()); useScaling=true;}
  
  if(varthreshold>=0.0 || varpoints>0) {
    varimg=localvariance(img,winsize);
  }
  
  if(padding) {
    DBG(20) << "Padding" << endl;
    paddingOffset=extractionSizes[extractionSizes.size()-1]+2;
    padded=zeropad(img,paddingOffset,paddingOffset);
  } else {
    DBG(20) << "Not padding" << endl;
    padded=img;
  }

  // stores a list of elements, where each element is a pair
  // containing a point to be extracted and a list of patch sizes for extraction
  vector< pair<Point, vector<uint> > > extractionPoints;

  /*------------------------------------------------------------
    Get extraction points
    ------------------------------------------------------------*/

  if(salientpoints>0) {
    DBG(20) << "Extracting SalientPoints" << endl;

    SalientPoints sp;
    if(forceGray) {
      sp=SalientPoints(img.layer(0));
    } else {
      if ((img.zsize() == 3) || (img.zsize() == 1)) {
        sp=SalientPoints(img);
      } else {
        ERR << "Was expecting image with 1 (gray) or 3 (color) layers, but got one with " << img.zsize() << " layers." << endl;
        exit(10);
      }
    }
    
    vector<Point> sPoints = sp.getSalientPoints(salientpoints);
    int i = 0;
    int selected = 0;
    while (true) {
      if ((i == (int) sPoints.size()) || (selected == (int) salientpoints)) {
        break;
      }
      if (hasAnnotation) {
        for (int a = 0; a < (int) annotation.nOfObjects(); a++) {
          PASObject po = annotation.object(a);
          if  ((sPoints[i].x >= (int) po.Xmin_ - 1) && (sPoints[i].x <= (int) po.Xmax_ + 1) &&
               (sPoints[i].y >= (int) po.Ymin_ - 1) && (sPoints[i].y <= (int) po.Ymax_ + 1)) {
            extractionPoints.push_back(pair<Point, vector<uint> >(sPoints[i], extractionSizes));
            selected++;
            break;
          }
        }
      } else {
        extractionPoints.push_back(pair<Point, vector<uint> >(sPoints[i], extractionSizes));
        selected++;
      }
      i++;
    }
    DBG(15) << "After SalientPoints: Having " << extractionPoints.size() << " feature extraction points." << endl;
  }
  
  if(varthreshold>=0.0) {
    DBG(20) << "Applying VarianceThreshold" << endl;
    for(uint x=0;x<varimg.xsize();++x) {
      for(uint y=0;y<varimg.ysize();++y) {
        if(varimg(x,y,0) > varthreshold) {
          p=Point(x,y);
          extractionPoints.push_back(pair<Point, vector<uint> >(p, extractionSizes));
        }
      }
    }
    DBG(15) << "After VarianceThreshold: Having " << extractionPoints.size() << " feature extraction points." << endl;
  }

  if(varpoints>0) {
    DBG(20) << "Extracting VariancePoints" << endl;
    vector< pair<double,Point> > variances;
    for(uint x=0;x<varimg.xsize();++x) {
      for(uint y=0;y<varimg.ysize();++y) {
        variances.push_back(pair<double,Point>(varimg(x,y,0),Point(x,y)));
      }
    }
    sort(variances.begin(), variances.end());
    for(uint i=0;i<varpoints;++i) {
      extractionPoints.push_back(pair<Point, vector<uint> >(variances[i].second, extractionSizes));
    }
    DBG(15) << "After VariancePoints: Having " << extractionPoints.size() << " feature extraction points." << endl;
  }
  
  if(uniformx >0 and uniformy>0) {
    DBG(20) << "Applying UniformGrid" << endl;
    uint stepx=max(uint(1),img.xsize()/(uniformx+1));
    uint stepy=max(uint(1),img.ysize()/(uniformy+1));
    if(uniformx>0 || uniformy>0) {
      for(uint x=1;x<=uniformx;++x) {
        for(uint y=1;y<=uniformy;++y) {
          p=Point(x*stepx,y*stepy);
          extractionPoints.push_back(pair<Point, vector<uint> >(p, extractionSizes));
        }
      }
    }
    DBG(15) << "After UniformGrid: Having " << extractionPoints.size() << " feature extraction points." << endl;
  }

  if (dogpoints > 0) {
    DBG(20) << "Applying Difference-of-Gaussian" << endl;
    DifferenceOfGaussian sift(img);
    vector<InterestPoint> interestPoints = sift.getInterestPoints(dogpoints);
    for (vector<InterestPoint>::const_iterator ipIterator = interestPoints.begin(); ipIterator != interestPoints.end(); ipIterator++) {
      Point p; p.x = ipIterator->x; p.y = ipIterator->y;
      if (ipIterator->scale >= 5) {
        extractionPoints.push_back(pair<Point, vector<uint> >(p, vector<uint>(1, ipIterator->scale / 2)));
      }
    }
    DBG(15) << "After Difference-of-Gaussian: Having " << extractionPoints.size() << " feature extraction points." << endl;
  }

  if (randompoints > 0) {
    srand(0); // always use the same random seed to achieve deterministic results
    int p = 0;
    while (p < (int) randompoints) {
      int x = rand() % img.xsize();
      int y = rand() % img.ysize();
      if (hasAnnotation) {
        for (int a = 0; a < (int) annotation.nOfObjects(); a++) {
          PASObject po = annotation.object(a);
          if  ((x >= (int) po.Xmin_ - 1) && (x <= (int) po.Xmax_ + 1) &&
               (y >= (int) po.Ymin_ - 1) && (y <= (int) po.Ymax_ + 1)) {
            extractionPoints.push_back(pair<Point, vector<uint> >(Point(x, y), extractionSizes));
            p++;
            break;
          }
        }
      } else {
        extractionPoints.push_back(pair<Point, vector<uint> >(Point(x, y), extractionSizes));
        p++;
      }
    }
    DBG(15) << "After Random: Having " << extractionPoints.size() << " feature extraction points." << endl;
  }

  if (alignedgrid) {
    for (uint sizes = 0; sizes < extractionSizes.size(); sizes++) {
      uint size = extractionSizes[sizes];
      uint x = size;
      uint y = size;
      while (y + size < img.ysize()) {
        while (x + size < img.xsize()) {
          extractionPoints.push_back(pair<Point, vector<uint> >(Point(x, y), vector<uint>(1, size)));
          if(overlapping and (x+2*size)<img.xsize() and (y+2*size)< img.ysize()) {
            extractionPoints.push_back(pair<Point, vector<uint> >(Point(x+size, y+size), vector<uint>(1, size)));
          }
          x += 2 * size + 1;
        }
        x = size;
        y += 2 * size + 1;
      }
    }
    DBG(15) << "After aligned grid: Having " << extractionPoints.size() << " feature extraction points." << endl;
  }



  /*------------------------------------------------------------
    Extraction
    ------------------------------------------------------------*/
  DBG(15) << "Starting to extract from " << extractionPoints.size() << " feature extraction points." << endl;

  uint savesize=winsize*2+1;
  ImageFeature patch;
  ImageFeature padded_sobelv;
  ImageFeature padded_sobelh;
  if (useSobel) {
    padded_sobelv = padded;
    padded_sobelh = padded;
    convolve(padded_sobelv, sobelFilterVertical);
    convolve(padded_sobelh, sobelFilterHorizontal);
    normalize(padded_sobelh);
    normalize(padded_sobelv);
  }

  for(uint i=0;i<extractionPoints.size();++i) {
    pair<Point, vector<uint> > extractionElement = extractionPoints[i];
    p=extractionElement.first;
    DBG(50) << "Extracting image patch " << i << " at (" << p.x << ", " << p.y << ") in size";
    for(uint j=0;j<extractionElement.second.size();++j) {
      int w=extractionElement.second[j];
      BLINK(50) << " " << w;

      if(padding) {
        vector<double> patchVector;
        extractPatch(padded, p.x + paddingOffset, p.y + paddingOffset, w, mvnormalize, winsize, savesize, patchVector);
        if (useSobel) {
          extractPatch(padded_sobelh, p.x + paddingOffset, p.y + paddingOffset, w, mvnormalize, winsize, savesize, patchVector);
          extractPatch(padded_sobelv, p.x + paddingOffset, p.y + paddingOffset, w, mvnormalize, winsize, savesize, patchVector);
        }
        lf.addLocalFeature(patchVector,p,w);
      } else {
        int largeMargin = (spatial? 3*w+1: w);
        int normalMargin = w;
        if ((p.x - largeMargin >= 0) and (p.y - largeMargin >=0) and (p.x + normalMargin < int(padded.xsize())) and 
            (p.y + normalMargin < int(padded.ysize()))) {
          vector<double> patchVector;
          extractPatch(padded, p.x, p.y, w, mvnormalize, winsize, savesize, patchVector);
          if (useSobel) {
            extractPatch(padded_sobelh, p.x, p.y, w, mvnormalize, winsize, savesize, patchVector);
            extractPatch(padded_sobelv, p.x, p.y, w, mvnormalize, winsize, savesize, patchVector);
          }
          lf.addLocalFeature(patchVector,p,w);
          if (spatial) {
            patchVector.clear();
            extractPatch(padded, p.x - 2 * w - 1, p.y, w, mvnormalize, winsize, savesize, patchVector);
            lf.addLocalFeature(patchVector,pair<uint,uint>(p.x - 2 * w - 1,p.y),w);
            patchVector.clear();
            extractPatch(padded, p.x, p.y - 2 * w - 1, w, mvnormalize, winsize, savesize, patchVector);
            lf.addLocalFeature(patchVector,pair<uint,uint>(p.x,p.y - 2 * w - 1),w);
          }
        }
      }
    }
    BLINK(50) << endl;
  }

  /*------------------------------------------------------------
    Marking of extraction points
    ------------------------------------------------------------*/
  if(markPoints) {
    vector<double> color(3,1.0);
    for(uint i=0;i<extractionPoints.size();++i) {
      pair<Point, vector<uint> > extractionElement = extractionPoints[i];
      p=extractionElement.first;
      // cross(padded,p.x+paddingOffset,p.y+paddingOffset,color,3);
      for(uint j=0;j<extractionElement.second.size();++j) {
        uint w=extractionElement.second[j];
        box(padded,p.x+paddingOffset, p.y+paddingOffset,color,w);
      }
    }
    // if there are annotations, visualize them, too
    if (hasAnnotation) {
      for (int a = 0; a < (int) annotation.nOfObjects(); a++) {
        PASObject po = annotation.object(a);
        rect(padded, po.Xmin_, po.Ymin_, po.Xmax_ - po.Xmin_, po.Ymax_ - po.Ymin_, color, true);
      }
    }
    DBG(10) << "Saving positionmarkimage to " << filename+".featureextractionpoints.png" << endl;
    padded.save(filename+".featureextractionpoints.png");
  }

  /*------------------------------------------------------------
    End
    ------------------------------------------------------------*/

  lf.numberOfFeatures()=lf.size();
  if (lf.size() > 0) {
    lf.dim()=lf[0].size();
  }
  DBG(10) << lf.size() << " patches extracted." << endl;
  return lf;
}
                             
int main(int argc, char** argv) {
  GetPot cl(argc,argv);
  if(cl.search(2,"--help","-h")) USAGE();
  if(cl.size()<2) USAGE();

  string TMPdir;

  bool useSobel = cl.search("--sobel");
  ImageFeature sobelFilterHorizontal;
  ImageFeature sobelFilterVertical;
  if (useSobel) {
    int sobelMatrixSize = cl.next(0);
    if ((sobelMatrixSize != 3) && (sobelMatrixSize != 5) && (sobelMatrixSize != 7)) {
      ERR << "sobel filter matrix must be either 3x3, 5x5 or 7x7" << endl;
      exit(1);
    } else {
      sobelFilterHorizontal = ImageFeature(sobelMatrixSize, sobelMatrixSize, 1);
      sobelFilterVertical = ImageFeature(sobelMatrixSize, sobelMatrixSize, 1);
      if (sobelMatrixSize == 3) {
      	sobelFilterHorizontal(0, 0, 0) = -1; sobelFilterHorizontal(1, 0, 0) = 0; sobelFilterHorizontal(2, 0, 0) = 1;
        sobelFilterHorizontal(0, 1, 0) = -2; sobelFilterHorizontal(1, 1, 0) = 0; sobelFilterHorizontal(2, 1, 0) = 2;
        sobelFilterHorizontal(0, 2, 0) = -1; sobelFilterHorizontal(1, 2, 0) = 0; sobelFilterHorizontal(2, 2, 0) = 1;
      } else if (sobelMatrixSize == 5) {
      	sobelFilterHorizontal(0, 0, 0) = -1; sobelFilterHorizontal(1, 0, 0) = -1; sobelFilterHorizontal(2, 0, 0) = 0; sobelFilterHorizontal(3, 0, 0) = 1; sobelFilterHorizontal(4, 0, 0) = 1;
        sobelFilterHorizontal(0, 1, 0) = -1; sobelFilterHorizontal(1, 1, 0) = -2; sobelFilterHorizontal(2, 1, 0) = 0; sobelFilterHorizontal(3, 1, 0) = 2; sobelFilterHorizontal(4, 1, 0) = 1;
        sobelFilterHorizontal(0, 2, 0) = -1; sobelFilterHorizontal(1, 2, 0) = -2; sobelFilterHorizontal(2, 2, 0) = 0; sobelFilterHorizontal(3, 2, 0) = 2; sobelFilterHorizontal(4, 2, 0) = 1;
        sobelFilterHorizontal(0, 3, 0) = -1; sobelFilterHorizontal(1, 3, 0) = -2; sobelFilterHorizontal(2, 3, 0) = 0; sobelFilterHorizontal(3, 3, 0) = 2; sobelFilterHorizontal(4, 3, 0) = 1;
        sobelFilterHorizontal(0, 4, 0) = -1; sobelFilterHorizontal(1, 4, 0) = -1; sobelFilterHorizontal(2, 4, 0) = 0; sobelFilterHorizontal(3, 4, 0) = 1; sobelFilterHorizontal(4, 4, 0) = 1;
      } else if (sobelMatrixSize == 7) {
      	sobelFilterHorizontal(0, 0, 0) = -1; sobelFilterHorizontal(1, 0, 0) = -1; sobelFilterHorizontal(2, 0, 0) = -1; sobelFilterHorizontal(3, 0, 0) = 0; sobelFilterHorizontal(4, 0, 0) = 1; sobelFilterHorizontal(5, 0, 0) = 1; sobelFilterHorizontal(6, 0, 0) = 1;
        sobelFilterHorizontal(0, 1, 0) = -1; sobelFilterHorizontal(1, 1, 0) = -2; sobelFilterHorizontal(2, 1, 0) = -2; sobelFilterHorizontal(3, 1, 0) = 0; sobelFilterHorizontal(4, 1, 0) = 2; sobelFilterHorizontal(5, 1, 0) = 2; sobelFilterHorizontal(6, 1, 0) = 1;
        sobelFilterHorizontal(0, 2, 0) = -1; sobelFilterHorizontal(1, 2, 0) = -2; sobelFilterHorizontal(2, 2, 0) = -3; sobelFilterHorizontal(3, 2, 0) = 0; sobelFilterHorizontal(4, 2, 0) = 3; sobelFilterHorizontal(5, 2, 0) = 2; sobelFilterHorizontal(6, 2, 0) = 1;
        sobelFilterHorizontal(0, 3, 0) = -1; sobelFilterHorizontal(1, 3, 0) = -2; sobelFilterHorizontal(2, 3, 0) = -3; sobelFilterHorizontal(3, 3, 0) = 0; sobelFilterHorizontal(4, 3, 0) = 3; sobelFilterHorizontal(5, 3, 0) = 2; sobelFilterHorizontal(6, 3, 0) = 1;
        sobelFilterHorizontal(0, 4, 0) = -1; sobelFilterHorizontal(1, 4, 0) = -2; sobelFilterHorizontal(2, 4, 0) = -3; sobelFilterHorizontal(3, 4, 0) = 0; sobelFilterHorizontal(4, 4, 0) = 3; sobelFilterHorizontal(5, 4, 0) = 2; sobelFilterHorizontal(6, 4, 0) = 1;
        sobelFilterHorizontal(0, 5, 0) = -1; sobelFilterHorizontal(1, 5, 0) = -2; sobelFilterHorizontal(2, 5, 0) = -2; sobelFilterHorizontal(3, 5, 0) = 0; sobelFilterHorizontal(4, 5, 0) = 2; sobelFilterHorizontal(5, 5, 0) = 2; sobelFilterHorizontal(6, 5, 0) = 1;
        sobelFilterHorizontal(0, 6, 0) = -1; sobelFilterHorizontal(1, 6, 0) = -1; sobelFilterHorizontal(2, 6, 0) = -1; sobelFilterHorizontal(3, 6, 0) = 0; sobelFilterHorizontal(4, 6, 0) = 1; sobelFilterHorizontal(5, 6, 0) = 1; sobelFilterHorizontal(6, 6, 0) = 1;
      }
      double componentSum = 0.0;
      for (int x = 0; x < sobelMatrixSize; x++) {
      	for (int y = 0; y < sobelMatrixSize; y++) {
          sobelFilterVertical(x, y, 0) = sobelFilterHorizontal(y, x, 0);
          componentSum += fabs(sobelFilterHorizontal(y, x, 0));
        }
      }
      for (int x = 0; x < sobelMatrixSize; x++) {
      	for (int y = 0; y < sobelMatrixSize; y++) {
          sobelFilterVertical(x, y, 0) = sobelFilterVertical(x, y, 0) * 2.0 / componentSum;
          sobelFilterHorizontal(x, y, 0) = sobelFilterHorizontal(x, y, 0) * 2.0 / componentSum;
        }
      }
      DBG(10) << "adding sobel filtered data, using " << sobelMatrixSize << "x" << sobelMatrixSize << " matrices" << endl;
    }
  }

  bool tmpforlf=cl.search("--tmpforlf");
  if(tmpforlf) TMPdir=cl.follow("/tmp/","--tmpforlf");


  string suffix=cl.follow("lf.gz","--suffix");
  string annotationDir = cl.follow("", "--annotation");
  if (annotationDir != "") {
    DBG(10) << "using annotation in directory " << annotationDir << endl;
  }

  uint winsize=cl.follow(5,"--winsize");
  bool padding=cl.search("--padding");
  bool normalizeFeatures=cl.search("--normalize");

  vector<uint> extractionSizes;
  if(cl.search("--extractionSize")) {
    int exS=cl.next(-1);
    while(exS!=-1) {
      extractionSizes.push_back(exS);
      exS=cl.next(-1);
    }
  }
  
  bool trans_pca=cl.search("--loadPCA");  
  bool calc_pca=cl.search("--pca")&&!trans_pca;
  bool markPoints=cl.search("--markPoints");
  
  uint zsize; bool forceGray=false;
  if(cl.search("--color")) { zsize=3; forceGray=false;}
  else if(cl.search("--gray")) { zsize=1; forceGray=true;}
  else { USAGE(); exit(20);}
  
  PCA pca((2*winsize+1)*(2*winsize+1)*zsize);
  PCA sobelvPCA((2*winsize+1)*(2*winsize+1)*zsize);
  PCA sobelhPCA((2*winsize+1)*(2*winsize+1)*zsize);
  
  // get the information where the features are extracted
  uint varpoints=cl.follow(0,"--varpoints");
  double varthreshold=cl.follow(-1.0,"--varthreshold");
  uint salientpoints=cl.follow(0,"--salientpoints");
  uint dogpoints=cl.follow(0, "--dogpoints");
  uint randompoints=cl.follow(0, "--randompoints");
  uint uniformx=cl.follow(0,"--uniformx");
  uint uniformy=cl.follow(0,"--uniformy");
  bool alignedgrid=cl.search("--alignedgrid");
  bool overlapping=cl.search("--overlapping");
  bool spatial = cl.search(1, "--spatial");
  bool histonorm = cl.search(1, "--histonorm");
  bool minmaxnorm = cl.search(1, "--minmaxnorm");

  if (useSobel && spatial) {
    ERR << "sobel filter and spatial derivate cannot be used in combination !" << endl;
    exit(1);
  }

  vector<string> images;
  if(cl.search("--images")) {
    string filename=cl.next(" ");
    while(filename!=" ") {
      images.push_back(filename);
      filename=cl.next(" ");
    }
  } else if (cl.search("--filelist")) {
    string filename="test";
    igzstream ifs; ifs.open(cl.follow("list","--filelist"));
    if(!ifs.good() || !ifs) {
      ERR << "Cannot open filelist " <<cl.follow("list","--filelist")  << ". Aborting." << endl;
      exit(20);
    }
    while(!ifs.eof() && filename!="") {
      getline(ifs,filename);
      if(filename!="") {
        images.push_back(filename);
      }
    }
    ifs.close();
    
  } else {
    USAGE();
    exit(20);
  }
  
  for(uint i=0;i<images.size();++i) {
    string filename=images[i];
    DBG(10) << "Processing '"<< filename << "'.(" << i << "/" << images.size()<< ")" <<endl;
    ImageFeature img; img.load(filename,forceGray);
    // apply histogram normalization, if desired
    if (histonorm) {
      histogramNormalization(img);
    }
    // apply mininum/maximum normalization, if desired
    if (minmaxnorm) {
      normalize(img);
    }
    
    PascalAnnotationFeature annotation;
    if (annotationDir != "") {
      string fname = filename;
      for (int i = 0; i < (int) fname.size(); i++) {
        if (fname[i] == '/') {
          fname[i] = '_';
        }
      }
      DBG(15) << "loading annotation from file " << annotationDir << "/" << fname << ".annotation" << endl;
      annotation.load(annotationDir + "/" + fname + ".annotation");
    }
    
    LocalFeatures lf;
    lf = extractPatches(img, winsize, padding, salientpoints, dogpoints, randompoints, varpoints, uniformx, uniformy, alignedgrid, overlapping, varthreshold, markPoints,
                        filename, extractionSizes, normalizeFeatures, forceGray, useSobel, spatial, sobelFilterHorizontal, sobelFilterVertical,
                        annotationDir != "", annotation);
    
    if(tmpforlf) {
      string fn=filename; 
      for(uint k=0;k<fn.size();++k) if (fn[k]=='/') fn[k]='_';
      lf.save(TMPdir+"/"+fn+"."+suffix);
    } else {
      lf.save(filename+"."+suffix);
    }
    DBG(10) << "Saved " << lf.numberOfVectors() << " local features to '" << filename << "." << suffix <<"'." << endl;
    
    if(calc_pca) {
      if (!useSobel) {
        for(uint i=0;i<lf.numberOfVectors();++i) {
          pca.putData(lf[i]);
        }
      } else {
        for(uint i = 0; i < lf.numberOfVectors(); ++i) {
          if ((lf.dim() % 3) != 0) {
            ERR << "dimension of local features should be dividable by 3, as it contains additional sobel data!" << endl;
            exit(1);
          }
          int partSize = lf.dim() / 3;
          vector<double> lfpart(partSize);
          for (int j = 0; j < partSize; j++) {
            lfpart[j] = lf[i][j];
          }
          pca.putData(lfpart);
          for (int j = 0; j < partSize; j++) {
            lfpart[j] = lf[i][j + partSize];
          }
          sobelhPCA.putData(lfpart);
          for (int j = 0; j < partSize; j++) {
            lfpart[j] = lf[i][j + 2 * partSize];
          }
          sobelvPCA.putData(lfpart);
        }
      }
    }
  }

  if(calc_pca||trans_pca) {
    if(calc_pca) {
      DBG(10) << "PCA received " << pca.counter() << " input vectors." << endl;
      pca.dataEnd();
      DBG(10) << "Starting to calculate PCA" << endl;
      pca.calcPCA();
      if(cl.search("--savePCA")) {
        pca.save(cl.follow("localfeatures.pca.gz","--savePCA"));
      }
      if (useSobel) {
        DBG(10) << "vertical sobel PCA received " << sobelvPCA.counter() << " input vectors." << endl;
        sobelvPCA.dataEnd();
        DBG(10) << "Starting to calculate vertical sobel PCA" << endl;
        sobelvPCA.calcPCA();
        if(cl.search("--savePCA")) {
          string infix = ".sobelv";
          string tmpfilename = cl.follow("localfeatures.pca.gz","--savePCA");
          if (tmpfilename.find(".") == tmpfilename.npos) {
            tmpfilename.append(infix);
          } else {
            tmpfilename.insert(tmpfilename.find("."), infix);
          }
          DBG(15) << "saving vertical sobel PCA to '" << tmpfilename << "'" << endl;
          sobelvPCA.save(tmpfilename);
        }
        DBG(10) << "horizontal sobel PCA received " << sobelhPCA.counter() << " input vectors." << endl;
        sobelhPCA.dataEnd();
        DBG(10) << "Starting to calculate horizontal sobel PCA" << endl;
        sobelhPCA.calcPCA();
        if(cl.search("--savePCA")) {
          string infix = ".sobelh";
          string tmpfilename = cl.follow("localfeatures.pca.gz","--savePCA");
          if (tmpfilename.find(".") == tmpfilename.npos) {
            tmpfilename.append(infix);
          } else {
            tmpfilename.insert(tmpfilename.find("."), infix);
          }
          DBG(15) << "saving horizontal sobel PCA to '" << tmpfilename << "'" << endl;
          sobelhPCA.save(tmpfilename);
        }
      }
    } else {
      pca.load(cl.follow("localfeatures.pca.gz","--loadPCA"));
      if (useSobel) {
        string infix = ".sobelv"; 
        string tmpfilename = cl.follow("localfeatures.pca.gz","--loadPCA");
        if (tmpfilename.find(".") == tmpfilename.npos) {
          tmpfilename.append(infix);
        } else {
          tmpfilename.insert(tmpfilename.find("."), infix);
        }
        DBG(15) << "loading vertical sobel PCA from '" << tmpfilename << "'" << endl;
        sobelvPCA.load(tmpfilename);
        infix = ".sobelh";
        tmpfilename = cl.follow("localfeatures.pca.gz","--loadPCA");
        if (tmpfilename.find(".") == tmpfilename.npos) {
          tmpfilename.append(infix);
        } else {
          tmpfilename.insert(tmpfilename.find("."), infix);
        }
        DBG(15) << "loading horizontal sobel PCA from '" << tmpfilename << "'" << endl;
        sobelhPCA.load(tmpfilename);
      }
    }
    
    
    for(uint i=0;i<images.size();++i) {

      if (!useSobel) {
        uint dim=cl.follow(40,"--pca");
        string filename=images[i];
        DBG(10) << "PCA transforming '" << filename  << "." << suffix <<"'." << endl;
        LocalFeatures lf1; 
        if(tmpforlf) {
          string fn=filename; for(uint k=0;k<fn.size();++k) if (fn[k]=='/') fn[k]='_';
          lf1.load(TMPdir+"/"+fn+"."+suffix);
        } else {
          lf1.load(filename+"."+suffix);
        }
        for(uint i=0;i<lf1.numberOfVectors();++i) {
          DBGI(30,{DBG(30) << "untrans:";for(uint j=0;j<lf1[i].size();++j) BLINK(30) << lf1[i][j] << " ";BLINK(30) << endl;});
          vector<double> transformed=vector<double>(pca.transform(lf1[i],dim));
          DBGI(30,{DBG(30) << "transformed:"; for(uint j=0;j<transformed.size();++j) BLINK(30) << transformed[j] << " "; BLINK(30) << endl;});
          lf1[i].resize(dim);
          for(uint j=0;j<dim;++j) { lf1[i][j]=transformed[j];}
        }
        if (spatial) {
          if ((lf1.numberOfFeatures() % 3) != 0) {
            ERR << "was expecting # of local features to be dividable by 3 !" << endl;
            exit(1);
          }
          uint numFeatures = lf1.numberOfFeatures() / 3;
          for (uint n = 0; n < numFeatures; n++) {
            if ((lf1.position(n).x <= lf1.position(n + 1).x) && (lf1.position(n).y != lf1.position(n + 1).y)) {
              ERR << "local features not aligned for spatial derivatives! (this should not happen)" << endl;
              exit(1);
            }
            if ((lf1.position(n).x != lf1.position(n + 2).x) && (lf1.position(n).y <= lf1.position(n + 2).y)) {
              ERR << "local features not aligned for spatial derivatives! (this should not happen)" << endl;
              exit(1);
            }
            lf1[n].resize(dim * 3);
            for (uint d = 0; d < dim; d++) {
              lf1[n][d + dim] = lf1[n][d] - lf1[n + 1][d];
              lf1[n][d + 2 * dim] = lf1[n][d] - lf1[n + 2][d];
            }
            lf1.removeLocalFeature(n + 2);
            lf1.removeLocalFeature(n + 1);
          }
          if (lf1.numberOfFeatures() != numFeatures) {
            ERR << "unexpected number of features, is " << lf1.numberOfFeatures() << ", should be " << numFeatures << endl;
            exit(1);
          }
          lf1.dim() = dim * 3;
        } else {
          lf1.dim() = dim;
        }
        lf1.save(filename+".pca."+suffix);
        filename=cl.next(" ");

      } else {
        uint dimTotal = 3 * cl.follow(40,"--pca");
        uint dimPixel = dimTotal / 3;
        string filename=images[i];
        DBG(10) << "PCA transforming '" << filename  << "." << suffix <<"'." << endl;
        LocalFeatures lfPixel; 
        if(tmpforlf) {
          string fn=filename; for(uint k=0;k<fn.size();++k) if (fn[k]=='/') fn[k]='_';
          lfPixel.load(TMPdir + "/" + fn + "." + suffix);
        } else {
          lfPixel.load(filename + "." + suffix);
        }
        uint partSize = lfPixel.dim() / 3;
        for(uint i = 0; i < lfPixel.numberOfVectors(); ++i) {
          vector<double> lfpart1(partSize), lfpart2(partSize), lfpart3(partSize);
          for (int j = 0; j < (int) partSize; j++) { 
            lfpart1[j] = lfPixel[i][j]; 
            lfpart2[j] = lfPixel[i][j + partSize];
            lfpart3[j] = lfPixel[i][j + 2 * partSize];
          }
          DBGI(30,{DBG(30) << "untrans:";for(uint j = 0; j < partSize; ++j) BLINK(30) << lfpart1[j] << " ";BLINK(30) << endl;});
          DBGI(30,{DBG(30) << "untrans:";for(uint j = 0; j < partSize; ++j) BLINK(30) << lfpart2[j] << " ";BLINK(30) << endl;});
          DBGI(30,{DBG(30) << "untrans:";for(uint j = 0; j < partSize; ++j) BLINK(30) << lfpart3[j] << " ";BLINK(30) << endl;});
          vector<double> transformed1 = vector<double>(pca.transform(lfpart1, dimPixel));
          vector<double> transformed2 = vector<double>(sobelhPCA.transform(lfpart2, dimPixel));
          vector<double> transformed3 = vector<double>(sobelvPCA.transform(lfpart3, dimPixel));
          DBGI(30,{DBG(30) << "transformed:"; for(uint j=0;j<transformed1.size();++j) BLINK(30) << transformed1[j] << " "; BLINK(30) << endl;});
          DBGI(30,{DBG(30) << "transformed:"; for(uint j=0;j<transformed2.size();++j) BLINK(30) << transformed2[j] << " "; BLINK(30) << endl;});
          DBGI(30,{DBG(30) << "transformed:"; for(uint j=0;j<transformed3.size();++j) BLINK(30) << transformed3[j] << " "; BLINK(30) << endl;});
          lfPixel[i].resize(dimTotal);
          for(uint j = 0;j<dimPixel;++j) { lfPixel[i][j] = transformed1[j];}
          for(uint j = dimPixel; j < 2 * dimPixel; ++j) { lfPixel[i][j] = transformed2[j - dimPixel];}
          for(uint j = 2 * dimPixel; j < 3 * dimPixel; ++j) { lfPixel[i][j] = transformed3[j - 2 * dimPixel]; }
        }
        lfPixel.dim() = dimTotal;
        lfPixel.save(filename + ".pca." + suffix);
        filename = cl.next(" ");
	
      } 
    }
  }

  DBG(10) << "cmdline was: "; printCmdline(argc,argv);
}
