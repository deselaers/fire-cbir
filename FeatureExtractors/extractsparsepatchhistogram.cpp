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
#include "pca.hpp"
#include "getpot.hpp"
#include "gzstream.hpp"
#include "diag.hpp"
#include "imagelib.hpp"
#include "sparsehistogramfeature.hpp"
#include "differenceofgaussian.hpp"
#include "pasannotation.hpp"

using namespace std;

struct less_than_string {
  bool operator()(const string s1, const string s2) {
    return strcmp(s1.c_str(), s2.c_str()) < 0;
  }
};

void USAGE() {
  cout << "USAGE:" << endl
       << "extractsparsepatchhistogram [options] (--gray|--color) (--images filename1 [filename2 ... ]|--filelist <filelist>)" << endl
       << "   Options:" << endl
       << "    ~~~~~~~~~~~~~~~~~ general options (always applicable) ~~~~~~~~~~~~~~~~" << endl
       << "    -h,--help          shows this help" << endl
       << "    -su,--suffix=S     uses S as the suffix of output files, default is 'patch.sparsehisto.gz'" << endl
       << "    -c,--commonhisto=S does not create a histogram for each input file but merges files with a" << endl
       << "                       common prefix into a single histogram file with the prefix's name." << endl
       << "                       the prefixes are read from a file named S" << endl
       << "    -st,--steps=N      uses N steps for each dimension in the histogram, default value is 4" << endl
       << "    -pr,--prune=N      prunes bins that contain less than N% items" << endl
       << "                       N must be in the range of [0.0, 100.0], default: 0.0 (no pruning)" << endl
       << "    -sc, -scale=N      scales the images to a fixed height of N pixels and a" << endl
       << "                       variable width with constant aspect ratio, default: 0 (no scaling)" << endl
       << "                       if 0.0 < N < 1.0, then the image is scaled in both directions by the factor N" << endl
       << "    -pos,--position    treat x/y coordinates as additional dimensions" << endl
       << "    ~~~~~~~~~~~~~~~~~~~~~~~~ pca reduction options ~~~~~~~~~~~~~~~~~~~~~~~" << endl
       << "    -r,--reduction=N   reduces dimension of patches to N using PCA, default: 0 (no reduction)" << endl
       << "    -spca,--savePCA=S  saves the calculated pca and the histogram's" << endl
       << "                       step sizes to a file named S, this option requires -r to be set also" << endl
       << "    --noextraction     only calculate PCA and save it into file, do not create histograms" << endl
       << "    -lpca,--loadPCA=S  loads pca and histograms's step sizes, from a file named S" << endl
       << "                       requires -r to be set also, must no be used in combination with -mo" << endl
       << "                       because it is not necessary since the bin-boundary information is saved in a file beside the PCA" << endl
       << "    -mo,--minmaxmode=N calcuale histogram's stepsize base on:" << endl
       << "                       N=1: precomputed minimal/maximal values (default)" << endl
       << "                       N=2: actual minimal/maximal values" << endl
       << "                       N=3: actual mean/deviation (cut outliers to fit)" << endl
       << "                       requires -r to be set also, must not be used in combination with -lpca" << endl
       << "    -i, --interval=a   when using minmaxmode 3, use a to set the histogram's boundaries" << endl
       << "                       to [mean - a * sdev, mean + a * sdev], default is a=1.25" << endl
       << "    -l,--leaveoutfirst leaves out the first dimension. requires -r to be set also" << endl
       << "    -d,--decsteps=N    decrease steps with increasing dimension number" << endl
       << "                       from <steps> for the 1. dimension to N for the" << endl
       << "                       last dimension. requires -r to be set also" << endl
       << "    ~~~~~~~~~~~~~~~~~~~~~~~~~ patch size options ~~~~~~~~~~~~~~~~~~~~~~~~~" << endl
       << "    ~~~~~          either both -hp and -vp, both -mp and -ws         ~~~~~" << endl
       << "    ~~~~~            or both -dp and -ws must be specified           ~~~~~" << endl
       << "    -hp,--hpatch=N         uses N as horizontal patch size" << endl
       << "    -vp,--vpatch=N         uses N as vertical path size" << endl
       << "    -mp,--multipatch=sizes use multiple quadratic patches, whose sizes are given by 'sizes'," << endl
       << "                           a space separated list of integers" << endl
       << "    -dp,--DoGpatches       determine patch size by Differece-of-Gaussian method" << endl
       << "    -ws,--winsize=N        scale the patches to N, the side length of a quadratic patch" << endl
       << "    ~~~~~~~~~~~~~~~~~~~~~~~~~ annotations ~~~~~~~~~~~~~~~~~~~~~~~~~" << endl
       << "   --annotationlist read annotation files from this filelist in the same order as from the other filelist" << endl
       << "       this option does only have an effect for the extraction of histograms itself, not for the estimation of "<< endl
       << "       the PCA matrix neither for the estimation of the stepsizes" << endl
       << "   --featcls <suffix>    save featcls files" << endl
       << endl;
  exit(20);
}

// Extracts a patch of dimension hExtract, vExtract with the center x,y from the image, and, if necessary, scales it to
// the dimensions hPatch, vPatch. The extracted patch is stored in SHPoint
void getPatch(const ImageFeature& img, int x, int y, int hExtract, int vExtract, int hPatch, int vPatch, SHPoint& pt) {
  pt.clear();
  // extract the patch
  ImageFeature extractedPatch = getPatch(img, x - hExtract / 2, y - vExtract / 2, x + (hExtract - 1) / 2, y + (vExtract - 1) / 2);
  // scale the extracted patch, if necessary
  if ((hPatch != hExtract) || (vPatch != vExtract)) {
    ImageFeature scaledPatch = scale(extractedPatch, hPatch, vPatch);
    pt = toVector(scaledPatch);
  } else {
    pt = toVector(extractedPatch);
  }
  if ((int) pt.size() != (hPatch * vPatch * (int) img.zsize())) {
    ERR << "error in getPatch(...)!" << endl;
    ERR << "size of patch vector: " << pt.size() << endl;
    ERR << "hPatch=" << hPatch << ", vPatch=" << vPatch << endl;
    ERR << "img.zsize() = " << img.zsize() << endl;
    exit(1);
  }
}

void saveStepSizes(const string &filename, const vector<double> max, const vector<double> min) {
  ogzstream ofs; ofs.open(filename.c_str());
  if(!ofs.good()) {
    ERR << "Cannot write step sizes to '" << filename << "'."<< endl;
  } else {
    ofs << max.size() << endl;
    for (vector<double>::const_iterator it = max.begin(); it != max.end(); it++) {
      ofs << *it << " " << endl;
    }
    ofs << min.size() << endl;
    for (vector<double>::const_iterator it = min.begin(); it != min.end(); it++) {
      ofs << *it << " " << endl;
    }
    ofs << -1 << endl;
  }
}

void loadStepSizes(const string &filename, vector<double>& max, vector<double>& min) {
  igzstream ifs; ifs.open(filename.c_str());
  int dimmax, dimmin, posy;

  if (!ifs.good()) {
    ERR << "Cannot open '" << filename << "' to read step sizes" << endl;
  } else {
    ifs >> dimmax;
    max = vector<double>(dimmax);
    for (int i = 0; i < dimmax; i++) {
      ifs >> max[i];
    }
    ifs >> dimmin;
    if (dimmax != dimmin) {
      ERR << "Strange: max and min vector have different sizes!" << endl;
      exit(20);
    }
    min = vector<double>(dimmin);
    for (int i = 0; i < dimmin; i++) {
      ifs >> min[i];
    }
    ifs >> posy;
    if (posy != -1) {
      ERR << "Reading step sizes was strange. Did not end with -1." << endl;
      exit(20);
    }
  }
}

void printSettings(string suffix, bool commonHistograms, int steps, double threshold,
                   double scaleSize, bool position, int reduction, bool savePCA, bool loadPCA,
                   int minMaxMode, double interval,  bool leaveoutfirst, int decsteps, vector<int> extrSizes,
                   int winsizeX, int winsizeY, bool DoG) {
  DBG(10) << "settings:" << endl;
  DBG(10) << "===========================" << endl;
  DBG(10) << "suffix=" << suffix << endl;
  DBG(10) << "common histograms? " << (commonHistograms? "yes": "no") << endl;
  DBG(10) << "steps=" << steps << endl;
  DBG(10) << "prune threshold=" << threshold << endl;
  DBG(10) << "scale=" << scaleSize << endl;
  DBG(10) << "position? " << (position? "yes": "no") << endl;
  DBG(10) << "dimensionality reduction=" << reduction << endl;
  if (reduction > 0) {
    DBG(10) << "save PCA? " << savePCA << endl;
    DBG(10)  << "load PCA? " << loadPCA << endl;
    DBG(10) << "discretisation mode=" << minMaxMode << endl;
    if (minMaxMode == 3) {
      DBG(10) << "boundaries: [mean - " << interval << "*sdev, mean + " << interval << "*sdev]" << endl;
    }
    DBG(10) << "leaveoutfirst? " << (leaveoutfirst? "yes": "no") << endl;
    DBG(10) << "decrease stepsizes=" << decsteps << endl;
  }
  if (DoG) {
    DBG(10) << "extraction sizes with DoG method" << endl;
  } else {
    DBG(10) << "extraction sizes: ";
    for (int i = 0; i < (int) extrSizes.size(); i++) {
      BLINK(10) << extrSizes[i] << " ";
    }
    BLINK(10) << endl;
  }
  DBG(10) << "winsize: " << winsizeX << "," << winsizeY << endl;
}

int main(int argc, char** argv) {
  GetPot cl(argc,argv);

  const int STEP_SIZE_DEFAULT = 4,
    HPATCH_SIZE_DEFAULT = 3,
    VPATCH_SIZE_DEFAULT = 3,
    REDUCTION_DEFAULT = 0,
    MINMAXMODE_DEFAULT = 1;
  const uint MIN_SCALE = 5;
  const string SUFFIX_DEFAULT = "patch.sparsehisto.gz";
  const double PRUNE_DEFAULT = 0.0,
    SCALE_DEFAULT = 0.0,
    INTERVAL_DEFAULT = 1.25;

  int steps = STEP_SIZE_DEFAULT,
    reduction = REDUCTION_DEFAULT,
    hPatch = HPATCH_SIZE_DEFAULT,
    vPatch = VPATCH_SIZE_DEFAULT,
    decSteps = STEP_SIZE_DEFAULT,
    minMaxMode = MINMAXMODE_DEFAULT;
  string suffix = SUFFIX_DEFAULT,
    pcaFile,
    stepSizesFile;
  double prune = PRUNE_DEFAULT,
    interval = INTERVAL_DEFAULT,
    scaleSize = SCALE_DEFAULT;

  bool colorImage = false, 
    commonHistograms = false, 
    withReduction = false, 
    hasPca = false, 
    savePca = false,
    leaveOutFirstDimension = false,
    spatialInformation = false,
    noextraction = false;

  map<string, SparseHistogramFeature*, less_than_string> prefixMap;
  vector<string> prefixes;
  vector<int> extractionSizes; // stores extraction sizes in the format x1, y1, x2, y2, ...
  vector<string> infiles;  // for the list of files to be processed

  PCA pca; // for dimensionality reduction, either calculated or loaded from file
  vector<double> max, min, premax, premin; // for histogram step size, either calculated or loaded from file

  //command line parsing
  if (cl.search(2, "--help", "-h")) {
    USAGE();
  }

  if (cl.search(2, "--suffix", "-su")) {
    suffix = cl.next(" ");
  }

  // read prefixes for common histograms
  if (cl.search(2, "-c", "--commonhisto")) {
    commonHistograms = true;
    string prefixFile = cl.next(" ");
    igzstream ifs; ifs.open(prefixFile.c_str());
    if(!ifs.good() || !ifs) {
      ERR << "Cannot open prefix file " << prefixFile  << ". Aborting." << endl;
      exit(20);
    }
    string prefix = " ";
    DBG(10) << "Create common histograms for files with the following prefixes: " << endl;
    while(!ifs.eof() && (prefix != "")) {
      getline(ifs, prefix);
      if (prefix.size() > 0) {
        prefixes.push_back(prefix);
        DBG(10) << "\t" << prefix << endl;
      }
    }
    ifs.close();
  }

  if (cl.search(2, "-st", "--stepsize")) {
    steps = cl.next(STEP_SIZE_DEFAULT);
    decSteps = steps;
  }

  if (cl.search(2, "-pos", "--position")) {
    DBG(10) << "adding spatial information to histogram" << endl;
    spatialInformation = true;
  }

  // read size of patch(es) to extract
  bool DoG = false;
  if (cl.search(2, "-vp", "--vpatch") && cl.search(2, "-hp", "--hpatch")) {
    // single patch
    cl.search(2, "-vp", "--vpatch");
    vPatch = cl.next(VPATCH_SIZE_DEFAULT);
    cl.search(2, "-hp", "--hpatch");
    hPatch = cl.next(HPATCH_SIZE_DEFAULT);
    extractionSizes.push_back(hPatch);
    extractionSizes.push_back(vPatch);
  } else if (cl.search(2, "-mp", "--multipatch") && cl.search(2, "-ws", "--winsize")) {
    // multiple patches
    cl.search(2, "-ws", "--winsize");
    int winsize = cl.next(0); // ignore the 0
    hPatch = winsize;
    vPatch = winsize;
    cl.search(2, "-mp", "--multipatch");
    int exS=cl.next(-1);
    while(exS!=-1) {
      extractionSizes.push_back(exS);
      extractionSizes.push_back(exS); // push twice for horizontal and vertical size of patch
      exS=cl.next(-1);
    }
    if (extractionSizes.empty()) {
      USAGE();
    }
  } else if (cl.search(2, "-dp", "--DoGpatch") && cl.search(2, "-ws", "--winsize")) {
    DoG = true;
    DBG(15) << "using Difference-of-Gaussian for patch extraction" << endl;
    cl.search(2, "-ws", "--winsize");
    int winsize = cl.next(0); // ignore the 0
    hPatch = winsize;
    vPatch = winsize;
  } else {
    USAGE();
  }
  // sanity check for patch size
  if ((hPatch < 1) || (vPatch < 1)) {
    ERR << "Invalid patch size, must be >= 1 !" << endl;
    exit(20);
  }

  if (cl.search(2, "-pr", "--prune")) {
    prune = cl.next(PRUNE_DEFAULT);
  }

  // process pca options
  if (cl.search(2, "-r", "--reduction")) {
    reduction = cl.next(REDUCTION_DEFAULT);
    withReduction = true;
  }
  if (cl.search(2, "-mo", "--minmaxmode")) {
    if (!withReduction) {
      USAGE();
    }
    minMaxMode = cl.next(MINMAXMODE_DEFAULT);
    if (cl.search(2, "-lpca", "--loadpca")) {
      USAGE();
    }
    if (cl.search(2, "-i", "--interval")) {
      interval = cl.next(INTERVAL_DEFAULT);
    }
  } else {
    if (cl.search(2, "-i", "--interval")) {
      USAGE();
    }
    if (cl.search(2, "-lpca", "--loadPCA")) {
      if (!withReduction) {
        USAGE();
      }
      hasPca = true;
      // load pca and step sizes
      string loadFilePca = cl.next(" ");
      string loadFileStepSizes = loadFilePca;
      loadFilePca.append("-pca.gz");
      loadFileStepSizes.append("-stepSizes.gz");
      pca.load(loadFilePca);
      loadStepSizes(loadFileStepSizes, max, min);
    }
  }
  if (cl.search(2, "-spca", "--savePCA")) {
    if (!withReduction) {
      USAGE();
    }
    savePca = true;
    pcaFile = cl.next(" ");
    stepSizesFile = pcaFile;
    pcaFile.append("-pca.gz");
    stepSizesFile.append("-stepSizes.gz");
    noextraction = cl.search(1, "--noextraction");
  } else {
    if (cl.search(1, "--noextraction")) {
      ERR << "option --noextraction required -spca/--savePCA to be set also!" << endl;
      USAGE();
    }
  }
  if (cl.search(2, "-l", "--leaveoutfirst")) {
    if (!withReduction) {
      USAGE();
    }
    leaveOutFirstDimension = true;
  }
  if (cl.search(2, "-d", "--decsteps")) {
    if (!withReduction) {
      USAGE();
    }
    decSteps = cl.next(steps);
  }

  if (cl.search(2, "-sc", "--scale")) {
    scaleSize = cl.next(SCALE_DEFAULT);
  }

  if(cl.search("--color")) {
    colorImage = true;
  } else if (cl.search("--gray")) {
    colorImage = false;
  } else {
    USAGE();
  }

  bool featcls=cl.search("--featcls");
  string featclssuffix="featcls";
  if(featcls) {
    featclssuffix=cl.follow("featcls","--featcls");
  }

  bool useAnnotation=false;
  std::vector< std::string> annotationfiles;
  if(cl.search("--annotationlist")) {
    DBG(10) << "Using annotation" << endl;
    useAnnotation=true;
    string annotationlist=cl.follow("annolist","--annotationlist");
    igzstream ifs; ifs.open(annotationlist.c_str());
    std::string filename;
    while(ifs.good() && !ifs.eof()) {
      getline(ifs,filename);
      annotationfiles.push_back(filename);
    }
    ifs.close();
  }

  // read files to process
  if(cl.search("--images")) {
    string filename=cl.next(" ");;
    while(filename!=" ") {
      infiles.push_back(filename);
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
        infiles.push_back(filename);
      }
    }
    ifs.close();
  } else {
    USAGE();
  }

  // print settings
  printSettings(suffix, commonHistograms, steps, prune, scaleSize, spatialInformation, reduction, 
                savePca, hasPca, minMaxMode, interval, leaveOutFirstDimension, decSteps, extractionSizes,
                hPatch, vPatch, DoG);

  // end of command line parsing

  int patchNrDimensions = hPatch * vPatch * (colorImage? 3: 1);
  int histoNrDimensions = patchNrDimensions;
  if (reduction > patchNrDimensions) {
    cout << "cannot reduce to " << reduction << " dimensions, patches contain only " << patchNrDimensions << " dimensions" << endl;
    exit(20);
  } else if (reduction > 0) {
    histoNrDimensions = reduction;
    DBG(20) << "reducing patch to " << histoNrDimensions << " dimensions" << endl;
    if (leaveOutFirstDimension) {
      DBG(20) << "leaving out first dimension (brightness normalization)" << endl;
    }
  }

  // used only if patch sizes are determined by Difference-of-Gaussian
  vector<short int> sizesDoG;

  if (!hasPca) {
    // no pca yet loaded
    max = vector<double>(histoNrDimensions + (spatialInformation? 2: 0));
    min = vector<double>(histoNrDimensions + (spatialInformation? 2: 0));
    premax = vector<double>(histoNrDimensions + (spatialInformation? 2: 0));
    premin = vector<double>(histoNrDimensions + (spatialInformation? 2: 0));
    if (reduction > 0) {
      // calculated it
      // for PCA reduction, run over all files to prepare PCA
      pca = PCA(patchNrDimensions);
      DBG(10) << "collection patches for PCA" << endl;
      for(uint i = 0; i < infiles.size(); i++) {
        string filename = infiles[i];
        DBG(10) << "processing " << filename << "(" << i <<"/" << infiles.size() << ")" << endl;
        ImageFeature img;
        if (colorImage) {
          img.load(filename);
        } else {
          img.load(filename, true);
        }
        if (scaleSize > 1.0) {
          img = scale(img, (img.xsize() * (int) scaleSize) / img.ysize(), (int) scaleSize);
        } else if (scaleSize > 0.0) {
          img = scale(img, (int) ((double) img.xsize() * scaleSize), (int) ((double) img.ysize() * scaleSize));
        }
        int patchCount = 0;
        if (DoG) {
          // extract patches for patch size calculated by Difference-of-Gaussian
          DifferenceOfGaussian sift(img);
          vector<uint> allScales = sift.getAllScales();
          //double meanScale = 0.0;
          //double devScale = 0.0;
          int pixelCount = 0;
          for (int x = 0; x < (int) img.xsize(); x++) {
            for (int y = 0; y < (int) img.ysize(); y++) {
              uint curScale = allScales[pixelCount];
              //meanScale += double(curScale);
              //devScale += double(curScale * curScale);
              //cout << "scale for (" << x << "," << y << ") is " << curScale << endl;
              if (curScale >= MIN_SCALE) {
                SHPoint pt;
                getPatch(img, x, y, curScale / 2, curScale / 2, hPatch, vPatch, pt);
                pca.putData(pt);
                patchCount++;
              }
              pixelCount++;
            }
          }
          //exit(1);
          //meanScale /= double(count);
          //devScale /= double(count);
          //devScale -= (meanScale * meanScale);
          //DBG(15) << "scale of extracted DoG patches: mean=" << meanScale << " var=" << devScale << endl;
        } else {
          // extract patches for all sizes
          for (vector<int>::iterator it = extractionSizes.begin(); it != extractionSizes.end(); it++) {
            int xSize = *it;
            it++;
            int ySize = *it;
            for (int x = xSize / 2; x < (int) img.xsize() - ((xSize - 1) / 2); x++) {
              for (int y = ySize / 2; y < (int) img.ysize() - ((ySize - 1) / 2); y++) {
                patchCount++;
                SHPoint pt; 
                getPatch(img, x, y, xSize, ySize, hPatch, vPatch, pt);
                pca.putData(pt);
              }
            }
          }
        }
        BLINK(15) << " (" << patchCount << ")" << endl;
      }
      pca.dataEnd();
      DBG(10) << "calculating PCA" << endl;
      pca.calcPCA();
      // now that we have calculated the PCA, we need to determine the histogram boundaries
      // this is controlled by the 'minMaxMode' option
      if (minMaxMode != 1) {
        // determine the boundaries by either the actual minimal and maximal value
        // or determine the interval based on medium value and variance
        // both methods required another run over the data :-(
        DBG(10) << "calculating histogram boundaries" << endl;
        vector<double> mean(histoNrDimensions), var(histoNrDimensions);
        for (int d = 0; d < reduction; d++) {
          mean[d] = 0.0;
          var[d] = 0.0;
        }
        bool first = true;
        int patchCounter = 0;
        for(uint i = 0; i < infiles.size(); i++) {
          string filename = infiles[i];
          ImageFeature img;
          if (colorImage) {
            img.load(filename);
          } else {
            img.load(filename, true);
          }
          if (scaleSize > 1.0) {
            img = scale(img, (img.xsize() * (int) scaleSize) / img.ysize(), (int) scaleSize);
          } else if (scaleSize > 0.0) {
            img = scale(img, (int) ((double) img.xsize() * scaleSize), (int) ((double) img.ysize() * scaleSize));
          }
          if (DoG) {
            // extract patches for patch size calculated by Difference-of-Gaussian
            DifferenceOfGaussian sift(img);
            vector<uint> allScales = sift.getAllScales();
            int count = 0;
            for (int x = 0; x < (int) img.xsize(); x++) {
              for (int y = 0; y < (int) img.ysize(); y++) {
                uint curScale = allScales[count];
                if (curScale >= MIN_SCALE) {
                  SHPoint pt;
                  getPatch(img, x, y, curScale / 2, curScale / 2, hPatch, vPatch, pt);
                  SHPoint reducedPt(reduction);
                  reducedPt = pca.transform(pt, reduction);
                  for (int d = 0; d < reduction; d++) {
                    mean[d] += reducedPt[d];
                    var[d] += reducedPt[d] * reducedPt[d];
                    if ((reducedPt[d] > max[d]) || first) {
                      max[d] = reducedPt[d];
                    }
                    if ((reducedPt[d] < min[d]) || first) {
                      min[d] = reducedPt[d];
                    }
                  }
                  first = false;
                  patchCounter++;
                }
                count++;
              }
            }
          } else {
            for (vector<int>::iterator it = extractionSizes.begin(); it != extractionSizes.end(); it++) {
              int xSize = *it;
              it++;
              int ySize = *it;
              for (int x = xSize / 2; x < (int) img.xsize() - ((xSize - 1) / 2); x++) {
                for (int y = ySize / 2; y < (int) img.ysize() - ((ySize - 1) / 2); y++) {
                  SHPoint pt; 
                  getPatch(img, x, y, xSize, ySize, hPatch, vPatch, pt);
                  SHPoint reducedPt(reduction);
                  reducedPt = pca.transform(pt, reduction);
                  for (int d = 0; d < reduction; d++) {
                    mean[d] += reducedPt[d];
                    var[d] += reducedPt[d] * reducedPt[d];
                    if (first || (reducedPt[d] > max[d])) {
                      max[d] = reducedPt[d];
                    }
                    if (first || (reducedPt[d] < min[d])) {
                      min[d] = reducedPt[d];
                    }
                  }
                  first = false;
                  patchCounter++;
                }
              }
            }
          }
        }
        for (int d = 0; d < reduction; d++) {
          mean[d] /= double(patchCounter);
          var[d] /= double(patchCounter);
          var[d] -= mean[d] * mean[d];
          DBG(15) << "maximal value for dimension " << d << " is " << max[d] << endl;
          DBG(15) << "minimal value for dimension " << d << " is " << min[d] << endl;
          DBG(15) << "mean value for dimension " << d << " is " << mean[d] << endl;
          DBG(15) << "variance for dimension " << d << " is " << var[d] << endl;
          if (minMaxMode == 3) {
            max[d] = mean[d] + interval * sqrt(var[d]);
            min[d] = mean[d] - interval * sqrt(var[d]);
          }
        }
      } else if (minMaxMode == 1) {
        // in this case, we use information provided by pca to calculate the boundaries for the histogram
        // this is quick and safe, but somehow inaccurate
        for (int d = 0; d < histoNrDimensions; d++) {
          premax[d] = 0.0;
          premin[d] = 0.0;
          for (int c = 0; c < patchNrDimensions; c++) {
            double inMax = (pca.mean()[d] > 0.5)? pca.mean()[d]: 1.0 - pca.mean()[d];
            premax[d] += fabs(pca.covariance()[d][c]) * inMax;
            premin[d] -= fabs(pca.covariance()[d][c]) * inMax;
          }
          DBG(10) << "precomputed maximal value for dimension " << d << " is " << premax[d] << endl;
          DBG(10) << "precomputed minimal value for dimension " << d << " is " << premin[d] << endl;
          // set precomputed maximal/minimal values
          max[d] = premax[d];
          min[d] = premin[d];
        }
      }
      // the boundaries for position dimensions do not required precalculation, can be set to [0.0, 1.0]
      if (spatialInformation) {
        max[max.size() - 2] = 1.0;
        min[min.size() - 2] = 0.0;
        max[max.size() - 1] = 1.0;
        min[min.size() - 1] = 0.0;
      }
      // save pca and max/min values, if desired
      if (savePca) {
        DBG(10) << "saving pca to " << pcaFile << endl;
        pca.save(pcaFile);
        DBG(10) << "saving step sizes to " << stepSizesFile << endl;
        saveStepSizes(stepSizesFile, max, min);
        if (noextraction) {
          DBG(10) << "--noextraction was set, therefore exiting" << endl;
          exit(0);
        }
      }
    } else {
      // no pca at all
      for (int d = 0; d < (int) max.size(); d++) {
        max[d] = 1.0;
        min[d] = 0.0;
      }
    }
  }

  // if the first dimension should be left out, adjust max/min vectors
  if (leaveOutFirstDimension) {
    max.erase(max.begin());
    min.erase(min.begin());
  }

  // set step size vector for sparse histogram
  // implement step size reduction for increasing dimension number, if enabled
  vector<uint> stepVector;
  for (int i = (leaveOutFirstDimension? 1: 0); i < histoNrDimensions; i++) {
    if ((reduction > 0) && (decSteps > 0)) {
      stepVector.push_back(steps - (int) (double((steps - decSteps) * i) / double(histoNrDimensions - 1)));
    } else {
      stepVector.push_back(steps);
    }
  }
  if (spatialInformation) { 
    // two additional dimensions for spatial information
    stepVector.push_back(steps);
    stepVector.push_back(steps);
  }

  // if all output goes into a single file create the sparse histogram here
  // otherwise do it later for each file individually
  SparseHistogramFeature* shf = NULL;
  if (commonHistograms) {
    for (vector<string>::iterator it = prefixes.begin(); it != prefixes.end(); it++) {
      shf = new SparseHistogramFeature(stepVector);
      shf->max() = max;
      shf->min() = min;
      shf->initStepsize();
      prefixMap[*it] = shf;
    }
  }

  // processing the files
  DBG(10) << "creating histograms" << endl;
  for(uint i=0;i<infiles.size();++i) {
    string filename=infiles[i];
    DBG(10) << "Processing '" << filename << "' (" << i+1<< "/" << infiles.size() << ")." << endl;

    // load the image and scale it, if necessary
    ImageFeature img;
    if (colorImage) {
      img.load(filename);
    } else {
      img.load(filename, true);
    }
    
    PascalAnnotationFeature pasanno;
    if(useAnnotation) { pasanno.load(annotationfiles[i]); }
    
    if (scaleSize > 1.0) {
      int height = img.ysize();
      int width = img.xsize();
      img = scale(img, (width * (int) scaleSize) / height, (int) scaleSize);
      DBG(20) << "scaled image to " << img.xsize() << "x" << img.ysize();
      if(useAnnotation) {ERR << "Scaling annotation to fixed width not supported" << endl;}
    } else if (scaleSize > 0.0) {
      int height = img.ysize();
      int width = img.xsize();
      img = scale(img, (int) ((double) width * scaleSize), (int) ((double) height * scaleSize));
      DBG(20) << "scaled image to " << img.xsize() << "x" << img.ysize();
      if(useAnnotation) {pasanno.resize(scaleSize);}
    }

    // if we use common histograms, see if the current file matches any of the mentioned prefixes
    bool foundPrefix = false;
    shf = NULL;
    if (commonHistograms) {
      for (vector<string>::iterator it = prefixes.begin(); it != prefixes.end(); it++) {
        if (strncmp(it->c_str(), filename.c_str(), it->size()) == 0) {
          shf = prefixMap[*it];
          foundPrefix = true;
          break;
        }
      }
    }

    std::vector<SparseHistogramFeature> histograms;//(pasanno.nOfObjects()+1,SparseHistogramFeature(stepVector,min,max));
    if (!foundPrefix) {
      // create a new sparse histogram, if necessary, and initialize it
      shf = new SparseHistogramFeature(stepVector);
      shf->max() = max;
      shf->min() = min;
      shf->initStepsize();
    }
    
    // fill the sparse histogram

    if (DoG) { 
      // extract patches for patch size calculated by Difference-of-Gaussian
      DifferenceOfGaussian sift(img);
      vector<uint> allScales = sift.getAllScales();
      int count = 0;
      for (int x = 0; x < (int) img.xsize(); x++) {
        for (int y = 0; y < (int) img.ysize(); y++) {
          uint curScale = allScales[count];
          if (curScale >= MIN_SCALE) {
            SHPoint pt;
            getPatch(img, x, y, curScale / 2, curScale / 2, hPatch, vPatch, pt);
            if (reduction > 0) {
              SHPoint reducedPt(reduction);
              reducedPt = pca.transform(pt, reduction);
              if (leaveOutFirstDimension > 0) {
                // remove the first dimension, if desired
                reducedPt.erase(reducedPt.begin());
              }
              shf->feed(reducedPt);
            } else {
              shf->feed(pt);
            }
          }
          count++;
        }
      }
    } else { // !DoG
      if(!useAnnotation) {
        Position pos; ofstream featclsOs;
        if(featcls) {
          string fn=filename+"."+featclssuffix;
          featclsOs.open(fn.c_str());
        }
        for (vector<int>::iterator it = extractionSizes.begin(); it != extractionSizes.end(); it++) {
          int xSize = *it;
          it++;
          int ySize = *it;
          for (int x = xSize / 2; x < (int) img.xsize() - ((xSize - 1) / 2); x++) {
            for (int y = ySize / 2; y < (int) img.ysize() - ((ySize - 1) / 2); y++) {
              SHPoint pt; pt.clear();
              getPatch(img, x, y, xSize, ySize, hPatch, vPatch, pt);
              if (reduction > 0) {
                SHPoint reducedPt(reduction);
                reducedPt = pca.transform(pt, reduction);
                if (leaveOutFirstDimension > 0) {
                  // remove the first dimension, if desired
                  reducedPt.erase(reducedPt.begin());
                }
                if (spatialInformation) {
                  // add x and y coordinates, if desired
                  reducedPt.push_back(double(x) / double(img.xsize()));
                  reducedPt.push_back(double(y) / double(img.ysize()));
                }
                pos=shf->feed(reducedPt);
                if(featcls) {
                  featclsOs << x << " " << y << " " << xSize << " " << pos << endl;
                }
              } else {
                if (spatialInformation) {
                  // add x and y coordinates, if desired
                  pt.push_back(double(x) / double(img.xsize()));
                  pt.push_back(double(y) / double(img.ysize()));
                }
                shf->feed(pt);
              }
            }
          }
        }
        featclsOs.close();
      } else { // useAnnotation
        histograms=vector<SparseHistogramFeature>(pasanno.nOfObjects()+1,SparseHistogramFeature(stepVector,min,max));
        for (vector<int>::iterator it = extractionSizes.begin(); it != extractionSizes.end(); it++) {
          int xSize = *it;
          it++;
          int ySize = *it;
          for (int x = xSize / 2; x < (int) img.xsize() - ((xSize - 1) / 2); x++) {
            for (int y = ySize / 2; y < (int) img.ysize() - ((ySize - 1) / 2); y++) {
              SHPoint pt; pt.clear();
              getPatch(img, x, y, xSize, ySize, hPatch, vPatch, pt);
              if (reduction > 0) {
                SHPoint reducedPt(reduction);
                reducedPt = pca.transform(pt, reduction);
                if (leaveOutFirstDimension > 0) {
                  // remove the first dimension, if desired
                  reducedPt.erase(reducedPt.begin());
                }
                if (spatialInformation) {
                  // add x and y coordinates, if desired
                  reducedPt.push_back(double(x) / double(img.xsize()));
                  reducedPt.push_back(double(y) / double(img.ysize()));
                }

                bool inside=false;
                for(uint o=0;o<pasanno.nOfObjects();++o) {
                  if(pasanno.object(o).inside(x,y)) {
                    histograms[o].feed(reducedPt);
                    inside=true;
                  }
                }
                if(! inside) {
                  histograms[pasanno.nOfObjects()].feed(reducedPt);
                }
              }
            }
          }
        }
      }
    }
    if(useAnnotation) {
      for(uint i=0;i<pasanno.nOfObjects();++i) {
        histograms[i].calcRelativeFrequency();
        if (prune > 0) {
          histograms[i].prunePoorBins(prune / 100.0);
        }

        ostringstream fns; fns<< filename << "-obj" << i <<"-label-"<< pasanno.object(i).label_ <<"."<<suffix;
        histograms[i].save(fns.str());
      }

      int bgno=pasanno.nOfObjects();
      histograms[bgno].calcRelativeFrequency();
      if (prune > 0) {
        histograms[bgno].prunePoorBins(prune / 100.0);
      }
     
      histograms[bgno].save(filename+".bg."+suffix);
      


    } else {
      shf->calcRelativeFrequency();
      // write the histogram to a file if it is not a common histogram
      if (!foundPrefix) {
        // if desired, prune poorly filled bins before writing
        if (prune > 0) {

           shf->prunePoorBins(prune / 100.0);
        }
        string output = filename + "." + suffix;
        DBG(20) << "saving histogram to file '" << output  << "'..." << endl;
        shf->save(output);
        delete shf;
      }
    }
    DBG(20) << "Finished with '" << filename << "'." << endl;
  }

  // if we have common sparse histogram, this is the time to save them
  if (commonHistograms) {
    for (vector<string>::iterator it = prefixes.begin(); it != prefixes.end(); it++) {
      shf = prefixMap[*it];
      // if desired, prune poorly filled bins
      if (prune > 0) {
        shf->prunePoorBins(prune / 100.0);
      }
      // write the sparse histogram to a file
      string output = (*it) + "." + suffix;
      shf->save(output);
    }
  }

  DBG(1) << "cmdline was: "; printCmdline(argc,argv);
}
