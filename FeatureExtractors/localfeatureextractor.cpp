
#include "localfeatureextractor.hpp"
#include "diag.hpp"
#include "imagelib.hpp"
#include "differenceofgaussian.hpp"
#include "salientpoints.hpp"
#include "histogramfeature.hpp"

#include <vector>
#include <sstream>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

string LocalFeatureExtractorSettings::suffix() const {
  ostringstream oss;

  if(patches) oss << "patches";
  else if(sift) oss << "sift";
  else if(histo) oss << "histo";
  else {oss << "unclear"; ERR << "Undefined suffix" << endl;}

  if(waveletPoints>0) oss << "-wave" << waveletPoints;
  if(dogPoints>0) oss << "-dog" << dogPoints;
  if(randomPoints>0) oss << "-rp" << randomPoints;
  if(gridX>0) oss << "-ux" << gridX;
  if(gridY>0) oss << "-uy" << gridY;
  oss << "-ws" << winsize;
  if(alignedGrid) oss << "-ag";
  if(alignGrid) oss << "-ga";
  if(forceGray) oss << "-gray";
  if(padding>0) oss << "-pad" << padding;
  if(pca) oss << "-pca" << pcadim;
  if(addonsuffix != "") oss << "-" << addonsuffix;

  oss << ".lf.gz";
  return oss.str();
}

int LocalFeatureExtractorSettings::dim() const {
  int dim=0;
  if(patches) {
    dim=2*winsize+1;
    dim*=dim;

    if(not forceGray) {
      dim*=3;
    }
  }
  else if(histo)
  {
	if(forceGray)
		dim = histosteps;
	else
		dim = histosteps*histosteps*histosteps;
  }
  else if(sift)
  {
    ERR << "NYI" << endl;
  }


  DBG(35) << VAR(dim) << endl;
  return dim;

}


void LocalFeatureExtractor::extractFromDatabase(const Database &db) {
  string suffix=settings_.suffix();

  DBG(10) << "Settings for extraction:" << endl
          << settings_ << endl;

  PCA pca(settings_.dim());
  if(settings_.loadpca) {
    pca.load(settings_.pcafile);
  }

  string tmppath=settings_.tmppath+"/"+getenv("USER");
  mkdir(tmppath.c_str(), 0xFFFFFF);
  tmppath+="/extractlf";
  mkdir(tmppath.c_str(), 0xFFFFFF);

  for(uint n=0;n<db.size();++n) {
    string imagefilename=db.path()+"/"+db.filename(n);
    string imagebasename=db.filename(n).substr(db.filename(n).find_last_of("/")+1);

    LocalFeatures lf;
    ImageFeature img;
    DBG(10) << "Extracting features from '" <<imagefilename << "': " << n+1 << "/" << db.size() << endl;
    img.load(imagefilename);
    lf=extractFromImage(img,lf);
    if(settings_.pcaonly){
    	lf.save(tmppath+"/"+imagebasename+"."+suffix);
    } else {
    	lf.save(imagefilename+"."+suffix);
    }

    if(settings_.pca and not settings_.loadpca) {
      for(uint l=0;l<lf.data_.size();++l) {
        pca.putData(lf[l]);
      }
    } else if(settings_.loadpca) {
      pcaTransform(pca,lf);
      lf.save(imagefilename+".pca."+suffix);
    }
    DBG(10) << "Extracted " << lf.numberOfFeatures_ << " features." << endl;
  }

  if(settings_.pca and not settings_.loadpca) {
    DBG(10) << "Estimating PCA transformation matrix." << endl;
    pca.dataEnd();
    pca.calcPCA();
    pca.save(settings_.pcafile+"."+suffix+".pca.gz");

    for(uint n=0;n<db.size();++n) {
      string imagefilename=db.path()+"/"+db.filename(n);
      string imagebasename=db.filename(n).substr(db.filename(n).find_last_of("/")+1);
      DBG(10) << "PCA transforming features for '" << imagefilename << "': " << n+1 << "/" << db.size() << endl;

      LocalFeatures lf;
      if(settings_.pcaonly){
    	  lf.load(tmppath+"/"+imagebasename+"."+suffix);
    	  remove((tmppath+"/"+imagebasename+"."+suffix).c_str());
      } else {
    	  lf.load(imagefilename+"."+suffix);
      }

      pcaTransform(pca,lf);
      lf.save(imagefilename+".pca."+suffix);
    }
  }

  rmdir(tmppath.c_str());
}


LocalFeatures& LocalFeatureExtractor::extractFromImage(const ImageFeature &image, LocalFeatures& lf) {
  // get a copy of the image
  ImageFeature img(image);
  vector<FeatureExtractionPosition> positions;

  preprocess(img);
  getExtractionPositions(img,positions);

  if(settings_.patches) { extractPatches(img, positions, lf);}
  else if(settings_.sift) { extractSIFT(img, positions, lf);}
  else if(settings_.histo) { extractHisto(img, positions, lf);}
  else { ERR << "Nothing to extract?" << endl; }

  return lf;
}

void LocalFeatureExtractor::preprocess(ImageFeature &img) {
  if(settings_.forceGray) {
    img=makeGray(img,1);
  }

  // padding is done in the extractPatches routine
}

void LocalFeatureExtractor::getExtractionPositions(const ImageFeature &img, ::std::vector<FeatureExtractionPosition> &extractionPositions) {

  int currentPosition=0;

  if(settings_.dogPoints>0) {
    DBG(15) << "Extracting DoG salient points" << endl;

    DifferenceOfGaussian dog(img);
    vector<InterestPoint> interestPoints = dog.getInterestPoints(settings_.dogPoints);

    FeatureExtractionPosition f;
    extractionPositions.resize(extractionPositions.size()+interestPoints.size());

    for(uint i=0;i<interestPoints.size();++i) {
      f.x=interestPoints[i].x;
      f.y=interestPoints[i].y;
      f.s=max(uint(3),interestPoints[i].scale/2);
      extractionPositions[currentPosition++]=f;
    }
    DBG(15) << "Now having " << extractionPositions.size() << " extraction positions." << endl;
  }

  if (settings_.waveletPoints>0) {
    DBG(15) << "Extracting wavelet-based salient points" << endl;
    SalientPoints sp(img);
    vector<Point> sPoints=sp.getSalientPoints(settings_.waveletPoints);

    extractionPositions.resize(extractionPositions.size()+sPoints.size()*settings_.extractionSizes.size());

    FeatureExtractionPosition f;

    for(uint i=0;i<sPoints.size();++i) {
      for(uint s=0;s<settings_.extractionSizes.size();++s) {
        f.x=sPoints[i].x;
        f.y=sPoints[i].y;
        f.s=settings_.extractionSizes[s];

        extractionPositions[currentPosition++]=f;
      }
    }
    DBG(15) << "Now having " << extractionPositions.size() << " extraction positions." << endl;
  }


  if(settings_.gridX>0 and settings_.gridY>0 and not settings_.alignGrid) {
    DBG(15) << "Determining grid positions" << endl;
    int curSize=extractionPositions.size();
    extractionPositions.resize(curSize+(settings_.gridX*settings_.gridY*settings_.extractionSizes.size()));

    uint stepx=max(uint(1),img.xsize()/(settings_.gridX+1));
    uint stepy=max(uint(1),img.ysize()/(settings_.gridY+1));

    FeatureExtractionPosition f;

    for(uint x=1;x<=settings_.gridX;++x) {
      for(uint y=1;y<=settings_.gridY;++y) {
        for(uint s=0;s<settings_.extractionSizes.size();++s) {
          f.x=x*stepx;
          f.y=y*stepy;
          f.s=settings_.extractionSizes[s];
          extractionPositions[currentPosition++]=f;
        }
      }
    }
    DBG(15) << "Now having " << extractionPositions.size() << " extraction positions." << endl;
  }

  if(settings_.gridX>0 and settings_.gridY>0 and settings_.alignGrid) {
    DBG(15) << "Determining aligned grid with given grid resolution" << endl;

    uint exs=img.xsize()/settings_.gridX;
    if(img.ysize()/settings_.gridY>exs) {exs=img.ysize()/settings_.gridY;}

    int curSize=extractionPositions.size();

    extractionPositions.resize(curSize+(settings_.gridX*settings_.gridY*1));

    exs/=2;

    uint stepx=max(uint(1),img.xsize()/(settings_.gridX+1));
    uint stepy=max(uint(1),img.ysize()/(settings_.gridY+1));

    FeatureExtractionPosition f;

    for(uint x=1;x<=settings_.gridX;++x) {
      for(uint y=1;y<=settings_.gridY;++y) {
        f.x=x*stepx;
        f.y=y*stepy;
        f.s=exs;
        extractionPositions[currentPosition++]=f;
      }
    }
    DBG(15) << "Now having " << extractionPositions.size() << " extraction positions." << endl;
  }




  if(settings_.randomPoints>0) {
    DBG(15) << "Determining random positions" << endl;

    srand(0);
    extractionPositions.resize(extractionPositions.size()+(settings_.randomPoints*settings_.extractionSizes.size()));

    FeatureExtractionPosition f;
    for(uint i=0;i<settings_.randomPoints;++i) {
      for(uint s=0;s<settings_.extractionSizes.size();++s) {
        f.x=rand() % img.xsize();
        f.y=rand() % img.ysize();
        f.s=settings_.extractionSizes[s];
        extractionPositions[currentPosition++]=f;
      }
    }
    DBG(15) << "Now having " << extractionPositions.size() << " extraction positions." << endl;
  }


  if(settings_.alignedGrid) {
    DBG(15) << "Determining positions on an aligned grid." << endl;

    FeatureExtractionPosition f;

    for(uint s=0;s<settings_.extractionSizes.size();++s) {
      uint size=2*settings_.extractionSizes[s]+1;

      uint y=size/2+1;

      while(y+size/2+1<img.ysize()) {
        uint x=size/2+1;
        while(x+size/2+1<img.xsize()) {
          f.x=x;
          f.y=y;
          f.s=size;

          DBG(55) << VAR(f) << endl;

          extractionPositions.push_back(f);
          x+=size;
        }
        y+=size;
      }
    }
    DBG(15) << "Now having " << extractionPositions.size() << " extraction positions." << endl;
    currentPosition=extractionPositions.size()-1;
  }

  if(extractionPositions.size()==0) {
    DBG(10) << "NO EXTRACTION POSITION FOUND! ADDING IMAGE CENTER FOR STABILITY!" << endl;
    FeatureExtractionPosition f;
    f.x=img.xsize()/2;
    f.y=img.ysize()/2;
    for(uint s=0;s<settings_.extractionSizes.size();++s) {
      uint size=2*settings_.extractionSizes[s]+1;
      f.s=size;
      extractionPositions.push_back(f);
    }
  }
  

}

void LocalFeatureExtractor::extractPatches(ImageFeature &img, const ::std::vector<FeatureExtractionPosition> &positions, LocalFeatures &lf) {
  if(settings_.padding>0) {
    img=zeropad(img,settings_.padding, settings_.padding);
  }

  lf.winsize_=settings_.winsize;
  lf.padding_=settings_.padding;
  lf.numberOfFeatures_=positions.size();
  lf.zsize_=img.zsize();

  int windiameter=2*settings_.winsize+1;

  lf.dim_=windiameter*windiameter*img.zsize();
  lf.imageSizeX_=img.xsize();
  lf.imageSizeY_=img.ysize();
  lf.filename_=img.filename();

  int savesize=settings_.winsize*2+1;

  lf.positions_=positions;
  lf.data_.resize(positions.size());

  for(uint i=0;i<positions.size();++i) {
    const FeatureExtractionPosition &pos=positions[i];

    ImageFeature p=getPatch(img,pos.x,pos.y,pos.s);
    if(uint(pos.s)!=settings_.winsize) {
      p=scale(p,savesize,savesize);
    }
    lf.data_[i]=toVector(p);
  }
}

void LocalFeatureExtractor::extractSIFT(ImageFeature &img, const ::std::vector<FeatureExtractionPosition> &positions, LocalFeatures &lf) {
  if(settings_.padding>0) {
    img=zeropad(img,settings_.padding, settings_.padding);
  }
}

void LocalFeatureExtractor::extractHisto(ImageFeature &img, const ::std::vector<FeatureExtractionPosition> &positions, LocalFeatures &lf) {

	uint n_colors;
	if(settings_.forceGray == true)
		n_colors = 1;
	else
		n_colors = 3;

	// Extract patches
	if(settings_.padding>0) {
		img=zeropad(img,settings_.padding, settings_.padding);
	}

	lf.winsize_=settings_.winsize;
	lf.padding_=settings_.padding;
	lf.numberOfFeatures_=positions.size();
	lf.zsize_=img.zsize();

	int windiameter=2*settings_.winsize+1;

	if(n_colors == 1)
		lf.dim_=settings_.histosteps;
	else if(n_colors == 3)
		lf.dim_=settings_.histosteps*settings_.histosteps*settings_.histosteps;
	lf.imageSizeX_=img.xsize();
	lf.imageSizeY_=img.ysize();
	lf.filename_=img.filename();

	int savesize=settings_.winsize*2+1;

	lf.positions_=positions;
	lf.data_.resize(positions.size());

	uint steps=settings_.histosteps;

	for(uint i=0;i<positions.size();++i) {
		const FeatureExtractionPosition &pos=positions[i];

		ImageFeature p=getPatch(img,pos.x,pos.y,pos.s);
		if(uint(pos.s)!=settings_.winsize) {
			p=scale(p,savesize,savesize);
		}

		// Histogramize patch
		HistogramFeature histo(vector<uint>(n_colors,steps));
		histo.min()=vector<double>(n_colors,0.0);
		histo.max()=vector<double>(n_colors,1.0);
		histo.initStepsize();

		vector<double> tofeed(n_colors);

		for(uint x=0;x<p.xsize();++x) {
			for(uint y=0;y<p.ysize();++y) {
				for(uint c=0;c<n_colors;++c) {
					tofeed[c]=p(x,y,c);
				}
				histo.feed(tofeed);
			}
		}
		lf.data_[i]=histo.data();
	}
}

void LocalFeatureExtractor::pcaTransform(const PCA &pca, LocalFeatures &lf) {
  for(uint l=0;l<lf.numberOfFeatures_;++l) {
    lf[l]=pca.transform(lf[l],settings_.pcadim);
  }
  lf.dim_=settings_.pcadim;
}

bool LocalFeatureExtractorSettings::read(::std::istream& is) const {
  ERR << "NYI" << endl;
  return false;
}
bool LocalFeatureExtractorSettings::write(::std::ostream &os) const {
  os << "waveletPoints "     << waveletPoints      << endl
     << "dogPoints "         << dogPoints          << endl
     << "randomPoints "      << randomPoints       << endl
     << "gridX "             << gridX              << endl
     << "gridY "             << gridY              << endl
     << "alignedGrid "       << alignedGrid        << endl
     << "forceGray "         << forceGray          << endl
     << "padding "           << padding            << endl
     << "patches "           << patches            << endl
     << "nOfextractionSizes "<< extractionSizes.size() << endl
     << "extractionSizes";
  for(uint i=0;i<extractionSizes.size();++i) os << " " <<extractionSizes[i];
  os << endl
     << "winsize "           << winsize            << endl
     << "sift "              << sift               << endl
     << "histo "             << histo     	       << endl
     << "pca "               << pca                << endl
     << "loadpca "           << loadpca            << endl
     << "pcadim "            << pcadim             << endl
     << "pcafile "           << pcafile            << endl;

  os << "Results in: " << endl
     << "suffix " << suffix() << endl
     << "dim " << dim() << endl;

  return true;
}
