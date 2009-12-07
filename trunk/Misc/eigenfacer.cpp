#include <iostream>
#include <string>
#include "pca.hpp"
#include "imagefeature.hpp"
#include "imagelib.hpp"
#include <fstream>
#include <sstream>
#include "getpot.hpp"
#include <limits>

#ifdef HAVE_FFT_LIBRARY
extern "C" {
  #include FFTW_INCLUDE
}
#endif


using namespace std;

void USAGE() {
  cout << "USAGE: " <<endl
       << "eigenfacer <command> <options>" << endl
       << "   commands/options:" << endl
       << "    -E <filelist> Estimate PCA from images in filelist"<< endl
       << "        -c <filename> save covariance to file <filename>" << endl
       << "        -t <filename> save transformation to file <filename>" << endl
       << "    -T <filelist> Transform images from filelist using PCA" << endl
       << "        -t <filename> load transformation from file <filename>" << endl
       << "        -d <dim> set dimensionality for transformation"<< endl
       << "    -B <filelist> backtransform images from filelist using PCA" << endl
       << "         -t <filename> load transformation from file <filename>" << endl
       << "         -x <widht>, -y <height> set height and width of output images" << endl
       << "    -M <image> create face probability map for image" << endl
       << "         -t <filename> load transformation from file <filename>" << endl
       << "         -x <widht>, -y <height> set height and width of output images" << endl
       << "         -d <dim> set dimensionality for transformation"<< endl
       << "         -s <scaleingfactor> factor to scale image with (default: 0.83333=1/1.2)"<< endl
       << "    -F <image> create face probability map for image using Fourier implementation" << endl
       << "         -t <filename> load transformation from file <filename>" << endl
       << "         -x <widht>, -y <height> set height and width of output images" << endl
       << "         -d <dim> set dimensionality for transformation"<< endl
       << "         -s <scaleingfactor> factor to scale image with (default: 0.83333=1/1.2)"<< endl
       << endl;
}


double getEnergy(const vector<double>& vec) {
  double sum=0.0;
  for(uint i=0;i<vec.size();++i) {
    sum+=vec[i]*vec[i];
  }
  return sum;
}

#ifdef HAVE_FFT_LIBRARY
double energTrans(fftw_complex **F,const uint i,const uint D) {
  double sum=.0;
  for(uint d=0;d<D;++d) {
    sum+=F[d][i].re*F[d][i].re;
  }
  return sum;
}
#endif

double energImg(const ImageFeature &img, const uint x, const uint y,const uint w,const uint h,const vector<double>& mean) {
  double sum=.0;
  double im;
  for(uint i=x;i<x+w;++i) {
    for(uint j=y;j<y+h;++j) {
      if(j<img.ysize() && i<img.xsize()) im=img(i,j,0); 
      else im=0.0;
      double tmp=im-mean[(j-y)*w+(i-x)];
      sum+=tmp*tmp;
    }
  }
  return sum;
}



ImageFeature detect(const ImageFeature& img, const PCA& pca,uint dim, uint w, uint h) {
  uint imgdim=img.xsize()*img.ysize();
  uint padx=img.xsize()+w; uint pady=img.ysize()+h;
  uint paddim=padx*pady;
  
#ifdef HAVE_FFT_LIBRARY

  // get memory for calculations
  fftw_complex *FIMG=NULL, **FFILTER=NULL, **FMULT=NULL;
  FIMG=new fftw_complex[paddim]; 
  for(uint i=0;i<paddim;++i) {FIMG[i].re=0.0;FIMG[i].im=0.0;}
  FFILTER=new fftw_complex*[dim];
  FMULT=new fftw_complex*[dim];
  for(uint i=0;i<dim;++i) {
    FFILTER[i]=new fftw_complex[paddim];
    for(uint j=0;j<imgdim;++j) {FFILTER[i][j].re=0.0;FFILTER[i][j].im=0.0;}
    
    FMULT[i]=new fftw_complex[paddim];
    for(uint j=0;j<imgdim;++j) {FMULT[i][j].re=0.0;FMULT[i][j].im=0.0;}
  }
  
  DBG(10) << "Memory allocated" << endl;
  
  vector<double> meanTrans=pca.transform(pca.mean(),dim);
  
  // create strategy for fft
  fftwnd_plan plan = fftw2d_create_plan(padx,pady, FFTW_FORWARD, FFTW_ESTIMATE  | FFTW_IN_PLACE); 
  fftwnd_plan planb = fftw2d_create_plan(padx,pady, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE); 
  
  DBG(10) << "Strategies for FFT created" << endl;
  
  //copy image into fourier transform data structure
  for(uint x=0;x<img.xsize();++x) { for(uint y=0;y<img.ysize();++y) { FIMG[y*padx+x].re=img(x,y,0); } }
  
  // fourier transform the image
  fftwnd_one(plan,FIMG,NULL);
  
  DBG(10) << "Image Transformed" << endl;
  
  
  // fourier transform the filters
  for(uint d=0;d<dim;++d) {
    for(uint x=0;x<w;++x) {
      for(uint y=0;y<h;++y) {
        uint i=y*padx+x;
        FFILTER[d][i].re=pca.eigenvector(d)[y*w+x];
      }
    }
    fftwnd_one(plan,FFILTER[d],NULL);
    DBG(10) << "Filter " << d << " transformed." << endl;
  }
  
  // multiplication in fourier domain
  for(uint d=0;d<dim;++d) {
    for(uint i=0;i<paddim;++i) {
      FMULT[d][i].re=FIMG[i].re*FFILTER[d][i].re-FIMG[i].im*FFILTER[d][i].im;
      FMULT[d][i].im=FIMG[i].re*FFILTER[d][i].im+FIMG[i].im*FFILTER[d][i].re;
    }
    DBG(10) << "Filter " << d << " applied." ;
    //fourier back transform
    fftwnd_one(planb,FMULT[d],NULL);
    BLINK(10) << "... backtransformed.." ;
    
    
    // subtract transformed mean
    for(uint i=0;i<paddim;++i) {
      FMULT[d][i].re=FMULT[d][i].re/paddim-meanTrans[d];
    }
    BLINK(10) << ". Mean subtracted." << endl;
  }
  
  
  double energyTrans, energyImg;
  
  ImageFeature fpm(img.xsize(),img.ysize(),1);
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      uint i=y*padx+x;
      energyTrans=energTrans(FMULT,i,dim);
      energyImg=energImg(img,x,y,w,h,pca.mean());
      DBG(25) << VAR(x) << " " << VAR(y) << " " << VAR(energyTrans) << " " << VAR(energyImg) << endl;
      fpm(x,y,0)=energyImg-energyTrans;
    }
  }

  DBG(10) << "fpm generation." << endl;
  return fpm;
  #endif

}


int main(int argc , char **argv) {
  
  GetPot cl(argc,argv);
  
  if(cl.search("-h")) {USAGE(); exit(0);}

  string  line;
  if(cl.search("-E")) {
    cout << "Estimating PCA" << endl;
    
    ifstream filelist(cl.follow("filelist","-E"));
    ImageFeature img;    
    getline(filelist,line);

    img.load(line,true);
    cout << line << endl;
    
    PCA pca(img.size());
    cout << img.size() << endl;
    
    pca.putData(img.layerVector(0));
    
    while(getline(filelist,line)) {
      cout << line << endl;
      img.load(line,true);
      pca.putData(img.layerVector(0));
    }
    pca.dataEnd();
    pca.save(cl.follow("covariance.pca","-c"));
    pca.calcPCA();
    pca.save(cl.follow("transformation.pca","-t"));
    filelist.close();
    
  } else if(cl.search("-T")) {
    PCA pca;
    pca.load(cl.follow("transformation.pca","-t"));
    int dim=cl.follow(20,"-d");
    
    ifstream filelist(cl.follow("filelist","-T"));
    vector<double> tmp;
    ImageFeature img;
    while(getline(filelist,line)) {
      cout << line << endl;
      img.load(line,true);
      tmp=pca.transform(img.layerVector(0),img.size());
      tmp.resize(dim);
      VectorFeature tmpvec(tmp);
      tmpvec.save(line+".pca.vec.gz");
    }
    filelist.close();
  } else if(cl.search("-B")) {
    PCA pca;
    pca.load(cl.follow("transformation.pca","-t"));
    int x=cl.follow(16,"-x");
    int y=cl.follow(16,"-y");
    ifstream filelist(cl.follow("filelist","-B"));
    vector<double> tmp;
    VectorFeature tmpvec;
    while(getline(filelist,line)) {
      cout << line << endl;
      tmpvec.load(line);
      vector<double> tmp;
      tmp=pca.backTransform(tmpvec.data());
      ImageFeature img(tmp,x,y);
      cutoff(img);
      img.save(line+".backpca.png");
    }
    filelist.close();
  } else if(cl.search("-M")) {
    ImageFeature img; img.load(cl.follow("image.png","-M"),true);
    PCA pca;
    vector<double> backproj,vec;
    pca.load(cl.follow("transformation.pca","-t"));
    int w=cl.follow(16,"-x");
    int h=cl.follow(16,"-y");
    double scalefac=cl.follow(0.8333,"-s");
    int dim=cl.follow(20,"-d");
    
    ImageFeature scimg(img);
    ImageFeature patch;
    uint minX=100000, minY=100000; 
    uint maxX=100000, maxY=100000; 
   
    while(int(scimg.xsize())>=w and int(scimg.ysize())>=h) {
      DBG(10) << VAR(scimg.xsize()) << " x " << VAR(scimg.ysize()) << endl;
      ImageFeature faceprobmap(scimg.xsize(),scimg.ysize(),1);
      double maxDist=0.0;
      double minDist=numeric_limits<double>::max();
      
      vector<double> tmpvec;

      for(uint x=0;x<scimg.xsize();++x) {
        DBG(10) << VAR(x) << endl;
        for(uint y=0;y<scimg.ysize();++y) {
          patch=getPatch(scimg,x,y,x+w,y+h);
          vec=patch.layerVector(0);
          
          vector<double> imgMinMean=vec;
          for(uint i=0;i<imgMinMean.size();++i) {
            imgMinMean[i]-=pca.mean()[i];
          }
          double energyImg=getEnergy(imgMinMean);
          
          tmpvec=pca.transform(vec,dim);
          double energyTrans=getEnergy(tmpvec);
          backproj=pca.backTransform(tmpvec);
          
          double energyBack=getEnergy(backproj);

          double d=0;
          double tmp;
          for(uint i=0;i<backproj.size();++i) {
            tmp=backproj[i]-vec[i];
            d+=tmp*tmp;
          }
          faceprobmap(x,y,0)=d;
          
          DBG(10) << VAR(energyImg) << " "
                  << VAR(energyTrans) << " " 
                  << VAR(energyBack) << " "
                  << VAR(energyImg-energyTrans) << " "
                  << VAR(d) << endl;


          if(minDist>d) {
            minDist=d;
            minX=x; minY=y;
          }
          if(maxDist<d) {
            maxDist=d;
            maxX=x; maxY=y;
          }
        }
      }
      
      DBG(10) << VAR(scimg.xsize()) << " " << VAR(scimg.ysize()) << endl;
      DBG(10) << VAR(minDist) << " " << VAR(minX) << " " << VAR(minY) << endl;
      DBG(10) << VAR(maxDist) << " " << VAR(maxX) << " " << VAR(maxY) << endl << endl;

      normalize(faceprobmap);
      ostringstream filenamestream;
      filenamestream << cl.follow("image.png","-M") << ".fpm." << (scimg.xsize()) <<".png";
      faceprobmap.save(filenamestream.str());
      
      uint newW=int(scimg.xsize()*scalefac);
      uint newH=int(scimg.ysize()*scalefac);
      scimg=scale(scimg,newW,newH);
    }
    
  } else if(cl.search("-F")) {
#ifdef HAVE_FFT_LIBRARY
    ImageFeature img; img.load(cl.follow("image.png","-F"),true);
    PCA pca; pca.load(cl.follow("transformation.pca","-t"));
    int w=cl.follow(16,"-x");
    int h=cl.follow(16,"-y");
    double scalefac=cl.follow(0.8333,"-s");
    uint dim=cl.follow(20,"-d");

    DBG(10) << "Eigenfaces loaded" << endl;

    
    ImageFeature scimg=img;
    while(int(scimg.xsize())>w and int(scimg.ysize())>h) {
      ImageFeature fpm=detect(scimg,pca,dim,w,h);

      pair<uint,uint> p;

      p=argmax(fpm);
      
      DBG(10) << scimg.xsize() << "x" << scimg.ysize() << " (" <<p.first<<", "<< p.second << ") " << VAR(maximum(fpm)) ;

      p=argmin(fpm);
      BLINK(10) << " (" <<p.first<<", "<< p.second << ") " << VAR(minimum(fpm)) << endl;

      normalize(fpm);
      
      ostringstream filenamestream;
      filenamestream << cl.follow("image.png","-F") << ".fpm." << (scimg.xsize()) <<".png";
      fpm.save(filenamestream.str());
      
      scimg=scale(scimg,int(scimg.xsize()*scalefac),int(scimg.ysize()*scalefac));
    }
  
    

#else
    DBG(10) << "compiled without FFT lib. this does not work. use -M option" << endl;
#endif
    
  } else {
    USAGE();
    exit(20);
  }
}
