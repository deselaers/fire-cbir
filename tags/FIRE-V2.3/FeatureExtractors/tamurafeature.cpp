#include <vector>
#include <algorithm>
#include "tamurafeature.hpp"
#include "imagelib.hpp"
#include <cmath>
#include <limits>

using namespace std;

ImageFeature calculate(const ImageFeature &img) {
  ImageFeature result(img.xsize(), img.ysize(), 3);

  ImageFeature tmp;
  // process the layers
  for(uint z=0;z<img.zsize();++z) {
    tmp=directionality(img,z);
    add(result,tmp,0);
    
    tmp=coarseness(img,z);
    add(result,tmp,1);

    tmp=contrast(img,z);
    add(result,tmp,2);
  }
  
  divide(result,img.zsize());
//  normalize(result);
//  result.save("result.png");
  return result;
}


ImageFeature contrast(const ImageFeature &image, uint z) {
  int xdim=image.xsize();
  int ydim=image.ysize();
  
  ImageFeature result(xdim,ydim,1);
  
  double min=numeric_limits<double>::max();
  double max=0;
  double tmp;
  DBG(20) << "Starting to get contrast for image." << endl;

  for(int x=0;x<xdim;++x) {
    for(int y=0;y<ydim;++y) {
      tmp=getLocalContrast(image,x,y,z);
      result(x,y,0)=tmp;
      min=::std::min(min,tmp);
      max=::std::max(max,tmp);
    }
  }
  
  DBG(20) << "minContrast: " << min << " maxContrast: " << max << endl;
  return result;
}

double getLocalContrast(const ImageFeature &image, int xpos, int ypos, uint z) {
  int xdim=image.xsize();  
  int ydim=image.ysize();
  int ystart=::std::max(0,ypos-5);
  int xstart=::std::max(0,xpos-5);
  int ystop=::std::min(ydim,ypos+6);
  int xstop=::std::min(xdim,xpos+6);
  
  int size=(ystop-ystart)*(xstop-xstart);
  
  double mean=0.0, sigma=0.0, kurtosis=0.0, tmp;

  for(int y=ystart;y<ystop;++y) {
    for(int x=xstart;x<xstop;++x) {
      tmp=image(x,y,z);
      mean+=tmp;
      sigma+=tmp*tmp;
    }
  }
  mean/=size;
  sigma/=size;
  sigma-=mean*mean;
  
  for(int y=ystart;y<ystop;++y) {
    for(int x=xstart;x<xstop;++x) {
      tmp=image(x,y,z)-mean;
      tmp*=tmp;
      tmp*=tmp;
      kurtosis+=tmp;
    }
  }
  kurtosis/=size;
  //double alpha4=kurtosis/(sigma*sigma);
  //double contrast=sqrt(sigma)/sqrt(sqrt(alpha4));
  //return contrast;
  
  /// if we don't have this exception, there are numeric problems!
  /// if the region is homogeneous: kurtosis and sigma are numerically very close to zero
  if(kurtosis<numeric_limits<double>::epsilon()) {
    return 0.0;
  }

  return sigma/sqrt(sqrt(kurtosis));
}


ImageFeature directionality(const ImageFeature &image, uint z) {
  DBG(20) << "Starting to get directionality" << endl;  
  
  //init
  uint xdim=image.xsize();
  uint ydim=image.ysize();  
  
  ImageFeature deltaH(xdim, ydim, 1);
  ImageFeature deltaV(xdim, ydim, 1);
  
  for(uint x=0;x<xdim;++x) {
    for(uint y=0;y<ydim;++y) {
      deltaH(x,y,0)=image(x,y,z);
      deltaV(x,y,0)=image(x,y,z);
    }
  }
  
  //step1
  ImageFeature matrixH(3,3,1), matrixV(3,3,1);
  matrixH(0,0,0)=-1;  matrixH(0,1,0)=-2;  matrixH(0,2,0)=-1;
  matrixH(2,0,0)=1;   matrixH(2,1,0)=2;   matrixH(2,2,0)=1;
  
  matrixV(0,0,0)=1;  matrixV(1,0,0)=2;  matrixV(2,0,0)=1;
  matrixV(0,2,0)=-1; matrixV(1,2,0)=-2; matrixV(2,2,0)=-1;
  
  convolve(deltaH,matrixH);
  convolve(deltaV,matrixV);

  //step2
  //ImageFeature deltaG(xdim,ydim,1);
  ImageFeature phi(xdim,ydim,1);
  for(uint y=0;y<ydim;++y) {
    for(uint x=0;x<xdim;++x) {
      //deltaG(x,y,0)=fabs(deltaH(x,y,0))+fabs(deltaV(x,y,0));
      //deltaG(x,y,0)*=0.5;

      if(deltaH(x,y,0)!=0.0) {
        phi(x,y,0)=atan(deltaV(x,y,0)/deltaH(x,y,0))+(M_PI/2.0+0.001); //+0.001 because otherwise sometimes getting -6.12574e-17
      }
    }
  }
  
  return phi;
}


double efficientLocalMean(const int x,const int y,const int k,const ImageFeature &laufendeSumme) {
  int k2=k/2;
  
  int dimx=laufendeSumme.xsize();
  int dimy=laufendeSumme.ysize();
  
  //wanting average over area: (y-k2,x-k2) ... (y+k2-1, x+k2-1)
  int starty=y-k2;
  int startx=x-k2;
  int stopy=y+k2-1;
  int stopx=x+k2-1;
  
  if(starty<0) starty=0;
  if(startx<0) startx=0;
  if(stopx>dimx-1) stopx=dimx-1;
  if(stopy>dimy-1) stopy=dimy-1;
  
  double unten, links, oben, obenlinks;
  
  if(startx-1<0) links=0; 
  else links=laufendeSumme(startx-1,stopy,0);
  
  if(starty-1<0) oben=0;
  else oben=laufendeSumme(stopx,starty-1,0);
  
  if((starty-1 < 0) || (startx-1 <0)) obenlinks=0;
  else obenlinks=laufendeSumme(startx-1,starty-1,0);
  
  unten=laufendeSumme(stopx,stopy,0);
  
//   cout << "obenlinks=" << obenlinks << " oben=" << oben << " links=" << links << " unten=" <<unten << endl;
  int counter=(stopy-starty+1)*(stopx-startx+1);
  return (unten-links-oben+obenlinks)/counter;
}


ImageFeature coarseness(const ImageFeature &image,const uint z) {
  const uint yDim=image.ysize();
  const uint xDim=image.xsize();

  ImageFeature laufendeSumme(xDim,yDim,1);
  
  // initialize for running sum calculation
  double links, oben, obenlinks;
  for(int y=0;y<int(yDim);++y) {
    for(int x=0;x<int(xDim);++x) {
      if(x<1) links=0;
      else links=laufendeSumme(x-1,y,0);
      
      if(y<1) oben=0;
      else oben=laufendeSumme(x,y-1,0);
      
      if(y<1 || x<1) obenlinks=0;
      else obenlinks=laufendeSumme(x-1,y-1,0);
      
      laufendeSumme(x,y,0)=image(x,y,z)+links+oben-obenlinks;
    }
  }
  
  ImageFeature Ak(xDim,yDim,5);
  ImageFeature Ekh(xDim,yDim,5);
  ImageFeature Ekv(xDim,yDim,5);

  ImageFeature Sbest(xDim,yDim,1);
  
  DBG(25) << "  ... step 1 ... " << endl;
  //step 1
  int lenOfk=1;
  for(int k=1;k<=5;++k) {
    lenOfk*=2;
    for(int y=0;y<int(yDim);++y) {
      for(int x=0;x<int(xDim);++x) {
        Ak(x,y,k-1)=efficientLocalMean(x,y,lenOfk,laufendeSumme);
      }
    }
  }
  
  DBG(25) << "  ... step 2 ... " << endl;
  //step 2
  lenOfk=1;
  for(int k=1;k<=5;++k) {
    int k2=lenOfk;
    lenOfk*=2;
    for(int y=0;y<int(yDim);++y) {
      for(int x=0;x<int(xDim);++x) {
        
        int posx1=x+k2;
        int posx2=x-k2;
        
        int posy1=y+k2;
        int posy2=y-k2;
        if(posx1<int(xDim) && posx2>=0) {
          Ekh(x,y,k-1)=fabs(Ak(posx1,y,k-1)-Ak(posx2,y,k-1));
        } else {
          Ekh(x,y,k-1)=0;
        }
        if(posy1<int(yDim) && posy2>=0) {
          Ekv(x,y,k-1)=fabs(Ak(x,posy1,k-1)-Ak(x,posy2,k-1));
        } else {
          Ekv(x,y,k-1)=0;
        }
      }
    }
  }
  double sum=0.0;
  DBG(25) << "  ... step 3 ... " << endl;  
  //step3
  for(int y=0;y<int(yDim);++y) {
    for(int x=0;x<int(xDim);++x) {
      double maxE=0;
      int maxk=0;
      for(int k=1;k<=5;++k) {
        if(Ekh(x,y,k-1)>maxE) {
          maxE=Ekh(x,y,k-1);
          maxk=k;
        }
        if(Ekv(x,y,k-1)>maxE) {
          maxE=Ekv(x,y,k-1);
          maxk=k;
        }
      }
      Sbest(x,y,0)=maxk;
      sum+=maxk;
    }
  }
 
  sum/=((yDim-32)*(xDim-32));
  DBG(25) << "Average Coarseness: " << sum << endl;
  return Sbest;
}

vector<uint> histogramImageBinImage(const ImageFeature &image) {
  vector<uint> result(image.xsize()*image.ysize());
  vector<double> t(3,0);
  
  vector<uint> steps(3,8);
  HistogramFeature tmp(steps);
  tmp.min()[0]=0;    tmp.max()[0]=1.0;
  tmp.min()[1]=0;    tmp.max()[1]=1.0;
  tmp.min()[2]=0;    tmp.max()[2]=1.0;
  
  tmp.initStepsize();
  
  for(uint x=0;x<image.xsize();++x) {
    for(uint y=0;y<image.ysize();++y) {
      for(uint z=0;z<3;++z) {t[z]=image(x,y,z);}
      result[y*image.xsize()+x]=tmp.posToBin(tmp.pointToPos(t));
    }
  }
  
  return result;
}


HistogramFeature histogramize(const ImageFeature image,
                              const vector<uint> &bins, 
                              const uint left, 
                              const uint top, 
                              const uint right, 
                              const uint bottom) {
  vector<uint> steps(3,8);
  HistogramFeature result(steps);
  result.min()[0]=0;    result.max()[0]=1.0;
  result.min()[1]=0;    result.max()[1]=1.0;
  result.min()[2]=0;    result.max()[2]=1.0;
  result.initStepsize();
  
  uint xsize=image.xsize();
  
  for(uint x=left;x<right;++x) {
    for(uint y=top;y<bottom;++y) {
      //      DBG(10) << y << " " << x << " "<< y*xsize+x << " " << bins[y*xsize+x] << endl;
      result.feedbin(bins[y*xsize+x]);
    }
  }
  
  return result;
}


HistogramFeature histogramize(const ImageFeature &image) {
  // min -> max: 0,0,0 -> pi,5,128
  // steps: 8,8,8
  
  vector<uint> steps(3,8);
  HistogramFeature result(steps);
  result.min()[0]=0;    result.max()[0]=1.0;
  result.min()[1]=0;    result.max()[1]=1.0;
  result.min()[2]=0;    result.max()[2]=1.0;
  result.initStepsize();

  vector<double> tmp(3,0.0);

  for(uint x=0;x<image.xsize();++x) {
    for(uint y=0;y<image.ysize();++y) {
      for(uint z=0;z<3;++z) {tmp[z]=image(x,y,z);}
      result.feed(tmp);
    }
  }
  return result;
}
