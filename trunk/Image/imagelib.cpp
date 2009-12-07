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
#include <cmath>
#include <math.h>
#include <limits>
#include "imagelib.hpp"
#include "diag.hpp"
#include "imagefeature.hpp"
#include "interpolatingimage.hpp"
#include <map>
#ifdef HAVE_FFT_LIBRARY
extern "C" {
  #include FFTW_INCLUDE
}
#endif
#include <map>

using namespace std;

void sobelv(ImageFeature &img) {
  ImageFeature filter(3,3,1);
  filter(0,0,0)=-1; filter(1,0,0)=0; filter(2,0,0)=1;
  filter(0,1,0)=-2; filter(1,1,0)=0; filter(2,1,0)=2;
  filter(0,2,0)=-1; filter(1,2,0)=0; filter(2,2,0)=1;
  convolve(img,filter);
}
void sobelh(ImageFeature &img) {
  ImageFeature filter(3,3,1);
  filter(0,0,0)=-1; filter(1,0,0)=-2; filter(2,0,0)=-1;
  filter(0,1,0)=0; filter(1,1,0)=0; filter(2,1,0)=0;
  filter(0,2,0)=1; filter(1,2,0)=2; filter(2,2,0)=1;
  convolve(img,filter);
}
void sobeldiagonal1(ImageFeature &img) {
  ImageFeature filter(3,3,1);
  filter(0,0,0)=-2; filter(1,0,0)=-1; filter(2,0,0)=0;
  filter(0,1,0)=-1; filter(1,1,0)=0; filter(2,1,0)=1;
  filter(0,2,0)=0; filter(1,2,0)=1; filter(2,2,0)=2;
  convolve(img,filter);
}
void sobeldiagonal2(ImageFeature &img) {
  ImageFeature filter(3,3,1);
  filter(0,0,0)=0; filter(1,0,0)=-1; filter(2,0,0)=-2;
  filter(0,1,0)=1; filter(1,1,0)=0; filter(2,1,0)=-1;
  filter(0,2,0)=2; filter(1,2,0)=1; filter(2,2,0)=0;
  convolve(img,filter);
}

void laplace(ImageFeature &img) {
  ImageFeature filter(3,3,1);
  filter(0,0,0)=1; filter(1,0,0)=-2; filter(2,0,0)=1;
  filter(0,1,0)=-2; filter(1,1,0)=4; filter(2,1,0)=-2;
  filter(0,2,0)=1; filter(1,2,0)=-2; filter(2,2,0)=1; 
  convolve(img, filter);
}


void fftconvolve(ImageFeature &img, const ImageFeature &filter) {
#ifdef HAVE_FFT_LIBRARY

  // get size and smallest power of 2 that is larger
  uint width=img.xsize(); uint height=img.ysize();
  uint dim=width; if(height>dim) dim=height;
  uint maxnn=1; while(maxnn<dim) maxnn*=2;
  uint padded=dim;
  dim*=dim;
  
  // variables for image and filter in transformed
  fftw_complex *FIMG=NULL, *FFILTER=NULL;
  FIMG=new fftw_complex[dim]; FFILTER=new fftw_complex[dim];
  
  //copy image to temp variable
  uint xoffset=padded/2-width/2; uint yoffset=padded/2-height/2;
  for(uint y=0;y<height;++y) {
    for(uint x=0;x<width;++x) {
      uint idx=((y+yoffset)*padded)+x+xoffset;
      FIMG[idx].re=img(x,y,0);
    }
  }
  
  //copy filter to temp variable
  xoffset=padded/2-filter.xsize()/2; yoffset=padded/2-filter.ysize()/2;
  for(uint y=0;y<filter.ysize();++y) {
    for(uint x=0;x<filter.xsize();++x) {
      FFILTER[y*padded+x].re=filter(x,y,0);
    }
  }
  
  //fourier transform
  fftwnd_plan plan = fftw2d_create_plan(padded, padded, FFTW_FORWARD, FFTW_ESTIMATE  | FFTW_IN_PLACE); 
  fftwnd_one(plan,FIMG,NULL);
  fftwnd_one(plan,FFILTER,NULL);
  fftwnd_destroy_plan(plan);
  
  //multiplication in fourier domain
  for(uint x=0;x<dim;++x) {
    double re=FIMG[x].re*FFILTER[x].re-FIMG[x].im*FFILTER[x].im;
    double im=FIMG[x].re*FFILTER[x].im+FIMG[x].im*FFILTER[x].re;
    
    FIMG[x].re=re;
    FIMG[x].im=im;
  }
  
  //fourier transform backwards
  fftwnd_plan planb = fftw2d_create_plan(padded, padded, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE); 
  fftwnd_one(planb,FIMG,NULL);
  fftwnd_destroy_plan(planb);

  //copy back
  xoffset=padded/2-width/2;  yoffset=padded/2-height/2;
  for(uint y=0;y<height;++y) {
    for(uint x=0;x<width;++x) {
      uint idx=((y+yoffset)*padded)+x+xoffset;
      img(x,y,0)=FIMG[idx].re;
    }
  }
  delete[] FFILTER;
  delete[] FIMG;
  


#else
  DBG(10) << "Compiled without fftw library. Thus using normal convolution." << endl;
  convolve(img,filter);
#endif  
}


void convolve(ImageFeature &img, const ImageFeature &filter) {
  ImageFeature copy=img;
  
  int height2=filter.ysize()/2;
  int width2=filter.xsize()/2;

  double tmp;
  for(uint c=0;c<img.zsize();++c) {
    for(uint x=0;x<img.xsize();++x) {
      for(uint y=0;y<img.ysize();++y) {
        tmp=0.0;
        for(int i=-width2;i<=width2;++i) {
          int xx=x+i;
          if(xx<int(img.xsize()) && int(xx) >= 0) {
            for(int j=-height2;j<=height2;++j) {
              int yy=y+j; 
              if(int(yy)>=0 && yy < int(img.ysize())) {
                tmp+=filter(i+width2,j+height2,0)*copy(xx,yy,c);
              }
            }
          }
        }
        img(x,y,c)=tmp;
      }
    }
  }
}

void multiply(ImageFeature &img1, const ImageFeature& img2) {
  if ((img1.xsize() != img2.xsize()) ||
      (img1.ysize() != img2.ysize()) ||
      (img1.zsize() != img2.zsize())) {
    DBG(10) << "image features have different sizes, cannot be multiplied.";
  }
  for (uint z = 0; z < img1.zsize(); z++)
    for (uint x = 0; x < img1.xsize(); x++)
      for (uint y = 0; y < img1.ysize(); y++) {
	img1(x, y, z) = img1(x, y, z) * img2(x, y, z);
      }
}

void multiply(ImageFeature &img, const double s) {
  for(uint z=0;z<img.zsize();++z) {
    for(uint x=0;x<img.xsize();++x) {
      for(uint y=0;y<img.ysize();++y)  {
        img(x,y,z)*=s;
      }
    }
  }
}


void square(ImageFeature &img1) {
  multiply(img1, img1);
}

void power(ImageFeature &img1, int p) {
  for (uint z = 0; z < img1.zsize(); z++)
    for (uint x = 0; x < img1.xsize(); x++)
      for (uint y = 0; y < img1.ysize(); y++) {
	double tmp = 1.0;
	double val = img1(x, y, z);
	for (int t = 0; t < p; t++) {
	  tmp *= val;
	}
	img1(x, y, z) = tmp;
      }
}

void signedPower(ImageFeature &img1, int p) {
  for (uint z = 0; z < img1.zsize(); z++)
    for (uint x = 0; x < img1.xsize(); x++)
      for (uint y = 0; y < img1.ysize(); y++) {
	double tmp = 1.0;
	for (int t = 0; t < p; t++) {
	  tmp *= img1(x, y, z);
	}
	if ((img1(x, y, z) < 0.0) && ((p % 2) == 0)) {
	  tmp *= -1.0;
	}
	img1(x, y, z) = tmp;
      }
}

void absolute(ImageFeature &img1) {
  for (uint z = 0; z < img1.zsize(); z++)
    for (uint x = 0; x < img1.xsize(); x++)
      for (uint y = 0; y < img1.ysize(); y++) {
	if (img1(x, y, z) < 0.0) {
	  img1(x, y, z) = -img1(x, y, z);
	} 
      }
}

void normalize(ImageFeature &img) {
  double min, max;
  
  for(uint c=0;c<img.zsize();++c) {
    min=minimum(img,c);
    max=maximum(img,c);
    if (max - min != 0.0) {
      double d=(1.0/(max-min));
      for(uint x=0;x<img.xsize();++x) {
        for(uint y=0;y<img.ysize();++y) {
          img(x,y,c)=(img(x,y,c)-min)*d;
        }
      }
    }
  }
}

void gammaCorrection(ImageFeature &img, double gamma) {
  for(uint c=0;c<img.zsize();++c) {
    for(uint x=0;x<img.xsize();++x) {
      for(uint y=0;y<img.ysize();++y) {
        img(x,y,c)=exp(gamma*log(img(x,y,c)));
      }
    }
  }
}

double minimum(const ImageFeature &img, uint c) {
  double min=numeric_limits<double>::max();
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      min=::std::min(min,img(x,y,c));
    }
  }
  return min;
}

double maximum(const ImageFeature &img, uint c) {
  double max=-1.0;
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      max=::std::max(max,img(x,y,c));
    }
  }
  return max;
}


::std::pair<uint,uint> argmin(const ImageFeature &img,uint c) {
  pair<uint,uint> res;
  double min=numeric_limits<double>::max();
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      if(min>img(x,y,c)) {
        min=img(x,y,c); res.first=x; res.second=y;
      }
    }
  }
  return res;
}

::std::pair<uint,uint> argmax(const ImageFeature &img,uint c) {
  pair<uint,uint> res;
  double max=0.0;
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      if(max<img(x,y,c)) {
        max=img(x,y,c); res.first=x; res.second=y;
      }
    }
  }
  return res;
}



void add(ImageFeature &a, const ImageFeature &b, const uint z) {
  for(uint x=0;x<a.xsize();++x) {
    for(uint y=0;y<a.ysize();++y) {
      a(x,y,z)+=b(x,y,0);
    }
  }
}


void divide(ImageFeature &a, const double &d) {
  for(uint x=0;x<a.xsize();++x) {
    for(uint y=0;y<a.ysize();++y) {
      for(uint z=0;z<a.zsize();++z) {
        a(x,y,z)/=d;
      }
    }
  }
}

// default for type is set in imagelib.hpp
// type = 0: maximum (as in HSV)
// type = 1: mean 
// type = 2: luminance (use only for RGB images!)
ImageFeature makeGray(const ImageFeature &img, const uint type) {
  ImageFeature result(img.xsize(), img.ysize(),1) ;
  double value;
  
  if( (type == 2) && (img.zsize() != 3) )
  {
    ERR << "Image has invalid color resolution!" << endl;
    return result;
  }
  
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      value=0.0;
      switch(type){
      case 0: //max
	for(uint c=0;c<img.zsize();++c) {
	  value=::std::max(value,img(x,y,c));
	}
	break;
      case 1: //mean
        for(uint c=0;c<img.zsize();++c) {
          value+=img(x,y,c);
        }
        value/=img.zsize();
        break;
      case 2: // luminance
        value = 0.3 * img(x,y,0) + 0.59 * img(x,y,1) + 0.11 * img(x,y,2);
        break;
      }
      result(x,y,0)=value;
    }
  }
  return result;
}

ImageFeature rotate90(const ImageFeature &img) {
  uint X=img.xsize();
  ImageFeature result(img.ysize(),img.xsize(),img.zsize());
    for(uint x=0;x<img.ysize();++x) {
      for(uint y=0;y<img.xsize();++y) {
        for(uint z=0;z<img.zsize();++z) {
          result(x,y,z)=img(y,X-x-1,z);
        }
      }
    }
  return result;
}

ImageFeature flip(const ImageFeature &img) {
  ImageFeature result(img.xsize(),img.ysize(),img.zsize());
  
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      for(uint z=0;z<img.zsize();++z) {
        result(img.xsize()-x-1,y,z)=img(x,y,z);
      }
    }
  }
  return result;
}

void contrast(ImageFeature &img,const double min, const double max) {
  double m=max-min;
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      for(uint z=0;z<img.zsize();++z) {
        img(x,y,z)=img(x,y,z)*m+min;
      }
    }
  }
}

void saveLayer(uint idx, const ImageFeature &img, const string & filename) {
  ImageFeature image(img);
  normalize(image);
  image.save(filename,idx,idx,idx);
}



void calc_histo(int x1, int x2, int y1, int y2, float histo[256], const ImageFeature& image,int c) {
  int i,j;

  for(j=0; j < 256; j++) 
    histo[j] = 0.0;
  
  for(j=y1; j < y2; j++) {
    for(i=x1; i < x2; i++) {
      histo[(int)rint(image(i,j,c)*255.9-0.45)]++;
    }
  }
  
  const float t=1.0/((y2-y1)*(x2-x1));
  for(j=0; j < 256; j++) histo[j]*=t;
}


void eq_histo(int x1, int x2, int y1, int y2, const ImageFeature &image, int c, ImageFeature& acu, ImageFeature& divisor) {
  float histo[256];
  int i,j;
  
  calc_histo(x1,x2,y1,y2,histo,image,c);
  // cumulative histogram
  for(j=1; j < 256; j++){
    histo[j] += histo[j-1];
  }

  for(j=y1; j < y2; j++) {
    for(i=x1; i < x2; i++) {
      acu(i,j,0) += histo[(int)rint(image(i,j,c)*255.9-0.45)];
      // new grayvalue = (cum. histo.)(old grayvalue) 
      // example: 30 % of the grayvalues are below 120
      // -> new grayvalue (120) = 0.3
      divisor(i,j,0)+=1.0;
    }
  }
}

//default values for fact_eq and offset are set in imagelib.hpp.... by paredes
void equalize(ImageFeature &image, uint wndsize, double bnd, uint c, double fact_eq, uint offset) {
  double fact_orig=1.0-fact_eq;
  int xdim=image.xsize();
  int ydim=image.ysize();
  
  int i,j,xini,xend,yini,yend;
  ImageFeature acu(xdim,ydim,1);  
  ImageFeature divisor(xdim,ydim,1);
  
  for(i=0; i < int(ydim-wndsize); i+=offset) {
    yini = i;
    yend = i+wndsize;
    for(j=0; j < int(xdim-wndsize); j+=offset) {
      xini = j;
      xend = j+wndsize;
      eq_histo(xini,xend,yini,yend,image,c,acu,divisor);
    }
  }
  
  for(j=0; j < ydim; j++) {
    for(i=0; i < xdim; i++) {
      if ( (image(i,j,c) > bnd) && (divisor(i,j,0) > 0.0) ) {
        image(i,j,c)=fact_eq*(acu(i,j,0)/divisor(i,j,0)) + fact_orig*image(i,j,c);  // weighted average between original and equalized image
      }
    }
  }
}


ImageFeature zeropad(const ImageFeature& image, uint xpad, uint ypad) {
  ImageFeature result(image.xsize()+2*xpad, image.ysize()+2*ypad, image.zsize());
  
  for(uint x=xpad;x<image.xsize()+xpad;++x) {
    for(uint y=ypad;y<image.ysize()+ypad;++y) {
      for(uint c=0;c<image.zsize();++c) {
        result(x,y,c)=image(x-xpad,y-ypad,c);
      }
    }
  }
  return result;
}


ImageFeature affineTransformation(const ImageFeature &image, const AffineTransformation& at) {

  InterpolatingImage interpolated(image);

  ImageFeature result(image.xsize(), image.ysize(), image.zsize());
  
  float x,y;
  double value;
  double inverted[2][2];
  
  inverted[0][0]=at.A[1][1]/(at.A[0][0]*at.A[1][1]-at.A[0][1]*at.A[1][0]);
  inverted[0][1]=at.A[0][1]/(at.A[0][0]*at.A[1][1]-at.A[0][1]*at.A[1][0]);
  inverted[1][0]=at.A[1][0]/(at.A[0][0]*at.A[1][1]-at.A[0][1]*at.A[1][0]);
  inverted[1][1]=at.A[0][0]/(at.A[0][0]*at.A[1][1]-at.A[0][1]*at.A[1][0]);
  
  DBG(20) << "The inverted matrix is:" << endl
          << "A=/ " << inverted[0][0] << " " << inverted[0][1] << " \\ "<< endl
          << "  \\ "<< inverted[1][0] << " " << inverted[1][1] << " /  "<< endl;
  
  
  for(uint i=0;i<image.ysize();++i) {
    for(uint j=0;j<image.xsize();++j) {
      y=inverted[0][0]*(i-at.b[1])+inverted[0][1]*(j-at.b[0]);
      x=inverted[1][0]*(i-at.b[1])+inverted[1][1]*(j-at.b[0]);
      
      for(uint z=0;z<image.zsize();++z) {
        value=interpolated(x,y,z);
        if(value > 1.0) value=1.0;
        if(value < 0.0) value=0.0;
        result(j,i,z)=value;
      }
    }
  }
  
  return result;
}

ImageFeature localvariance(const ImageFeature &image, const uint winsize) {
  ImageFeature result(image.xsize(),image.ysize(),1);

  for(uint x=winsize;x<image.xsize()-winsize;++x) {
    for(uint y=winsize;y<image.ysize()-winsize;++y) {
      result(x,y,0)=localvariance(image, winsize, x,y);
    }
  }
  return result;
}

double localvariance(const ImageFeature &image, const uint winsize, const uint xpos, const uint ypos) {

  double mean=0.0, var=0.0; uint counter=0;
  for(uint x=xpos-winsize;x<=xpos+winsize;++x) {
    for(uint y=ypos-winsize;y<=ypos+winsize;++y) {
      for(uint z=0;z<image.zsize();++z) {
        mean+=image(x,y,z);
        var+=image(x,y,z)*image(x,y,z);
        ++counter;
      }
    }
  }

  mean/=double(counter);  var/=double(counter);
  
  var-=mean*mean;
  return var;
}

ImageFeature getPatch(const ImageFeature &image, const uint xpos, const uint ypos, const uint winsize) {
  uint dim=2*winsize+1;

  ImageFeature result(dim,dim,image.zsize());
  
  for(int x=(int)xpos-(int)winsize;x<=(int)xpos+(int)winsize;++x) {
    for(int y=(int)ypos-(int)winsize;y<=(int)ypos+(int)winsize;++y) {
      for(uint z=0;z<image.zsize();++z) {
        if(inImage(image,x,y,z)) {
          result(x+winsize-xpos,y+winsize-ypos,z)=image(x,y,z);
        } else {
          result(x+winsize-xpos,y+winsize-ypos,z)=0.0;
        }
      }
    }
  }
  return result;
}

vector<double> toVector(const ImageFeature &image) {
  vector<double> result(image.size());
  for(uint i=0;i<image.size();++i) {
    result[i]=image[i];
  }
  return result;
}



ImageFeature localEntropy(const ImageFeature &image, const uint winsize) {
  ImageFeature result(image.xsize(), image.ysize(), image.zsize());
  
  for(uint x=0;x<image.xsize();++x) {
    for(uint y=0;y<image.ysize();++y) {
      for(uint z=0; z<image.zsize();++z) {
        result(x,y,z)=localEntropy(image, x,y,z,winsize);
      }
    }
  }
  
  return result;
}

double localEntropy(const ImageFeature &image, const uint xpos, const uint ypos, const uint zpos, const uint winsize) {
  double result=0.0;
  typedef map<const double, uint> MapType;
  MapType el;

  double count=(2*winsize+1); count*=count;
  double val;  
  for(int x=int(xpos)-int(winsize);x<=int(xpos+winsize);++x) {
    for(int y=int(ypos)-int(winsize);y<=int(ypos+winsize);++y) {
      if(x>=0 && x<int(image.xsize()) && y>=0 && y<int(image.ysize())) val=image(x,y,zpos);
      else val=0.0;
      el[val]++;
    }
  }

  for(MapType::const_iterator i=el.begin();i!=el.end();++i) {
    double p=double((*i).second)/count;
    result+= p*log2(p);
  }
  return -result;
}

ImageFeature threshold(const ImageFeature &image, const double threshold) {
  ImageFeature result(image.xsize(), image.ysize(), image.zsize());
  for(uint x=0;x<image.xsize();++x) {
    for(uint y=0;y<image.ysize();++y) {
      for(uint z=0;z<image.zsize();++z) {
        if(image(x,y,z)>=threshold) {
          result(x,y,z)=1.0;
        } else {
          result(x,y,z)=0.0;
        }
      }
    }
  }
  return result;
}

ImageFeature getPatch(const ImageFeature &image, const uint left, const uint top, const uint right, const uint bottom) {
  DBG(50) << "get image patch from (" << left << "," << top << ") to (" << right << "," << bottom << ") " << endl;
  ImageFeature result(right-left+1, bottom-top+1,image.zsize());
  
  for(uint x=left;x<=right;++x) {
    for(uint y=top;y<=bottom;++y) {
      for(uint z=0;z<image.zsize();++z) {
        if(x<image.xsize() && y<image.ysize() && x>0 && y>0) {
          result(x-left,y-top,z)=image(x,y,z);
        }
      }
    }
  }
  return result;
}


void setPatch(ImageFeature& image, const uint xpos, const uint ypos, const ImageFeature &patch) {
  DBG(50) << "set image patch from (" << xpos << "," << ypos << ") to (" << xpos+patch.xsize() << "," << ypos+patch.ysize() << ") " << endl;
  for(uint x=0;x<patch.xsize();++x) {
    for(uint y=0;y<patch.ysize();++y) {
      for(uint z=0;z<patch.zsize();++z) {
        if(x+xpos<image.xsize() && y+ypos<image.ysize()) {
          image(x+xpos,y+ypos,z)=patch(x,y,z);
        }
      }
    }
  }
}


/**
 * fftshift operation as in octave
 *
 *  1 1 2 2         3 3 4 4 
 *  1 1 2 2         3 3 4 4
 *  4 4 3 3   -->   2 2 1 1
 *  4 4 3 3         2 2 1 1
 * 
 * 
 *  1 1 1 2 2       4 4 3 3 3 
 *  1 1 1 2 2       4 4 3 3 3 
 *  1 1 1 2 2 -->   2 2 1 1 1 
 *  4 4 4 3 3       2 2 1 1 1 
 *  4 4 4 3 3       2 2 1 1 1 
 */
void fftshift(ImageFeature &img){
  ImageFeature tmp(img);

  uint w=img.xsize();
  uint h=img.ysize();

  uint h2=(h+1)/2;
  uint w2=(w+1)/2;
  
  uint dh2=h/2;
  uint dw2=w/2;

  //1 quadrant -> 3
  for(uint x=0;x<w2;++x) {
    for(uint y=0;y<h2;++y) {
      img(dw2+x,dh2+y,0)=tmp(x,y,0);
    }
  }
  
  // 2 quandrant -> 4
  for(uint x=0;x<dw2;++x) {
    for(uint y=0;y<h2;++y) {
      img(x,dh2+y,0)=tmp(w2+x,y,0);
    }
  }
  
  // 3 quadrant -> 1
  for(uint x=0;x<dw2;++x) {
    for(uint y=0;y<dh2;++y) {
      img(x,y,0)=tmp(x+w2,y+h2,0);
    }
  }
  
  //4 quadrant -> 2
  for(uint x=0;x<w2;++x) {
    for(uint y=0;y<dh2;++y) {
      img(x+dw2,y,0)=tmp(x,y+h2,0);
    }
  }
}

ImageFeature fft(const ImageFeature &img) {
  ImageFeature result(img.xsize(),img.ysize(),2);
#ifdef HAVE_FFT_LIBRARY
  fftw_complex FIMG[img.xsize()][img.ysize()];

  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      FIMG[x][y].re=img(x,y,0);
      FIMG[x][y].im=0;
    }
  }
  
  fftwnd_plan plan = fftw2d_create_plan(img.xsize(),img.ysize(), FFTW_FORWARD, FFTW_ESTIMATE  | FFTW_IN_PLACE); 
  fftwnd_one(plan,&FIMG[0][0],NULL);
  fftwnd_destroy_plan(plan);
  
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      result(x,y,0)=FIMG[x][y].re;
      result(x,y,1)=FIMG[x][y].im;
    }
  }
#else
#warning "FFT is used without FFT library enabled."
  ERR << "FFT not available, thus returning empty image instead of fourier transformed version. All that comes now is probably bogus." << endl;
#endif
  return result;
}

ImageFeature fftscale(const ImageFeature& image, const uint newWidth, const uint newHeight) {
#ifdef HAVE_FFT_LIBRARY

  uint imgheight=image.ysize();
  uint imgwidth=image.xsize();

  // variables for image and filter in transformed
  fftw_complex *FIMG=new fftw_complex[image.xsize()*image.ysize()];
  fftw_complex *FRESULT=new fftw_complex[newWidth*newHeight];
  for(uint i=0;i<newWidth*newHeight;++i) { FRESULT[i].re=0.0; FRESULT[i].im=0.0;}
  
  //copy image to temp variable
  for(uint x=0;x<image.xsize();++x) {
    for(uint y=0;y<image.ysize();++y) {
      FIMG[y*imgwidth+x].re=image(x,y,0);
      FIMG[y*imgwidth+x].im=0.0;
    }
  }
  
  //fourier transform
  fftwnd_plan plan = fftw2d_create_plan(imgheight,imgwidth, FFTW_FORWARD, FFTW_ESTIMATE  | FFTW_IN_PLACE); 
  fftwnd_one(plan,FIMG,NULL);
  fftwnd_destroy_plan(plan);

  // copy into destination with zeropadding
  
  uint w2=imgwidth/2;
  uint h2=imgheight/2;
  
  uint gw2=(imgwidth+1)/2;
  uint gh2=(imgheight+1)/2;

  uint nh2=newHeight/2;
  uint nw2=newWidth/2;
  uint gnh2=(newHeight+1)/2;
  uint gnw2=(newWidth+1)/2;
  
  
  if(newHeight>imgheight) {
    nh2=h2; gnh2=gh2;
  }
  if(newWidth>imgwidth) {
    nw2=w2; gnw2=gw2;
  }

  //oben links
  for(uint y=0;y<gh2 and y<gnh2;++y) {
    for(uint x=0;x<gw2 and x<gnw2;++x) {
      FRESULT[y*newWidth+x]=FIMG[y*imgwidth+x];
    }
  }
  

  //oben rechts
  for(uint y=0;y<gh2 and y<gnh2;++y) {
    for(uint x=0;x<w2 and x<nw2;++x) {
      FRESULT[y*newWidth+newWidth-nw2+x]=FIMG[y*imgwidth+imgwidth-nw2+x];
    }
  }  

  //unten rechts
  for(uint y=0;y<h2 and y<nh2;++y) {
    for(uint x=0;x<w2 and x<nw2;++x) {
      //FRESULT[newWidth-nw2+x][newHeight-nh2+y]=FIMG[image.xsize()-nw2+x][image.ysize()-nh2+y];
      FRESULT[(newHeight-nh2+y)*newWidth+newWidth-nw2+x]=FIMG[(image.ysize()-nh2+y)*imgwidth+image.xsize()-nw2+x];
    }
  }

  //unten links
  for(uint y=0;y<h2 and y<nh2;++y) {
    for(uint x=0;x<gw2 and x<gnw2;++x) {
      //FRESULT[x][newHeight-nh2+y]=FIMG[x][image.ysize()-nh2+y];
      FRESULT[(newHeight-nh2+y)*newWidth+x]=FIMG[(image.ysize()-nh2+y)*imgwidth+x];
    } 
  }

  
  //fourier transform backwards
  fftwnd_plan planb = fftw2d_create_plan(newHeight,newWidth, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE); 
  fftwnd_one(planb,FRESULT,NULL);
  fftwnd_destroy_plan(planb);
  

  int divider=(imgwidth*imgheight);
  
  ImageFeature res(newWidth, newHeight,1);
  for(uint x=0;x<newWidth;++x) {
    for(uint y=0;y<newHeight;++y) {
      res(x,y,0)=sqrt((FRESULT[y*newWidth+x].re*FRESULT[y*newWidth+x].re)+(FRESULT[y*newWidth+x].im*FRESULT[y*newWidth+x].im))/divider;
    }
  }
  
  cutoff(res);
  return res;

#else 
#warning "FFTscale is used without fftw enabled. Replacing by normal (Bresenham) scaling."
  ERR << "FFT lib is not available, thus using normal (Bresenham) scaling." << endl;
  return scale(image,  newWidth,  newHeight);
#endif
}


// simple Bresenham scale
// works well if image gets smaller, otherwise interpolation / affine transform gives better results
ImageFeature scale(const ImageFeature& image, const uint newWidth, const uint newHeight) {
  uint oldHeight=image.ysize();
  uint oldWidth=image.xsize();
  
  ImageFeature scaled(newWidth, newHeight, image.zsize());
  
  double x_weight =  double(newWidth)  /double(oldWidth) ;
  double y_weight =  double(newHeight)/double(oldHeight);
  
  for(uint z=0;z<image.zsize();++z) {
    ImageFeature work(newWidth,oldHeight,1);
    uint  index_orig, index_scaled;
    double akku;
  
    for(uint y=0; y< oldHeight ; y++) {
      // scale line y from oldWidth to newWidth
      index_orig=0;
      index_scaled=0;
      akku = 0.0;

      while(index_orig < oldWidth ) {
        
        if( (akku+x_weight) >= 1.0) { 
          // target pixel left
        
          work(index_scaled,y,0) += image(index_orig,y,z) * (1.0 - akku) ;
          akku += x_weight - 1.0;
          index_scaled++;
          
          //target pixel centered (if enlarging)
          while( akku >= 1.0 ) {
            work(index_scaled,y,0) = image(index_orig,y,z);
            index_scaled++;
            akku -= 1.0;
          }

          if( index_scaled < newWidth ) {
            work(index_scaled,y,0)=image(index_orig,y,z) * (akku);
          }
        } else {
          // border to new target pixel not crossed
          work(index_scaled,y,0) += image(index_orig,y,z) * x_weight;
          akku += x_weight;
        }
        index_orig++;
      } // while
    } // for(y...
    
    // now scale in y-direction
    // -> column from oldHeight -> newHeight

    for(uint x=0; x < newWidth; x++) {
      index_orig = 0;
      index_scaled = 0;
      akku = 0.0;
      scaled(x,0,z)=0.0;
      
      while( index_orig < oldHeight ) {
        
        if ( (akku + y_weight) >= 1.0 ) {  //going to new pixel?
          //target pixel at the top (ready)
          
          scaled(x,index_scaled,z) += (1.0 - akku ) * work(x,index_orig,0);
          
          akku += y_weight - 1.0;
          index_scaled ++;
          
          //target pixel vertically in the middle (only for enlargment)
          while ( akku >= 1.0 ) {
            scaled(x,index_scaled,z) = work(x,index_orig,0);
            
            index_scaled ++;
            akku -= 1.0;
          }
          
          // target pixel on bottom (new started)
          if( index_scaled < newHeight ) {
            scaled(x,index_scaled,z) = work(x,index_orig,0) * akku;
          } //if
          
        } // if
        else {
          // no border crossing to new target pixel
          scaled(x,index_scaled,z) += work(x,index_orig,0) * y_weight;
          akku += y_weight;
        } // else
        index_orig ++;
      } // while
    } // for(x
  } // for z
  return scaled;
} //scale
 

bool inImage(const ImageFeature &img, const int x, const int y, const int z) {
  return (x<int(img.xsize()) && y<int(img.ysize()) && y>=0 && x>=0 && z>=0 && z<int(img.zsize()));
}

void cross(ImageFeature &img, const uint x, const uint y,const vector<double>& color,const uint len) {
  //mark center
  if (inImage(img,x,y)) {for(uint c=0;c<img.zsize();++c) {img(x,y,c)=color[c];}}
  //mark outer pars
  for(uint i=1;i<=len;++i) {
    if(inImage(img,x+i,y)) {for(uint c=0;c<img.zsize();++c) {img(x+i,y,c)=color[c];}}
    if(inImage(img,x-i,y)) {for(uint c=0;c<img.zsize();++c) {img(x-i,y,c)=color[c];}}
    if(inImage(img,x,y+i)) {for(uint c=0;c<img.zsize();++c) {img(x,y+i,c)=color[c];}}
    if(inImage(img,x,y-i)) {for(uint c=0;c<img.zsize();++c) {img(x,y-i,c)=color[c];}}
  }
}

void diagcross(ImageFeature &img, const uint x, const uint y,const vector<double>& color,const uint len) {
  //mark center
  if (inImage(img,x,y)) {for(uint c=0;c<img.zsize();++c) {img(x,y,c)=color[c];}}
  //mark outer pars
  for(uint i=1;i<=len;++i) {
    if(inImage(img,x+i,y+i)) {for(uint c=0;c<img.zsize();++c) {img(x+i,y+i,c)=color[c];}}
    if(inImage(img,x-i,y-i)) {for(uint c=0;c<img.zsize();++c) {img(x-i,y-i,c)=color[c];}}
    if(inImage(img,x-i,y+i)) {for(uint c=0;c<img.zsize();++c) {img(x-i,y+i,c)=color[c];}}
    if(inImage(img,x+i,y-i)) {for(uint c=0;c<img.zsize();++c) {img(x+i,y-i,c)=color[c];}}
  }
}

void box(ImageFeature &img, const uint centerx, const uint centery, const vector<double> &color, const uint winsize) {
  for(uint x=centerx-winsize;x<=centerx+winsize;++x) {
    if(inImage(img,x,centery+winsize)) {for(uint c=0;c<img.zsize();++c) {img(x,centery+winsize,c)=color[c];}}
    if(inImage(img,x,centery-winsize)) {for(uint c=0;c<img.zsize();++c) {img(x,centery-winsize,c)=color[c];}}
  }
  
  for(uint y=centery-winsize;y<=centery+winsize;++y) {
    if(inImage(img,centerx+winsize,y)) {for(uint c=0;c<img.zsize();++c) {img(centerx+winsize,y,c)=color[c];}}
    if(inImage(img,centerx-winsize,y)) {for(uint c=0;c<img.zsize();++c) {img(centerx-winsize,y,c)=color[c];}}
  }
}

void meanandvariance(const ImageFeature &img, double &mean, double &variance, const uint layer) {
  mean=0.0; variance=0.0;
  for(uint y=0;y<img.ysize();++y) {
    for(uint x=0;x<img.xsize();++x) {
      mean+=img(x,y,layer);
      variance+=img(x,y,layer)*img(x,y,layer);
    }
  }

  uint size=img.xsize()*img.ysize();

  mean/=double(size);
  variance/=double(size);
  variance-=mean*mean;
}


void meanAndVarianceNormalization(ImageFeature &img, const double mu, const double sigma) {
  double mean, variance,stddev;
  for(uint c=0;c<img.zsize();++c) {
    meanandvariance(img,mean,variance,c);
    mean+=mu;
    stddev=sqrt(variance);
    stddev/=sqrt(sigma);
    for(uint x=0;x<img.xsize();++x) {
      for(uint y=0;y<img.ysize();++y) {
        img(x,y,c)=(img(x,y,c)-mean)/stddev;
      }
    }
  }
}


void cutoff(ImageFeature &img, const double minimum, const double maximum) {
  for(uint c=0;c<img.zsize();++c) {
    for(uint x=0;x<img.xsize();++x) {
      for(uint y=0;y<img.ysize();++y) {
        img(x,y,c)=max(minimum,img(x,y,c));
        img(x,y,c)=min(maximum,img(x,y,c));
      }
    }
  }
}


void shift(ImageFeature &img, const double offset) {
  for(uint c=0;c<img.zsize();++c) {
    for(uint x=0;x<img.xsize();++x) {
      for(uint y=0;y<img.ysize();++y) {
        img(x,y,c)+=offset;
      }}}
}


void histogramNormalization(ImageFeature &img) {
  const uint bins=256;
  vector<int> H(bins);
  vector<double> T(bins);

  for(uint c=0;c<img.zsize();++c) {
    //DBG(10) << "c=" << c << endl;
    
    for(uint x=0;x<img.xsize();++x) {
      for(uint y=0;y<img.ysize();++y) {
        ++H[int(img(x,y,c)*255)];
      }
    }
    
    for(uint p=1;p<bins;++p) {
      H[p]+=H[p-1];
    }
    
    for(uint p=0;p<bins;++p) {
      T[p]=(1.0/double(img.size()))*double(H[p]);
    }
    
    for(uint x=0;x<img.xsize();++x) {
      for(uint y=0;y<img.ysize();++y) {
        img(x,y,c)=T[int(img(x,y,c)*255)];
      }
    }
    
  }
}
void appendVertically(ImageFeature &imga, const ImageFeature &imgb) {
  if(imga.xsize()!=imgb.xsize() || imga.zsize() != imgb.zsize()) {
    ERR << "Images have to have same width and depth" << endl;
  } else {
    ImageFeature tmp(imga);
    imga.resize(imga.xsize(),imga.ysize()+imgb.ysize(),imga.zsize());

    for(uint x=0;x<tmp.xsize();++x) {
      for(uint z=0;z<tmp.zsize();++z) {
        for(uint y=0;y<tmp.ysize();++y) {
          imga(x,y,z)=tmp(x,y,z);
        }
        for(uint y=tmp.ysize();y<tmp.ysize()+imgb.ysize();++y) {
          imga(x,y,z)=imgb(x,y-tmp.ysize(),z);
        }
      }
    }
  }
}

void appendHorizontally(ImageFeature &imga, const ImageFeature &imgb) {
  if(imga.ysize()!=imgb.ysize() || imga.zsize() != imgb.zsize()) {
    ERR << "Images have to have same width and depth" << endl;
  } else {
    ImageFeature tmp(imga);
    imga.resize(imga.xsize()+imgb.xsize(),imga.ysize(),imga.zsize());
    
    for(uint y=0;y<tmp.ysize();++y) {
      for(uint z=0;z<tmp.zsize();++z) {      
        for(uint x=0;x<tmp.xsize();++x) {
          imga(x,y,z)=tmp(x,y,z);
        }
        for(uint x=tmp.xsize();x<tmp.xsize()+imgb.xsize();++x) {
          imga(x,y,z)=imgb(x-tmp.xsize(),y,z);
        }
      }
    }
  }
}


void shave(ImageFeature &imga, const uint shaver) {
  ImageFeature tmp(imga);
  imga.resize(imga.xsize()-2*shaver, imga.ysize()-2*shaver, imga.zsize());
  for(uint x=0;x<tmp.xsize();++x) {
    if(x>=shaver and x<tmp.xsize()-shaver) {
      for(uint y=0;y<tmp.ysize();++y) {
        if(y>=shaver and y<tmp.ysize()-shaver) {
          for(uint z=0;z<tmp.zsize();++z) {
            imga(x-shaver,y-shaver,z)=tmp(x,y,z);
          }
        }
      }
    }
  }
}

void rect(ImageFeature &img, const uint x, const uint y, const uint w, const uint h, const ::std::vector<double> &color, bool dash) {
  uint dashCounterMax = 3;
  uint dashCounter = dashCounterMax;
  uint dashMode = 1;
  // upper line
  for(uint i = x; i < x + w; ++i) { 
    dashCounter--;
    if (dashCounter == 0) {
      dashCounter = dashCounterMax;
      dashMode = 1 - dashMode;
    }
    for(uint c = 0; c < img.zsize(); ++c) {
      if (!dash || (dashMode == 1)) {
	img(i, y, c) = color[c];
      }
    }
  }
  // right line
  for(uint i = y; i < y + h; ++i) { 
    dashCounter--;
    if (dashCounter == 0) {
      dashCounter = dashCounterMax;
      dashMode = 1 - dashMode;
    }
    for(uint c = 0; c < img.zsize(); ++c) {
      if (!dash || (dashMode == 1)) {
	img(x + w - 1, i, c) = color[c];
      }
    }
  }
  // lower line
  for(int i = x + w - 1; i >= (int) x; --i) { 
    dashCounter--;
    if (dashCounter == 0) {
      dashCounter = dashCounterMax;
      dashMode = 1 - dashMode;
    }
    for(uint c = 0; c < img.zsize(); ++c) {
      if (!dash || (dashMode == 1)) {
	img(i, y + h - 1, c)=color[c];
      }
    }
  }
  // left line
  for(int i = y + h - 1; i >= (int) y; --i) { 
    dashCounter--;
    if (dashCounter == 0) {
      dashCounter = dashCounterMax;
      dashMode = 1 - dashMode;
    }
    for(uint c = 0; c < img.zsize(); ++c) {
      if (!dash || (dashMode == 1)) {
	img(x, i, c) = color[c];
      }
    }
  }
}

void RGBtoHSV(ImageFeature &img) {
  HSVPixel hsv;
  for(uint x=0;x<img.xsize();++x) {
    for(uint y=0;y<img.ysize();++y) {
      hsv=img.hsvPixel(x,y);
      img(x,y,0)=hsv.h;
      img(x,y,1)=hsv.s;
      img(x,y,2)=hsv.v;
    }
  }
}

void HSVtoRGB(ImageFeature &) {
#warning not implemented: HSVtoRGB 
}


// helper for DCTs
inline double alpha(uint u, uint N) {
  if (u == 0)
    return sqrt(1/double(N));
  else if ((u > 0) && (u < N))
    return sqrt(2/double(N));
  else
    return 0.0;
  
  // Silence strange "control may reach end of non-void function" warning
  return 0.0;
}


void DCT(ImageFeature &src, ImageFeature &dest) {
  dest.resize(src.xsize(),src.ysize(),src.zsize()); //NB Not efficient
  
  for(uint c=0;c<src.zsize();++c) {

    uint i,j,k;
    double sum;
    // Transform in x direction
    ImageFeature horizontal(src.xsize(),src.ysize(),1);
    for (j = 0; j < src.ysize(); j++) {
      for (i = 0; i < src.xsize(); i++) {
        sum = 0.0;
        for (k = 0; k < src.xsize(); k++) {
          sum += src(k,j,c) * cos(double(2*k+1)*M_PI*double(i)/(double(2*src.xsize())));
        }
        horizontal(i,j,0) = sum;
      }
    }
    // Transform in y direction
    for (i = 0; i < src.xsize(); i++) {
      for (j = 0; j < src.ysize(); j++) {
        sum = 0.0;
        for (k = 0; k < src.ysize(); k++) {
          sum += horizontal(i,k,0) * cos(double(2*k+1)*M_PI*double(j)/(double(2*src.ysize()))); 
        }
        dest(i,j,c) = alpha(i,src.xsize())*alpha(j,src.ysize())*sum; 
      }
    }
  }
}


void iDCT(ImageFeature &src, ImageFeature &dest) {
  dest.resize(src.xsize(),src.ysize(),src.zsize()); //NB Not efficient
  
  for(uint c=0;c<src.zsize();++c) {
  
    uint i,j,k;
    double sum;
    // Transform in x direction
    ImageFeature horizontal(src.xsize(), src.ysize(),1);
    for (j = 0; j < src.ysize(); j++) {
      for (i = 0; i < src.xsize(); i++) {
        sum = 0.0;
        for (k = 0; k < src.xsize(); k++) {
          sum += alpha(k,src.xsize())*src(k,j,0)* cos(double(2*i+1)*M_PI*double(k)/(double(2*src.xsize())));  
        }
        horizontal(i,j,0) = sum;
      }
    }

    // Transform in y direction
    for (i = 0; i < src.xsize(); i++) {
      for (j = 0; j < src.ysize(); j++) {
        sum = 0.0;
        for (k = 0; k < src.ysize(); k++) {
          sum += alpha(k, src.xsize())*horizontal(i,k,0)* cos(double(2*j+1)*M_PI*double(k)/(double(2*src.ysize()))); 
        }
        dest(i,j,c) = sum;
      }
    }
  }
}
