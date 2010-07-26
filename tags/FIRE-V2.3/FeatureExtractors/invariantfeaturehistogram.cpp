#include <string>
#include <sstream>
#include "invariantfeaturehistogram.hpp"
#include "kernelfunction.hpp"
#include "point.hpp"
#include "kernelfunctionmaker.hpp"

using namespace std;


// calculates the histogram  
// by using the pixels of the given positions
void getInvariantFeatureHistogramByList(HistogramFeature &histo,
                                        const InterpolatingImage &iimg, 
                                        KernelFunction *f,  
                                        const vector<Point> &points,
                                        const int R) {
																			
  // iterating over the list of positions
  for (uint i = 0; i < points.size(); ++i) {
    //cout << "sample no " << i << endl;
    vector<double> pixel(iimg.zsize());
    DBG(45) << "no of layer " << iimg.zsize() << endl;
    // obtaining value for each layer (color)
    for (uint z = 0; z < iimg.zsize(); ++z) {
      DBG(45) << "layer no " << z << endl;
      pixel[z] = .0;
      for (int r = 0; r < R; ++r) {
        pixel[z] += f->calc(iimg, points.at(i).x, points.at(i).y, r, R, z);
      }
      pixel[z] /= (double) R;	
    }
    DBG(45) << "pixel = " << pixel[0] << " " << pixel[1] << " " << pixel[2] << " " << endl;
    // feeding resulting pixel into the histogram
    histo.feed(pixel);
  }
}


HistogramFeature getInvariantFeatureHistogramPart(const ImageFeature &image, 
                                              const ::std::string &kernelFunction, 
                                              const int samples,
                                              const int steps,
                                              const int R, 
                                              int startx=-1,
                                              int stopx=-1,
                                              int starty=-1,
                                              int stopy=-1) {
                                              
  DBG(30) << "startx: " << startx << " stopx: " << stopx << " starty: " << starty << " stopy: " << stopy << endl;                                              
                                              
  // create the appropriate kernel function
  KernelFunction *f=getNewKernelFunction(kernelFunction);
  // obtaining image with interpolation capabilities
  InterpolatingImage iimg(image); 
         
  // Image that contains the pixels that are fed into the histogram  
  // (only used when high debug-level)                                                
  ImageFeature histoImage(image.xsize(), image.ysize(), image.zsize());    
	  
  // initializing histogram
  vector<uint> stepsVec(image.zsize(),steps); 
  HistogramFeature result(stepsVec);  
  for(uint i=0;i<image.zsize();++i) {
    result.min()[i]=0.0; 
    result.max()[i]=1.0;
  }
  result.initStepsize();

  // check if the given parameters are within range  
  if ((startx >= stopx) || (starty >= stopy)) {
    DBG(30) << "start value greater than stop value" << endl;
    return (result);
  }
  if ((startx < 0) || (stopx > (int) image.xsize()) || (starty < 0) || (stopy > (int) image.ysize())) {
    DBG(30) << "start or stop value out of range" << endl;
    return (result);
  }
  
  // check whether we want to use monte carlo method or not
  if(samples < 0) {
    // not using monte carlo, iterate over all pixels (within the start/stop-range):
    vector<double> pixel(image.zsize());
    for (uint x = startx; x < (uint) stopx; ++x) {
      for (uint y = starty; y < (uint) stopy; ++y) {
        DBG(30) << "x: " << x << " y: " << y << endl;
        // retrieving value for each layer (color)
        for (uint z = 0; z < image.zsize(); ++z) {
          pixel[z] = .0;
          for (int r = 0; r < R; ++r) {
            pixel[z] += f->calc(iimg, x, y, r, R, z);
          }
          pixel[z] /= (double) R;	
        }
        // checking if the pixels color is < 1.0
        DBGI(45,{for (uint z = 0; z < image.zsize(); ++z) {
                 if (pixel[z] > 1.0) {
                   DBG(25) << "Pixel value too big (whiter than white)!" << endl;
                   pixel[z] = 1.0;
                 }
               }});
        // feeding resulting pixel into the histogram
        result.feed(pixel);
        DBGI(45,{if (isnan(pixel[0]) || isnan(pixel[1]) || isnan(pixel[2])) {
                 histoImage(x,y,0) = 1;
                 histoImage(x,y,1) = 0;
                 histoImage(x,y,2) = 0;
               } else {
                   histoImage(x,y,0) = pixel[0];
                   histoImage(x,y,1) = pixel[1];
                   histoImage(x,y,2) = pixel[2];
                 }
             });
      }    
    }
     
  } else {
    // using monte carlo: use samples - many random pixels
		
    // initializing list of points
    vector<Point> points;
				
    // making a list of random positions
    srand((unsigned) time (0));
    for (int j = 0; j < samples; ++j) {
      Point point(rand()%image.xsize(),rand()%image.xsize());
      points.push_back(point);
    }
    
    // fill the histogramm by using the random positions
    getInvariantFeatureHistogramByList(result, iimg, f, points, R);
  }
  
  DBGI(45, histoImage.save("histoimage.png"));
  
  return result; 
}                                              

                                              


HistogramFeature getInvariantFeatureHistogram(const ImageFeature &image, 
                                              const string &kernelFunction, 
                                              const int samples, 
                                              const int steps,
                                              const int R) {
                                              
  return (getInvariantFeatureHistogramPart(image, kernelFunction, samples, steps, R, 0, image.xsize(), 0, image.ysize()));                                              

}
