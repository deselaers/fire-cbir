// THIS SOFTWARE IS BASED ON THE PAPER
// "Wavelet-based Salient Points for Image Retrieval"
// by E.Loupias, N.Sebe (Leiden Institute of Advanced Computer Science, Leiden University, The Netherlands, nicu@wi.leidenuniv.nl
// and S.Bres, J.-M.Jolion (Laboratoire Reconnaisannce de Formes et Vision, INSA Lyon, {loupias,jolion}@rfv.insa-lyon.fr
// AND SOFTWARE WRITTEN BY E.LOPIAS

#ifndef __SALIENTPOINTS_HPP__
#define __SALIENTPOINTS_HPP__

#include <vector>
#include <string>
#include "imagefeature.hpp"
#include "point.hpp"

::std::vector< ::std::vector<float> > array1Dto2D( ::std::vector<float>& in, uint height, uint width );

class SalientPoints  {
private:
  // members declaration
  ::std::vector<int> Pixels;
  
  int nb_niveaux, nivX, nivY, BORD;
  int limX[32], limY[32];
  bool pixelsSet, imageSet;
  
  // methods declaration
  int grad( int x, int y );
  void getMultires();
  Point getNiveaux(int l, int c);
  
  ::std::vector<Point> dernierCoeff2pixels(const ::std::string& wavelet, int la, int ca, int niv_l, int niv_c);
  
  void coordEnfant(int la, int ca, int niv_l, int niv_c, Point& coordEnfAbs);

  ::std::vector<Point> getBlocEnfant(const ::std::string& wavelet, int la, int ca, int niv_l, int niv_c);
  Point followPt(const ::std::string& wavelet, int la_parent, int ca_parent, float& som_multires,  ::std::vector< ::std::vector<int> >& gradient);
  ::std::vector<Point> sPoints(const ::std::string & wavelet,::std::vector<float>& floatPixels, int bord);

public:

  ::std::vector<float> floatPixels;
  int width, height;
  bool floatSet;

  SalientPoints(){}
  SalientPoints(const ImageFeature& img);
  ~SalientPoints() {}
  ::std::vector<Point> getSalientPoints(int nrPoints);
  ::std::vector<Point> getSalientPoints();
  
  int getWidth() { return width; }
  int getHeight() { return height; }
  ::std::vector<int>& getPixels();
};

#endif
