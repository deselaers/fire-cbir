#include <string>
#include "diag.hpp"
#include "imagefeature.hpp"
#include "histogramfeature.hpp"


ImageFeature calculate(const ImageFeature &img);

ImageFeature contrast(const ImageFeature &image, uint z);
ImageFeature directionality(const ImageFeature &image, uint z);
ImageFeature coarseness(const ImageFeature &image, uint z);

HistogramFeature histogramize(const ImageFeature &image);

double getLocalContrast(const ImageFeature &image, int xpos, int ypos, uint z);
double efficientLocalMean(const uint y,const uint x,const uint k,const ImageFeature &laufendeSumme);

::std::vector<uint> histogramImageBinImage(const ImageFeature &image);
HistogramFeature histogramize(const ImageFeature image, const ::std::vector<uint> &bins, const uint left, const uint top, const uint right, const uint bottom);

