// THIS SOFTWARE IS BASED ON THE PAPER
// "Wavelet-based Salient Points for Image Retrieval"
// by E.Loupias, N.Sebe (Leiden Institute of Advanced Computer Science, Leiden University, The Netherlands, nicu@wi.leidenuniv.nl
// and S.Bres, J.-M.Jolion (Laboratoire Reconnaisannce de Formes et Vision, INSA Lyon, {loupias,jolion}@rfv.insa-lyon.fr
// AND SOFTWARE WRITTEN BY E.LOPIAS

#ifndef __WAVELET_HPP__
#define __WAVELET_HPP__

#include "salientpoints.hpp"
#include <vector>
#include <string>

::std::vector<int> linearisation(int dim,  ::std::vector<int>& in);
::std::vector<int> recadrage(int dim,  ::std::vector<float>& in, int min, int max);
::std::vector<int> recadrage(int dim,  ::std::vector<int>& in, int min, int max);

class Wavelet {

private:

  ::std::string wavelet;
  bool std;

  ::std::vector<float> transformee(const ::std::string& wave, bool std, ::std::vector<float>& floatPixels, int DimX, int DimY, int isign);
  ::std::vector<float> daub4NS(::std::vector <float>& FloatImage, int DimX, int DimY, int isign );
  void daub4_step(::std::vector<float>& data, int nn, int isign);
  ::std::vector<float> daub4Std(::std::vector<float>& FloatImage, int DimX, int DimY, int isign);
  void Daub4_1D(::std::vector<float> &data, int nn, int isign);

public:
  Wavelet(const ::std::string &m_Wavelet, bool std_, int isign, SalientPoints& m_Exemple, ::std::vector<int>& transfoAffichable);
  ~Wavelet() {  }
  ::std::vector<float> transformee(const ::std::string& wave, bool std, SalientPoints& m_Exemple, int isign);
};


#endif
