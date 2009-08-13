#ifndef _SVD_HPP_
#define _SVD_HPP_

#ifdef HAVE_LAPACK
#include <vector>

typedef int u16;
typedef ::std::vector< ::std::vector<double> > Matrix;
typedef ::std::vector<double> Vector;

#include "Lapack.hh"

const double svdThreshold=1.0e-12;
const long SMLSIZ=100;

void svd(Matrix &u_, Vector &w_,Matrix& v_);
#else
#warning "compiling SVD without lapack library. This won't work"
#endif

#endif
