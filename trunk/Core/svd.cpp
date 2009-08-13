#ifdef HAVE_LAPACK

#include "Lapack.hh"
#include "svd.hpp"
#include <iostream>


using namespace Lapack;

void svd(Matrix &u_, Vector &w_,Matrix& v_){
  
  char jobz = 'A';
  u16 dim_=u_.size();
  
  //TD changes here
  //   long m=u_.nRows();
  //   long n=u_.nColumns();
  long m=u_.size();
  long n=u_[0].size();
  //changes until here

  double *a= new double[m*n];
  long lda=n;
  double *s= new double[n];
  double *u = new double[n*n];
  long ldu=n;
  double *vt = new double[n*n];
  long ldvt=n;
  long lwork = ::std::max(14*::std::min(m,n)+4, 10*::std::min(m,n)+2+ SMLSIZ*(SMLSIZ+8))+ ::std::max(m,m);
  double* work= new double[lwork];
  long * iwork=new long[8*::std::min(m,n)];
  long info;	
  
  
  for(u16 rowIdx = 0; rowIdx < dim_; rowIdx++)
    for(u16 columnIdx = 0; columnIdx < dim_; columnIdx++) {
      a[columnIdx*dim_ + rowIdx] = u_[rowIdx][columnIdx];
    }
  
  
  lwork=-1;
  //::std::cout << "lwork=" << lwork << ::std::endl;
  dgesdd( &jobz, &m, &n,  a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info );

  lwork=long(work[0]);
  //::std::cout << "lwork=" << lwork << ::std::endl;
  delete[] work;
  work=new double[lwork];
  
  dgesdd( &jobz, &m, &n,  a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info );
  
  for(u16 rowIdx = 0; rowIdx < dim_; rowIdx++){
    for(u16 columnIdx = 0; columnIdx < dim_; columnIdx++) {
      u_[rowIdx][columnIdx] = u[columnIdx*dim_ + rowIdx];
      v_[rowIdx][columnIdx] = vt[columnIdx*dim_ + rowIdx];
    }	    
    w_[rowIdx]=s[rowIdx];
  }	
  delete[] a;
  delete[] s;
  delete[] u;
  delete[] vt;
  delete[] work;
  delete[] iwork;	
}

#else
#warning "compiling SVD without lapack library. This won't work"
#endif
