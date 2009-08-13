#include "Lapack.hh"

#define F77NAME(x) x##_

// Type dependend Lapack definitions:

namespace Lapack{

  extern "C" {

    void F77NAME(dgelsd)(long* m, long* n, long* nrhs, double* a, long* lda, double* b, long* ldb, 
                         double* s, double* rcond, long* rank, double* work, long* lwork, long* iwork,
                         long* info);

			 
  
    void F77NAME(dggevx)(char* balanc, char* jobvl, char* jobvr, char* sense,
                         long* n, double* a, long* lda, double* b, long* ldb,
                         double* alphar, double* alphai, double* beta, double* vl, long* ldvl, double* vr, long* ldvr,
                         long* ilo, long* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, 
                         double* rconde, double* rcondv, double* work, long* lwork, long* iwork, bool* bwork, long* info);
  
  
    void F77NAME(sggev)(char *jobvl,char *jobvr,int *n,double *a,int *lda,double *b, int *ld,
                        double *alphar,double *alphai,double *beta,double *vl,
                        int *ldvl,double *vr,int *ldvr,double *work,int *lwork,
                        int *info) ;
    void F77NAME(dggev)(char *jobvl,char *jobvr,int *n,double *a,int *lda,double *b, int *ld,
                        double *alphar,double *alphai,double *beta,double *vl,
                        int *ldvl,double *vr,int *ldvr,double *work,int *lwork,
                        int *info) ;
  
  
    void F77NAME(dsygvd)(long* itype, char* jobz, char* uplo, long* n, double* a, long* lda, double* b, long* ldb,
                         double* w, double* work, long* lwork, long* iwork, long* liwork, long* info);

    void F77NAME(dgels)(char* trans, long* m, long* n, long* nrhs, double* A, long* lda, 
                        double* B, long* ldb, double* work, long* lwork, long* info);

    void F77NAME(sgels)(char* trans, long* m, long* n, long* nrhs, float* A, long* lda, 
                        float* B, long* ldb, float* work, long* lwork, long* info);

    void F77NAME(dgelss)(long* m, long* n, long* nrhs, double* A, long* lda, double* B, long* ldb, 
                         double* s, double* rcond, long* rank, double* work, long* lwork, long* info);

    void F77NAME(sgelss)(long* m, long* n, long* nrhs, float* A, long* lda, float* B, long* ldb, 
                         float* s, float* rcond, long* rank, float* work, long* lwork, long* info);

    void F77NAME(dgebrd)(long* m, long *n, double* a, long* lda, double* d, double* e, 
                         double* tauq,double* taup, double* work, long* lwork,long* info);

    void F77NAME(dbdsdc) (char* uplo, char* compq, long* n, double* d, double* e, double* u, 
                          long* ldu, double* vt, long* ldvt, double* q, long* iq, double* work, 
                          long* iwork, long* info );

    void F77NAME(dgesdd) (char* jobz, long* m, long* n, double* a, long* lda, double* s, 
                          double* u, long *ldu, double* vt, long* ldvt, double* work, 
                          long* lwork, long* iwork, long* info );

  }

  void dgelsd(long* m, long* n, long* nrhs, double* a, long* lda, double* b, long* ldb, 
              double* s, double* rcond, long* rank, double* work, long* lwork, long* iwork, long* info){
#ifdef HAVE_LAPACK
    F77NAME(dgelsd)( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work,lwork, iwork, info);
#endif
  }
    


  void dggevx(char* balanc, char* jobvl, char* jobvr, char* sense,
              long* n, double* a, long* lda, double* b, long* ldb,
              double* alphar, double* alphai, double* beta, double* vl, long* ldvl, double* vr, long* ldvr,
              long* ilo, long* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, 
              double* rconde, double* rcondv, double* work, long* lwork, long* iwork, bool* bwork, long* info)
  {
#ifdef HAVE_LAPACK
    F77NAME(dggevx)(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb,
                    alphar, alphai, beta, vl, ldvl, vr, ldvr,
                    ilo, ihi, lscale, rscale, abnrm, bbnrm, 
                    rconde, rcondv, work, lwork, iwork, bwork, info);
#endif
  }


  void sggev(char *jobvl,char *jobvr,int *n,double *a,int *lda,double *b, int *ld,
             double *alphar,double *alphai,double *beta,double *vl,
             int *ldvl,double *vr,int *ldvr,double *work,int *lwork,
             int *info) {
#ifdef HAVE_LAPACK
    F77NAME(sggev)(jobvl,jobvr,n,a,lda,b,ld,alphar,alphai,beta,vl,ldvl,vr,ldvr,work,lwork,info);
#endif
  }
  void dggev(char *jobvl,char *jobvr,int *n,double *a,int *lda,double *b, int *ld,
             double *alphar,double *alphai,double *beta,double *vl,
             int *ldvl,double *vr,int *ldvr,double *work,int *lwork,
             int *info) {
#ifdef HAVE_LAPACK
    F77NAME(dggev)(jobvl,jobvr,n,a,lda,b,ld,alphar,alphai,beta,vl,ldvl,vr,ldvr,work,lwork,info);
#endif
  }


  void dsygvd(long* itype, char* jobz, char* uplo, long* n, double* a, long* lda, double* b, long* ldb,
              double* w, double* work, long* lwork, long* iwork, long* liwork, long* info)
  {
#ifdef HAVE_LAPACK
    F77NAME(dsygvd)(itype, jobz, uplo, n, a, lda, b, ldb,
                    w, work, lwork, iwork, liwork, info);
#endif
  }


  void gels(char* trans, long* m, long* n, long* nrhs, double* A, long* lda, 
            double* B, long* ldb, double* work, long* lwork, long* info)
  {
#ifdef HAVE_LAPACK
    F77NAME(dgels)(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
#endif
  }

  void gels(char* trans, long* m, long* n, long* nrhs, float* A, long* lda, 
            float* B, long* ldb, float* work, long* lwork, long* info)
  {
#ifdef HAVE_LAPACK
    F77NAME(sgels)(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
#endif
  }

  void gelss(long* m, long* n, long* nrhs, double* A, long* lda, double* B, long* ldb, 
             double* s, double* rcond, long* rank, double* work, long* lwork, long* info)
  {
#ifdef HAVE_LAPACK
    F77NAME(dgelss)(m, n, nrhs, A, lda, B, ldb, s, rcond, rank, work, lwork, info);
#endif
  }

  void gelss(long* m, long* n, long* nrhs, float* A, long* lda, float* B, long* ldb, 
             float* s, float* rcond, long* rank, float* work, long* lwork, long* info)
  {
#ifdef HAVE_LAPACK
    F77NAME(sgelss)(m, n, nrhs, A, lda, B, ldb, s, rcond, rank, work, lwork, info);
#endif
  }  


  void dgebrd(long* m, long *n, double* a, long* lda, double* d, double* e, double* tauq, 
              double* taup, double* work, long* lwork,long* info){

#ifdef HAVE_LAPACK
    F77NAME(dgebrd)( m, n, a,  lda, d, e, tauq, taup, work, lwork, info);
#endif
  }

  void dbdsdc ( char* uplo, char* compq, long* n, double* d, double* e, double* u, long* ldu, 
                double* vt, long* ldvt, double* q, long* iq, double* work, long* iwork, long* info ){
#ifdef HAVE_LAPACK
    F77NAME(dbdsdc) (uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info );
#endif
  }

  void dgesdd( char* jobz, long* m, long* n, double* a, long* lda, double* s, double* u, 
               long *ldu, double* vt, long* ldvt, double* work, long* lwork, long* iwork, 
               long* info ){

#ifdef HAVE_LAPACK
    F77NAME(dgesdd)(jobz,  m,  n,  a,  lda,  s,  u, ldu,  vt,  ldvt,  work,  lwork,  iwork, info );
#endif
  }

} //namespace
