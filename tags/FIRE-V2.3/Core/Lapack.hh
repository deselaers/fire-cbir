#ifndef _LAPACK_HH
#define _LAPACK_HH

namespace Lapack{

  void dgelsd(long* m, long* n, long* nrhs, double* a, long* lda, double* b, long* ldb, 
              double* s, double* rcond, long* rank, double* work, long* lwork, long* iwork, long* info);
  
  void dggevx(char* balanc, char* jobvl, char* jobvr, char* sense,
              long* n, double* a, long* lda, double* b, long* ldb,
              double* alphar, double* alphai, double* beta, double* vl, long* ldvl, double* vr, long* ldvr,
              long* ilo, long* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, 
              double* rconde, double* rcondv, double* work, long* lwork, long* iwork, bool* bwork, long* info); 
  
  /**  routine is deprecated and has been replaced by routine SGGEV 
       void dgegv_(char *jobvl,char *jobvr,int *n,double *a,int *lda,double *b, int *ld
                   double *alphar,double *alphai,double *beta,double *vl,
                   int *ldvl,double *vr,int *ldvr,double *work,int *lwork,
                   int *info) ;
  */

  void sggev(char *jobvl,char *jobvr,int *n,double *a,int *lda,double *b, int *ld,
             double *alphar,double *alphai,double *beta,double *vl,
             int *ldvl,double *vr,int *ldvr,double *work,int *lwork,
             int *info) ;

  void dggev(char *jobvl,char *jobvr,int *n,double *a,int *lda,double *b, int *ld,
             double *alphar,double *alphai,double *beta,double *vl,
             int *ldvl,double *vr,int *ldvr,double *work,int *lwork,
             int *info) ;

  
  void dsygvd(long* itype, char* jobz, char* uplo, long* n, double* a, long* lda, double* b, long* ldb,
              double* w, double* work, long* lwork, long* iwork, long* liwork, long* info);
  
  void gels(char* trans, long* m, long* n, long* nrhs, double* A, long* lda, 
            double* B, long* ldb, double* work, long* lwork, long* info);
  
  void gels(char* trans, long* m, long* n, long* nrhs, float* A, long* lda, 
            float* B, long* ldb, float* work, long* lwork, long* info);
    
  void gelss(long* m, long* n, long* nrhs, double* A, long* lda, double* B, long* ldb, 
             double* s, double* rcond, long* rank, double* work, long* lwork, long* info);
    
  void gelss(long* m, long* n, long* nrhs, float* A, long* lda, float* B, long* ldb, 
             float* s, float* rcond, long* rank, float* work, long* lwork, long* info);
    
  void dgebrd(long* m, long *n, double* a, long* lda, double* d, double* e, double* tauq, 
              double* taup, double* work, long* lwork,long* info);

  void dbdsdc ( char* uplo, char* compq, long* n, double* d, double* e, double* u, long* ldu, 
                double* vt, long* ldvt, double* q, long* iq, double* work, long* iwork, long* info );

  void dgesdd( char* jobz, long* m, long* n, double* a, long* lda, double* s, double* u, 
               long *ldu, double* vt, long* ldvt, double* work, long* lwork, long* iwork, 
               long* info );



}

#endif // _LAPACK_HH
