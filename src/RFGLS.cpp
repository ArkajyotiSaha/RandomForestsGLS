#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <stdio.h>
#include <limits>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define Ind(a,b) (((a)==(b))?(1.0):(0.0))
#define swapInt(a, b) ((a ^= b), (b ^= a), (a ^= b))
#define NODE_TERMINAL -1
#define NODE_TOSPLIT  -2
#define NODE_INTERIOR -3
#define a 1.0e-10
#define a_jump 1
#define inf_temp 1/1e-22

//Description: update B and F.
void updateBF_org(double *B, double *F, double *c, double *C, double *D, double *d, int *nnIndxLU, int *CIndx, int n, double *theta, int covModel, int nThreads){
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';

  double nu = 0;
  //check if the model is 'matern'
  if (covModel == 2) {
    nu = theta[2];
  }

  double *bk = (double *) calloc(nThreads*(static_cast<int>(1.0+5.0)), sizeof(double));


  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(5.0));
  int threadID = 0;

#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID)
#endif
  for(i = 0; i < n; i++){
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    if(i > 0){
      for(k = 0; k < nnIndxLU[n+i]; k++){
        c[nnIndxLU[i]+k] = spCor(d[nnIndxLU[i]+k], theta[1], nu, covModel, &bk[threadID*nb]);
        for(l = 0; l <= k; l++){
          C[CIndx[i]+l*nnIndxLU[n+i]+k] = spCor(D[CIndx[i]+l*nnIndxLU[n+i]+k], theta[1], nu, covModel, &bk[threadID*nb]);
          if(l == k){
            C[CIndx[i]+l*nnIndxLU[n+i]+k] += theta[0];
          }
        }
      }
      F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotri failed\n");}
      F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[CIndx[i]], &nnIndxLU[n+i], &c[nnIndxLU[i]], &inc, &zero, &B[nnIndxLU[i]], &inc);
      F[i] = 1 - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[nnIndxLU[i]], &inc) + theta[0];
    }else{
      B[i] = 0;
      F[i] = 1 + theta[0];
    }
  }
  free(bk);
}

extern "C" {

  double pinv_dgelsy_rss_cpp(double *A, double *b, int nrowA, int ncolA){


    int inc = 1;
    int nlengthb=nrowA;
    int nlengthA=nrowA * ncolA;

    double *y0 = (double *) malloc(nlengthb * sizeof(double));
    F77_NAME(dcopy)(&nlengthb, b, &inc, y0, &inc);

    double *X0 = (double *) malloc(nlengthA * sizeof(double));
    F77_NAME(dcopy)(&nlengthA, A, &inc, X0, &inc);

    int    m = nrowA, n = ncolA, nrhs = 1, lda = nrowA, ldb = nrowA, info, rank;
    double rcond = -1.0;
    /* Local arrays */
    int *jvpt;
    jvpt = new int [ncolA];
    for(int zp = 0; zp < ncolA; zp++){
      jvpt[zp] = 0;
    }
    int lwork = -1;

    int* iwork = NULL;
    double* work = NULL;

    double work_query;
    F77_NAME(dgelsy)(&m, &n, &nrhs, X0, &lda, y0, &ldb, jvpt, &rcond, &rank, &work_query, &lwork, &info);

    lwork = (int)work_query;

    work = (double *) malloc(lwork * sizeof(double));
    F77_NAME(dgelsy)(&m, &n, &nrhs, X0, &lda, y0, &ldb, jvpt, &rcond, &rank, work, &lwork, &info);



    //create Xbeta
    char const *ntran = "N";
    const double one = 1.0;
    const double zero = 0.0;
    double *tmp_n = (double *) malloc(m * sizeof(double));
    F77_NAME(dgemv)(ntran, &m, &n, &one, A, &m, y0, &inc, &zero, tmp_n, &inc);


    //create residual
    const double negOne = -1.0;
    F77_NAME(daxpy)(&m, &negOne, b, &inc, tmp_n, &inc);

    double rss = 0.0;
    for(int ilp = 0; ilp < m; ilp++){
      rss = rss + pow(tmp_n[ilp], 2.0);
    }

    free(y0);
    free(X0);
    free(iwork);
    free(work);
    free(jvpt);
    free(tmp_n);
    return(rss);

  }

  double* pinv_dgelsy_beta_cpp(double *A, double *b, int nrowA, int ncolA){


    int inc = 1;
    int nlengthb=nrowA;
    int nlengthA=nrowA * ncolA;

    double *y0 = (double *) malloc(nlengthb * sizeof(double));
    F77_NAME(dcopy)(&nlengthb, b, &inc, y0, &inc);

    double *X0 = (double *) malloc(nlengthA * sizeof(double));
    F77_NAME(dcopy)(&nlengthA, A, &inc, X0, &inc);

    int    m = nrowA, n = ncolA, nrhs = 1, lda = nrowA, ldb = nrowA, info, rank;
    double rcond = -1.0;
    /* Local arrays */
    int *jvpt;
    jvpt = new int [ncolA];
    for(int zp = 0; zp < ncolA; zp++){
      jvpt[zp] = 0;
    }
    int lwork = -1;


    double* work = NULL;

    double work_query;
    F77_NAME(dgelsy)(&m, &n, &nrhs, X0, &lda, y0, &ldb, jvpt, &rcond, &rank, &work_query, &lwork, &info);

    lwork = (int)work_query;

    work = (double *) malloc(lwork * sizeof(double));
    F77_NAME(dgelsy)(&m, &n, &nrhs, X0, &lda, y0, &ldb, jvpt, &rcond, &rank, work, &lwork, &info);

    free(X0);
    free(work);
    free(jvpt);
    return(y0);

  }

  double pinv_dgelsd_rss_cpp(double *A, double *b, int nrowA, int ncolA){

    int inc = 1;
    int nlengthb=nrowA;
    int nlengthA=nrowA * ncolA;

    double *y0 = (double *) malloc(nlengthb * sizeof(double));
    F77_NAME(dcopy)(&nlengthb, b, &inc, y0, &inc);

    double *X0 = (double *) malloc(nlengthA * sizeof(double));
    F77_NAME(dcopy)(&nlengthA, A, &inc, X0, &inc);

    int    m = nrowA, n = ncolA, nrhs = 1, lda = nrowA, ldb = nrowA, info, rank;
    double rcond = -1.0;
    /* Local arrays */
    double *s;
    s = new double [nrowA];
    //double s[nrowA];
    int lwork = -1;
    int liwork;
    int* iwork = NULL;
    double* work = NULL;
    int iwork_query;
    double work_query;
    F77_NAME(dgelsd)(&m, &n, &nrhs, X0, &lda, y0, &ldb, s, &rcond, &rank, &work_query, &lwork, &iwork_query, &info);
    liwork = (int)iwork_query;
    lwork = (int)work_query;
    iwork = (int *) malloc(liwork * sizeof(int));
    work = (double *) malloc(lwork * sizeof(double));
    F77_NAME(dgelsd)(&m, &n, &nrhs, X0, &lda, y0, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);

    //calculate rss
    double rss = 0.0;
    for(int lip = 0; lip < (nrowA - ncolA -1); lip++){
      rss = rss + pow(y0[ncolA+1+lip], 2.0);
    }

    free(y0);
    free(X0);
    free(iwork);
    free(work);
    free(s);
    return(rss);
  }

  double* pinv_dgelsd_beta_cpp(double *A, double *b, int nrowA, int ncolA){

    int inc = 1;
    int nlengthb=nrowA;
    int nlengthA=nrowA * ncolA;

    double *y0 = (double *) malloc(nlengthb * sizeof(double));
    F77_NAME(dcopy)(&nlengthb, b, &inc, y0, &inc);

    double *X0 = (double *) malloc(nlengthA * sizeof(double));
    F77_NAME(dcopy)(&nlengthA, A, &inc, X0, &inc);

    int    m = nrowA, n = ncolA, nrhs = 1, lda = nrowA, ldb = nrowA, info, rank;
    double rcond = -1.0;
    /* Local arrays */
    double *s;
    s = new double [nrowA];
    //double s[nrowA];
    int lwork = -1;
    int liwork;
    int* iwork = NULL;
    double* work = NULL;
    int iwork_query;
    double work_query;
    F77_NAME(dgelsd)(&m, &n, &nrhs, X0, &lda, y0, &ldb, s, &rcond, &rank, &work_query, &lwork, &iwork_query, &info);
    liwork = (int)iwork_query;
    lwork = (int)work_query;
    iwork = (int *) malloc(liwork * sizeof(int));
    work = (double *) malloc(lwork * sizeof(double));
    F77_NAME(dgelsd)(&m, &n, &nrhs, X0, &lda, y0, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);

    double *beta = (double *) malloc(ncolA * sizeof(double));
    //calculate beta
    for(int lip = 0; lip < ncolA; lip++){
      beta[lip] = y0[lip];
    }

    free(y0);
    free(X0);
    free(iwork);
    free(work);
    free(s);
    return(beta);
  }


  SEXP RFGLS_BFcpp(SEXP n_r, SEXP m_r, SEXP coords_r, SEXP covModel_r, SEXP alphaSqStarting_r, SEXP phiStarting_r, SEXP nuStarting_r,
                              SEXP sType_r, SEXP nThreads_r, SEXP verbose_r){

    int i, j, k, l, nProtect=0;

    //get args
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    double *coords = REAL(coords_r);

    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);

    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];



#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Calculation covariance with %i locations.\n\n", n);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    }

    //parameters
    int nTheta;

    if(corName != "matern"){
      nTheta = 2;//tau^2 = 0, phi = 1
    }else{
      nTheta = 3;//tau^2 = 0, phi = 1, nu = 2;
    }
    //starting
    double *theta = (double *) R_alloc (nTheta, sizeof(double));

    theta[0] = pow(REAL(alphaSqStarting_r)[0], 2.0);
    theta[1] = pow(REAL(phiStarting_r)[0], 2.0);

    if(corName == "matern"){
      theta[2] = pow(REAL(nuStarting_r)[0], 2.0);
    }

    //allocated for the nearest neighbor index vector (note, first location has no neighbors).
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
    SEXP nnIndx_r; PROTECT(nnIndx_r = allocVector(INTSXP, nIndx)); nProtect++; int *nnIndx = INTEGER(nnIndx_r);
    SEXP d_r; PROTECT(d_r = allocVector(REALSXP, nIndx)); nProtect++; double *d = REAL(d_r);
    SEXP nnIndxLU_r; PROTECT(nnIndxLU_r = allocVector(INTSXP, 2*n)); nProtect++; int *nnIndxLU = INTEGER(nnIndxLU_r); //first column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but will simlifying some parallelization).
    //make the neighbor index
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tBuilding neighbor index\n");
#ifdef Win32
      R_FlushConsole();
#endif
    }

    if(INTEGER(sType_r)[0] == 0){
      mkNNIndx(n, m, coords, nnIndx, d, nnIndxLU);
    }else{
      mkNNIndxTree0(n, m, coords, nnIndx, d, nnIndxLU);
    }


    SEXP CIndx_r; PROTECT(CIndx_r = allocVector(INTSXP, 2*n)); nProtect++; int *CIndx = INTEGER(CIndx_r); //index for D and C.
    for(i = 0, j = 0; i < n; i++){//zero should never be accessed
      j += nnIndxLU[n+i]*nnIndxLU[n+i];
      if(i == 0){
        CIndx[n+i] = 0;
        CIndx[i] = 0;
      }else{
        CIndx[n+i] = nnIndxLU[n+i]*nnIndxLU[n+i];
        CIndx[i] = CIndx[n+i-1] + CIndx[i-1];
      }
    }

    SEXP j_r; PROTECT(j_r = allocVector(INTSXP, 1)); nProtect++; INTEGER(j_r)[0] = j;

    SEXP D_r; PROTECT(D_r = allocVector(REALSXP, j)); nProtect++; double *D = REAL(D_r);

    for(i = 0; i < n; i++){
      for(k = 0; k < nnIndxLU[n+i]; k++){
        for(l = 0; l <= k; l++){
          D[CIndx[i]+l*nnIndxLU[n+i]+k] = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);
        }
      }
    }

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tCalculationg the approximate Cholesky Decomposition\n");
#ifdef Win32
      R_FlushConsole();
#endif
    }

    SEXP B_r; PROTECT(B_r = allocVector(REALSXP, nIndx)); nProtect++; double *B = REAL(B_r);

    SEXP F_r; PROTECT(F_r = allocVector(REALSXP, n)); nProtect++; double *F = REAL(F_r);

    double *c =(double *) R_alloc(nIndx, sizeof(double));
    double *C = (double *) R_alloc(j, sizeof(double)); zeros(C, j);


    updateBF_org(B, F, c, C, D, d, nnIndxLU, CIndx, n, theta, covModel, nThreads);

    //return stuff
    SEXP result_r, resultName_r;
    int nResultListObjs = 8;



    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, B_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("B"));

    SET_VECTOR_ELT(result_r, 1, F_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("F"));

    SET_VECTOR_ELT(result_r, 2, nnIndxLU_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("nnIndxLU"));

    SET_VECTOR_ELT(result_r, 3, CIndx_r);
    SET_VECTOR_ELT(resultName_r, 3, mkChar("CIndx"));

    SET_VECTOR_ELT(result_r, 4, D_r);
    SET_VECTOR_ELT(resultName_r, 4, mkChar("D"));

    SET_VECTOR_ELT(result_r, 5, d_r);
    SET_VECTOR_ELT(resultName_r, 5, mkChar("d"));

    SET_VECTOR_ELT(result_r, 6, nnIndx_r);
    SET_VECTOR_ELT(resultName_r, 6, mkChar("nnIndx"));

    SET_VECTOR_ELT(result_r, 7, j_r);
    SET_VECTOR_ELT(resultName_r, 7, mkChar("Length.D"));

    namesgets(result_r, resultName_r);

    //unprotect
    UNPROTECT(nProtect);


    return(result_r);
  }
  //invP_val = contains inverse of P. is of nsample length
  //invP_loc[i,i+1] = starting and ending location of invP(i). is of n+1 length.
  //invZ_val = contains inverse of nn mapping for every i.
  //invZ_loc[i,i+1] = starting and ending location of invZ(i). is of n+1 length.
  //double *PQZ is a n x t_n length vector (a matrix, entered as vector).

  void PQZ_update(int *P, int *Z, int *invP_val, int *invP_loc, double *B, double *F, int *nnIndx, int *nnIndxLU, int n, int rc, double *PQZ){
    int i, j, l, t, temp_t;
    double tt_1, tt_2;
    for(i = 0; i < n; i++){
      temp_t = invP_loc[i+1] - invP_loc[i];
      if(temp_t > 0){
        for(l = 0; l < rc; l++){
          tt_1 = 0;
          for(j = 0; j < nnIndxLU[n+i]; j++){
            tt_1 += B[nnIndxLU[i]+j]*Ind(Z[nnIndx[nnIndxLU[i]+j]],l);
          }
          tt_2 = (Ind(Z[i],l) - tt_1)/sqrt(F[i]);
          for(t = 0; t < temp_t; t++){
            PQZ[n*l + invP_val[invP_loc[i]+t] ] = tt_2;
          }
        }
      }
    }
  }

  void PQy_update(int *P, double *y, int *invP_val, int *invP_loc, double *B, double *F, int *nnIndx, int *nnIndxLU, int n, double *PQy){
    int i, j, t, temp_t;
    double tt_1, tt_2;
    for(i = 0; i < n; i++){
      temp_t = invP_loc[i+1] - invP_loc[i];
      if(temp_t > 0){
        tt_1 = 0;
        for(j = 0; j < nnIndxLU[n+i]; j++){
          tt_1 += B[nnIndxLU[i]+j]*y[nnIndx[nnIndxLU[i]+j]];
        }
        tt_2 = (y[i] - tt_1)/sqrt(F[i]);
        for(t = 0; t < temp_t; t++){
          PQy[invP_val[invP_loc[i] + t]] = tt_2;
        }
      }
    }
  }

  SEXP RFGLS_invZcpp(SEXP n_r, SEXP nnIndx_r, SEXP nnIndxLU_r, SEXP invZ_freq_r, SEXP invZ_val_r, SEXP invZ_loc_r, SEXP i_loc_r){

    int n = INTEGER(n_r)[0];
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *invZ_val = INTEGER(invZ_val_r);
    int *invZ_loc = INTEGER(invZ_loc_r);
    int *invZ_freq = INTEGER(invZ_freq_r);
    int i, j;
    int *i_loc = INTEGER(i_loc_r);


    for(i = 0; i < n; i++){
      for(j = 0; j < nnIndxLU[n+i]; j++){
        invZ_freq[nnIndx[nnIndxLU[i]+j]] += 1;
      }
    }

    int cumul_count = 0;
    invZ_loc[0] = 0;
    for(i = 1; i < n; i++){
      cumul_count += invZ_freq[i-1];
      invZ_loc[i] = cumul_count;
    }
    invZ_loc[n] = cumul_count;

    for(i = 0; i < n; i++){
      i_loc[i] = 0;
    }

    for(i = 0; i < n; i++){
      for(j = 0; j < nnIndxLU[n+i]; j++){
        invZ_val[invZ_loc[nnIndx[nnIndxLU[i]+j]] + i_loc[nnIndx[nnIndxLU[i]+j]]] = i;
        i_loc[nnIndx[nnIndxLU[i]+j]] = i_loc[nnIndx[nnIndxLU[i]+j]] + 1;
      }
    }


    SEXP result_r, resultName_r;
    int nResultListObjs = 2, nProtect=0;


    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, invZ_val_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("invZ_val"));

    SET_VECTOR_ELT(result_r, 1, invZ_loc_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("invZ_loc"));

    namesgets(result_r, resultName_r);

    //unprotect
    UNPROTECT(nProtect);


    return(result_r);
  }

  void findBestSplit(double *X, int n, int *jdex, int mdim, int nsample,
                     int ndstart, int ndend, int &msplit, double &decsplit,
                     double &ubest, int &ndendl, int &jstat, int mtry,
                     int nodecnt, int *Z_index, int rc, int *invP_val,
                     int *invP_loc, int *invZ_val, int *invZ_loc, double *PQZ,
                     double *PQy, double *B, double *F, int *nnIndx,
                     int *nnIndxLU, int pinv_choice) {
    int last, nl, nr;
    int nwhole_l, nwhole_r;
    int i, j, kv, tieVal, tieVar;
    double ubestt;
    double crit, critmax, critvar;

    double *ut = (double *) calloc (n , sizeof(double));
    double *xt = (double *) calloc (n , sizeof(double));
    double *v = (double *) calloc (n , sizeof(double));
    int *yl_index = (int *) calloc (n , sizeof(int));
    int *mind = (int *) calloc (mdim , sizeof(int));
    int *ncase = (int *) calloc (n , sizeof(int));
    int *Z_index_local = (int *) calloc (n, sizeof(int));
    double *PQy_local = (double *) calloc (nsample, sizeof(double));
    int PQZ_length = nsample * (rc + 1);
    double *PQZ_local = (double *) calloc (PQZ_length, sizeof(double));

    int temp_t, temp_i, l, t, ji;
    double tt_1;

    int PQy_length = nsample;


    /* START BIG LOOP */
    msplit = -1;
    decsplit = 0.0;
    //critmax = arma::datum::inf; for convenience, use a very large number
    critmax = 10000000000;
    ubestt = 0.0;

    for (i=0; i < mdim; ++i) mind[i] = i;

    last = mdim - 1;
    tieVar = 1;

    for(i = 0; i < mtry; ++i){
      //critvar = std::numeric_limits<double>::infinity(); for convenience, use a very large number
      critvar = 10000000000;
      j = (int) (unif_rand() * (last+1));
      kv = mind[j];
      swapInt(mind[j], mind[last]);
      last--;

      nwhole_l = 0;
      nwhole_r = 0;
      int npop_tot = 0;
      int inc = 1;

      F77_NAME(dcopy)(&PQZ_length, PQZ, &inc, PQZ_local, &inc);
      F77_NAME(dcopy)(&PQy_length, PQy, &inc, PQy_local, &inc);


      //copy the x data in this node
      for (j = ndstart; j <= ndend; ++j) {
        xt[j] = X[kv + (jdex[j] - 1) * mdim];
        yl_index[j] = jdex[j] - 1;
        npop_tot = npop_tot + (invP_loc[jdex[j]] - invP_loc[jdex[j] - 1]);
      }

      for (j = ndstart; j <= ndend; ++j) v[j] = xt[j];
      for (j = 1; j <= n; ++j) ncase[j - 1] = j;

      R_qsort_I(v, ncase, ndstart + 1, ndend + 1);

      if (v[ndstart] >= v[ndend]) continue;
      //crit = std::numeric_limits<double>::infinity(); for convenience, use a very large number
      crit = 10000000000;
      nwhole_l = 0;
      nwhole_r = npop_tot;
      tieVal = 1;

      for(j = ndstart; j <= ndend - 1; ++j) {
        temp_i = yl_index[ncase[j] - 1];
        Z_index_local[temp_i] = rc;
        temp_t = invP_loc[temp_i+1] - invP_loc[temp_i];
        if(temp_t > 0){
          nwhole_l = nwhole_l + temp_t;
          nwhole_r = nwhole_r - temp_t;

          l = rc;
          tt_1 = 0;
          for(ji = 0; ji < nnIndxLU[n+temp_i]; ji++){
            tt_1 += B[nnIndxLU[temp_i]+ji]*Ind(Z_index_local[nnIndx[nnIndxLU[temp_i]+ji]],l);
          }
          for(t = 0; t < temp_t; t++){
            PQZ_local[n*l + invP_val[invP_loc[temp_i]+t] ] = (1 - tt_1)/sqrt(F[temp_i]);
          }

          l = Z_index[temp_i];
          tt_1 = 0;
          for(ji = 0; ji < nnIndxLU[n+temp_i]; ji++){
            tt_1 += B[nnIndxLU[temp_i]+ji]*Ind(Z_index_local[nnIndx[nnIndxLU[temp_i]+ji]],l);
          }
          for(t = 0; t < temp_t; t++){
            PQZ_local[n*l + invP_val[invP_loc[temp_i]+t] ] = (0 - tt_1)/sqrt(F[temp_i]);
          }
        }

        if (v[j] < v[j+1]) {
          if(nwhole_l == 0){
            crit = 10000000000;
          }
          if(nwhole_r == 0){
            crit = 10000000000;
          }
          if(nwhole_l > 0 && nwhole_r > 0){
            if(pinv_choice == 1){
              crit = pinv_dgelsy_rss_cpp(PQZ_local, PQy_local, nsample, rc+1)/nsample;
            }
            if(pinv_choice == 0){
              crit = pinv_dgelsd_rss_cpp(PQZ_local, PQy_local, nsample, rc+1)/nsample;
            }
          }

          if (crit < critvar) {
            ubestt = (v[j] + v[j+1]) / 2.0;
            critvar = crit;
            tieVal = 1;
          }
          if (crit == critvar) {
            tieVal++;
            if(unif_rand() < 1.0/tieVal){
              ubestt =  (v[j] + v[j+1]) / 2.0;
              critvar = crit;
            }
          }
        }
      }


      if (critvar < critmax) {
        ubest = ubestt;
        msplit = kv + 1;
        critmax = critvar;
        for (j = ndstart; j <= ndend; ++j) {
          ut[j] = xt[j];
        }
        tieVar = 1;
      }
      if (critvar == critmax) {
        tieVar++;
        if (unif_rand() < 1.0 / tieVar) {
          ubest = ubestt;
          msplit = kv + 1;
          critmax = critvar;
          for (j = ndstart; j <= ndend; ++j) {
            ut[j] = xt[j];
          }
        }
      }


    }
    decsplit = critmax;

    /* If best split can not be found, set to terminal node and return. */
    if (msplit != -1) {
      nl = ndstart;
      for (j = ndstart; j <= ndend; ++j) {
        if (ut[j] <= ubest) {
          nl++;
          ncase[nl-1] = jdex[j];
        }
      }
      ndendl = imax2(nl - 1, ndstart);
      nr = ndendl + 1;
      for (j = ndstart; j <= ndend; ++j) {
        if (ut[j] > ubest) {
          if (nr >= n) break;
          nr++;
          ncase[nr - 1] = jdex[j];
          Z_index[jdex[j] - 1] = rc;
        }
      }
      if (ndendl >= ndend) ndendl = ndend - 1;
      for (j = ndstart; j <= ndend; ++j) jdex[j] = ncase[j];

    }
    else jstat = 1;

    free(ut);
    free(xt);
    free(v);
    free(yl_index);
    free(mind);
    free(ncase);
    free(Z_index_local);
    free(PQy_local);
    free(PQZ_local);
  }

  void predictRegTree(double *x, int nsample, int mdim,
                      int *lDaughter, int *rDaughter, int *nodestatus,
                      double *ypred, double *upper, double *avnode,
                      int *mbest) {
    int i, k, m;

    for (i = 0; i < nsample; ++i) {
      k = 0;
      while (nodestatus[k] != NODE_TERMINAL) { /* go down the tree */
    m = mbest[k] - 1;
        k = (x[m + i*mdim] <= upper[k]) ?
        lDaughter[k] - 1 : rDaughter[k] - 1;
      }
      /* terminal node: assign prediction and move on to next */
      ypred[i] = avnode[k];
    }
  }

  //basic structure of RFtree

  SEXP RFGLStree_cpp(SEXP X_r, SEXP y_r, SEXP B_r, SEXP F_r, SEXP nnIndx_r, SEXP nnIndxLU_r, SEXP invZ_val_r, SEXP invZ_loc_r,
                         SEXP mtry_r, SEXP n_r, SEXP p_r, SEXP nsample_r, SEXP nthsize_r, SEXP nrnodes_r, SEXP treeSize_r, SEXP pinv_choice_r, SEXP Xtest_r, SEXP ntest_r, SEXP nThreads_r, SEXP q_r){
    double *X = REAL(X_r);
    double *y = REAL(y_r);
    double *B = REAL(B_r);
    double *F = REAL(F_r);
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *invZ_val = INTEGER(invZ_val_r);
    int *invZ_loc = INTEGER(invZ_loc_r);
    int mtry = INTEGER(mtry_r)[0];
    int n = INTEGER(n_r)[0];
    int p = INTEGER(p_r)[0];
    int nsample = INTEGER(nsample_r)[0];
    int nthsize = INTEGER(nthsize_r)[0];
    int nrnodes = INTEGER(nrnodes_r)[0];
    int treeSize = INTEGER(treeSize_r)[0];
    int pinv_choice = INTEGER(pinv_choice_r)[0];
    double *Xtest = REAL(Xtest_r);
    int ntest = INTEGER(ntest_r)[0];

    int q = INTEGER(q_r)[0];


    //assign memories

    int nProtect = 0;


    SEXP P_index_r; PROTECT(P_index_r = allocVector(INTSXP, nsample)); nProtect++; int *P_index = INTEGER(P_index_r);
    int *invP_freq = (int *) R_alloc (n, sizeof(int));

    int *invP_loc = (int *) R_alloc ((n + 1), sizeof(int));

    int *invP_val = (int *) R_alloc (n, sizeof(int));

    int *i_loc = (int *) R_alloc (n, sizeof(int));

    int i, k, l;
    double xrand;

    SEXP lDaughter_r; PROTECT(lDaughter_r = allocVector(INTSXP, nrnodes)); nProtect++; int *lDaughter = INTEGER(lDaughter_r);

    SEXP rDaughter_r; PROTECT(rDaughter_r = allocVector(INTSXP, nrnodes)); nProtect++; int *rDaughter = INTEGER(rDaughter_r);

    SEXP nodestatus_r; PROTECT(nodestatus_r = allocVector(INTSXP, nrnodes)); nProtect++; int *nodestatus = INTEGER(nodestatus_r);

    int *nodestart = (int *) R_alloc (nrnodes, sizeof(int));
    int *nodepop = (int *) R_alloc (nrnodes, sizeof(int));

    SEXP mbest_r; PROTECT(mbest_r = allocVector(INTSXP, nrnodes)); nProtect++; int *mbest = INTEGER(mbest_r);
    SEXP upper_r; PROTECT(upper_r = allocVector(REALSXP, nrnodes)); nProtect++; double *upper = REAL(upper_r);

    int *avnode_number = (int *) R_alloc (nrnodes, sizeof(int));
    for(int avcount = 0; avcount < nrnodes; avcount++){
      avnode_number[avcount] = -99;
    }
    SEXP avnode_r; PROTECT(avnode_r = allocVector(REALSXP, nrnodes)); nProtect++; double *avnode = REAL(avnode_r);


    for(l = 0; l < nrnodes; l++){
      lDaughter[l]  = 0;
      rDaughter[l] = 0;
      upper[l]  = 0;
      avnode[l] = 0;
      nodestatus[l] = 0;
      mbest[l] = 0;
    }

    for(i = 0; i < nsample; i++){
      invP_freq[i] = 0;
    }

    //Sampling with replacement
    GetRNGstate();
    for (i = 0; i < nsample; i++) {
      xrand = unif_rand();
      k = xrand * (n - q) + q;
      P_index[i] = k;
      invP_freq[k] = invP_freq[k] + 1;
    }
    int cumul_count = 0;
    invP_loc[0] = 0;
    for(i = 1; i < n; i++){
      cumul_count += invP_freq[i-1];
      invP_loc[i] = cumul_count;
    }
    invP_loc[n] = cumul_count + invP_freq[n-1];

    for(i = 0; i < n; i++){
      i_loc[i] = 0;
    }

    for(i = 0; i < nsample; i++){
      invP_val[invP_loc[P_index[i]] + i_loc[P_index[i]]] = i;
      i_loc[P_index[i]] = i_loc[P_index[i]] + 1;
    }


    int *jdex = (int *) R_alloc (n, sizeof(int));

    int *Z_index = (int *) R_alloc (n, sizeof(int));
    double *PQy = (double *) R_alloc (nsample, sizeof(double));


    PQy_update(P_index, y, invP_val, invP_loc, B, F, nnIndx, nnIndxLU, n, PQy);


    int ncur, ndstart, ndend, ndendl, nodecnt, jstat, msplit = 0;
    double decsplit, ubest;


    for (i = 1; i <= n; ++i){
      jdex[i-1] = i;
      Z_index[i-1] = 0;
    }

    ncur = 0;
    nodestart[0] = 0;
    nodepop[0] = n;
    nodestatus[0] = NODE_TOSPLIT;
    int PQZ_length;

    int rc = 1;

    double *PQZ;

    for (k = 0; k < nrnodes - 2; ++k){
      if (k > ncur || ncur >= nrnodes - 2) break;
      /* skip if the node is not to be split */
      if (nodestatus[k] != NODE_TOSPLIT) continue;

      /* initialize for next call to findbestsplit */

      ndstart = nodestart[k];
      ndend = ndstart + nodepop[k] - 1;
      nodecnt = nodepop[k];
      jstat = 0;
      decsplit = 0.0;
      PQZ_length = nsample * (rc + 1);
      PQZ = (double *) malloc (PQZ_length * sizeof(double));

      PQZ_update(P_index, Z_index, invP_val, invP_loc, B, F, nnIndx, nnIndxLU, n, (rc+1), PQZ);

      findBestSplit(X, n, jdex, p, nsample, ndstart, ndend, msplit, decsplit, ubest, ndendl, jstat, mtry, nodecnt, Z_index, rc, invP_val,
                    invP_loc, invZ_val, invZ_loc, PQZ, PQy, B, F, nnIndx, nnIndxLU, pinv_choice);



      if (jstat == 1) {
        /* Node is terminal: Mark it as such and move on to the next. */
        nodestatus[k] = NODE_TERMINAL;
        avnode_number[k] = Z_index[jdex[nodestart[k]] - 1];
        continue;
      }

      mbest[k] = msplit;
      upper[k] = ubest;
      nodestatus[k] = NODE_INTERIOR;

      nodepop[ncur + 1] = ndendl - ndstart + 1;
      nodepop[ncur + 2] = ndend - ndendl;
      nodestart[ncur + 1] = ndstart;
      nodestart[ncur + 2] = ndendl + 1;

      nodestatus[ncur + 1] = NODE_TOSPLIT;
      avnode_number[ncur + 1] = Z_index[jdex[nodestart[ncur + 1]] - 1];
      if (nodepop[ncur + 1] <= nthsize) {
        nodestatus[ncur + 1] = NODE_TERMINAL;
      }

      nodestatus[ncur + 2] = NODE_TOSPLIT;
      avnode_number[ncur + 2] = Z_index[jdex[nodestart[ncur + 2]] - 1];
      if (nodepop[ncur + 2] <= nthsize) {
        nodestatus[ncur + 2] = NODE_TERMINAL;
      }

      lDaughter[k] = ncur + 1 + 1;
      rDaughter[k] = ncur + 2 + 1;
      /* Augment the tree by two nodes. */
      ncur += 2;
      rc = rc+1;
      free(PQZ);
    }
    treeSize = nrnodes;
    for (k = nrnodes - 1; k >= 0; --k) {
      if (nodestatus[k] == 0) (treeSize)--;
      if (nodestatus[k] == NODE_TOSPLIT) {
        nodestatus[k] = NODE_TERMINAL;
      }
    }
    PQZ_length = nsample * (rc);
    PQZ = (double *) malloc (PQZ_length * sizeof(double));
    PQZ_update(P_index, Z_index, invP_val, invP_loc, B, F, nnIndx, nnIndxLU, n, rc, PQZ);

    double *beta = (double *) R_alloc (rc, sizeof(double));

    if(pinv_choice == 0){
      beta = pinv_dgelsd_beta_cpp(PQZ, PQy, nsample, rc);
    }
    if(pinv_choice == 1){
      beta = pinv_dgelsy_beta_cpp(PQZ, PQy, nsample, rc);
    }

    free(PQZ);

    for(k = 0; k < treeSize; ++k){
      if(nodestatus[k] == NODE_TERMINAL){
        avnode[k] = beta[avnode_number[k]];
      }
    }

    PutRNGstate();
    SEXP ytest_r; PROTECT(ytest_r = allocVector(REALSXP, ntest)); nProtect++; double *ytest = REAL(ytest_r);

    predictRegTree(Xtest, ntest, p, lDaughter, rDaughter, nodestatus, ytest, upper, avnode, mbest);


    SEXP result_r, resultName_r;
    int nResultListObjs = 8;


    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, P_index_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("P_index"));

    SET_VECTOR_ELT(result_r, 1, ytest_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("ytest"));

    SET_VECTOR_ELT(result_r, 2, lDaughter_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("lDaughter"));

    SET_VECTOR_ELT(result_r, 3, rDaughter_r);
    SET_VECTOR_ELT(resultName_r, 3, mkChar("rDaughter"));

    SET_VECTOR_ELT(result_r, 4, nodestatus_r);
    SET_VECTOR_ELT(resultName_r, 4, mkChar("nodestatus"));

    SET_VECTOR_ELT(result_r, 5, upper_r);
    SET_VECTOR_ELT(resultName_r, 5, mkChar("upper"));

    SET_VECTOR_ELT(result_r, 6, avnode_r);
    SET_VECTOR_ELT(resultName_r, 6, mkChar("avnode"));

    SET_VECTOR_ELT(result_r, 7, mbest_r);
    SET_VECTOR_ELT(resultName_r, 7, mkChar("mbest"));

    namesgets(result_r, resultName_r);

    //unprotect
    UNPROTECT(nProtect);


    return(result_r);
  }

  SEXP RFGLSpredicttree_cpp(SEXP Xtest_r, SEXP ntest_r, SEXP mdim_r,
                 SEXP lDaughter_r, SEXP rDaughter_r, SEXP nodestatus_r,
                 SEXP upper_r, SEXP avnode_r,
                 SEXP mbest_r){

    double *Xtest = REAL(Xtest_r);
    int ntest = INTEGER(ntest_r)[0];
    int p = INTEGER(mdim_r)[0];
    int *lDaughter = INTEGER(lDaughter_r);
    int *rDaughter = INTEGER(rDaughter_r);
    int *nodestatus = INTEGER(nodestatus_r);
    double *upper = REAL(upper_r);
    double *avnode = REAL(avnode_r);
    int *mbest = INTEGER(mbest_r);

    int nProtect = 0;
    SEXP ytest_r; PROTECT(ytest_r = allocVector(REALSXP, ntest)); nProtect++; double *ytest = REAL(ytest_r);

    predictRegTree(Xtest, ntest, p, lDaughter, rDaughter, nodestatus, ytest, upper, avnode, mbest);

    SEXP result_r, resultName_r;
    int nResultListObjs = 1;


    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, ytest_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("ytest"));

    namesgets(result_r, resultName_r);

    //unprotect
    UNPROTECT(nProtect);


    return(result_r);
  }
}



