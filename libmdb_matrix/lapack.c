#include <assert.h>
#include <stdio.h>

#include "lapack.h"


int eposv(const int order, const char uplo,
	   const int N, const int NRHS,
	   elem *A, const int lda,
	   elem *B, const int ldb) {
#ifdef DOUBLE_ELEM
  return LAPACKE_dposv(order, uplo, N, NRHS, A, lda, B, ldb);
#elif defined FLOAT_ELEM
  return LAPACKE_sposv(order, uplo, N, NRHS, A, lda, B, ldb);
#else
#error ?
#endif
}


int epotrf(const int order, const char uplo,
	   const int N, elem *A, const int lda) {
#ifdef DOUBLE_ELEM
  return LAPACKE_dpotrf(order, uplo, N, A, lda);
#elif defined FLOAT_ELEM
  return LAPACKE_spotrf(order, uplo, N, A, lda);
#else
#error ?
#endif
}


int epotrs(const int order, const char uplo,
	   const int N, const int NRHS,
	   const elem *A, const int lda, elem *B, const int ldb) {
#ifdef DOUBLE_ELEM
  return LAPACKE_dpotrs(order, uplo, N, NRHS, A, lda, B, ldb);
#elif defined FLOAT_ELEM
  return LAPACKE_spotrs(order, uplo, N, NRHS, A, lda, B, ldb);
#else
#error ?
#endif
}


int eppsv(const int order, const char uplo,
	  const int N, const int NRHS,
	  elem *A, elem *B, const int ldb) {
#ifdef DOUBLE_ELEM
    return LAPACKE_dppsv(order, uplo, N, NRHS, A, B, ldb);
#elif FLOAT_ELEM
    return LAPACKE_sppsv(order, uplo, N, NRHS, A, B, ldb);
#else
#error ?
#endif
}


int epptrf(const int order, const char uplo,
	   const int N, elem *A) {
#ifdef DOUBLE_ELEM
    return LAPACKE_dpptrf(order, uplo, N, A);
#elif FLOAT_ELEM
    return LAPACKE_spptrf(order, uplo, N, A);
#else
#error ?
#endif
}


int epptrs(const int order, const char uplo,
	   const int N, const int NRHS,
	   const elem *A, elem *B, const int ldb) {
#ifdef DOUBLE_ELEM
    return LAPACKE_dpptrs(order, uplo, N, NRHS, A, B, ldb);
#elif FLOAT_ELEM
    return LAPACKE_spptrs(order, uplo, N, NRHS, A, B, ldb);
#else
#error ?
#endif
}
