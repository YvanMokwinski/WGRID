#ifndef WMESH_BLAS_H
#define WMESH_BLAS_H

#include "wmesh-types.h"

#if 0

extern wmesh_int_t sgeev_( const char* jobvl, const char* jobvr,
                          const wmesh_int_t* n, float* a, const wmesh_int_t* lda, float* wr,
                          float* wi, float* vl, const wmesh_int_t *ldvl, float* vr,
                          const wmesh_int_t *ldvr );
extern wmesh_int_t dgeev_( const char *jobvl, const char* jobvr,
                          const wmesh_int_t *n, double* a, const  wmesh_int_t* lda, double* wr,
                          double* wi, double* vl, const wmesh_int_t* ldvl, double* vr,
                          const wmesh_int_t* ldvr );

#endif

extern "C" void sgeev_( const char* jobvl, const char* jobvr, wmesh_int_t* n, float* a,
		    wmesh_int_t* lda, float* wr, float* wi, float* vl,
		    wmesh_int_t* ldvl, float* vr, wmesh_int_t* ldvr, float* work,
		    wmesh_int_t* lwork, wmesh_int_t *info );
extern "C" void dgeev_( const char* jobvl, const char* jobvr, wmesh_int_t* n, double* a,
                   wmesh_int_t* lda, double* wr, double* wi, double* vl,
                   wmesh_int_t* ldvl, double* vr, wmesh_int_t* ldvr, double* work,
                   wmesh_int_t* lwork, wmesh_int_t *info );

#ifdef WMESH_OPEN_BLAS


#include </usr/lib/openblas/include/f77blas.h>
// #include </usr/lib/openblas/include/lapacke.h>


#define BLAS_dscal 	dscal_
#define BLAS_dcopy 	dcopy_
#define BLAS_dgemm 	dgemm_
#define BLAS_ddot 	ddot_

#define LAPACK_dgesv 	dgesv_
#define LAPACK_dgeev 	dgeev_
#define LAPACK_sgeev 	sgeev_

#else

#ifdef WMESH_MKL_BLAS

#if WMESH_ILP64
#define MKL_ILP64 1
#endif

#include "mkl.h"
#define BLAS_dscal 	dscal
#define BLAS_dcopy 	dcopy
#define BLAS_dgemm 	dgemm
#define BLAS_ddot  	ddot
#define LAPACK_dgesv 	dgesv
#define LAPACK_dgeev 	dgeev
#define LAPACK_sgeev 	sgeev

#else
#error choose WMESH_MKL_BLAS or WMESH_OPEN_BLAS
#endif

#endif

#endif
