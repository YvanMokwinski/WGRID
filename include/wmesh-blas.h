#ifndef WMESH_BLAS_H
#define WMESH_BLAS_H

#include "wmesh-types.h"

#ifdef WMESH_OPEN_BLAS
#include </usr/lib/openblas/include/f77blas.h>
#define BLAS_dcopy 	dcopy_
#define BLAS_dgemm 	dgemm_
#define BLAS_ddot 	ddot_
#define LAPACK_dgesv 	dgesv_
#else

#ifdef WMESH_MKL_BLAS

#if WMESH_ILP64
#define MKL_ILP64 1
#endif

#include "mkl.h"

#define BLAS_dcopy 	dcopy
#define BLAS_dgemm 	dgemm
#define BLAS_ddot  	ddot
#define LAPACK_dgesv 	dgesv

#else
#error choose WMESH_MKL_BLAS or WMESH_OPEN_BLAS
#endif

#endif

#endif
