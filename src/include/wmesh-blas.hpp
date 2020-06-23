#pragma once
#include "wmesh-blas.h"



template <typename T>
inline void xgemm(const char*transa,const char*transb,const_wmesh_int_p m,const_wmesh_int_p n,const_wmesh_int_p k,const T* alpha, const T* a,const_wmesh_int_p lda,const T *b,const_wmesh_int_p ldb, const T* beta,  T *c,const_wmesh_int_p ldc);

template <>
inline void xgemm<float>(const char*transa,
			 const char*transb,
			 const_wmesh_int_p m,
			 const_wmesh_int_p n,
			 const_wmesh_int_p k,
			 const float* alpha,
			 const float* a,
			 const_wmesh_int_p lda,
			 const float *b,
			 const_wmesh_int_p ldb,
			 const float* beta,
			 float *c,
			 const_wmesh_int_p ldc)
{
  sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}


template <>
inline void xgemm<double>(const char*transa,
			 const char*transb,
			 const_wmesh_int_p m,
			 const_wmesh_int_p n,
			 const_wmesh_int_p k,
			 const double* alpha,
			 const double* a,
			 const_wmesh_int_p lda,
			 const double *b,
			 const_wmesh_int_p ldb,
			 const double* beta,
			 double *c,
			 const_wmesh_int_p ldc)
{
  dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}




template <typename T>
inline void xger(const_wmesh_int_p m,const_wmesh_int_p n,const T* alpha,const T* x,const_wmesh_int_p xinc,const T *y,const_wmesh_int_p yinc,T *a,const_wmesh_int_p lda);

template <>
inline void xger<float>(const_wmesh_int_p m,const_wmesh_int_p n,const float* alpha,const float* x,const_wmesh_int_p xinc,const float *y,const_wmesh_int_p yinc,float *a,const_wmesh_int_p lda)
{
  sger(m,n,alpha,x,xinc,y,yinc,a,lda);
}
template <>
inline void xger<double>(const_wmesh_int_p m,const_wmesh_int_p n,const double* alpha,const double* x,const_wmesh_int_p xinc,const double *y,const_wmesh_int_p yinc,double *a,const_wmesh_int_p lda)
{
  dger(m,n,alpha,x,xinc,y,yinc,a,lda);
}

template <typename T>
inline void xcopy(const_wmesh_int_p n,const T* x,const_wmesh_int_p xinc,T *y,const_wmesh_int_p yinc);
template <>
inline void xcopy<double>(const_wmesh_int_p n,const double *x,const_wmesh_int_p xinc,double *y,const_wmesh_int_p yinc)
{
  dcopy(n,x,xinc,y,yinc);
}
template <>
inline void xcopy<float>(const_wmesh_int_p n,const float *x,const_wmesh_int_p xinc,float *y,const_wmesh_int_p yinc)
{
  scopy(n,x,xinc,y,yinc);
}


template <typename T>
inline void xscal(const_wmesh_int_p n,const T* x,T *y,const_wmesh_int_p yinc);
template <>
inline void xscal<double>(const_wmesh_int_p n,const double *x,double *y,const_wmesh_int_p yinc)
{
  dscal(n,x,y,yinc);
}
template <>
inline void xscal<float>(const_wmesh_int_p n,const float *x,float *y,const_wmesh_int_p yinc)
{
  sscal(n,x,y,yinc);
}





template <typename T>
inline void xgesv(const_wmesh_int_p n,const_wmesh_int_p nrhs,T* a,const_wmesh_int_p a_ld,wmesh_int_p perm,T *b,const_wmesh_int_p b_ld_,wmesh_int_p info_lapack_);
template <>
inline void xgesv<float>(const_wmesh_int_p n,const_wmesh_int_p nrhs,float* a,const_wmesh_int_p a_ld,wmesh_int_p perm,float *b,const_wmesh_int_p b_ld,wmesh_int_p info_lapack)
{

  
  sgesv(n,
	nrhs,
	a,
	a_ld,
	perm,
	b,
	b_ld,
	info_lapack);

}
template <>
inline void xgesv<double>(const_wmesh_int_p n,const_wmesh_int_p nrhs,double* a,const_wmesh_int_p a_ld,wmesh_int_p perm,double *b,const_wmesh_int_p b_ld,wmesh_int_p info_lapack)
{

  
  dgesv(n,
	nrhs,
	a,
	a_ld,
	perm,
	b,
	b_ld,
	info_lapack);

}
