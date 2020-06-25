#pragma once
#include "wmesh-blas.h"

template <typename T>
inline void xgemv(const char*transa,
			 const_wmesh_int_p m,
			 const_wmesh_int_p n,
			 const T* alpha,
			 const T* a,
			 const_wmesh_int_p lda,
			 const T *b,
			 const_wmesh_int_p incb,
			 const T* beta,
			 T *c,
			  const_wmesh_int_p incc);


template <>
inline void xgemv<float>(const char*transa,
			 const_wmesh_int_p m,
			 const_wmesh_int_p n,
			 const float* alpha,
			 const float* a,
			 const_wmesh_int_p lda,
			 const float *b,
			 const_wmesh_int_p incb,
			 const float* beta,
			 float *c,
			 const_wmesh_int_p incc)
{
  sgemv(transa,m,n,alpha,a,lda,b,incb,beta,c,incc);
}

template <>
inline void xgemv<double>(const char*transa,
			 const_wmesh_int_p m,
			 const_wmesh_int_p n,
			 const double* alpha,
			 const double* a,
			 const_wmesh_int_p lda,
			 const double *b,
			 const_wmesh_int_p incb,
			 const double* beta,
			 double *c,
			 const_wmesh_int_p incc)
{
  dgemv(transa,m,n,alpha,a,lda,b,incb,beta,c,incc);
}

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








template<typename T>
inline void wmesh_mat_gemm(T 				alpha_,
			   const wmesh_mat_t<T>&	a_,
			   const wmesh_mat_t<T>&	b_,
			   T 				beta_,
			   wmesh_mat_t<T>&		c_)
{
  xgemm("N",
	"N",
	&c_.m,
	&c_.n,
	&a_.n,
	&alpha_,
	a_.v,
	&a_.ld,
	b_.v,
	&b_.ld,
	&beta_,
	c_.v,
	&c_.ld);  
}

template<typename T>
inline wmesh_status_t wmesh_mat_gemm(const char * 		transa_,
				     const char * 		transb_,
				     T 				alpha_,
				     const wmesh_mat_t<T>&	a_,
				     const wmesh_mat_t<T>&	b_,
				     T 				beta_,
				     wmesh_mat_t<T>&		c_)
{
#ifndef NDEBUG
  wmesh_int_t am = (transa_[0] == 'N') ? a_.m : a_.n;
  wmesh_int_t an = (transa_[0] == 'N') ? a_.n : a_.m;
  wmesh_int_t bm = (transb_[0] == 'N') ? b_.m : b_.n;
  wmesh_int_t bn = (transb_[0] == 'N') ? b_.n : b_.m;
  wmesh_int_t cm = c_.m;
  wmesh_int_t cn = c_.n;
  WMESH_CHECK(an == bm);
  WMESH_CHECK(am == cm);
  WMESH_CHECK(bn == cn);
#endif
  xgemm(transa_,
	transb_,
	&c_.m,
	&c_.n,
	(transa_[0] == 'N') ? &a_.n : &a_.m,
	&alpha_,
	a_.v,
	&a_.ld,
	b_.v,
	&b_.ld,
	&beta_,
	c_.v,
	&c_.ld);
  return WMESH_STATUS_SUCCESS;
}




