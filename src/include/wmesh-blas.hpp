#pragma once
#include "wmesh-blas.h"

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
