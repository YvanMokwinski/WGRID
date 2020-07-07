//
// treat 2d
//
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-math.hpp"
#include "wmesh-status.h"
#include "wmesh.hpp"
#include "wmesh_utils.hpp"
#include "wmesh-blas.h"
#include "bms.h"

#include <chrono>
#include <iostream>
#include <math.h>
#include "bms_templates.hpp"

#if 0
template <wmesh_int_t 	ALPHA_,
	  wmesh_int_t 	BETA_,
	  wmesh_int_t 	N_>
struct AA
{
  wmesh_status_t bms_template_jacobip(wmesh_int_t 			x_n_,
				      const double * __restrict__  	x_,
				      wmesh_int_t  			x_ld_,
				      double *  __restrict__ 		y_,
				      wmesh_int_t  			y_ld_,
				      wmesh_int_t 			work_n_,			   
				      double *  __restrict__ 		work_)
  {
    WMESH_CHECK_POSITIVE(x_n_);
    WMESH_CHECK(work_n_  >= 2*x_n_);
    WMESH_CHECK_POINTER(x_);
    WMESH_CHECK_POINTER(y_);
    WMESH_CHECK_POINTER(work_);
    
    static constexpr double
      r2 = double(2.0),
      r3 = double(3.0),
      r1 = double(1.0),
      r0 = double(0.0);
  
    static constexpr const double
      ab  = ALPHA_ + BETA_,
      ab1 = ALPHA_+ BETA_ + 1,
      a1  = ALPHA_ + 1,
      b1  = BETA_ + 1;

    static constexpr const double gamma0 = Pow2<double>(ab1)*Gamma<double>(a1)*Gamma<double>(b1)/Factorial<double>(ab1);
    
    double
      aold   = r0,
      anew   = r0,
      bnew   = r0,
      h1     = r0,
      gamma1 = r0;
  
    // Initial values P_0(x) and P_1(x)
    const double y0 = r1 / wmesh_math<double>::xsqrt(gamma0);
    for (wmesh_int_t i=0;i<x_n_;++i)
      {
	y_[i * y_ld_] = y0;
      }
  
    wmesh_int_t n1=1;
    double * __restrict__ pii 	= work_;
    double * __restrict__  pi 	= work_ + x_n_;
    if (N_>0)
      {
	BLAS_dcopy(&x_n_,y_,&y_ld_,pii,&n1);
	//      auto pii = y;
	gamma1 = (a1)*(b1)/(ab+3.0)*gamma0;
	for (wmesh_int_t i=0;i<x_n_;++i)
	  {      
	    y_[i*y_ld_] = ((ab+r2)*x_[i]/r2 + (ALPHA_-BETA_)/r2) / wmesh_math<double>::xsqrt(gamma1);
	  }
	if (N_>1)
	  {

	    //
	    // pi = y
	    //
	    BLAS_dcopy(&x_n_,y_,&y_ld_,pi,&n1);

	    // Repeat value in recurrence.
	    aold = r2 / (r2+ab) * wmesh_math<double>::xsqrt((a1)*(b1)/(ab+r3));	  
	    // Forward recurrence using the symmetry of the recurrence.
	    for (int i=1; i<=(N_-1); ++i)
	      {
		h1 = r2*i+ab;
		double ri = static_cast<double>(i);
		anew = r2/(h1+r2)*wmesh_math<double>::xsqrt( (ri+1.0)*(ri+ab1)*(ri+a1)*(ri+b1)/(h1+r1)/(h1+r3));
		bnew = -(ALPHA_*ALPHA_-BETA_*BETA_) / ( h1*(h1+r2) );
	      
		for (wmesh_int_t i=0;i<x_n_;++i)
		  {
		    y_[i*y_ld_] = (x_[i]-bnew) * pi[i] - aold*pii[i];
		  }
		for (wmesh_int_t i=0;i<x_n_;++i)
		  {
		    y_[i*y_ld_] *= r1/anew;
		  }
	      
		//	      y = (x-bnew) * pi - aold*pii;
		//	      BLAS_daxpy(&x_n_,&s,y_,&n1);
		//	      y *= r1/anew;
		//	      double s = r1 / anew;
		//	      BLAS_dscal(&x_n_,&s,y_,&n1);
		aold = anew;
		BLAS_dcopy(&x_n_,pi,&n1,pii,&n1);
		BLAS_dcopy(&x_n_,y_,&y_ld_,pi,&n1);
		//  pii = pi;
		//  pi = y;
	      }	  
	  }     
      }
    return WMESH_STATUS_SUCCESS;

  }
};


template <wmesh_int_t  ALPHA_,wmesh_int_t BETA_,wmesh_int_t N_>
struct AA
{
  wmesh_status_t bms_template_jacobip(double 				scal_,
				      wmesh_int_t 			x_n_,
				      const double * __restrict__  	x_,
				      wmesh_int_t  			x_inc_,
				      double *  __restrict__ 		y_,
				      wmesh_int_t  			y_inc_,
				      wmesh_int_t 			work_n_,			   
				      double *  __restrict__ 		work_)
  {

    double * __restrict__ y1 	= work_;
    double * __restrict__ y2 	= work_ + x_n_;
    double * __restrict__ tmp 	= work_ + x_n_ * 2;
    work_ 			= work_ + x_n_ * 3;
    work_n_ 			-= x_n_ * 3;

    static constexpr wmesh_int_t c	= ( 2*N_*(N_ + ALPHA_ + BETA_)*(2*N_ + ALPHA_ + BETA_ - 2) );
    static constexpr wmesh_int_t c10 	= ( 2 * N_ + ALPHA_ + BETA_ - 1);
    static constexpr wmesh_int_t c11	= ( 2 * N_ + ALPHA_ + BETA_)*( 2 * N_ + ALPHA_ + BETA_ - 2);
    static constexpr wmesh_int_t c12 	= ( ALPHA_*ALPHA_-BETA_*BETA_);
    static constexpr wmesh_int_t c2     = 2 * (N_ + ALPHA_ - 1)*(N_ + BETA_ - 1)*(2*N_ + ALPHA_ + BETA_);

    static constexpr wmesh_int_t scal0  = c2 / c;
    static constexpr wmesh_int_t scal1  = (c10 * c11) / c;
    static constexpr wmesh_int_t scal2  = (c10 * c12) / c;

    

    for (wmesh_int_t i=0;i<x_n_;++i)
      {
	tmp[i] = scal2;
      }    
    double sc = scal1;
    daxpy(&sc,x_,&x_inc_,work_,&n1);

    AA<ALPHA_,BETA_,0>::bms_template_jacobi(x_n_,
					    x_,
					    x_inc_,
					    y0,
					    1,
					    work_n_,
					    work_);
    
    AA<ALPHA_,BETA_,1>::bms_template_jacobi(x_n_,
					    x_,
					    x_inc_,
					    y1,
					    1,
					    work_n_,
					    work_);

    for (wmesh_int_t k=2;k<=N_;++k)
      {
	for (wmesh_int_t i=0;i<x_n_;++i)
	  {
	    y_[y_inc_* i] = (a1 * ( a2 * tmp[i] + a3) * y1[i] + a4 * y0[i]) / a0;
	  }
        
	BLAS_dcopy(&x_n_,y1,&n1,y0,&n1);
	BLAS_dcopy(&x_n_,y_,&n1,y1,&n1);
      }
    
    
    AA<ALPHA_,BETA_,N_-1>::bms_template_jacobi(x_n_,
					       x_,
					       x_inc_,
					       y_,
					       y_inc_,
					       work_n_,
					       work_);
    
    
    

    dscal(&x_n_,&scal0,y0,&n1);



    
    WMESH_CHECK_POSITIVE(x_n_);
    WMESH_CHECK_POINTER(x_);
    WMESH_CHECK_POINTER(y_);
    
    static constexpr double
      r1 = double(1.0);
    
    for (wmesh_int_t i=0;i<x_n_;++i)
      {
	y_[i * y_ld_] = r1;
      }
    return 0;
  }

};

template <wmesh_int_t  ALPHA_,wmesh_int_t BETA_>
struct AA<ALPHA_,BETA_,0>
{
  wmesh_status_t bms_template_jacobip(double 				scal_,
				      wmesh_int_t 			x_n_,
				      const double * __restrict__  	x_,
				      wmesh_int_t  			x_ld_,
				      double *  __restrict__ 		y_,
				      wmesh_int_t  			y_ld_,
				      wmesh_int_t 			work_n_,			   
				      double *  __restrict__ 		work_)
  {
    
    WMESH_CHECK_POSITIVE(x_n_);
    WMESH_CHECK_POINTER(x_);
    WMESH_CHECK_POINTER(y_);
    
    for (wmesh_int_t i=0;i<x_n_;++i)
      {
	y_[i * y_ld_] = scal_;
      }
    return 0;
  }

};


template <wmesh_int_t  ALPHA_,wmesh_int_t BETA_>
struct AA<ALPHA_,BETA_,1>
{
  wmesh_status_t bms_template_jacobip(double scal_,
				      wmesh_int_t 			x_n_,
				      const double * __restrict__  	x_,
				      wmesh_int_t  			x_inc_,
				      double *  __restrict__ 		y_,
				      wmesh_int_t  			y_inc_,
				      wmesh_int_t 			work_n_,			   
				      double *  __restrict__ 		work_)
  {
    
    WMESH_CHECK_POSITIVE(x_n_);
    WMESH_CHECK_POINTER(x_);
    WMESH_CHECK_POINTER(y_);
    
    static constexpr double
      r1 = double(1.0);
    static constexpr double
      r2 = double(2.0);

    for (wmesh_int_t i=0;i<x_n_;++i)
      {
	y_[y_inc_ * i] = scal_ * ( (ALPHA_+1) + (ALPHA_+BETA_+2) * (x_[xinc_ * i] - r1 ) / r2 );
      }
    return 0;
  }
};



#if 0
template <wmesh_int_t 	ALPHA_,
	  wmesh_int_t 	BETA_,
	  wmesh_int_t 	N_>
wmesh_status_t bms_template_jacobip(wmesh_int_t 			x_n_,
				    const double * __restrict__  	x_,
				    wmesh_int_t  			x_ld_,
				    double *  __restrict__ 		y_,
				    wmesh_int_t  			y_ld_,
				    wmesh_int_t 			work_n_,			   
				    double *  __restrict__ 		work_)
{
};

#endif
#endif
  


#if 0
void bms_basis_monomial_triangle(cst_pI 	degree,
				 cst_pI 	n,
				 pR 		r,
				 cst_pI 	roff_,
				 cst_pR 	p,
				 cst_pI 	poff_,
				 pR 		rwork,
				 cst_pI 	rwork_n,
				 pI 		err_)
{
  I i,q,k,j;
  err_[0] = (I)0;
  if (degree[0]==0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (R)1.0;
    }
  else if (degree[0]==1)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)1.0;
	  r[i*roff_[0]+1] = x;
	  r[i*roff_[0]+2] = y;
	}
    }
  else if (degree[0]==2)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)1.0;
	  r[i*roff_[0]+1] = x;
	  r[i*roff_[0]+2] = y;
	  r[i*roff_[0]+3] = x*x;
	  r[i*roff_[0]+4] = y*x;
	  r[i*roff_[0]+5] = y*y;
	}
    }
  else if (degree[0]==3)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x = p[i];
	  const R y = p[poff_[0]+i];
	  r[i*roff_[0]] = (R)1.0;
	  r[i*roff_[0]+1] = x;
	  r[i*roff_[0]+2] = y;
	  r[i*roff_[0]+3] = x*x;
	  r[i*roff_[0]+4] = y*x;
	  r[i*roff_[0]+5] = y*y;
	  r[i*roff_[0]+6] = r[i*roff_[0]+3]*x;
	  r[i*roff_[0]+7] = r[i*roff_[0]+3]*y;
	  r[i*roff_[0]+8] = r[i*roff_[0]+5]*x;
	  r[i*roff_[0]+9] = r[i*roff_[0]+5]*y;
	}
    }
  else if (degree[0]==4)
    {
      for (i=0;i<n[0];++i)
	{
	  const R x 	= p[i];
	  const R y 	= p[poff_[0]+i];
	  r[i*roff_[0]] 	= (R)1.0;
	  r[i*roff_[0]+1] 	= x;
	  r[i*roff_[0]+2] 	= y;
	  r[i*roff_[0]+3] 	= x*x;
	  r[i*roff_[0]+4] 	= y*x;
	  r[i*roff_[0]+5] 	= y*y;
	  r[i*roff_[0]+6] 	= r[i*roff_[0]+3]*x;
	  r[i*roff_[0]+7] 	= r[i*roff_[0]+3]*y;
	  r[i*roff_[0]+8] 	= r[i*roff_[0]+5]*x;
	  r[i*roff_[0]+9] 	= r[i*roff_[0]+5]*y;
	  r[i*roff_[0]+10] 	= r[i*roff_[0]+6]*x;
	  r[i*roff_[0]+11] 	= r[i*roff_[0]+6]*y;
	  r[i*roff_[0]+12] 	= r[i*roff_[0]+3]*r[i*roff_[0]+5];
	  r[i*roff_[0]+13] 	= r[i*roff_[0]+9]*x;
	  r[i*roff_[0]+14] 	= r[i*roff_[0]+9]*y;
	}
    }
  else
    {
      for (i=0;i<n[0];++i)
	{
	  const R x 		= p[i];
	  const R y 		= p[poff_[0]+i];
	  r[i*roff_[0]] 	= (R)1.0;
	  r[i*roff_[0]+1] 	= x;
	  r[i*roff_[0]+2] 	= y;
	  r[i*roff_[0]+3] 	= x*x;
	  r[i*roff_[0]+4] 	= y*x;
	  r[i*roff_[0]+5] 	= y*y;
	  r[i*roff_[0]+6] 	= r[i*roff_[0]+3]*x;
	  r[i*roff_[0]+7] 	= r[i*roff_[0]+3]*y;
	  r[i*roff_[0]+8] 	= r[i*roff_[0]+5]*x;
	  r[i*roff_[0]+9] 	= r[i*roff_[0]+5]*y;
	  r[i*roff_[0]+10] 	= r[i*roff_[0]+6]*x;
	  r[i*roff_[0]+11] 	= r[i*roff_[0]+6]*y;
	  r[i*roff_[0]+12] 	= r[i*roff_[0]+3]*r[i*roff_[0]+5];
	  r[i*roff_[0]+13] 	= r[i*roff_[0]+9]*x;
	  r[i*roff_[0]+14] 	= r[i*roff_[0]+9]*y;
	}
      for (k=5;k<=degree[0];++k)
	for (q = (k*(k+1))/((I)2),j=0;j<=k;++j)
	  for (i=0;i<n[0];++i)
	    r[i*roff_[0]+q+j] = nsPOW(p[i],((R)k-j)) * nsPOW(p[poff_[0]+i],((R)j));    
    }
}



void mkS_canonic_tetra(cst_mpsint 	degree,
			     cst_mpsint 	n,
			     mpsreal 		r,
			     cst_mpsint 	roff_,
			     cst_mpsreal 	p,
			     cst_mpsint 	poff_,
			     mpsreal 		rwork,
			     cst_mpsint 	rwork_n,
			     mpsint 		err_)
{
  nsINT i;
  err_[0] = (nsINT)0;
  if (degree[0]==0)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)1.0;
    }
  else if (degree[0]==1)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)1.0;
      nsblas_dcopy(n,(nsREAL*)&p[0],&__vmps_blas_negal1,&r[1],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[poff_[0]],&__vmps_blas_negal1,&r[2],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[2*poff_[0]],&__vmps_blas_negal1,&r[3],roff_);      
    }
  else if (degree[0]==2)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)1.0;
      nsblas_dcopy(n,(nsREAL*)&p[0],&__vmps_blas_negal1,&r[1],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[poff_[0]],&__vmps_blas_negal1,&r[2],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[2*poff_[0]],&__vmps_blas_negal1,&r[3],roff_);            
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+4] = p[i]*p[i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+5] = p[i]*p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+6] = p[i]*p[2*poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+7] = p[poff_[0]+i]*p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+8] = p[poff_[0]+i]*p[2*poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+9] = p[2*poff_[0]+i]*p[2*poff_[0]+i];
    }
  else  if (degree[0]==3)
    {
      for (i=0;i<n[0];++i)
	r[i*roff_[0]] = (nsREAL)1.0;
      nsblas_dcopy(n,(nsREAL*)&p[0],&__vmps_blas_negal1,&r[1],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[poff_[0]],&__vmps_blas_negal1,&r[2],roff_);
      nsblas_dcopy(n,(nsREAL*)&p[2*poff_[0]],&__vmps_blas_negal1,&r[3],roff_);            

      for (i=0;i<n[0];++i)
	r[i*roff_[0]+4] = p[i]*p[i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+5] = p[i]*p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+6] = p[i]*p[2*poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+7] = p[poff_[0]+i]*p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+8] = p[poff_[0]+i]*p[2*poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+9] = p[2*poff_[0]+i]*p[2*poff_[0]+i];

      for (i=0;i<n[0];++i)
	{
	  double r = p[i];
	  double s = p[poff_[0]+i];
	  double t = p[2*poff_[0]+i];
	  
	  r[i*roff_[0]+4] = r*r;
	  r[i*roff_[0]+5] = r*s;
	  r[i*roff_[0]+6] = r*t;
	  r[i*roff_[0]+7] = s*s;
	  r[i*roff_[0]+8] = s*t;
	  r[i*roff_[0]+9] = t*t;
	  
	  for (i=0;i<n[0];++i)
	    r[i*roff_[0]+10] = p[i]*p[i]*p[i];
	}
      
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+10] = p[i]*p[i]*p[i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+11] = p[i]*p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+12] = p[i]*p[2*poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+13] = p[poff_[0]+i]*p[poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+14] = p[poff_[0]+i]*p[2*poff_[0]+i];
      for (i=0;i<n[0];++i)
	r[i*roff_[0]+15] = p[2*poff_[0]+i]*p[2*poff_[0]+i];

      
    }
}


//
// TEMPLATE DEFINITION
//
template<typename T>
wmesh_status_t
wfe_shape_eval_edge_canonic
(

 wmesh_int_t 		rst_storage_,
 wmesh_int_t 		rst_m_,
 wmesh_int_t 		rst_n_,
 const T * 		rst_v_,
 wmesh_int_t 		rst_ld_,
 
 wmesh_int_t storagee_,
 wmesh_int_t me_,
 wmesh_int_t ne_,
 T* e_,
 wmesh_int_t lde_,
 size_t work_n_,
 void* work_)
{
  
  //  if (WFE::storage_t::is_invalid(storagec_)) return WFE::status_t::invalid_enum;
  //  if (WFE::storage_t::is_invalid(storagee_)) return WFE::status_t::invalid_enum;
  
  WMESH_CHECK_POSITIVE(rst_m_);
  WMESH_CHECK_POSITIVE(rst_n_);
  WMESH_CHECK( rst_ld_ >= rst_m_ );
  WMESH_CHECK_POINTER( rst_v_);

  WMESH_CHECK_POSITIVE(me_);
  WMESH_CHECK_POSITIVE(ne_);
  WMESH_CHECK( lde_ >= me_ );
  WMESH_CHECK_POINTER( e_);

  bool e_is_block = (storagee_ == WCOMMON_STORAGE_BLOCK);
  bool c_is_block = (storagec_ == WCOMMON_STORAGE_BLOCK);
  
  wmesh_int_t mc 	= c_is_block ? rst_m_ : rst_n_;
  wmesh_int_t nc 	= c_is_block ? rst_n_ : rst_m_;
  
  wmesh_int_t me 	= e_is_block ? me_ : ne_;
  wmesh_int_t ne 	= e_is_block ? ne_ : me_;

  wmesh_int_t ndofs =  mc;
  
  WMESH_CHECK(nc == 1);
  WMESH_CHECK(mc == me);
  // WMESH_CHECK(ne == ndofs);

  if (e_is_block)
    {
      for (wmesh_int_t i=0;i<mc;++i)
	{
	  e_[i] = 1.0;
	}
    }
  else
    {
      for (wmesh_int_t i=0;i<mc;++i)
	{
	  e_[i * lde_] = 1.0;
	}	
    }
  
  wmesh_int_t option =
    (c_is_block && e_is_block)    ? 1 :
    (c_is_block && !e_is_block)   ? 2 :
    (!c_is_block && e_is_block)   ? 3 :
    (!c_is_block && !e_is_block)  ? 4 : 0;
  
  if (1 == ndofs)
    {
      return WFE::status_t::success;
    }
  
  switch(option)
    {
    case 1: 
      {    
	for (wmesh_int_t i=0;i<mc;++i)
	  {
	    for (wmesh_int_t j=0;j<nc;++j)
	      {
		e_[ (1 + j) * lde_ + i] = rst_v_[j * rst_ld_ + i];
	      }
	  }
	break;
      }
    case 2: 
      {
	for (wmesh_int_t i=0;i<mc;++i)
	  {
	    for (wmesh_int_t j=0;j<nc;++j)
	      {
		e_[i * lde_ + 1 + j] = = rst_v_[j * rst_ld_ + i];
	      }
	  }
	break;
      }
    case 3: 
      {    
	for (wmesh_int_t i=0;i<mc;++i)
	  {
	    for (wmesh_int_t j=0;j<nc;++j)
	      {
		e_[ (1 + j) * lde_ + i] = rst_v_[i * ldc_ + j];
	      }
	  }
	break;
      }
    case 4: 
      {
	for (wmesh_int_t i=0;i<mc;++i)
	  {
	    for (wmesh_int_t j=0;j<nc;++j)
	      {
		e_[i* lde_ + 1 + j] = rst_v_[i * ldc_ + j];
	      }
	  }
	break;
      }
    }

  if (dim == 1)
    {
      switch(option)
	{
	case 1:
	case 3:
	  {
	
	    for (wmesh_int_t d=2;d < ndofs;++d)
	      {
		for (wmesh_int_t i=0;i<mc;++i)
		  {
		    e_[lde_ * d + i] = e_[lde_ + i] * e_[lde_ * (d - 1) + i];
		  }
	      }
	  
	    break;
	  }
	
	case 2:
	case 4:
	  {	  
	    for (wmesh_int_t d=2;d < ndofs;++d)
	      {
		for (wmesh_int_t i=0;i<mc;++i)
		  {
		    e_[lde_ * i + d] = e_[lde_ * i + 1] * e_[lde_ * i + (d - 1)];
		  }
	      }
	    break;
	  }
	}
    }
  else if (dim == 2)
    {
      for (wmesh_int_t i=0;i<n;++i)
	{
	}
      
      //
      // 1 x y x*x y*x y*y
      //
      // 3 0
      // 2 1
      // 1 2
      // 0 1
      //
      // 1 x y z x*x x*y x*z y*y y*z z*z      x*x*x x*x*y x*x*z x*y*y x*y*z x*z*z
      //      
      for (int i=3;i>=0;--i)
	{
	  for (int j=0;j<3-i;++j)
	    {
	      
	    }
	}
      
      // 3 0 0
      // 2 1 0
      // 2 0 1
      // 1 2 0
      // 1 1 1      
      // 1 0 2      
      // 0 3 0
      // 0 2 1
      // 0 1 2
      // 0 0 3
      //
      switch(option)
	{
	case 1:
	case 3:
	  {	
	    for (wmesh_int_t d=2;d < ndofs;++d)
	      {
		for (wmesh_int_t i=0;i<mc;++i)
		  {
		    e_[lde_ * d + i] = e_[lde_ + i] * e_[lde_ * (d - 1) + i];
		  }
	      }
	  
	    break;
	  }
	
	case 2:
	case 4:
	  {	  
	    for (wmesh_int_t d=2;d < ndofs;++d)
	      {
		for (wmesh_int_t i=0;i<mc;++i)
		  {
		    e_[lde_ * i + d] = e_[lde_ * i + 1] * e_[lde_ * i + (d - 1)];
		  }
	      }
	    break;
	  }
	}
    }
  else  if (dim == 3)
    {
      switch(option)
	{
	case 1:
	case 3:
	  {
	
	    for (wmesh_int_t d=2;d < ndofs;++d)
	      {
		for (wmesh_int_t i=0;i<mc;++i)
		  {
		    e_[lde_ * d + i] = e_[lde_ + i] * e_[lde_ * (d - 1) + i];
		  }
	      }
	  
	    break;
	  }
	
	case 2:
	case 4:
	  {	  
	    for (wmesh_int_t d=2;d < ndofs;++d)
	      {
		for (wmesh_int_t i=0;i<mc;++i)
		  {
		    e_[lde_ * i + d] = e_[lde_ * i + 1] * e_[lde_ * i + (d - 1)];
		  }
	      }
	    break;
	  }
	}

    }
    
  return WFE::status_t::success;
}




template<typename T>
wmesh_status_t  bms_shape_nodal_calculate(wmesh_int_t 			degree_,
					  wmesh_int_t			rst_storage_,
					  wmesh_int_t 			rst_m_,
					  wmesh_int_t			rst_n_,
					  const T * __restrict__ 	rst_x_,
					  wmesh_int_t			rst_ld_)
{
  
  WMESH_CHECK(WCOMMON_STORAGE_UNKNOWN != rst_storage_);
  WMESH_CHECK(rst_m_ > 0);
  WMESH_CHECK(rst_n_ > 0);
  WMESH_CHECK(rst_ld_ >= rst_m_);
  WMESH_CHECK_POINTER(rst_x_);

  wmesh_int_t dim   = (WCOMMON_STORAGE_INTERLEAVE == rst_storage_) ? rst_m_ : rst_n_;
  wmesh_int_t ndofs = (WCOMMON_STORAGE_INTERLEAVE == rst_storage_) ? rst_n_ : rst_m_;

  size_t required_memory = wfe_shape_eval_memsize_edge_lagrange(self_,neval);
  if (work_n_ < required_memory)
    {
      return WFE::status_t::invalid_work_size;
    }

  //
  // Get the local dofs coordinates.
  //
  
  T * __restrict__ vandermonde 	= (T*__restrict__)work_;
  T * __restrict__ canonic	= vandermonde + ndofs * ndofs;  
  wmesh_int_p ipiv		= (wmesh_int_p)(canonic + ndofs * neval);

  

  
  const wmesh_int_t nshapes = ndofs;
  WFE::status_t status =  wfe_shape_eval_edge_canonic(&shape_eval_edge_canonic,
						      WFE::storage_t::block,
						      nshapes,
						      1,
						      lc,
						      nshapes,
						      WFE::storage_t::interleave,
						      nshapes,
						      nshapes,
						      vandermonde,
						      nshapes,
						      work_n,
						      work);
  

  status =  wfe_shape_eval_edge_canonic(&shape_eval_edge_canonic,
					storagec_,
					mc_,
					nc_,
					c_,
					ldc_,
					WFE::storage_t::interleave,
					nshapes,
					neval,
					canonic,
					nshapes,
					work_n,
					work);
  

  wmesh_int_t info_lapack;
  gesv(&nshapes,
       &neval,
       vandermonde,
       &nshapes,
       ipiv,
       canonic,
       &nshapes,						     
       &info_lapack);
  if (info_lapack!=0)
    {
      printf("info lapack error\n");
    }
  
  //
  // The result is interleave
  //
  if (storagee_ == WFE::storage_t::block)
    {
      //
      // transpose
      //
      for (wmesh_int_t j=0;j<neval;++j)
	{
	  for (wmesh_int_t i=0;i<ndofs;++i)
	    {
	      e_[i*lde_ + j] = canonic[j * ndofs + i];
	    }	    
	}	
    }
  else
    {
      //
      // Copy
      //
      for (wmesh_int_t j=0;j<neval;++j)
	{
	  for (wmesh_int_t i=0;i<ndofs;++i)
	    {
	      e_[j*lde_ + i] = canonic[j * ndofs + i];
	    }	    
	}
    }
  
  return WMESH_STATUS_SUCCESS;

};

template <typename T>
wmesh_status_t bms_basis_nodal(wmesh_int_t     		cell_type_,

			       wmesh_int_t		rst_storage_,
			       wmesh_int_t 		rst_m_,
			       wmesh_int_t		rst_n_,
			       const T * __restrict__ 	rst_x_,
			       wmesh_int_t		rst_ld_,
			       
			       wmesh_int_t		eval_storage_,
			       wmesh_int_t		eval_m_,
			       wmesh_int_t 		eval_n_,
			       T *__restrict__ 		eval_,
			       wmesh_int_t 		eval_ld_)
{     
  static constexpr  T s_zero 	= ((T)0.0);
  static constexpr  T s_one 	= ((T)1.0);
  if (cell_type_==0)
    {
      for (wmesh_int_t k=0;k<coo_n_;++k)
	{
	  const auto
	    r = coo_x_[coo_ld_*k+0],
	    s = coo_x_[coo_ld_*k+1],
	    t = coo_x_[coo_ld_*k+2];
	  
	  x_[k*x_ld_+0] = s_one - (r+s+t);
	  x_[k*x_ld_+1] = r;
	  x_[k*x_ld_+2] = s;
	  x_[k*x_ld_+3] = t;
	}		  
    }
  else if (cell_type_==1)
    {
      for (wmesh_int_t k=0;k<coo_n_;++k)
	{
	  const auto
	    r = coo_x_[coo_ld_*k+0],
	    s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	    
	    const auto
	      tscale = (t < s_one) ? s_one / (s_one - t) : s_zero,
	      lts = s_one - t - s,
	      ltr = s_one - t - r;
	    
	    x_[k*x_ld_+0] = ltr  * lts * tscale;
	    x_[k*x_ld_+1] = r * lts * tscale;
	    x_[k*x_ld_+2] = r * s * tscale;
	    x_[k*x_ld_+3] = ltr * s  * tscale;
	    x_[k*x_ld_+4] = t;
	  }		  
      }
    else if (cell_type_==2)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	   
	    const auto
	      l = one - (r+s),
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = l * lt;
	    x_[k*x_ld_+1] = r * lt;
	    x_[k*x_ld_+2] = s * lt;
	    x_[k*x_ld_+3] = l * t;
	    x_[k*x_ld_+4] = r * t;
	    x_[k*x_ld_+5] = s * t;
	  }		  
      }
    else if (cell_type_==3)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];

	    const auto
	      lr = one - r,
	      ls = one - s,
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = lr* ls * lt;
	    x_[k*x_ld_+1] = r * ls * lt;
	    x_[k*x_ld_+2] = r *  s * lt;
	    x_[k*x_ld_+3] = lr*  s * lt;
	    x_[k*x_ld_+4] = lr* ls * t;
	    x_[k*x_ld_+5] =  r* ls * t;
	    x_[k*x_ld_+6] =  r*  s * t;
	    x_[k*x_ld_+7] = lr*  s * t;
	  }
      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
  return WMESH_STATUS_SUCCESS;
}
#endif
using namespace std::chrono;
extern "C"
{

  


  wmesh_status_t bms_djacobip(wmesh_int_t 	alpha_,
			      wmesh_int_t 	beta_,
			      wmesh_int_t 	N_,
			      wmesh_int_t 	x_n_,
			      const double * 	x_,
			      wmesh_int_t  	x_ld_,
			      double * 		y_,
			      wmesh_int_t  	y_ld_,
			      wmesh_int_t 	work_n_,			   
			      double * 		work_)
  
  {

    
    return bms_jacobip(alpha_,
		       beta_,
		       N_,
		       x_n_,
		       x_,
		       x_ld_,
		       y_,
		       y_ld_,
		       work_n_,			   
		       work_);
  }
  
#if 0
  wmesh_status_t bms_basis_eval(wmesh_int_t     cell_type_,

				wmesh_int_t	rst_storage_,
				wmesh_int_t 	rst_m_,
				wmesh_int_t 	rst_n_,
				double * 	rst_x_,
				wmesh_int_t 	rst_ld_,
				
				wmesh_int_t	eval_storage_,
				wmesh_int_t 	eval_m_,
				wmesh_int_t 	eval_n_,
				double * 	eval_,
				wmesh_int_t 	eval_ld_)
  {     
    static const double s_zero 	= ((double)0.0);
    static const double s_one 	= ((double)1.0);
    static constexpr double one = 1.0;
    if (cell_type_==0)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	    
	    x_[k*x_ld_+0] = s_one - (r+s+t);
	    x_[k*x_ld_+1] = r;
	    x_[k*x_ld_+2] = s;
	    x_[k*x_ld_+3] = t;
	  }		  
      }
    else if (cell_type_==1)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	    
	    const auto
	      tscale = (t < s_one) ? s_one / (s_one - t) : s_zero,
	      lts = s_one - t - s,
	      ltr = s_one - t - r;
	    
	    x_[k*x_ld_+0] = ltr  * lts * tscale;
	    x_[k*x_ld_+1] = r * lts * tscale;
	    x_[k*x_ld_+2] = r * s * tscale;
	    x_[k*x_ld_+3] = ltr * s  * tscale;
	    x_[k*x_ld_+4] = t;
	  }		  
      }
    else if (cell_type_==2)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	   
	    const auto
	      l = one - (r+s),
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = l * lt;
	    x_[k*x_ld_+1] = r * lt;
	    x_[k*x_ld_+2] = s * lt;
	    x_[k*x_ld_+3] = l * t;
	    x_[k*x_ld_+4] = r * t;
	    x_[k*x_ld_+5] = s * t;
	  }		  
      }
    else if (cell_type_==3)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];

	    const auto
	      lr = one - r,
	      ls = one - s,
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = lr* ls * lt;
	    x_[k*x_ld_+1] = r * ls * lt;
	    x_[k*x_ld_+2] = r *  s * lt;
	    x_[k*x_ld_+3] = lr*  s * lt;
	    x_[k*x_ld_+4] = lr* ls * t;
	    x_[k*x_ld_+5] =  r* ls * t;
	    x_[k*x_ld_+6] =  r*  s * t;
	    x_[k*x_ld_+7] = lr*  s * t;
	  }
      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
    return WMESH_STATUS_SUCCESS;
  }
#endif


  
  wmesh_status_t wmeshspace_generate_coodofs(wmeshspace_t * 	self_,
					     wmesh_int_t 	coo_storage_,
					     wmesh_int_t 	coo_m_,
					     wmesh_int_t 	coo_n_,
					     double * 		coo_,
					     wmesh_int_t 	coo_ld_)
  {    
    static constexpr double r0 = 0.0;
    static constexpr double r1 = 1.0;
    wmesh_status_t 	status;
    const wmesh_t * 	mesh 	= self_->get_mesh();
    wmesh_int_t         topodim = mesh->m_topology_dimension;
    wmesh_int_t 	coo_m  	= mesh->m_coo_m;
    
    wmesh_int_t 	num_types;
    wmesh_int_t 	elements[4];
    //    double * 		refevals[4] {};
    double 		cell_xyz[32];
    wmesh_int_t 	cell_xyz_ld = coo_m;
    wmesh_mat_t<double> eval_basis[4];

    status = bms_topodim2elements(topodim,
				  &num_types,
				  elements);
    WMESH_STATUS_CHECK(status);        
    for (wmesh_int_t l=0;l<num_types;++l)
      {
	if (mesh->m_c2n.m_n[l]>0)
	  {
	    wmesh_int_t element = (topodim==3) ? (4+l) : ((topodim==2)? (2+l) : 1+l);
	    wmesh_int_t element_num_nodes;

	    status = bms_elements_num_nodes(1,&element,&element_num_nodes);
	    WMESH_STATUS_CHECK(status);        
	    
	    wmesh_t* 		rmacro			= self_->get_refinement_pattern(l);
	    double * 	rmacro_coo 		= wmesh_get_coo(rmacro);
	    wmesh_int_t 	rmacro_coo_ld 		= rmacro->m_coo_ld;
	    wmesh_int_t 	rmacro_num_nodes 	= rmacro->m_num_nodes;

	    const wmesh_int_t  	mat_rmacro_coo_storage = WMESH_STORAGE_INTERLEAVE;
	    wmesh_mat_t<double> mat_rmacro_coo;
	    wmesh_mat_t<double>::define(&mat_rmacro_coo,topodim,rmacro_num_nodes,rmacro_coo,rmacro_coo_ld);

	    const wmesh_int_t  	eval_basis_storage = WMESH_STORAGE_INTERLEAVE;
	    wmesh_mat_t<double>::alloc(&eval_basis[l], element_num_nodes, rmacro_num_nodes);
	    wmesh_shape_t shape_element;
	    wmesh_shape_def(&shape_element,element,WMESH_SHAPE_FAMILY_LAGRANGE,1);
	    wmesh_shape_calculate_eval(shape_element,
				       mat_rmacro_coo_storage,
				       mat_rmacro_coo,
				       eval_basis_storage,
				       eval_basis[l]);

#if 0
	    if ((topodim==0)&&(l==0))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*rmacro_num_nodes*4);
		for (wmesh_int_t k=0;k<rmacro_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    double lr = ( ((double)1.0)-r );
		    double ls = ( ((double)1.0)-s );
		    b[k*4+0] = lr* ls;
		    b[k*4+1] = r * ls;
		    b[k*4+2] = r *  s;
		    b[k*4+3] = lr * s;
		  }		  
		refevals[l] = b;		  
	      }
	    
	    if ((topodim==2)&&(l==0))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*rmacro_num_nodes*3);
		for (wmesh_int_t k=0;k<rmacro_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    
		    b[k*3+0] = ((double)1.0)-(r+s);
		    b[k*3+1] = r;
		    b[k*3+2] = s;
		  }		  
		refevals[l] = b;		  
	      }
	    
	    if ((topodim==2)&&(l==1))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*rmacro_num_nodes*4);
		for (wmesh_int_t k=0;k<rmacro_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    double lr = ((double)1.0)-r;
		    double ls = ((double)1.0)-s;
		    b[k*4+0] = lr* ls;
		    b[k*4+1] = r * ls;
		    b[k*4+2] = r *  s;
		    b[k*4+3] = lr * s;
		  }		  
		refevals[l] = b;		  
	      }

	    if ((topodim==3)&&(l==0))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*rmacro_num_nodes*4);
		for (wmesh_int_t k=0;k<rmacro_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    double t = rmacro_coo[rmacro_coo_ld*k+2];
		    
		    b[k*4+0] = ((double)1.0)-(r+s+t);
		    b[k*4+1] = r;
		    b[k*4+2] = s;
		    b[k*4+3] = t;
		  }		  
		refevals[l] = b;		  
	      }

	    if ((topodim==3)&&(l==1))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->get_refinement_pattern(l)->m_num_nodes*5);
		for (wmesh_int_t k=0;k<self_->get_refinement_pattern(l)->m_num_nodes;++k)
		  {
		    double r = rmacro_coo[rmacro_coo_ld*k+0];
		    double s = rmacro_coo[rmacro_coo_ld*k+1];
		    double t = rmacro_coo[rmacro_coo_ld*k+2];
		    // std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;
		    if (t<1.0)
		      {
			  
			double h = 1.0-t;
			  
			double phir0 = (h-r)/h;
			double phir1 = r/h;
		      
			double phis0 = (h-s)/h;
			double phis1 = s/h;
			  
			b[k*5+0] = phir0 * phis0 * (1.0-t);
			b[k*5+1] = phir1 * phis0 * (1.0-t);
			b[k*5+2] = phir1 * phis1 * (1.0-t);
			b[k*5+3] = phir0 * phis1 * (1.0-t);
			b[k*5+4] = t;

			b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			b[k*5+2] = (r * s) / (1.0-t);
			b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			b[k*5+4] = t;
			  
#if 0
			r = 2.0*r-1.0;
			s = 2.0*s-1.0;
			t = 2.0*t-1.0;
			double c = (1.0-t)/2.0;
			b[k*5+0] = ( (c-t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+1] = ( (c+t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+2] = ( (c+t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+3] = ( (c-t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+4] = (t+1.0)/2.0;
			  

#else
			b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			b[k*5+2] = (r * s) / (1.0-t);
			b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			b[k*5+4] = t;
#endif
		      }
		    else
		      {
			b[k*5+0] = 0.0;
			b[k*5+1] = 0.0;
			b[k*5+2] = 0.0;
			b[k*5+3] = 0.0;
			b[k*5+4] = 1.0;			  
		      }
		      
		    //
		    // Split into tetrahedra.
		    //
		    // 4
		    //
		    // 3 2
		    // 0 1
		    //
#if 0
		    if (r > s)
		      {
			b[k*5+0] = (1.0-r-s-t);
			b[k*5+1] = r;
			b[k*5+2] = s;
			b[k*5+3] = 0.0;
			b[k*5+4] = t;
		      }
		    else
		      {
			b[k*5+0] = (1.0-r-s-t);
			b[k*5+1] = 0.0;
			b[k*5+2] = r;
			b[k*5+3] = s;
			b[k*5+4] = t;
		      }
#endif
		      
		      
#if 0
		    b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+4] = (1.0+t)/2.0;
		    r=hr;
		    s=hs;
		    t=ht;
		    b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+4] = (1.0+t)/2.0;
		    b[k*5+0] = (1.0-r) * (1.0-s) * (1.0-t);
		    b[k*5+1] = (r) * (1.0-s) * (1.0-t);
		    b[k*5+2] = (r) * s * (1.0-t);
		    b[k*5+3] = (1.0-r) * s * (1.0-t);
		    b[k*5+4] = t;
#endif

#if 0
		    b[k*5+5] = ( one - r -t  )* ( one - s -t ) / (4.0 * (1-t));
		    b[k*5+1] = 
		      b[k*5+2] = ( ( r )* ( s ) * ( one - t ) );
		    b[k*5+3] = ( one - r )* ( s ) * ( one - t );
		    b[k*5+4] = ( one - r )* ( one - s ) * ( t );
#endif
		  }

		refevals[l] = b;		  
#if 0
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->get_refinement_pattern(l)->m_num_nodes*6);
		for (wmesh_int_t k=0;k<self_->get_refinement_pattern(l)->m_num_nodes;++k)
		  {
		    double r = self_->get_refinement_pattern(l)->m_coo[3*k+0];
		    double s = self_->get_refinement_pattern(l)->m_coo[3*k+1];
		    double t = self_->get_refinement_pattern(l)->m_coo[3*k+2];
		    double one=1.0;
		    b[k*8+0] = ( one - r - s) * ( one - t );
		    b[k*8+1] = r * ( one - t );
		    b[k*8+2] = s * ( one - t );
		    b[k*8+4] = ( one - r - s) * ( t );
		    b[k*8+5] = r * ( t );
		  }		  
		refevals[l] = b;
#endif
	      }
	      
	    if ((topodim==3)&&(l==2))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->get_refinement_pattern(l)->m_num_nodes*6);
		for (wmesh_int_t k=0;k<self_->get_refinement_pattern(l)->m_num_nodes;++k)
		  {
		    double r = self_->get_refinement_pattern(l)->m_coo[3*k+0];
		    double s = self_->get_refinement_pattern(l)->m_coo[3*k+1];
		    double t = self_->get_refinement_pattern(l)->m_coo[3*k+2];
		    double one=1.0;
		    b[k*6+0] = ( one - r - s) * ( one - t );
		    b[k*6+1] = r * ( one - t );
		    b[k*6+2] = s * ( one - t );
		    b[k*6+3] = ( one - r - s) * ( t );
		    b[k*6+4] = r * ( t );
		    b[k*6+5] = s * ( t );
		  }		  
		refevals[l] = b;		  
	      }
	    if ((topodim==3)&&(l==3))
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->get_refinement_pattern(l)->m_num_nodes*8);
		//		  std::cout << " " << std::endl;
		for (wmesh_int_t k=0;k<self_->get_refinement_pattern(l)->m_num_nodes;++k)
		  {
		    double r = self_->get_refinement_pattern(l)->m_coo[3*k+0];
		    double s = self_->get_refinement_pattern(l)->m_coo[3*k+1];
		    double t = self_->get_refinement_pattern(l)->m_coo[3*k+2];
		    //		      std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;
		    double one=1.0;
		    b[k*8+0] = ( one - r )* ( one - s ) * ( one - t );
		    b[k*8+1] = ( ( r )* ( one - s ) * ( one - t ) );
		    b[k*8+2] = ( ( r )* ( s ) * ( one - t ) );
		    b[k*8+3] = ( one - r )* ( s ) * ( one - t );
		    b[k*8+4] = ( one - r )* ( one - s ) * ( t );
		    b[k*8+5] = ( ( r )* ( one - s ) * ( t ) );
		    b[k*8+6] = ( ( r )* ( s ) * ( t ) );
		    b[k*8+7] = ( one - r )* ( s ) * ( t );
		  }

		refevals[l] = b;		  
	      }
#endif	      
	  }

      }
    
    wmesh_int_t rw_n = 0;
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
      {
	wmesh_int_t k = self_->m_c2d.m_m[l]*topodim;
	rw_n = (rw_n < k) ? k : rw_n;
      }



    
    double * rw = (double*)malloc(sizeof(double)*rw_n);
    
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
      {
	//	double * 	refeval = refevals[l];
	
	//
	// Local c2d.
	//	    
	// wmesh_int_t c2d_n	= self_->m_c2d.m_n[l];
	wmesh_int_t c2d_m 		= self_->m_c2d.m_m[l];
	wmesh_int_t c2d_ld 		= self_->m_c2d.m_ld[l];
	wmesh_int_p c2d_v 		= self_->m_c2d.m_data + self_->m_c2d.m_ptr[l];
	
	//
	// Local c2n.
	//
	wmesh_int_t c2n_n 		= mesh->m_c2n.m_n[l];
	wmesh_int_t c2n_m 		= mesh->m_c2n.m_m[l];
	wmesh_int_t c2n_ld		= mesh->m_c2n.m_ld[l];
	wmesh_int_p c2n_v 		= mesh->m_c2n.m_data + mesh->m_c2n.m_ptr[l];
	
	for (wmesh_int_t j=0;j<c2n_n;++j)
	  {
	    
	    //
	    // Get the coordinates of the cell.
	    //		
	    for (wmesh_int_t i=0;i<c2n_m;++i)
	      {
		wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
		for (wmesh_int_t k=0;k<mesh->m_coo_m;++k)
		  {
		    cell_xyz[cell_xyz_ld * i + k] = mesh->m_coo[mesh->m_coo_ld * idx + k];
		  }
	      }

#if 0
	    //
	    // Now 
	    //
	    if (mat_rmacro_coo_storage == WMESH_STORAGE_INTERLEAVE)
	      {
		wmesh_mat_gemm(static_cast<T>(1),ref_c,eval_basis,static_cast<T>(0),mat_rmacro_coo);
	      }
	    else
	      {
		
	      }
#endif

	    //	    wmesh_mat_gemm(static_cast<T>(1),cell_xyz,eval_basis[l],static_cast<T>(0),physical_coordinates);
	    
	    dgemm("N",
	      "N",
	      &coo_m ,
	      &c2d_m,
	      &c2n_m ,
	      &r1,
	      cell_xyz,
	      &cell_xyz_ld,
	      eval_basis[l].v,
	      &eval_basis[l].ld,
		  &r0,
		  rw,
		  &coo_m);

	    //
	    // Copy.
	    //
	    if (WMESH_STORAGE_INTERLEAVE == coo_storage_)
	      {
		for (wmesh_int_t i=0;i<c2d_m;++i)
		  {
		    wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		    for (wmesh_int_t k = 0;k<coo_m;++k)
		      {
			coo_[coo_ld_ * idx + k] = rw[coo_m * i + k];
		      }
		  }
	      }
	    else
	      {
		for (wmesh_int_t i=0;i<c2d_m;++i)
		  {
		    wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		    for (wmesh_int_t k = 0;k<coo_m;++k)
		      {
			coo_[coo_ld_ * k + idx] = rw[coo_m * i + k];
		      }
		  }
	      }	    
	  }
      }
    return WMESH_STATUS_SUCCESS;
  };










  wmesh_status_t wmeshspace_sublinearmesh(wmeshspace_t * 	self_,
					  wmesh_t ** 		mesh__)
  {    
    wmesh_status_t 	status;
    const wmesh_t * 	mesh 	= self_->get_mesh();
    wmesh_int_t         topodim = mesh->m_topology_dimension;
    wmesh_int_t 	coo_m  	= mesh->m_coo_m;
    
    wmesh_int_t 	num_types;
    wmesh_int_t 	elements[4];


    status = bms_topodim2elements(topodim,
				  &num_types,
				  elements);
    WMESH_STATUS_CHECK(status);    
    
    wmesh_int_t coo_dofs_m  	= coo_m;
    wmesh_int_t coo_dofs_n  	= self_->get_ndofs();
    wmesh_int_t coo_dofs_ld 	= coo_dofs_m;
    double * 	coo_dofs 	= (double*)malloc(sizeof(double) * coo_dofs_n * coo_dofs_ld);
    
    status =  wmeshspace_generate_coodofs(self_,
					  WMESH_STORAGE_INTERLEAVE,
					  coo_dofs_m,
					  coo_dofs_n,
					  coo_dofs,
					  coo_dofs_ld);
    WMESH_STATUS_CHECK(status);

    //
    // GENERATE SUBLINEAR CONNECTIVITY
    //
    {
      
      wmesh_int_t c2n_size = num_types;
      wmesh_int_t c2n_ptr[5];
      wmesh_int_t c2n_m[4]{};
      wmesh_int_t c2n_n[4]{};	
      wmesh_int_t c2n_ld[4]{};
      
      status = bms_elements_num_nodes(num_types,
				      elements,
				      c2n_m);
      WMESH_STATUS_CHECK(status);
      
      for (int i=0;i<num_types;++i)
	{
	  const wmesh_t * pattern = self_->get_refinement_pattern(i);
	  if (pattern!=nullptr)
	    {
	      if (i==1 && topodim == 3)
		{
		  c2n_n[0] += mesh->m_c2n.m_n[i] * pattern->m_c2n.m_n[0];
		  c2n_n[1] += mesh->m_c2n.m_n[i] * pattern->m_c2n.m_n[1];
		}
	      else
		{
		  c2n_n[i] = mesh->m_c2n.m_n[i] * pattern->m_num_cells;
		}
	    }
	}

      status = wmesh_int_sparsemat_init(c2n_size,
					c2n_ptr,
					c2n_m,
					c2n_n,
					c2n_ld);
      WMESH_STATUS_CHECK(status);      
      wmesh_int_p c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[c2n_size]);
      if (!c2n_v)
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);      
	}

      wmesh_int_t mxdofs = 0;
      for (int i=0;i<num_types;++i)
	{
	  const wmesh_t * pattern = self_->get_refinement_pattern(i);
	  //	  if (c2n_n[i] > 0) mxdofs = (self_->m_patterns[i]->m_num_nodes > mxdofs) ? self_->m_patterns[i]->m_num_nodes : mxdofs;
	  if (nullptr != pattern)
	    {
	      mxdofs = (pattern->m_num_nodes > mxdofs) ? pattern->m_num_nodes : mxdofs;
	    }
	}
      
      wmesh_int_p dofs = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);
      wmesh_int_p lidx = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);
      
      printf("generate connectivity.\n");
      wmesh_int_t idx[4]{0,0,0,0};
      for (wmesh_int_t l=0;l<num_types;++l)
	{
	  auto ref_c2n = &self_->get_refinement_pattern(l)->m_c2n;
	  
	  wmesh_int_t ncells = self_->m_c2d.m_n[l];
	  for (wmesh_int_t j=0;j<ncells;++j)
	    {
	      
	      //
	      // extract dofs.
	      //
	      for (wmesh_int_t i=0;i<self_->m_c2d.m_m[l];++i)
		{
		  dofs[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[l] + j * self_->m_c2d.m_ld[l] + i];
		}

		
	      //
	      // For each.
	      //
	      for (wmesh_int_t ref_cell_type=0;ref_cell_type<ref_c2n->m_size;++ref_cell_type)
		{
		  wmesh_int_t ref_ncells = ref_c2n->m_n[ref_cell_type];
		  wmesh_int_t ref_nnodes_per_cell = ref_c2n->m_m[ref_cell_type];
		  for (wmesh_int_t sj=0;sj<ref_ncells;++sj)
		    {
			
		      for (wmesh_int_t si=0;si<ref_nnodes_per_cell;++si)
			{
			  lidx[si] = ref_c2n->m_data[ref_c2n->m_ptr[ref_cell_type] + sj * ref_c2n->m_ld[ref_cell_type]  + si ] - 1;
			}
			
		      for (wmesh_int_t si=0;si<ref_c2n->m_m[ref_cell_type];++si)
			{
			  c2n_v[ c2n_ptr[ref_cell_type] + c2n_ld[ref_cell_type] * idx[ref_cell_type] + si] = dofs[lidx[si]];
			}
		      ++idx[ref_cell_type];		       
		    }
		}
	    }
	}
      
      //
      // Define the mesh.
      //
      status =  wmesh_def(mesh__,
			  mesh->m_topology_dimension,				 
			  c2n_size,
			  c2n_ptr,
			  c2n_m,
			  c2n_n,
			  c2n_v,
			  c2n_ld,
			  coo_dofs_m,
			  coo_dofs_n,
			  coo_dofs,
			  coo_dofs_ld);

      //
      // Copy the dof codes.
      // It's a bit tricky, need to be changed.
      //
      const wmesh_int_t ndofs = self_->get_ndofs();
      for (wmesh_int_t i=0;i<ndofs;++i)
	{
	  mesh__[0]->m_n_c.v[i] = self_->m_dof_codes[i];
	}
      
      WMESH_STATUS_CHECK(status);
    }
    return WMESH_STATUS_SUCCESS;
  }

};




