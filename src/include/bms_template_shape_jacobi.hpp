#pragma once
#include "bms.h"
#include "bms_templates.hpp"

#if 0
template <wmesh_int_t  ALPHA_,wmesh_int_t BETA_,wmesh_int_t N_, typename T>
struct bms_template_jacobi_t
{
  static wmesh_status_t eval(double 				scal_,
			     wmesh_int_t 			x_n_,
			     const T * __restrict__  		x_,
			     wmesh_int_t  			x_inc_,
			     T *  __restrict__ 			y_,
			     wmesh_int_t  			y_inc_,
			     wmesh_int_t 			work_n_,			   
			     T *  __restrict__ 			work_)
  {
    T * __restrict__ y1 	= work_;
    T * __restrict__ y2 	= work_ + x_n_;
    T * __restrict__ tmp 	= work_ + x_n_ * 2;
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
    static constexpr wmesh_int_t n1     = 1;

    for (wmesh_int_t i=0;i<x_n_;++i)
      {
	tmp[i] = scal2;
      }
    
    T sc = scal1;
    daxpy(&sc,x_,&x_inc_,work_,&n1);

    bms_template_jacobi_t<ALPHA_,BETA_,0,T>::eval(x_n_,
						x_,
						x_inc_,
						y0,
						1,
						work_n_,
						work_);

    bms_template_jacobi_t<ALPHA_,BETA_,1,T>::eval(x_n_,
						x_,
						x_inc_,
						y1,
						1,
						work_n_,
						work_);
    
    for (wmesh_int_t k = 2;k <= N_;++k)
      {
	
	for (wmesh_int_t i = 0;i < x_n_;++i)
	  {
	    y_[y_inc_* i] = (a1 * ( a2 * tmp[i] + a3) * y1[i] + a4 * y0[i]) / a0;
	  }
        
	BLAS_dcopy(&x_n_,y1,&n1,y0,&n1);
	BLAS_dcopy(&x_n_,y_,&n1,y1,&n1);
      }

    bms_template_jacobi_t<ALPHA_,BETA_,N_-1,T>::eval(x_n_,
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
    
    static constexpr T
      r1 = static_cast<T>(1);
    
    for (wmesh_int_t i=0;i<x_n_;++i)
      {
	y_[i * y_ld_] = r1;
      }
    return WMESH_STATUS_SUCCESS;
  }

  
};

template <wmesh_int_t  ALPHA_,wmesh_int_t BETA_,typename T>
struct bms_template_jacobi_t<ALPHA_,BETA_,0,T>
{
  wmesh_status_t eval(double 				scal_,
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
    return WMESH_STATUS_SUCCESS;
  }

};


template <wmesh_int_t  ALPHA_,wmesh_int_t BETA_,typename T>
struct bms_template_jacobi_t<ALPHA_,BETA_,1, T>
{
  static wmesh_status_t eval(double 				scal_,
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
    
    static constexpr T
      r1 = T(1.0);
    static constexpr T
      r2 = T(2.0);

    for (wmesh_int_t i=0;i<x_n_;++i)
      {
	y_[y_inc_ * i] = scal_ * ( (ALPHA_+1) + (ALPHA_+BETA_+2) * (x_[xinc_ * i] - r1 ) / r2 );
      }
    
    return WMESH_STATUS_SUCCESS;
  }
};


#endif

template<wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_template_shape_jacobi
{  
  static wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 	c_,
				     wmesh_int_t 	c_ld_,
		      
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,

			     
				     wmesh_int_t   	iw_n_,
				     wmesh_int_p   	iw_,
				     wmesh_int_t   	rw_n_,
				     T* 		rw_);
};

template<typename 	T>
struct bms_template_shape_jacobi<WMESH_ELEMENT_EDGE,T>
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 		c_,
				     wmesh_int_t 	c_ld_,
		      
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,
		      
				     wmesh_int_t   	iw_n_,
				     wmesh_int_p   	iw_,
				     wmesh_int_t   	rw_n_,
				     T* 		rw_)
  {
    wmesh_status_t status;
    const wmesh_int_t n 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
    const wmesh_int_t c_inc 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
    if (diff_[0] > 0)
      {
	if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
	  {
	    for (wmesh_int_t i=0;i<=degree_;++i)
	      {
		if (i>0)
		  {
		    status = bms_jacobip(1,
					  1,
					  i-1,		
					  n,
					  c_,
					  c_inc,
					  b_+i,
					  b_ld_,
					  rw_n_,
					  rw_);
		    WMESH_STATUS_CHECK(status);
#if 0
		    for (wmesh_int_t l=0;l<n;++l)
		      {
			b_[l*b_ld_ + i] *= static_cast<T>(i+1);
		      }
#endif
		  }
		else
		  {
		    for (wmesh_int_t l=0;l<n;++l)
		      {
			b_[l*b_ld_ + 0] = static_cast<T>(0);
		      }
		  }
	      }
	  }
	else
	  {
	    for (wmesh_int_t i=0;i<=degree_;++i)
	      {
		if (i>0)
		  {
		    status = bms_jacobip(1,
					  1,
					  i-1,		
					  n,
					  c_,
					  c_inc,
					  b_+b_ld_*i,
					  1,
					  rw_n_,
					  rw_);
		    WMESH_STATUS_CHECK(status);
#if 0
		    for (wmesh_int_t l=0;l<n;++l)
		      {
			b_[i*b_ld_ + l] *= static_cast<T>(i+1);
		      }
#endif
		  }
		else
		  {
		    for (wmesh_int_t l=0;l<n;++l)
		      {
			b_[0*b_ld_ + l] = static_cast<T>(0);
		      }
		  }
	      }
	  }
      }
    else
      {
	if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
	  {
	    for (wmesh_int_t i=0;i<=degree_;++i)
	      {
		status = bms_jacobip(0,
				      0,
				      i,		
				      n,
				      c_,
				      c_inc,
				      b_+i,
				      b_ld_,
				      rw_n_,
				      rw_);
		WMESH_STATUS_CHECK(status);
	      }
	  }
	else
	  {
	    for (wmesh_int_t i=0;i<=degree_;++i)
	      {
		status = bms_jacobip(0,
				      0,
				      i,		
				      n,
				      c_,
				      c_inc,
				      b_ + b_ld_ * i,
				      1,
				      rw_n_,
				      rw_);
		WMESH_STATUS_CHECK(status);
	      }
	  }
      }
    
    return WMESH_STATUS_SUCCESS;    
  }
};





template<typename 	T>
struct bms_template_shape_jacobi<WMESH_ELEMENT_QUADRILATERAL,T>
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 	c_,
				     wmesh_int_t 	c_ld_,
		      
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,
		      
				     wmesh_int_t   	iw_n_,
				     wmesh_int_p   	iw_,
				     wmesh_int_t   	rw_n_,
				     T* 		rw_)
  {
    wmesh_status_t status;
    const wmesh_int_t n 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
    const wmesh_int_t c_inc 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
    const T * cr 		= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? (c_ + 0) : (c_ + c_ld_ * 0);
    const T * cs 		= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? (c_ + 1) : (c_ + c_ld_ * 1);
    // phi(r,s) = p(r) q(s)
    // dr = p'(r) q(s)
    // ds = p(s) q'(s)
    T * tmpi = rw_;
    T * tmpj = rw_+n;
    rw_ += n*2;
    rw_n_ -= n*2;
    std::cout << "rw_n_ " << rw_n_ << std::endl;
#if 0
    {
      wmesh_int_t N = 100;
      T ax[101];
      T ay[101];
      FILE * f = fopen("info.0.txt","w");
      for (wmesh_int_t i=0;i<=N;++i)
	{
	  ax[i]=-1.0 + static_cast<T>(2 * i) / static_cast<T>(N);

	}

      status = bms_jacobip(1,
			   1,
			   2,		
			   N+1,
			   ax,
			   1,
			   ay,
			   1,
			   rw_n_,
			   rw_);
      
      for (wmesh_int_t i=0;i<=N;++i)
	{
	  fprintf(f,"%e %e\n",ax[i],ay[i]);
	}
      fclose(f);



      f = fopen("info.1.txt","w");
      for (wmesh_int_t i=0;i<=N;++i)
	{
	  ax[i]=-1.0 + static_cast<T>(2*i)/static_cast<T>(N);
	}
      status = bms_jacobip(2,
			   2,
			   1,		
			   N+1,
			   ax,
			   1,
			   ay,
			   1,
			   rw_n_,
			   rw_);
      for (wmesh_int_t i=0;i<=N;++i)
	{
	 
	  fprintf(f,"%e %e\n",ax[i],ay[i]);
	}
      fclose(f);

      //      exit(1);
    }
#endif
    if (diff_[0] == 0 && diff_[1] == 0)
      {
	//	fprintf(stdout,"EVSL F\n");
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {
	    status = bms_jacobip(0,
				 0,
				 i,		
				 n,
				 cr,
				 c_inc,
				 tmpi,
				 1,
				 rw_n_,
				 rw_);
	    WMESH_STATUS_CHECK(status);
	    
	    for (wmesh_int_t j=0;j<=degree_;++j)
	      {
		status = bms_jacobip(0,
				     0,
				     j,		
				     n,
				     cs,
				     c_inc,
				     tmpj,
				     1,
				     rw_n_,
				     rw_);
		WMESH_STATUS_CHECK(status);
		
		if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		  {
#if 0
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			fprintf(stdout,"%8.15e %8.15e %8.15e %8.15e\n",cr[c_inc*l],cs[c_inc*l],tmpi[l],tmpj[l]);
		      }
		    fprintf(stdout,"yyyyyyyyyyyyyy\n");
#endif
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*l + idx] = tmpi[l] * tmpj[l];
		      }
		  }
		else
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*idx + l] = tmpi[l] * tmpj[l];
		      }
		  }
		
		++idx;
	      }
	  }
      }
    else if (diff_[0] == 1)
      {
	//	fprintf(stdout,"EVSL DR\n");
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {
	    if (i>0)
	      {
		status = bms_jacobip(1,
				     1,
				     i-1,		
				     n,
				     cr,
				     c_inc,
				     tmpi,
				     1,
				     rw_n_,
				     rw_);
#if 1
		for (wmesh_int_t l =0;l<n;++l)
		  {
		    tmpi[l]  *= static_cast<T>(i+1) / static_cast<T>(2);
		  }
#endif
		WMESH_STATUS_CHECK(status);
	      }
	    else
	      {
		for (wmesh_int_t l =0;l<n;++l)
		  {
		    tmpi[l]  = static_cast<T>(0);
		  }
	      }
	    
	    for (wmesh_int_t j=0;j<=degree_;++j)
	      {
		status = bms_jacobip(0,
				      0,
				      j,		
				      n,
				      cs,
				      c_inc,
				      tmpj,
				      1,
				      rw_n_,
				      rw_);
	        WMESH_STATUS_CHECK(status);
		
		if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*l + idx] = tmpi[l] * tmpj[l];
		      }
		  }
		else
		  {

		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*idx + l] = tmpi[l] * tmpj[l];
		      }
		  }
		
		++idx;
	      }
	  }


      }
    else if (diff_[1] == 1)
      {
	// fprintf(stdout,"EVSL S\n");
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {

	    status = bms_jacobip(0,
				 0,
				 i,		
				 n,
				 cr,
				 c_inc,
				 tmpi,
				 1,
				 rw_n_,
				 rw_);
	    WMESH_STATUS_CHECK(status);
	  
	    
	    for (wmesh_int_t j=0;j<=degree_;++j)
	      {
		if (j>0)
		  {
		    status = bms_jacobip(1,
					 1,
					 j-1,		
					 n,
					 cs,
					 c_inc,
					 tmpj,
					 1,
					 rw_n_,
					 rw_);
#if 1
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			tmpj[l]  *= static_cast<T>(j+1) / static_cast<T>(2);
		      }
#endif
		  }
		else
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			tmpj[l]  = static_cast<T>(0);
		      }		    
		  }
		
		if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*l + idx] = tmpi[l] * tmpj[l];
		      }
		  }
		else
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*idx + l] = tmpi[l] * tmpj[l];
		      }
		  }
		++idx;
	      }
	  }

	
      }
    else      
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_SUCCESS);    
      }
    return WMESH_STATUS_SUCCESS;    
  }
};




template<typename 	T>
struct bms_template_shape_jacobi<WMESH_ELEMENT_HEXAHEDRON,T>
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 	c_,
				     wmesh_int_t 	c_ld_,
		      
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,
		      
				     wmesh_int_t   	iw_n_,
				     wmesh_int_p   	iw_,
				     wmesh_int_t   	rw_n_,
				     T* 		rw_)
  {
    wmesh_status_t status;
    const wmesh_int_t n 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
    const wmesh_int_t c_inc 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
    const T * cr 		= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? (c_ + 0) : (c_ + c_ld_ * 0);
    const T * cs 		= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? (c_ + 1) : (c_ + c_ld_ * 1);
    const T * ct 		= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? (c_ + 2) : (c_ + c_ld_ * 2);
    
    T * tmpi = rw_;
    T * tmpj = rw_+n;
    T * tmpk = rw_+n*2;
    rw_ += n*3;
    rw_n_ -= n*3;
    
    if (diff_[0] == 0 && diff_[1] == 0 && diff_[2] == 0)
      {
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {
	    status = bms_jacobip(0,
				 0,
				 i,		
				 n,
				 cr,
				 c_inc,
				 tmpi,
				 1,
				 rw_n_,
				 rw_);
	    
	    WMESH_STATUS_CHECK(status);
	    for (wmesh_int_t j=0;j<=degree_;++j)
	      {
		status = bms_jacobip(0,
				      0,
				      j,		
				      n,
				      cs,
				      c_inc,
				      tmpj,
				      1,
				      rw_n_,
				      rw_);
	    WMESH_STATUS_CHECK(status);

		for (wmesh_int_t k=0;k<=degree_;++k)
		  {
		    status = bms_jacobip(0,
					  0,
					  k,		
					  n,
					  ct,
					  c_inc,
					  tmpk,
					  1,
					  rw_n_,
					  rw_);
	    WMESH_STATUS_CHECK(status);
		    
		    if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		      {
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    b_[ b_ld_*idx + l] = tmpi[l] * tmpj[l] * tmpk[l];
			  }
		      }
		    else
		      {
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    b_[ b_ld_*l + idx] = tmpi[l] * tmpj[l]* tmpk[l];
			  }
		      }
		    
		    ++idx;
		  }
	      }
	  }
      }
    else if (diff_[0] == 1)
      {
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {
	    if (i>0)
	      {
		status = bms_jacobip(1,
				      1,
				      i-1,		
				      n,
				      cr,
				      c_inc,
				      tmpi,
				      1,
				     rw_n_,
				     rw_);
		WMESH_STATUS_CHECK(status);
		for (wmesh_int_t l =0;l<n;++l)
		      {
			tmpi[l]  *= (i+1) / 2;
		      }
	      }
	    else
	      {
		for (wmesh_int_t l =0;l<n;++l)
		  {
		    tmpi[l]  = static_cast<T>(0);
		  }
	      }
	    
	    for (wmesh_int_t j=0;j<=degree_;++j)
	      {
		status = bms_jacobip(0,
				      0,
				      j,		
				      n,
				      cs,
				      c_inc,
				      tmpj,
				      1,
				      rw_n_,
				      rw_);
	    WMESH_STATUS_CHECK(status);

		for (wmesh_int_t k=0;k<=degree_;++k)
		  {
		    status = bms_jacobip(0,
					  0,
					  k,		
					  n,
					  ct,
					  c_inc,
					  tmpk,
					  1,
					  rw_n_,
					  rw_);
	    WMESH_STATUS_CHECK(status);
		
		    if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		      {
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    b_[ b_ld_*idx + l] = tmpi[l] * tmpj[l] * tmpk[l];
			  }
		      }
		    else
		      {
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    b_[ b_ld_*l + idx] = tmpi[l] * tmpj[l]* tmpk[l];
			  }
		      }
		
		++idx;
		  }
	      }
	  }


      }
    else if (diff_[1] == 1)
      {
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {
	    status = bms_jacobip(0,
				  0,
				  i,		
				  n,
				  cr,
				  c_inc,
				  tmpi,
				  1,
				  rw_n_,
				  rw_);
	    WMESH_STATUS_CHECK(status);
	  
	    
	    for (wmesh_int_t j=0;j<=degree_;++j)
	      {
		if (j>0)
		  {
		    status = bms_jacobip(1,
					  1,
					  j-1,		
					  n,
					  cs,
					  c_inc,
					  tmpj,
					  1,
					  rw_n_,
					  rw_);
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			tmpj[l]  *= (j+1) / 2;
		      }
	    WMESH_STATUS_CHECK(status);
		  }
		else
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			tmpj[l]  = static_cast<T>(0);
		      }		    
		  }

		for (wmesh_int_t k=0;k<=degree_;++k)
		  {
		    status = bms_jacobip(0,
					  0,
					  k,		
					  n,
					  ct,
					  c_inc,
					  tmpk,
					  1,
					  rw_n_,
					  rw_);
	    WMESH_STATUS_CHECK(status);
		
		    if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		      {
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    b_[ b_ld_*idx + l] = tmpi[l] * tmpj[l] * tmpk[l];
			  }
		      }
		    else
		      {
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    b_[ b_ld_*l + idx] = tmpi[l] * tmpj[l]* tmpk[l];
			  }
		      }
		
		++idx;
		  }
	      }
	  }

      }
    else if (diff_[2] == 1)
      {
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {
	    status = bms_jacobip(0,
				 0,
				 i,		
				 n,
				 cr,
				 c_inc,
				 tmpi,
				 1,
				 rw_n_,
				 rw_);
	    WMESH_STATUS_CHECK(status);
	    
	    
	    for (wmesh_int_t j=0;j<=degree_;++j)
	      {
		status = bms_jacobip(0,
				     0,
				     j,		
				     n,
				     cs,
				     c_inc,
				     tmpj,
				     1,
				     rw_n_,
				     rw_);
		WMESH_STATUS_CHECK(status);
		
		for (wmesh_int_t k=0;k<=degree_;++k)
		  {
		    if (k>0)
		      {
			status = bms_jacobip(1,
					     1,
					     k-1,		
					     n,
					     cs,
					     c_inc,
					     tmpk,
					     1,
					     rw_n_,
					     rw_);
			WMESH_STATUS_CHECK(status);
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    tmpk[l]  *= (k+1) / 2;
			  }
		      }
		    else
		      {
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    tmpk[l]  = static_cast<T>(0);
			  }		    
		      }
		
		    if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		      {
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    b_[ b_ld_*idx + l] = tmpi[l] * tmpj[l] * tmpk[l];
			  }
		      }
		    else
		      {
			for (wmesh_int_t l =0;l<n;++l)
			  {
			    b_[ b_ld_*l + idx] = tmpi[l] * tmpj[l]* tmpk[l];
			  }
		      }
		    
		    ++idx;
		  }
	      }
	  }
	
      }
    else      
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_SUCCESS);    
      }
    return WMESH_STATUS_SUCCESS;    

  }
};


template<typename 	T>
struct bms_template_shape_jacobi<WMESH_ELEMENT_PYRAMID,T>
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 	c_,
				     wmesh_int_t 	c_ld_,
		      
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,
		      
				     wmesh_int_t   	iw_n_,
				     wmesh_int_p   	iw_,
				     wmesh_int_t   	rw_n_,
				     T* 		rw_)
  {
    return 0;
  }
};
  
template<typename 	T>
struct bms_template_shape_jacobi<WMESH_ELEMENT_WEDGE,T>
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 	c_,
				     wmesh_int_t 	c_ld_,
		      
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,
		      
				     wmesh_int_t   	iw_n_,
				     wmesh_int_p   	iw_,
				     wmesh_int_t   	rw_n_,
				     T* 		rw_)
  {
    return 0;
  }
};



template<typename 	T>
struct bms_template_shape_jacobi<WMESH_ELEMENT_TETRAHEDRON,T>
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 	c_,
				     wmesh_int_t 	c_ld_,
		      
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,
		      
				     wmesh_int_t   	iw_n_,
				     wmesh_int_p   	iw_,
				     wmesh_int_t   	rw_n_,
				     T* 		rw_)
  {
    return 0;
  }
};

//!
//! phi_ij = P_i^{0,0}( 2r / (1-s) - 1) * (1 - s)^i P_j^{2i+1,0}( 2s - 1 )
//!
//! f(g(x)) = f'(g())g'
//! A(r,s) = P_i^{0,0}( 2r / (1-s) - 1 )
//! B(r,s) = (1 - s)^i P_j^{2i+1,0}( 2s - 1 )
//!
//! dA/dr = ( (0 + 0 + i + 2) / 2 ) P_{i-1}^{0 + 1,0 + 1}(2r/(1-s) -1) * 2/(1-s) = (i+2) P_{i-1}^{1,1}( 2r/(1-s)-1 ) / (1-s)
//! dA/ds = ( (0 + 0 + i + 2) / 2 ) P_{i-1}^{0 + 1,0 + 1}(2r/(1-s) -1) * 2r/(1-s)/(1-s) = (i+2) P_{i-1}^{1,1}( 2r/(1-s)-1 ) r / (1-s) / (1-s)
//! dB/ds = (1-s)^{i-1} ( -P_j^{2i+1,0}( 2s - 1 ) + (1-s) (j + 2i + 3) P_j-1^{2i+2,1}( 2s - 1 ) )
//! dB/dr = 0
//!
//!
template<typename 	T>
struct bms_template_shape_jacobi<WMESH_ELEMENT_TRIANGLE,T>
{

  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 		c_,
				     wmesh_int_t 	c_ld_,
				     
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,
		      
				     wmesh_int_t   	iw_n_,
				     wmesh_int_p   	iw_,
				     wmesh_int_t   	rw_n_,
				     T* 		rw_)
  {


    wmesh_status_t status;
    const wmesh_int_t n 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_n_ : c_m_;
    const wmesh_int_t c_inc 	= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? c_ld_ : 1;
    const T * cr 		= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? (c_ + 0) : (c_ + c_ld_ * 0);
    const T * cs 		= (c_storage_ == WMESH_STORAGE_INTERLEAVE) ? (c_ + 1) : (c_ + c_ld_ * 1);
    for (wmesh_int_t i=0;i<n;++i)
      {
	std::cout << " " << cr[c_inc*i] << " " << cs[c_inc*i] << std::endl;
      }

    
    const wmesh_int_t cr_inc = c_inc;
    const wmesh_int_t cs_inc = c_inc;
    // phi(r,s) = p(r) q(s)
    // dr = p'(r) q(s)
    // ds = p(s) q'(s)
    T * tmpi = rw_;
    T * tmpj = rw_+n;
    T * tmpr = rw_ + 2*n;
    T * tmps = rw_ + 3*n;
    T * dtmpi = rw_ + 4*n;
    T * dtmpj = rw_ + 5*n;
    static constexpr T r1 = static_cast<T>(1);
    static constexpr T r2 = static_cast<T>(2);

    rw_ += n*6;
    rw_n_ -= n*6;
    std::cout << "rw_n_ " << rw_n_ << std::endl;
    
    for (wmesh_int_t l = 0;l < n;++l)
      {
	tmpr[l] = (cs[cs_inc*l] < r1) ? r2 * cr[cr_inc * l] / (r1 - cs[cs_inc*l]) - 1.0 : 1.0;
      }
    
    for (wmesh_int_t l = 0;l < n;++l)
      {
	tmps[l] = (r2 * cs[cr_inc * l] - r1);
      }
    
    //! phi_ij = P_i^{0,0}( 2r / (1-s) - 1) * (1 - s)^i P_j^{2i+1,0}( 2s - 1 )
    if (diff_[0] == 0 && diff_[1] == 0)
      {
	//	fprintf(stdout,"EVSL F\n");
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {
	    
	    
	    status = bms_jacobip(0,
				 0,
				 i,		
				 n,
				 tmpr,
				 1,
				 tmpi,
				 1,
				 rw_n_,
				 rw_);
	    WMESH_STATUS_CHECK(status);

	    for (wmesh_int_t j=0;j<=degree_-i;++j)
	      {
		status = bms_jacobip(2*i+1,
				     0,
				     j,		
				     n,
				     tmps,
				     1,
				     tmpj,
				     1,
				     rw_n_,
				     rw_);
		WMESH_STATUS_CHECK(status);
		for (wmesh_int_t k=1;k<=i;++k)
		  {
		    for (wmesh_int_t l = 0;l < n;++l)
		      {
			tmpj[l] *= (r1 - cs[cs_inc*l]);
		      }
		  }
		
		if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*l + idx] = tmpi[l] * tmpj[l];
		      }
		  }
		else
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*idx + l] = tmpi[l] * tmpj[l];
		      }
		  }
		
		++idx;
	      }
	  }
      }
    else if (diff_[0] == 1)
      {
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {
	    
	    if (i>0)
	      {	    
		status = bms_jacobip(1,
				     1,
				     i-1,		
				     n,
				     tmpr,
				     1,
				     tmpi,
				     1,
				     rw_n_,
				     rw_);
		WMESH_STATUS_CHECK(status);
		for (wmesh_int_t l =0;l<n;++l)
		  {
		    tmpi[l]  *= static_cast<T>(i+1);
		    std::cout << "cono " << tmpi[l] << std::endl;
		  }
		// / 2 * 2 / (r1 - cs[cs_inc*l])
	      }
	    else
	      {
		for (wmesh_int_t l =0;l<n;++l)
		  {
		    tmpi[l]  = static_cast<T>(0);
		  }
	      }
	    
	    for (wmesh_int_t j=0;j<=degree_-i;++j)
	      {
		status = bms_jacobip(2*i+1,
				     0,
				     j,		
				     n,
				     tmps,
				     1,
				     tmpj,
				     1,
				     rw_n_,
				     rw_);
		WMESH_STATUS_CHECK(status);
#if 0
		for (wmesh_int_t k=1;k<=i-1;++k) // i-1 rather than i because we didn't divide above.
		  {
		    for (wmesh_int_t l = 0;l < n;++l)
		      {
			tmpj[l] *= (r1 - cs[cs_inc*l]);
		      }
		  }

#endif		
		if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		  {
#if 0
		    if (i==1)
		      for (wmesh_int_t l =0;l<n;++l)
			{
			  std::cout << "roger " << dtmpi[l] << std::endl;
			}
#endif
		    for (wmesh_int_t l =0;l<n;++l)
		      {
#if 0
			if (i==0)
			  {
			    std::cout << "hello " << pow(r1-cs[cs_inc*l],i-1) << " " << pow(0.0,i-1) << " " << tmpi[l] << std::endl;
			  }
#endif
			b_[ b_ld_*l + idx] = tmpi[l] * pow(r1-cs[cs_inc*l],i-1) * tmpj[l];
			//			if (i==1)std::cout << " " << b_[ b_ld_*l + idx] << std::endl;
		      }
#if 0
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*l + idx] = tmpi[l] * tmpj[l];
		      }
#endif
		    //		    if (i==1)
		    //		    exit(1);
		  }
		else
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*idx + l] = tmpi[l] * tmpj[l];
		      }
		  }
		
		++idx;
	      }
	  }
      }
    else if (diff_[1] == 1)
      {
	wmesh_int_t idx = 0;
	for (wmesh_int_t i=0;i<=degree_;++i)
	  {
	    status = bms_jacobip(0,
				 0,
				 i,		
				 n,
				 tmpr,
				 1,
				 tmpi,
				 1,
				 rw_n_,
				 rw_);
	    if (i>0)
	      {	    
		status = bms_jacobip(1,
				     1,
				     i-1,		
				     n,
				     tmpr,
				     1,
				     dtmpi,
				     1,
				     rw_n_,
				     rw_);
		WMESH_STATUS_CHECK(status);
		for (wmesh_int_t l =0;l<n;++l)
		  {
		    dtmpi[l]  *= static_cast<T>(i+1) * cr[cr_inc*l]; // (r1-cs[cs_inc*l]);
		  }
	      }
	    else
	      {
		for (wmesh_int_t l =0;l<n;++l)
		  {
		    dtmpi[l]  = static_cast<T>(0);
		  }
	      }

	    
	    for (wmesh_int_t j=0;j<=degree_-i;++j)
	      {

		status = bms_jacobip(2*i+1,
				     0,
				     j,		
				     n,
				     tmps,
				     1,
				     tmpj,
				     1,
				     rw_n_,
				     rw_);
		if (j>0)
		  {	    
		    status = bms_jacobip(2*i+2,
					 1,
					 j-1,		
					 n,
					 tmps,
					 1,
					 dtmpj,
					 1,
					 rw_n_,
					 rw_);
		    WMESH_STATUS_CHECK(status);
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			dtmpj[l]  *= static_cast<T>( (2*i+1) + 0 + j + 1);
		      }
		  }
		else
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			dtmpj[l]  = static_cast<T>(0);
		      }
		  }
		//! phi_ij = P_i^{0,0}( 2r / (1-s) - 1) * (1 - s)^i P_j^{2i+1,0}( 2s - 1 )
		// phi(r,s) = a(r,s) * b(s) * q(s)
		// dr phi(r,s) = ( dr a(r,s) ) * b(s) * q(s)
		//
		// ds phi(r,s) = ( ds a(r,s) ) * b(s) * q(s) + a(r,s) * ( b'(s) * q(s) + b(s) * q'(s) )
		
		// p(r,s) * ( (1-s)^i  * q(s) )
		// p'(r,s)2*r/(1-s)^2 * ( (1-s)^i  * q(s) ) + p(r,s)     ( )
		
		if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		  {

		    for (wmesh_int_t l =0;l<n;++l)
		      {
			if (i<2)
			  {
			    std::cout << "rotrr " << dtmpi[l]  << " " << pow(r1-cs[cs_inc*l],i-2) << " " << tmpj[l] << std::endl;
			  }
			b_[ b_ld_*l + idx] = dtmpi[l] * pow(r1-cs[cs_inc*l],i-2) * tmpj[l] + tmpi[l] * ( -i * pow(r1-cs[cs_inc*l],i-1) * tmpj[l] + pow(r1-cs[cs_inc*l],i) * dtmpj[l]);
			//			if (i==1)std::cout << " " << b_[ b_ld_*l + idx] << std::endl;
		      }
		  }
		else
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*idx + l] = dtmpi[l] * pow(r1-cs[cs_inc*l],i-2) * tmpj[l] + tmpi[l] * ( -i * pow(r1-cs[cs_inc*l],i-1) * tmpj[l] + pow(r1-cs[cs_inc*l],i) * dtmpj[l]);
		      }
		  }

#if 0
		for (wmesh_int_t k=1;k<=i-2;++k) // i-1 because we didn't divide dtmpi a'(r,s) * b(s) * q(s)
		  {
		    for (wmesh_int_t l = 0;l < n;++l)
		      {
			dtmpi[l] *= r1 - cs[cs_inc*l];
		      }
		  }
		
		for (wmesh_int_t l = 0;l < n;++l)
		  {
		    dtmpi[l] *= tmpj[l];
		  }
		
		for (wmesh_int_t l = 0;l < n;++l)
		  {
		    tmpj[l] *= i;
		  }
		for (wmesh_int_t k=1;k<=i-1;++k) // q(s) * b'(s) 
		  {
		    for (wmesh_int_t l = 0;l < n;++l)
		      {
			tmpj[l] *= r1 - cs[cs_inc*l];
		      }
		  }
		for (wmesh_int_t k=1;k<=i;++k)  // q'(s) * b(s)
		  {
		    for (wmesh_int_t l = 0;l < n;++l)
		      {
			dtmpj[l] *= r1 - cs[cs_inc*l];
		      }
		  }
		for (wmesh_int_t l = 0;l < n;++l)
		  {
		    tmpi[l] *= tmpj[l] + dtmpj[l];
		  }

		if (b_storage_ == WMESH_STORAGE_INTERLEAVE)
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*l + idx] = dtmpi[l] *  tmpi[l];
		      }
		  }
		else
		  {
		    for (wmesh_int_t l =0;l<n;++l)
		      {
			b_[ b_ld_*idx + l] = dtmpi[l] * tmpi[l];
		      }
		  }
#endif		
		
		++idx;
	      }
	  }
	
      }
    else      
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);    
      }
    return WMESH_STATUS_SUCCESS;
  };
};


#if 0
template<wmesh_int_t 	ELEMENT_,
	 typename 	T>
struct bms_template_shape_jacobi
{  
  static inline  wmesh_status_t eval(wmesh_int_t 	degree_,
				     const_wmesh_int_p	diff_,				     
				     wmesh_int_t 	c_storage_,					  
				     wmesh_int_t 	c_m_,
				     wmesh_int_t 	c_n_,
				     const T * 	c_,
				     wmesh_int_t 	c_ld_,
		      
				     wmesh_int_t 	b_storage_,
				     wmesh_int_t 	b_m_,
				     wmesh_int_t 	b_n_,
				     T* 		b_,
				     wmesh_int_t 	b_ld_,
		      
				     wmesh_int_t   	iw_n_,
				     wmesh_int_p   	iw_,
				     wmesh_int_t   	rw_n_,
				     T* 		rw_)
  {
#ifdef FORWARD
#error FORWARD already defined
#endif
#define FORWARD					\
    diff_,					\
      c_storage_,				\
      c_m_,					\
      c_n_,					\
      c_,					\
      c_ld_,					\
						\
      b_storage_,				\
      b_m_,					\
      b_n_,					\
      b_,					\
      b_ld_,					\
    						\
      iw_n_,					\
      iw_,					\
      rw_n_,					\
      rw_

    switch(degree_)
      {
      case 0:
	{
	  return bms_template_shape_eval(FORWARD,
					 bms_template_shape_jacobi_splz<0,ELEMENT_,T>::basis);
	}
      case 1:
	{
	  return bms_template_shape_eval(FORWARD,
					 bms_template_shape_jacobi_splz<1,ELEMENT_,T>::basis);
	}
      default:
	{
	  return WMESH_STATUS_SUCCESS;
	}
      }

    return WMESH_STATUS_INVALID_ARGUMENT;
#undef FORWARD
    
  }

};
#endif
