
#include "bms.hpp"
#include "wmesh-math.hpp"
#include "wmesh-blas.h"
#include <iostream>
#ifndef NDEBUG
#include <iostream>
#endif
#define solve_params(_type) wmesh_int_t*n,_type*jacM,_type*wr,_type*wi,_type*vl,_type*vr,_type*work,wmesh_int_t*work_n_
template<typename T>
static void solve(solve_params(T));

template<>
void solve<double>(solve_params(double))
{
  static constexpr const  char trN[1] = {'N'};
  static constexpr const char trV[1] = {'V'};
  wmesh_int_t ldvl=1,info;
  LAPACK_dgeev(trN,
	trV,
	n,
	jacM,
	n, 
	wr, 
	wi, 
	vl, 
	&ldvl, 
	vr, 
	n, 
	work, 
	work_n_,
	&info);
}

template<>
void solve<float>(solve_params(float))
{
  static constexpr const  char trN[1] = {'N'};
  static constexpr const char trV[1] = {'V'};
  wmesh_int_t ldvl=1,info;
  LAPACK_sgeev(trN,
	trV,
	n,
	jacM,
	n, 
	wr, 
	wi, 
	vl, 
	&ldvl, 
	vr, 
	n, 
	work, 
	work_n_,
	&info);
}


    
static wmesh_status_t bms_cubature_legendre_buffer_size(wmesh_int_t nspl_,wmesh_int_p work_n_)
{
  work_n_[0] = (nspl_+2)*3 + 2*(nspl_+2)*(nspl_+2) + 4*(nspl_+1);
  return WMESH_STATUS_SUCCESS;
}

template<typename T>
wmesh_status_t bms_cubature_legendre(wmesh_int_t 		nspl_,
				     wmesh_int_t 		storagep_,
				     T*__restrict__ 		p_,
				     wmesh_int_t 		ldp_, 
				     T*__restrict__ 		w_,
				     wmesh_int_t 		incw_, 
				     wmesh_int_t 		work_n_,
				     T*__restrict__		work_)
{
  static constexpr T rhalf = static_cast<T>(0.5);
  static constexpr T r0 = static_cast<T>(0.0);
  static constexpr T r1 = static_cast<T>(1.0);
  static constexpr T r2 = static_cast<T>(2.0);
  static constexpr T r4 = static_cast<T>(4.0);
  
  //  if (WFE::storage_t::is_invalid(storagep_)) return WMESH_STATUS_INVALID_ENUM;
  
  if (nspl_ < 0)    	return WMESH_STATUS_INVALID_SIZE;
  if (!p_)      	return WMESH_STATUS_INVALID_POINTER;
  if (ldp_  < 1) 	return WMESH_STATUS_INVALID_SIZE;
  
  if (w_)
    {      
      if (incw_ < 1) 	return WMESH_STATUS_INVALID_SIZE;
    }
  wmesh_int_t required_work_n;
  wmesh_status_t status =  bms_cubature_legendre_buffer_size(nspl_,
							     &required_work_n);
  WMESH_STATUS_CHECK(status);
  if (required_work_n > work_n_)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }
  
  if (work_n_ < required_work_n ) return WMESH_STATUS_INVALID_SIZE;
  work_n_ -= (nspl_+2)*3 + 2*(nspl_+2)*(nspl_+2);
  
  bool p_is_block = (storagep_ == WMESH_STORAGE_BLOCK);
  if (p_is_block && ldp_ < nspl_) return WMESH_STATUS_INVALID_SIZE;  

  wmesh_int_t n    	= nspl_+1;
  wmesh_int_t nn 	= nspl_+2;
  T 
    * __restrict__ vl = NULL,
    * __restrict__ b    = &work_[0],
    * __restrict__ jacM = &work_[nn],
    * __restrict__ wr   = &work_[nn+nn*nn],
    * __restrict__ wi   = &work_[2*nn+nn*nn],
    * __restrict__ vr   = &work_[3*nn+nn*nn],
    * __restrict__ work  = &work_[3*nn+2*nn*nn];
  
#if 0
  T 
    * __restrict__ vl = NULL,
    * __restrict__ b    = &work_[nn],
    * __restrict__ jacM = &work_[2*nn],
    * __restrict__ wr   = &work_[2*nn+nn*nn],
    * __restrict__ wi   = &work_[2*nn+nn*nn+nn],
    * __restrict__ vr   = &work_[2*nn+nn*nn+2*nn],
    * __restrict__ work  = &work_[2*nn+2*nn*nn+2*nn];
#endif
  for (wmesh_int_t k=0;k<n;++k)
    {
      b[k]=(T)0.0;
    }
  
  b[0] = r2;
  b[1] = r0;
  for (wmesh_int_t k=2;k<n;++k)
    {
      b[k] = r1 / (r4 - r1/( ((T)(k-1))  * ((T)(k-1)) ) );
    }
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      for (wmesh_int_t j=0;j<n;++j)
	{
	  jacM[i+j*n]=r0;
	}
    }
  
  for (wmesh_int_t i=0;i<n;++i)
    {
      jacM[i+i*n]=(T)0.0;
    }
  for (wmesh_int_t i=0;i<n-1;++i)
    {
      jacM[i+i*n+1]=wmesh_math<T>::xsqrt(b[1+i]);
    }
  for (wmesh_int_t i=1;i<n;++i)
    {
      jacM[i+i*n-1]=wmesh_math<T>::xsqrt(b[i]);
    }

  //
  // Eigen values / vectors decompositi
  //
  solve(&n,
	jacM,
	wr,
	wi,
	vl,
	vr,
	work,
	&work_n_);
 
  if (p_is_block)
    {
      for (wmesh_int_t i=0;i<n-1;++i) 
	{
	  p_[i] = wr[i];
	}
    }
  else
    {
      for (wmesh_int_t i=0;i<n-1;++i) 
	{
	  p_[i * ldp_] = wr[i];
	}	  
    }

  if (w_)
    {
      for (wmesh_int_t i=0;i<n-1;++i) 
	{
	  w_[i * incw_] = vr[1+i*n]*vr[1+i*n] * r2;
	}
    }

  for (int i = 0; i < (n-1)-1; i++)
    {
      int min_idx = i;
      for (int j = i+1; j < (n-1); j++)
	{
	  if (p_[j] < p_[min_idx])
	    {
	      min_idx = j;
	    }
	}
      
      { auto tmp = p_[min_idx]; p_[min_idx] = p_[i]; p_[i] = tmp; }

      if (w_)
	{
	  { auto tmp = w_[min_idx]; w_[min_idx] = w_[i]; w_[i] = tmp; }
	}
	  
    }
  
  //
  // Force symmetry.
  //
  if (nspl_ % 2 > 0)
    {
      p_[ (nspl_/2) *ldp_] = 0.0;
    }

  if (p_is_block)
    {
      for (int i=0;i<nspl_/2;++i)
	{
	  T tmp = (p_[i ] - p_[(nspl_-1-i) ])* rhalf;
	  p_[(nspl_-1-i)] = -tmp;
	  p_[i] = tmp;
	}
    }
  else
    {
      for (int i=0;i<nspl_/2;++i)
	{
	  T tmp = (p_[i *ldp_] - p_[(nspl_-1-i) *ldp_])*0.5;
	  p_[(nspl_-1-i) *ldp_] = -tmp;
	  p_[i *ldp_] = tmp;
	}
    }

  return WMESH_STATUS_SUCCESS;
}



wmesh_status_t
bms_cubature_num_nodes(wmesh_int_t	element_,
		       wmesh_int_t	family_,
		       wmesh_int_t	degree_,
		       wmesh_int_p	num_nodes_)
{
  wmesh_int_t num_nodes = (degree_%2 > 0) ? ((degree_ + 1) / 2) : ((degree_ + 2) / 2);
  num_nodes_[0] = 0;
  switch(family_)
    {
    case WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE:
      {
	switch(element_)
	  {
	  case WMESH_ELEMENT_NODE:
	    {
	      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	      break;
	    }
	    
	  case WMESH_ELEMENT_EDGE:
	    {
	      num_nodes_[0] = num_nodes;
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_ELEMENT_TRIANGLE:
	  case WMESH_ELEMENT_QUADRILATERAL:
	    {
	      num_nodes_[0] = num_nodes*num_nodes;
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_ELEMENT_TETRAHEDRON:
	  case WMESH_ELEMENT_PYRAMID:
	  case WMESH_ELEMENT_WEDGE:
	  case WMESH_ELEMENT_HEXAHEDRON:
	    {
	      num_nodes_[0] = num_nodes*num_nodes*num_nodes;
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  }
	return WMESH_STATUS_INVALID_ENUM;
      }
      
    case WMESH_CUBATURE_FAMILY_GAUSSLOBATTO:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);	
      }

    }
  
  return WMESH_STATUS_INVALID_ENUM;  
}



template <typename T>
wmesh_status_t bms_template_cubature_transform(wmesh_int_t 	element_,
					       wmesh_int_t 	q_r_n_,
					       const T * 	q_r_v_,
					       wmesh_int_t 	q_r_inc_,
					       const T * 	q_w_v_,
					       wmesh_int_t 	q_w_inc_,
					       
					       wmesh_int_t	c_storage_,
					       wmesh_int_t	c_m_,
					       wmesh_int_t	c_n_,
					       T * 		c_v_,
					       wmesh_int_t	c_ld_,
					       
					       T * 		w_v_,
					       wmesh_int_t 	w_inc_)
#if 0
  wmesh_int_t		rst_storage_,
  wmesh_int_t		rst_m_,
  wmesh_int_t		rst_n_,
  const T*__restrict__ 	rst_v_,
  wmesh_int_t		rst_ld_;
#endif
#if 0  
T*  r_v = rst_v_ + (rst_storage_ == WMESH_STORAGE_BLOCK) ? rst_ld_ * 0 : 0;
T*  s_v = rst_v_ + (rst_storage_ == WMESH_STORAGE_BLOCK) ? rst_ld_ * 1 : 1;
T*  t_v = rst_v_ + (rst_storage_ == WMESH_STORAGE_BLOCK) ? rst_ld_ * 2 : 2;
#endif  

{

  WMESH_CHECK_POSITIVE(q_r_n_);
  WMESH_CHECK_POINTER(q_r_v_);
  WMESH_CHECK_POSITIVE(q_r_inc_);

  WMESH_CHECK_POINTER(q_w_v_);
  WMESH_CHECK_POSITIVE(q_w_inc_);

  WMESH_CHECK_POSITIVE(c_m_);
  WMESH_CHECK_POSITIVE(c_n_);
  WMESH_CHECK_POINTER(c_v_);
  WMESH_CHECK(c_ld_ >= c_m_);

  WMESH_CHECK_POINTER(w_v_);
  WMESH_CHECK_POSITIVE(w_inc_);
  
  static constexpr T s_one 	= static_cast<T>(1);
  static constexpr T s_two 	= static_cast<T>(2);
  static constexpr T s_three 	= static_cast<T>(3);
  static constexpr T s_four 	= static_cast<T>(4);

  const bool c_is_block = c_storage_ == WMESH_STORAGE_BLOCK;
  
  T*  r = c_v_ + ( c_is_block ? c_ld_ * 0 : 0);
  T*  s = c_v_ + ( c_is_block ? c_ld_ * 1 : 1);
  T*  t = c_v_ + ( c_is_block ? c_ld_ * 2 : 2);
  const wmesh_int_t inc = c_is_block  ? 1 : c_ld_;

  switch(element_)
    {
    case WMESH_ELEMENT_NODE:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	break;
      }
	    
    case WMESH_ELEMENT_EDGE:
      {		
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_TRIANGLE:
      {
	for (wmesh_int_t j=0;j<q_r_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<q_r_n_;++i)
	      {
		r[ ( j * q_r_n_ + i ) * inc] 	= ( s_one + q_r_v_[i*q_r_inc_] ) * 0.5;
		s[ ( j * q_r_n_ + i ) * inc] 	= ( s_one - q_r_v_[i*q_r_inc_] ) * ( s_one - q_r_v_[j*q_r_inc_] ) * 0.25;		
		w_v_[ ( j * q_r_n_ + i ) * w_inc_] 	= ( s_one - q_r_v_[i*q_r_inc_] ) * q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_] * 0.125;
	      }
	  }	
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_QUADRILATERAL:
      {
	for (wmesh_int_t j=0;j<q_r_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<q_r_n_;++i)
	      {
		r[ ( j * q_r_n_ + i ) * inc] 	= q_r_v_[i * q_r_inc_];
		s[ ( j * q_r_n_ + i ) * inc] 	= q_r_v_[j * q_r_inc_];
		w_v_[ ( j * q_r_n_ + i ) * w_inc_] 	= q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_];
	      }
	  }
	
	return WMESH_STATUS_SUCCESS;
      }


    case WMESH_ELEMENT_TETRAHEDRON:
      {
	for (wmesh_int_t k=0;k<q_r_n_;++k)
	  {
	    for (wmesh_int_t j=0;j<q_r_n_;++j)
	      {
		for (wmesh_int_t i=0;i<q_r_n_;++i)
		  {
		    T u  = ( s_one + q_r_v_[i] ) / s_two;
		    T v  = ( s_one + q_r_v_[j] ) / s_two;
		    T w  = ( s_one + q_r_v_[k] ) / s_two;
		    
		    r[ ( k * q_r_n_ * q_r_n_ + j * q_r_n_ + i ) * inc] 	= u * v * w;
		    s[ ( k * q_r_n_ * q_r_n_ + j * q_r_n_ + i ) * inc] 	= u * v * (s_one - w);
		    t[ ( k * q_r_n_ * q_r_n_ + j * q_r_n_ + i ) * inc] 	= u * (s_one - u);
		    w_v_[ ( k * q_r_n_ * q_r_n_ + j * q_r_n_ + i ) * w_inc_] 	= (q_w_v_[i*q_w_inc_] *q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_] * s_three) / s_four;
		  }
	      }
	  }		
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_ELEMENT_HEXAHEDRON:
      {
	for (wmesh_int_t k=0;k<q_r_n_;++k)
	  {
	    for (wmesh_int_t j=0;j<q_r_n_;++j)
	      {
		for (wmesh_int_t i=0;i<q_r_n_;++i)
		  {
		    r[ ( q_r_n_ * q_r_n_ * k +  j * q_r_n_ + i ) * inc] 	= q_r_v_[i * q_r_inc_];
		    s[ ( q_r_n_ * q_r_n_ * k +  j * q_r_n_ + i ) * inc] 	= q_r_v_[j * q_r_inc_];
		    t[ ( q_r_n_ * q_r_n_ * k +  j * q_r_n_ + i ) * inc] 	= q_r_v_[k * q_r_inc_];
		    w_v_[ ( q_r_n_ * q_r_n_ * k +  j * q_r_n_ + i ) * w_inc_] 	= q_w_v_[i*q_w_inc_] * q_w_v_[j*w_inc_] * q_w_v_[k*q_w_inc_];
		  }
	      }
	  }
	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_WEDGE:
      {

	for (wmesh_int_t k=0;k<q_r_n_;++k)
	  {	    
	    for (wmesh_int_t j=0;j<q_r_n_;++j)
	      {
		for (wmesh_int_t i=0;i<q_r_n_;++i)
		  {
		    r[ ( k * q_r_n_* q_r_n_+j * q_r_n_ + i ) * inc] 	= ( s_one + q_r_v_[i*q_r_inc_] ) / s_two;
		    s[ ( k * q_r_n_* q_r_n_+j * q_r_n_ + i ) * inc] 	= ( s_one - q_r_v_[i*q_r_inc_] ) * ( s_one - q_r_v_[j*q_r_inc_] ) / s_four;
		    t[ ( k * q_r_n_* q_r_n_+j * q_r_n_ + i ) * inc] 	= ( s_one + q_r_v_[k*q_r_inc_] ) / s_two;		
		    w_v_[ ( k * q_r_n_* q_r_n_ + j * q_r_n_ + i ) * w_inc_] = ( s_one - q_r_v_[i*q_r_inc_] ) * q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_] * 0.125 * q_w_v_[k *q_w_inc_];
		  }
	      }	    
	  }

	return WMESH_STATUS_SUCCESS;
      }

    case WMESH_ELEMENT_PYRAMID:
      {

	for (wmesh_int_t k=0;k<q_r_n_;++k)
	  {
	    T ww = q_r_v_[k*q_w_inc_];
	    T tt = (ww + s_one) / s_two;
	    
	    wmesh_int_t shiftk = k * q_r_n_* q_r_n_;	    
	    T scal = s_one - tt;
	    for (wmesh_int_t j=0;j<q_r_n_;++j)
	      {
		T vv = q_r_v_[j*q_w_inc_];
	    	wmesh_int_t shiftj = shiftk + j * q_r_n_;
		for (wmesh_int_t i=0;i<q_r_n_;++i)
		  {
		    T uu = q_r_v_[i*q_w_inc_];
		    wmesh_int_t at = shiftj + i;
		    
		    r[ at * inc] = uu * scal;
		    s[ at * inc] = vv * scal;
		    t[ at * inc] = tt;
		    
		    w_v_[ at * w_inc_] 	= q_w_v_[i*q_w_inc_] * q_w_v_[j*q_w_inc_] * q_w_v_[k *q_w_inc_] / s_two;
		  }
	      }	    
	  }
	return WMESH_STATUS_SUCCESS;
      }
      
    }
  
  return WMESH_STATUS_INVALID_ENUM;
}


template<typename T>
wmesh_status_t
bms_template_cubature(wmesh_int_t	element_,
		      wmesh_int_t	family_,
		      wmesh_int_t	n1d_,

		      wmesh_int_t	c_storage_,
		      wmesh_int_t	c_m_,
		      wmesh_int_t	c_n_,
		      T*__restrict__ 	c_v_,
		      wmesh_int_t	c_ld_,
	     
		      wmesh_int_t	w_n_,
		      T*__restrict__ 	w_v_,
		      wmesh_int_t	w_inc_,
	     
		      wmesh_int_t	rwork_n_,
		      T* __restrict__ 	rwork_)

{
  WMESH_CHECK_POSITIVE(c_m_);
  WMESH_CHECK_POSITIVE(c_n_);
  WMESH_CHECK_POINTER(c_v_);
  WMESH_CHECK(c_ld_ >= c_m_);

#ifndef NDEBUG
  std::cout << "// WMESH : bms_cubature : element = " << element_ << std::endl;
  std::cout << "// WMESH : bms_cubature : family  = " << family_ << std::endl;
  std::cout << "// WMESH : bms_cubature : n1d     = " << n1d_ << std::endl;
#endif  
  
  wmesh_status_t status;
  wmesh_int_t required_rwork_n;
  status = bms_cubature_buffer_size(element_,
				    family_,
				    n1d_,
				    &required_rwork_n);
  WMESH_STATUS_CHECK(status);

  if (required_rwork_n > rwork_n_)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }
  
  if (rwork_n_>0)
    {
      WMESH_CHECK_POINTER(rwork_);
    }
  
  switch(family_)
    {
    case WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE:
      {
	switch(element_)
	  {
	  case WMESH_ELEMENT_NODE:
	    {
	      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	      break;
	    }
	    
	  case WMESH_ELEMENT_EDGE:
	    {
	      status =  bms_cubature_legendre(n1d_,
					      c_storage_,
					      c_v_,
					      c_ld_, 
					      w_v_,					      
					      w_inc_,
					      rwork_n_,
					      rwork_);
	      WMESH_STATUS_CHECK(status);
	      return WMESH_STATUS_SUCCESS;
	    }
	    
	  case WMESH_ELEMENT_TRIANGLE:
	  case WMESH_ELEMENT_QUADRILATERAL:
	  case WMESH_ELEMENT_TETRAHEDRON:
	  case WMESH_ELEMENT_PYRAMID:
	  case WMESH_ELEMENT_WEDGE:
	  case WMESH_ELEMENT_HEXAHEDRON:
	    {
	      
	      T * p1d		=  rwork_;
	      T * w1d		=  rwork_ + n1d_;
	      rwork_n_ 	       -= n1d_*2;
	      rwork_ 	       += n1d_*2;
	      
	      status =  bms_cubature_legendre(n1d_,
					      WMESH_STORAGE_INTERLEAVE,
					      p1d,
					      1,
					      w1d,
					      1,
					      rwork_n_,
					      rwork_);
	      
	      WMESH_STATUS_CHECK(status);
	      status = bms_template_cubature_transform(element_,
						       n1d_,
						       p1d,
						       1,
						       w1d,
						       1,
						       c_storage_,
						       c_m_,
						       c_n_,
						       c_v_,
						       c_ld_,
						       w_v_,
						       w_inc_);
	      WMESH_STATUS_CHECK(status);
	      return WMESH_STATUS_SUCCESS;  
	    }
	  }	    
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_CUBATURE_FAMILY_GAUSSLOBATTO:
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);	
      }

    }
  
   return WMESH_STATUS_INVALID_ENUM;  
}

extern "C"
{

  
  wmesh_status_t bms_cubature_buffer_size(wmesh_int_t 	element_,
					  wmesh_int_t 	family_,
					  wmesh_int_t	n1d_,			
					  wmesh_int_p 	rwork_n_) 
  {
    WMESH_CHECK_POINTER(rwork_n_);

    rwork_n_[0] = 0;
    switch(family_)
      {
      case WMESH_NODES_FAMILY_GAUSSLOBATTO:
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);	
	}
      
      case WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE:
	{
	  wmesh_status_t status;	
	  status = bms_cubature_legendre_buffer_size(n1d_,
						     rwork_n_);
	  WMESH_STATUS_CHECK(status);
	
	  switch(element_)
	    {
	    
	    case WMESH_ELEMENT_NODE:
	      {
		WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	      }
	    
	    case WMESH_ELEMENT_EDGE:
	      {
		return WMESH_STATUS_SUCCESS;
	      }
	    
	    case WMESH_ELEMENT_TRIANGLE:
	    case WMESH_ELEMENT_QUADRILATERAL:
	      {
		rwork_n_[0] += n1d_*2;	      
		return WMESH_STATUS_SUCCESS;  
	      }
	    case WMESH_ELEMENT_TETRAHEDRON:
	    case WMESH_ELEMENT_PYRAMID:
	    case WMESH_ELEMENT_WEDGE:
	    case WMESH_ELEMENT_HEXAHEDRON:
	      {
		rwork_n_[0] += n1d_*2;	      
		return WMESH_STATUS_SUCCESS;  
	      }
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	

      }
  
    return WMESH_STATUS_INVALID_ENUM;  
  };

#if 1
  wmesh_status_t
  bms_cubature(wmesh_int_t 		element_,
	       wmesh_int_t 		family_,
	       wmesh_int_t		n1d_,
	       
	       wmesh_int_t		c_storage_,
	       wmesh_int_t		c_m_,
	       wmesh_int_t		c_n_,
	       double* 			c_v_,
	       wmesh_int_t		c_ld_,
	       
	       wmesh_int_t		w_n_,
	       double* 			w_v_,
	       wmesh_int_t		w_inc_,
	       
	       wmesh_int_t		rwork_n_,
	       double* __restrict__ 	rwork_)
  {
    return bms_template_cubature(element_,
				 family_,
				 n1d_,
				 c_storage_,
				 c_m_,
				 c_n_,
				 c_v_,
				 c_ld_,
				 w_n_,
				 w_v_,
				 w_inc_,
				 rwork_n_,
				 rwork_);    
  }
#endif
}

#if 0


//
//
// exact_degree = 2n-1;
// n mx_d
// 1 1
// 2 3
// 3 5
// 4 7
// 5 9
// 6 11

	
//
//
// Find min n to integrate d exactly
//
// if (d % 2 > 0) n = (d+1)/2 else (d+2)/2
//
//



  
#endif
