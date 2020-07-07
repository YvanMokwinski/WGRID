#include "bms.hpp"
#include "wmesh-math.hpp"
#include "wmesh-blas.hpp"
#include "wmesh-types.hpp"

#ifndef NDEBUG
#include <iostream>
#endif
#define solve_params(_type) wmesh_int_t*n,_type*jacM,_type*wr,_type*wi,_type*vl,_type*vr,_type*work,wmesh_int_t*work_n_
template<typename T>
void solve(solve_params(T));

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


static wmesh_status_t bms_nodes_legendre_buffer_size(wmesh_int_t nspl_,wmesh_int_p work_n_)
{
  work_n_[0] = (nspl_+2)*3 + 2*(nspl_+2)*(nspl_+2) + 4*(nspl_+1);
  return WMESH_STATUS_SUCCESS;
}
    
template<typename T>
wmesh_status_t bms_nodes_legendre(wmesh_int_t 		nspl_,
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
  wmesh_status_t status =  bms_nodes_legendre_buffer_size(nspl_,
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




template<typename T>
wmesh_status_t
bms_nodes(wmesh_int_t 		element_,
	  wmesh_int_t 		family_,
	  wmesh_int_t		degree_,

	  wmesh_int_t		b_storage_,
	  wmesh_int_t		b_m_,
	  wmesh_int_t		b_n_,
	  const_wmesh_int_p 	b_v_,
	  wmesh_int_t		b_ld_,

	  wmesh_int_t		c_storage_,
	  wmesh_int_t		c_m_,
	  wmesh_int_t		c_n_,
	  T*__restrict__ 	c_v_,
	  wmesh_int_t		c_ld_,

	  wmesh_int_t		iwork_n_,
	  wmesh_int_p 		iwork_,
	  wmesh_int_t		rwork_n_,
	  T* __restrict__ 	rwork_)  
{
  WMESH_CHECK_POSITIVE(c_m_);
  WMESH_CHECK_POSITIVE(c_n_);
  WMESH_CHECK_POINTER(c_v_);
  WMESH_CHECK(c_ld_ >= c_m_);

  WMESH_CHECK_POSITIVE(b_m_);
  WMESH_CHECK_POSITIVE(b_n_);
  WMESH_CHECK_POINTER(b_v_);
  WMESH_CHECK(b_ld_ >= b_m_);

#ifndef NDEBUG
  std::cout << "// WMESH : bms_nodes : element = " << element_ << std::endl;
  std::cout << "// WMESH : bms_nodes : family  = " << family_ << std::endl;
  std::cout << "// WMESH : bms_nodes : degree  = " << degree_ << std::endl;
#endif  

  wmesh_status_t status;
  wmesh_int_t required_iwork_n;
  wmesh_int_t required_rwork_n;
  status = bms_nodes_buffer_sizes(element_,
				  family_,
				  degree_,
				  &required_iwork_n,
				  &required_rwork_n);
  WMESH_STATUS_CHECK(status);
  if (required_iwork_n > iwork_n_)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }
  if (required_rwork_n > rwork_n_)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }

  if (rwork_n_>0)
    {
      WMESH_CHECK_POINTER(rwork_);
    }
  
  if (iwork_n_>0)
    {
      WMESH_CHECK_POINTER(iwork_);
    }
  
  wmesh_int_t num_nodes_P1;
  status = bms_elements_num_nodes(1,&element_,&num_nodes_P1);
  WMESH_STATUS_CHECK(status);
  wmesh_int_t num_nodes_Pk;
  status = bms_ndofs(element_,degree_,&num_nodes_Pk);
  WMESH_STATUS_CHECK(status);

  const wmesh_int_t topodim = (b_storage_ == WMESH_STORAGE_INTERLEAVE)  ? b_m_ : b_n_;
  const wmesh_int_t		c_storage = c_storage_;
  wmesh_mat_t<T> c;
  wmesh_mat_t<T>::define(&c,c_m_,c_n_,c_v_,c_ld_);
  
  
  const wmesh_int_t eval_storage 	= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T> eval;
  wmesh_mat_t<T>::alloc(&eval,
			num_nodes_P1,
			num_nodes_Pk);
  
  //  const wmesh_int_t ref_c_storage 	= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T> ref_c;
  wmesh_mat_t<T>::alloc(&ref_c,
			topodim,
			num_nodes_P1);
  switch(family_)
    {
    case WMESH_NODES_FAMILY_LAGRANGE:
      {
	if (degree_ > 0)
	  {
	    for (wmesh_int_t j=0;j<c_n_;++j)
	      {
		for (wmesh_int_t i=0;i<c_m_;++i)
		  {
		    c_v_[c_ld_*j+i] = ((T)b_v_[b_ld_*j+i]) / ((T)degree_);
		  }
	      }

	    //
	    // The ordering reference geometry is in the first nodes of c.
	    //

	    
	    //
	    // Compute the ordering linear shapes over the ordering reference cells.
	    //	    	   
	    status = bms_ordering_linear_shape(element_,
					       c_storage,
					       WMESH_MAT_FORWARD(c),
					       eval_storage,
					       WMESH_MAT_FORWARD(eval));
	    WMESH_STATUS_CHECK(status);
	    
	    //
	    // Get element geometry.
	    //
	    status = bms_element_geometry(element_,
					  ref_c.v);
	    WMESH_STATUS_CHECK(status);
	    
	    wmesh_mat_gemm(static_cast<T>(1),ref_c,eval,static_cast<T>(0),c);
	  }
	else
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	  }
	return WMESH_STATUS_SUCCESS;
      }
      
    case WMESH_NODES_FAMILY_GAUSSLOBATTO:
      {
	switch(element_)
	  {
	  case WMESH_ELEMENT_NODE:
	    {
	      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
	    }
	    
	  case WMESH_ELEMENT_EDGE:
	    {
	      c_v_[0]  = static_cast<T>(-1);
	      c_v_[c_ld_*1]  = static_cast<T>(1);
	      status =  bms_nodes_legendre(degree_-1,
					   c_storage_,
					   c_v_ + 2 * c_ld_,
					   c_ld_, 
					   (T*__restrict__)nullptr,
					   1,
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
	      //
	      // Flat the ordering.
	      //
	      wmesh_int_p perm = iwork_;	      
	      status =  bms_ordering_flat(degree_,
					  b_storage_,
					  b_m_,
					  b_n_,
					  b_v_,
					  b_ld_,
					  perm,
					  1);	      
	      WMESH_STATUS_CHECK(status);
	      
	      //
	      // Gauss-lobatto 1d.
	      //
	      wmesh_int_t n1d 	= (degree_+1);
	      T * p1d 		=  rwork_;
	      rwork_n_ 	       -= n1d;	      
	      rwork_ 	       += n1d;

	      //
	      // Be careful, this order -1 ... 1 is not the finite element order but we need it in this way for the transformation.
	      //
	      status =  bms_nodes_legendre(degree_-1,
					   WMESH_STORAGE_INTERLEAVE,
					   p1d + 1,
					   1,
					   (T*__restrict__)nullptr,
					   1,
					   rwork_n_,
					   rwork_);
	      WMESH_STATUS_CHECK(status);	      
	      p1d[0]  		= -1;
	      p1d[degree_]  	= 1;	    
	      
	      //
	      // The transformation is done on the ordering reference element.
	      //
	      status = bms_transform(element_,
				     n1d,
				     p1d,
				     1,
					 
				     c_storage_,
				     c_m_,
				     c_n_,
				     c_v_,					 
				     c_ld_,
				     
				     perm,
				     1);
	      WMESH_STATUS_CHECK(status);

	      //
	      // Compute the ordering linear shapes over the ordering reference cells.
	      //	    	   
	      status = bms_ordering_linear_shape(element_,
						 c_storage,
						 WMESH_MAT_FORWARD(c),
						 eval_storage,
						 WMESH_MAT_FORWARD(eval));
	      WMESH_STATUS_CHECK(status);
	      
	      //
	      // Get element geometry.
	      //
	      status = bms_element_geometry(element_,
					    ref_c.v);
	      WMESH_STATUS_CHECK(status);
	      
	      wmesh_mat_gemm(static_cast<T>(1),ref_c,eval,static_cast<T>(0),c);
	      
	      return WMESH_STATUS_SUCCESS;  
	    }
	    
	  }
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);	
      }

    }
  
  return WMESH_STATUS_INVALID_ENUM;  
}

template
wmesh_status_t
bms_nodes<double>(wmesh_int_t 		element_,
	  wmesh_int_t 		family_,
	  wmesh_int_t		degree_,

	  wmesh_int_t		b_storage_,
	  wmesh_int_t		b_m_,
	  wmesh_int_t		b_n_,
	  const_wmesh_int_p 	b_v_,
	  wmesh_int_t		b_ld_,

	  wmesh_int_t		c_storage_,
	  wmesh_int_t		c_m_,
	  wmesh_int_t		c_n_,
	double*__restrict__ 	c_v_,
	  wmesh_int_t		c_ld_,

	  wmesh_int_t		iwork_n_,
	  wmesh_int_p 		iwork_,
	  wmesh_int_t		rwork_n_,
		  double* __restrict__ 	rwork_)  ;

template
wmesh_status_t
bms_nodes<float>(wmesh_int_t 		element_,
	  wmesh_int_t 		family_,
	  wmesh_int_t		degree_,

	  wmesh_int_t		b_storage_,
	  wmesh_int_t		b_m_,
	  wmesh_int_t		b_n_,
	  const_wmesh_int_p 	b_v_,
	  wmesh_int_t		b_ld_,

	  wmesh_int_t		c_storage_,
	  wmesh_int_t		c_m_,
	  wmesh_int_t		c_n_,
	float*__restrict__ 	c_v_,
	  wmesh_int_t		c_ld_,

	  wmesh_int_t		iwork_n_,
	  wmesh_int_p 		iwork_,
	  wmesh_int_t		rwork_n_,
		  float* __restrict__ 	rwork_)  ;

  

extern "C"
{

  
  wmesh_status_t bms_nodes_buffer_sizes(wmesh_int_t 	element_,
				       wmesh_int_t 	family_,
				       wmesh_int_t	degree_,			
				       wmesh_int_p	iwork_n_,
				       wmesh_int_p 	rwork_n_) 
  {
    WMESH_CHECK_POINTER(iwork_n_);
    WMESH_CHECK_POINTER(rwork_n_);

    iwork_n_[0] = 0;
    rwork_n_[0] = 0;
    switch(family_)
      {
      case WMESH_NODES_FAMILY_LAGRANGE:
	{
	  return WMESH_STATUS_SUCCESS;
	}
      
      case WMESH_NODES_FAMILY_GAUSSLOBATTO:
	{
	  wmesh_status_t status;	
	  status = bms_nodes_legendre_buffer_size(degree_-1,
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
		rwork_n_[0] += degree_+1;	      
		iwork_n_[0] += (degree_+1)*(degree_+1);
		return WMESH_STATUS_SUCCESS;  
	      }
	    case WMESH_ELEMENT_TETRAHEDRON:
	    case WMESH_ELEMENT_PYRAMID:
	    case WMESH_ELEMENT_WEDGE:
	    case WMESH_ELEMENT_HEXAHEDRON:
	      {
		rwork_n_[0] += degree_+1;	      
		iwork_n_[0] += (degree_+1)*(degree_+1)*(degree_+1);
		return WMESH_STATUS_SUCCESS;  
	      }
	    }
	  return WMESH_STATUS_SUCCESS;
	}
	

      }
  
    return WMESH_STATUS_INVALID_ENUM;  
  };

  
  wmesh_status_t
  bms_snodes(wmesh_int_t 	element_,
	     wmesh_int_t 	family_,
	     wmesh_int_t	degree_,
	     
	     wmesh_int_t	b_storage_,
	     wmesh_int_t	b_m_,
	     wmesh_int_t	b_n_,
	     const_wmesh_int_p 	b_v_,
	     wmesh_int_t	b_ld_,

	     wmesh_int_t	c_storage_,
	     wmesh_int_t	c_m_,
	     wmesh_int_t	c_n_,
	     float* 		c_v_,
	     wmesh_int_t	c_ld_,
	     wmesh_int_t	iwork_n_,
	     wmesh_int_p 	iwork_,
	     wmesh_int_t	rwork_n_,
	     float* __restrict__ rwork_)
  {
    return bms_nodes(element_,
		     family_,
		     degree_,
		     b_storage_,
		     b_m_,
		     b_n_,
		     b_v_,
		     b_ld_,
		     c_storage_,
		     c_m_,
		     c_n_,
		     c_v_,
		     c_ld_,
		     iwork_n_,
		     iwork_,
		     rwork_n_,
		     rwork_);    
  }

  wmesh_status_t
  bms_dnodes(wmesh_int_t 	element_,
	     wmesh_int_t 	family_,
	     wmesh_int_t	degree_,
	     wmesh_int_t	b_storage_,
	     wmesh_int_t	b_m_,
	     wmesh_int_t	b_n_,
	     const_wmesh_int_p 	b_v_,
	     wmesh_int_t	b_ld_,
	     wmesh_int_t	c_storage_,
	     wmesh_int_t	c_m_,
	     wmesh_int_t	c_n_,
	     double* 		c_v_,
	     wmesh_int_t	c_ld_,
	     wmesh_int_t	iwork_n_,
	     wmesh_int_p 	iwork_,
	     wmesh_int_t	rwork_n_,
	     double* __restrict__ 	rwork_)
  {    
    return bms_nodes(element_,
		     family_,
		     degree_,
		     b_storage_,
		     b_m_,
		     b_n_,
		     b_v_,
		     b_ld_,
		     c_storage_,
		     c_m_,
		     c_n_,
		     c_v_,
		     c_ld_,
		     iwork_n_,
		     iwork_,
		     rwork_n_,
		     rwork_);    
  }

}
