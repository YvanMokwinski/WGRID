#include "bms.h"
#include "wmesh-utils.hpp"


extern "C"
{

  wmesh_status_t bms_sparse_buffer_size	(wmesh_int_t 		num_dofs_,
					 wmesh_int_t 		c2d_size_,
					 const_wmesh_int_p	c2d_m_,
					 const_wmesh_int_p	c2d_n_,
					 wmesh_int_p		iw_n_)
  {
    WMESH_CHECK_POINTER(c2d_m_);
    WMESH_CHECK_POINTER(c2d_n_);
    WMESH_CHECK_POINTER(iw_n_);
    wmesh_int_t num_table_coeffs 	= 0;    
    for (wmesh_int_t i=0;i<c2d_size_;++i)
      {
	num_table_coeffs += c2d_m_[i] * c2d_n_[i];
      }
    iw_n_[0] = ( ((num_dofs_ + 1) + num_coeffs_) + (2 * num_dofs_) );
    return WMESH_STATUS_SUCCESS;
  }
  
  wmesh_status_t bms_sparse	(wmesh_int_t 		num_dofs_,

				 wmesh_int_t 		c2d_size_,
				 const_wmesh_int_p	c2d_ptr_,
				 const_wmesh_int_p	c2d_m_,
				 const_wmesh_int_p	c2d_n_,
				 const_wmesh_int_p	c2d_v_,
				 const_wmesh_int_p	c2d_ld_,
				 
				 wmesh_int_p 		csr_ptr_,
				 wmesh_int_p 		csr_ind_,
				 wmesh_int_t		iw_n_,
				 wmesh_int_p		iw_)
  {
    wmesh_status_t 	status;

    wmesh_int_t required_iw_n;

    status = bms_sparse_buffer_size(num_dofs_,
				    c2d_size_,
				    c2d_m_,
				    c2d_n_,
				    &required_iw_n);    
    WMESH_STATUS_CHECK(status);
    if (iw_n_ < required_iw_n)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
      }
    
    wmesh_int_t 	n2d_m 	= num_dofs_;
    wmesh_int_p 	n2d_ptr = iw_;
    wmesh_int_p 	n2d_v 	= n2d_ptr + (num_dofs_+1);
    wmesh_int_p 	blank 	= n2d_v + num_table_coeffs;
    wmesh_int_p 	select 	= blank + num_dofs_;
    
    for (wmesh_int_t i=0;i<num_dofs_;++i)
      {
	blank[i] = 0;
      }
    
    status = bms_n2c(c2d_size_,
		     c2d_ptr_,
		     c2d_m_,
		     c2d_n_,
		     c2d_v_,
		     c2d_ld_,		     
		     n2d_ptr,
		     n2d_m,
		     n2d_v);
    
    WMESH_STATUS_CHECK(status);
    
    csr_ptr_[0] = 0;    
    for (wmesh_int_t idof = 0;idof < n2d_m;++idof)
      {
	wmesh_int_t select_n = 0;
	//
	// Union 
	//
	for (wmesh_int_t s = n2d_ptr[idof];s<n2d_ptr[idof+1];++s)
	  {
	    wmesh_int_t c = n2d_v[s];
	    wmesh_int_t cindex, ctype;
	    
	    status = bms_n2c_cindex(c,
				    &cindex);
	    WMESH_STATUS_CHECK(status);
	    
	    status = bms_n2c_ctype	(c,
					 &ctype);
	    WMESH_STATUS_CHECK(status);
	    
	    for (wmesh_int_t i=0;i<c2d_m_[ctype];++i)
	      {
		wmesh_int_t jdof = c2d_v_[c2d_ptr_[ctype] + cindex * c2d_ld_[ctype] + i];
		if (0 == blank[jdof])
		  {
		    select[select_n++] = jdof;
		    blank[jdof] = select_n;
		  }
	      }
	  }
	
	//
	// Reset.
	//
	for (wmesh_int_t s = 0;s<select_n;++s)
	  {
	    blank[select[s]] = 0;
	  }
	
	csr_ptr_[idof+1] = csr_ptr_[idof] + select_n;	
      }
    

    for (wmesh_int_t idof=0;idof<n2d_m;++idof)
      {
	wmesh_int_t select_n = 0;
	//
	// Union 
	//
	for (wmesh_int_t s = n2d_ptr[idof];s<n2d_ptr[idof+1];++s)
	  {
	    wmesh_int_t c = n2d_v[s];
	    wmesh_int_t cindex, ctype;
	    
	    status = bms_n2c_cindex	(c,
					 &cindex);
	    WMESH_STATUS_CHECK(status);

	    
	    status = bms_n2c_ctype	(c,
					 &ctype);
	    WMESH_STATUS_CHECK(status);
	    
	    for (wmesh_int_t i=0;i<c2d_m_[ctype];++i)
	      {
		wmesh_int_t jdof = c2d_m_[c2d_ptr_[ctype] + cindex * c2d_ld_[ctype] + i];
		if (0==blank[jdof])
		  {
		    select[select_n++] = jdof;
		    blank[jdof] = select_n;
		  }
	      }
	  }

	//
	// Reset.
	//
	for (wmesh_int_t s = 0;s<select_n;++s)
	  {
	    blank[select[s]] = 0;
	  }

	//
	// Sort.
	//
	qsort(select,
	      select_n,
	      sizeof(wmesh_int_t),
	      wmesh_qsort_increasing_predicate<wmesh_int_t>);

	//
	// Copy back.
	//
	for (wmesh_int_t s = 0;s<select_n;++s)
	  {
	    csr_ind_[csr_ptr_[idof] + s] = select[s];
	  }	
      }
    return WMESH_STATUS_SUCCESS;
  }

};
