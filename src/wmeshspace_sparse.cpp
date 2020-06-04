#include "wmesh.hpp"
#include "bms.h"
#include <iostream>

extern "C"
{
  wmesh_status_t wmeshspace_sparse	(const wmeshspace_t*__restrict__ 	self_,
					 wmesh_int_p 				csr_size_,
					 wmesh_int_p*__restrict__ 		csr_ptr_,
					 wmesh_int_p*__restrict__ 		csr_ind_)
  {
    
    WMESH_CHECK_POINTER(self_);
    WMESH_CHECK_POINTER(csr_size_);
    WMESH_CHECK_POINTER(csr_ptr_);
    WMESH_CHECK_POINTER(csr_ind_);
    
    wmesh_status_t status;    

    wmesh_int_t num_dofs 		= self_->m_ndofs;
    wmesh_int_t iw_n, num_table_coeffs;

    status = bms_sparse_buffer_size	(self_->m_ndofs,
					 self_->m_c2d.m_size,
					 self_->m_c2d.m_m,
					 self_->m_c2d.m_n,
					 &iw_n,
					 &num_table_coeffs);
    WMESH_STATUS_CHECK(status);
    
    
    wmesh_int_p iw  = (iw_n > 0) ? (wmesh_int_p)malloc(iw_n * sizeof(wmesh_int_t)) : nullptr;
    if (iw_n > 0 && !iw)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }

    wmesh_int_p csr_ptr = (wmesh_int_p)malloc( (num_dofs + 1) * sizeof(wmesh_int_t));
    if (!csr_ptr)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }

    status = bms_sparse_ptr(num_dofs,
			    WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2d),			 
			    csr_ptr,
			    iw_n,
			    iw);

    WMESH_STATUS_CHECK(status);


    wmesh_int_p csr_ind = (wmesh_int_p)malloc( csr_ptr[num_dofs] * sizeof(wmesh_int_t));
    if (!csr_ind)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    
    status = bms_sparse	(num_dofs,
			 WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2d),			 
			 csr_ptr,
			 csr_ind,
			 iw_n,
			 iw);
    WMESH_STATUS_CHECK(status);

    //
    // Free iw.
    //
    free(iw);
    iw = nullptr;

    csr_size_[0] 	= num_dofs;
    csr_ptr_[0] 	= csr_ptr;
    csr_ind_[0] 	= csr_ind;
    
    return WMESH_STATUS_SUCCESS;
  }
  
};
