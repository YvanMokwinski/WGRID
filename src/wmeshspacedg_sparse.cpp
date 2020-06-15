#include "wmesh.hpp"
#include "bms.h"
#include <iostream>

extern "C"
{
  
  
  wmesh_status_t wmeshspacedg_sparse	(const wmesh_t*__restrict__ 	self_,
					 wmesh_int_t 			degree_,
					 wmesh_int_p 			csr_size_,
					 wmesh_int_p*__restrict__ 	csr_ptr_,
					 wmesh_int_p*__restrict__ 	csr_ind_)
  {
    
    wmesh_int_t csr_n 		= 0;
    wmesh_int_t csr_m 		= 0;
    wmesh_int_t csr_nnz 	= 0;
    wmesh_int_p csr_ptr 	= nullptr;
    wmesh_int_p csr_ind 	= nullptr;
    wmesh_int_t size_blocks[4];
    
    wmesh_status_t status;

    // 1 1
    // 2 2
    // 4 4
    for (wmesh_int_t l=0;l<self_->m_c2c.m_size;++l)
      {
	wmesh_int_t ndofs;
	//
	// Trick.
	//
	wmesh_int_t element = self_->m_c2c.m_size + l;
	status = bms_ndofs(element,
			   degree_,
			   &ndofs);
	WMESH_STATUS_CHECK(status);
	size_blocks[l] = ndofs;
      }
    
    status = bms_sparse_dg_nnz(size_blocks,
			       WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2c),
			       &csr_n,
			       &csr_nnz);
    WMESH_STATUS_CHECK(status);
    csr_m = csr_n;

    csr_ptr = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*(csr_n+1));
    csr_ind = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*csr_nnz);    
    status = bms_sparse_dg(size_blocks,
			   WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2c),
			   csr_n,
			   csr_ptr,
			   csr_ind);
    WMESH_STATUS_CHECK(status);

    csr_size_[0] 	= csr_n;
    csr_ptr_[0] 	= csr_ptr;
    csr_ind_[0] 	= csr_ind;
    
    return WMESH_STATUS_SUCCESS;
  }
  
};
