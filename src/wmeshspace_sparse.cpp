#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh.hpp"
#include <chrono>
#include <iostream>
#include "bms.h"

#include "wmesh-utils.hpp"
using namespace std::chrono;
extern "C"
{
  
  wmesh_status_t wmeshspace_sparse	(const wmeshspace_t*__restrict__ 	self_,
					 wmesh_int_p 				csr_size_,
					 wmesh_int_p*__restrict__ 		csr_ptr_,
					 wmesh_int_p*__restrict__ 		csr_ind_)
  {
    wmesh_status_t status;    

    wmesh_int_t num_dofs 		= self_->m_ndofs;
    wmesh_int_t iw_n;
    status = bms_sparse_buffer_size	(self_->m_ndofs,
					 self_->m_c2d.m_size,
					 self_->m_c2d.m_m,
					 self_->m_c2d.m_n,
					 &iw_n);
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_p iw      = (iw_n > 0) ? (wmesh_int_p)malloc(iw_n * sizeof(wmesh_int_t)) : nullptr;
    wmesh_int_p csr_ptr = (wmesh_int_p)malloc( (num_dofs + 1)     * sizeof(wmesh_int_t));
    wmesh_int_p csr_ind = (wmesh_int_p)malloc( (num_table_coeffs) * sizeof(wmesh_int_t));
    
    status = bms_sparse	(wmesh_int_t 		num_dofs_,
			 
			 wmesh_int_t 		c2d_size_,
			 const_wmesh_int_p	c2d_ptr_,
			 const_wmesh_int_p	c2d_m_,
			 const_wmesh_int_p	c2d_n_,
			 const_wmesh_int_p	c2d_v_,
			 const_wmesh_int_p	c2d_ld_,
			 
			 wmesh_int_p 		csr_ptr_,
			 wmesh_int_p 		csr_ind_,
			 wmesh_int_t		iw_n_,
			 wmesh_int_p		iw_);
  


    

    
    wmesh_status_t status;
    
    wmesh_int_t n2c_m 	= mesh_->m_num_nodes;
    wmesh_int_p n2c_ptr = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*(n2c_m+1)); 
    wmesh_int_p n2c_v 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mesh_->m_c2n.m_ptr[mesh_->m_c2n.m_size]); 

    status = bms_n2c(WMESH_INT_SPARSEMAT_FORWARD(mesh_->m_c2n),
		     n2c_ptr,
		     n2c_m,
		     n2c_v);
    
    WMESH_STATUS_CHECK(status);

    wmesh_int_p blank 	= (wmesh_int_p)calloc(n2c_m,sizeof(wmesh_int_t)); 
    wmesh_int_p select 	= (wmesh_int_p)calloc(n2c_m,sizeof(wmesh_int_t));
    csr_size_[0] = n2c_m;
    csr_ptr_[0] = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*(n2c_m + 1));
    csr_ptr_[0][0] = 0;
    for (wmesh_int_t idof=0;idof<n2c_m;++idof)
      {
	wmesh_int_t select_n = 0;
	//
	// Union 
	//
	for (wmesh_int_t s = n2c_ptr[idof];s<n2c_ptr[idof+1];++s)
	  {
	    wmesh_int_t c = n2c_v[s];
	    wmesh_int_t cindex, ctype;
	    
	    status = bms_n2c_cindex(c,
				    &cindex);
	    WMESH_STATUS_CHECK(status);
	    
	    status = bms_n2c_ctype	(c,
					 &ctype);
	    WMESH_STATUS_CHECK(status);
	    
	    for (wmesh_int_t i=0;i<mesh_->m_c2n.m_m[ctype];++i)
	      {
		wmesh_int_t jdof = mesh_->m_c2n.m_data[mesh_->m_c2n.m_ptr[ctype] + cindex * mesh_->m_c2n.m_ld[ctype] + i];
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

	csr_ptr_[0][idof+1] = csr_ptr_[0][idof] + select_n;	
      }

    csr_ind_[0] = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*csr_ptr_[0][n2c_m]);
    for (wmesh_int_t idof=0;idof<n2c_m;++idof)
      {
	wmesh_int_t select_n = 0;
	//
	// Union 
	//
	for (wmesh_int_t s = n2c_ptr[idof];s<n2c_ptr[idof+1];++s)
	  {
	    wmesh_int_t c = n2c_v[s];
	    wmesh_int_t cindex, ctype;
	    
	    status = bms_n2c_cindex	(c,
					 &cindex);
	    WMESH_STATUS_CHECK(status);
	    
	    status = bms_n2c_ctype	(c,
					 &ctype);
	    WMESH_STATUS_CHECK(status);
	    
	    for (wmesh_int_t i=0;i<mesh_->m_c2n.m_m[ctype];++i)
	      {
		wmesh_int_t jdof = mesh_->m_c2n.m_data[mesh_->m_c2n.m_ptr[ctype] + cindex * mesh_->m_c2n.m_ld[ctype] + i];
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
	    csr_ind_[0][csr_ptr_[0][idof] + s] = select[s];
	  }
	
      }
    free(select);
    free(blank);
    free(n2c_ptr);
    free(n2c_v);
    return WMESH_STATUS_SUCCESS;
  }
  
};
