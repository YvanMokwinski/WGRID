#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh.hpp"
#include <chrono>
#include <iostream>
#include "bms.h"
#include "wmesh-utils.hpp"
extern "C"
{
  
  wmesh_status_t wmeshspace_laplace	(const wmeshspace_t*__restrict__ 	self_,
					 wmesh_int_t 				csr_size_,
					 const_wmesh_int_p			csr_ptr_,
					 const_wmesh_int_p			csr_ind_,
					 double * 				csr_val_,
					 double * 				rhs_)
  {
    return WMESH_STATUS_SUCCESS;
  }
};
    
#if 0


using namespace std::chrono;
extern "C"
{
  status = wmeshspace_laplace(meshspace,
			      csr_size,
			      csr_ptr,
			      csr_ind,
			      csr_val,
			      rhs);
  WMESH_STATUS_CHECK(status);
  
  wmesh_status_t wmeshspace_laplace	(const wmeshspace_t*__restrict__ 	self_,
					 wmesh_int_t 				csr_size_,
					 const_wmesh_int_p			csr_ptr_,
					 const_wmesh_int_p			csr_ind_,
					 double * 				csr_val_,
					 double * 				rhs_)
  {    
    WMESH_CHECK_POINTER(self_);
    WMESH_CHECK_POINTER(csr_size_);
    WMESH_CHECK_POINTER(csr_ptr_);
    WMESH_CHECK_POINTER(csr_ind_);
    
    wmesh_status_t status;    
    wmesh_int_p c2d = iw_;
    double cooelm[32];
    double c2d[32];
    
    //
    // Prepare data.
    //
    wmesh_int_t cubature_n[4]{};
    double * 	cubature_p[4]{};
    double * 	cubature_w[4]{};
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++i)
      {
	if (self_->m_c2d.m_n[l] > 0)
	  {
	    
	  }
      }
    
    //
    // Element shape.
    //
    double * 	eval_shape_element[4]{};
    double * 	eval_shape_element_dr[4]{};
    double * 	eval_shape_element_ds[4]{};
    double * 	eval_shape_element_dt[4]{};
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++i)
      {
	if (self_->m_c2d.m_n[l] > 0)
	  {
	    
	  }
      }

    //
    // Variable shape.
    //
    double * 	eval_shape_f[4]{};
    double * 	eval_shape_f_dr[4]{};
    double * 	eval_shape_f_ds[4]{};
    double * 	eval_shape_f_dt[4]{};
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++i)
      {
	if (self_->m_c2d.m_n[l] > 0)
	  {
	    
	  }
      }
    
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++i)
      {

	if (self_->m_c2d.m_n[l] > 0)
	  {
	    
	    //
	    // Need integrations points.
	    //

	    //
	    // Need to evaluate shape of the variable F on integrations points.
	    //

	    //
	    // Need to evaluate shape of the element on integrations points.
	    //

	    //
	    // Now loop over the element of this type.
	    //
	    for (wmesh_int_t j=0;j<self_->m_c2d.m_n[l];++j)
	      {

		//
		// Get the coordinates of the elements.
		//
		cooelm[0] = 0;
		
		//
		// Evaluate the derivatives of the transformation on  the cubature.
		//

		
		//
		// Element jacobian and determinant on the cubature.
		//


		//
		// Element jacobian and determinant on the cubature.
		//

		
		//
		// Load the dof indices.
		//
		for (wmesh_int_t i=0;i<self_->m_c2d.m_m[l];++i)
		  {
		    c2d[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[l] + j * self_->m_c2d.m_ld[l] + i] - 1;
		  }

		
		//
		// Assembly.
		//

	      }
	  }
      }

    
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

    std::cout << "ptr" << std::endl;
    status = bms_sparse_ptr(num_dofs,
			    WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2d),			 
			    csr_ptr,
			    iw_n,
			    iw);

    WMESH_STATUS_CHECK(status);
    std::cout << "ptr done" << std::endl;

    wmesh_int_p csr_ind = (wmesh_int_p)malloc( csr_ptr[num_dofs] * sizeof(wmesh_int_t));
    if (!csr_ind)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    std::cout << "sp .." << std::endl;
    
    status = bms_sparse	(num_dofs,
			 WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2d),			 
			 csr_ptr,
			 csr_ind,
			 iw_n,
			 iw);
    WMESH_STATUS_CHECK(status);
    std::cout << "sp done." << std::endl;

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
#endif
