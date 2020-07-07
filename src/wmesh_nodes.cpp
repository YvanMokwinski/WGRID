#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_t.hpp"
#include <chrono>
#include <iostream>
#include "bms.h"
#include "wmesh-utils.hpp"
#include "wmesh-blas.h"
#include "bms_templates.hpp"
#include "wmesh_nodes_info_t.hpp"
#include "wmesh_nodes_t.hpp"
extern "C"
{
  wmesh_status_t wmesh_nodes_info_def(wmesh_nodes_info_t*__restrict__ 	self_,
				      wmesh_int_t 			family_,
				      wmesh_int_t 			degree_)
  {
    WMESH_CHECK_POINTER(self_);
    memset(self_,0,sizeof(wmesh_nodes_info_t));
    self_->m_family 	= family_;
    self_->m_degree 	= degree_;
    return WMESH_STATUS_SUCCESS;
  };
};

template<typename T>
wmesh_status_t wmesh_nodes_def(wmesh_nodes_t<T>*__restrict__ self_,
			       wmesh_int_t 		element_,
			       wmesh_int_t 		family_,
			       wmesh_int_t 		degree_)
{
  memset(self_,0,sizeof(wmesh_nodes_t<T>));
  self_->m_element 	= element_;
  self_->m_family 	= family_;
  self_->m_degree 	= degree_;
  wmesh_status_t status;
  
  wmesh_int_t topodim;
  status = bms_element2topodim(element_,
			       &topodim);

  wmesh_int_t num_dofs;
  status = bms_ndofs(element_,
		     degree_,
		     &num_dofs);
  self_->m_c_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&self_->m_c, topodim, num_dofs);


  
  wmesh_int_t o_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<wmesh_int_t> o;
  wmesh_mat_t<wmesh_int_t>::alloc(&o, topodim, num_dofs);
  
  status = bms_ordering(element_,
			degree_,
			o_storage,
			WMESH_MAT_FORWARD(o));
  WMESH_STATUS_CHECK( status );
  
  wmesh_int_t iw_n, rw_n;
  status = bms_nodes_buffer_sizes(element_,
				  family_,
				  degree_,			
				  &iw_n,
				  &rw_n);
  WMESH_STATUS_CHECK( status );

  wmesh_int_p iw = (iw_n > 0) ? (wmesh_int_p)malloc(sizeof(wmesh_int_t)*iw_n) : nullptr;
  T * __restrict__ rw = (rw_n > 0) ? (T*__restrict__)malloc(sizeof(T)*rw_n) : nullptr;
  
  status = bms_nodes(element_,
		     family_,
		     degree_,
		     
		     o_storage,
		     WMESH_MAT_FORWARD(o),
		     
		     self_->m_c_storage,
		     WMESH_MAT_FORWARD(self_->m_c),

		     iw_n,
		     iw,
		     rw_n,
		     rw);
  WMESH_STATUS_CHECK( status );
  
  if (rw)
    {
      free(rw);
    }
  rw = nullptr;  
  if (iw)
    {
      free(iw);
    }
  iw = nullptr;
  
  return WMESH_STATUS_SUCCESS;
}

template
wmesh_status_t wmesh_nodes_def<float>(wmesh_nodes_t<float>*__restrict__ 	self_,
				      wmesh_int_t 				element_,
				      wmesh_int_t 				family_,
				      wmesh_int_t 				degree_);

template
wmesh_status_t wmesh_nodes_def<double>(wmesh_nodes_t<double>*__restrict__ 	self_,
				       wmesh_int_t 				element_,
				       wmesh_int_t 				family_,
				       wmesh_int_t 				degree_);
