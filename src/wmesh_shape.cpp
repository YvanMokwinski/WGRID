#include "wmesh-types.hpp"
#include "bms.h"
#include <string.h>

extern "C"
{
  wmesh_status_t wmesh_shape_def(wmesh_shape_t*__restrict__ self_,
				 wmesh_int_t 		element_,
				 wmesh_int_t 		family_,
				 wmesh_int_t 		degree_)
  {
    WMESH_CHECK_POINTER(self_);
    memset(self_,0,sizeof(wmesh_shape_t));
    self_->m_element 	= element_;
    self_->m_family 	= family_;
    self_->m_degree 	= degree_;

    wmesh_status_t status = bms_ndofs(element_,
				      degree_,
				      &self_->m_ndofs);
    WMESH_STATUS_CHECK(status);
    return WMESH_STATUS_SUCCESS;
  }
  
};
