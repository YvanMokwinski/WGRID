#include "wmesh_shape_info_t.hpp"
#include "bms.h"
#include <string.h>

extern "C"
{
  wmesh_status_t wmesh_shape_info_def(wmesh_shape_info_t**__restrict__ 	self__,
				      wmesh_int_t 			family_,
				      wmesh_int_t 			degree_)
  {
    self__[0] = new wmesh_shape_info_t(family_, degree_);
    return WMESH_STATUS_SUCCESS;
  };
  
  wmesh_status_t wmesh_shape_info_get_family	(const wmesh_shape_info_t*__restrict__ 	self_,
						 wmesh_int_p 				family_)
  {
    family_[0] = self_->get_family();
    return WMESH_STATUS_SUCCESS;
  };
  
  wmesh_status_t wmesh_shape_info_get_degree	(const wmesh_shape_info_t*__restrict__ 	self_,
						 wmesh_int_p 				degree_)
  {
    degree_[0] = self_->get_degree();
    return WMESH_STATUS_SUCCESS;
  };

}
