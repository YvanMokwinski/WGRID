#pragma once

#include "wmesh_shape_info_t.h"

extern "C"
{
  struct wmesh_shape_t;
  wmesh_status_t wmesh_shape_def(wmesh_shape_t**__restrict__ 	self_,
				 wmesh_int_t 			element_,
				 wmesh_int_t 			family_,
				 wmesh_int_t 			degree_);

};
