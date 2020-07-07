#ifndef WMESH_SHAPE_INFO_H
#define WMESH_SHAPE_INFO_H

#include "wmesh-types.h"

#ifdef __cplusplus
extern "C"
{
#endif

  struct wmesh_shape_info_t
  {
    wmesh_int_t		m_family;
    wmesh_int_t 	m_degree;
  };
  
  wmesh_status_t wmesh_shape_info_def(wmesh_shape_info_t*__restrict__ 	self_,
				      wmesh_int_t 			family_,
				      wmesh_int_t 			degree_);

#ifdef __cplusplus
}
#endif


#endif
