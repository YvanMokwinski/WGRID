#ifndef WMESH_CUBATURE_INFO_H
#define WMESH_CUBATURE_INFO_H

#include "wmesh-types.h"

#ifdef __cplusplus
extern "C"
{
#endif

  struct wmesh_cubature_info_t;
  
  wmesh_status_t wmesh_cubature_info_def	(wmesh_cubature_info_t**__restrict__ 		self_,
						 wmesh_int_t 					family_,
						 wmesh_int_t 					degree_);
  
  wmesh_status_t wmesh_cubature_info_get_family	(const wmesh_cubature_info_t*__restrict__ 	self_,
						 wmesh_int_p 					family_);
  
  wmesh_status_t wmesh_cubature_info_get_degree	(const wmesh_cubature_info_t*__restrict__ 	self_,
						 wmesh_int_p 					degree_);

#ifdef __cplusplus
}
#endif


#endif
