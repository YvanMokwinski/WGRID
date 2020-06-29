#pragma once
#include "wmesh-types.hpp"

extern "C"
{
  
  
  struct wmesh_cubature_info_t
  {
    wmesh_int_t		m_family;
    wmesh_int_t 	m_degree;
  };
  
  wmesh_status_t wmesh_cubature_info_def(wmesh_cubature_info_t*__restrict__ 	self_,
					 wmesh_int_t 				family_,
					 wmesh_int_t 				degree_);

};

