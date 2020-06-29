#pragma once

#include "wmesh_shape_info_t.hpp"


struct wmesh_shape_t
{
  wmesh_int_t		m_element;
  wmesh_int_t		m_family;
  wmesh_int_t		m_degree;
  wmesh_int_t		m_ndofs;
  wmesh_int_t		m_nodes_family;  
};
  
wmesh_status_t wmesh_shape_def(wmesh_shape_t*__restrict__ 	self_,
			       wmesh_int_t 			element_,
			       wmesh_int_t 			family_,
			       wmesh_int_t 			degree_);
