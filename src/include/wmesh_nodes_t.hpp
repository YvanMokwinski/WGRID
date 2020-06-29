#pragma once

#include "wmesh_nodes_info_t.hpp"

template<typename T>
struct wmesh_nodes_t
{
  wmesh_int_t 		m_element;
  wmesh_int_t 		m_family;
  wmesh_int_t 		m_degree;
  wmesh_int_t 		m_c_storage;
  wmesh_mat_t<T> 	m_c;
};

template<typename T>
wmesh_status_t wmesh_nodes_def(wmesh_nodes_t<T>*__restrict__ 	self_,
			       wmesh_int_t 			element_,
			       wmesh_int_t 			family_,
			       wmesh_int_t 			degree_);
