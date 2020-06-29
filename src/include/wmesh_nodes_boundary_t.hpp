#pragma once

#include "wmesh-types.hpp"

template<typename T>
struct wmesh_nodes_boundary_t
{

  wmesh_int_t 		m_element;
  wmesh_int_t 		m_nodes_family;
  wmesh_int_t 		m_nodes_degree;
  
  wmesh_int_t 		m_num_facets;  
  wmesh_int_t 		m_facets[6];
  wmesh_int_t 		m_facets_nodes_storage;
  wmesh_mat_t<T>	m_facets_nodes[6][8];
  
};

template<typename T>
wmesh_status_t wmesh_nodes_boundary_def(wmesh_nodes_boundary_t<T>*__restrict__ 	self_,
					wmesh_int_t 				element_,
					wmesh_int_t 				nodes_family_,
					wmesh_int_t 				nodes_degree_);
