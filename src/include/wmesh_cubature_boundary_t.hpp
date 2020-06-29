#if 0

#pragma once

#include "wmesh-types.hpp"
#include "wmesh_cubature_t.hpp"

template<typename T>
struct wmesh_cubature_boundary_t
{
  wmesh_int_t 		m_cubature_family;
  wmesh_int_t 		m_cubature_degree;
  wmesh_int_t 		m_element;
  wmesh_int_t           m_num_facets;
  wmesh_int_t 		m_facets[6];
  wmesh_cubature_t<T>	m_facets_cubature[6][8];
};

template<typename T>
wmesh_status_t wmesh_cubature_boundary_def(wmesh_cubature_boundary_t<T>*__restrict__ 	self_,
					   wmesh_int_t 				element_,
					   wmesh_int_t 				cubature_family_,
					   wmesh_int_t 				cubature_degree_);

#endif
