#pragma once

#include "wmesh-types.hpp"
#include "wmesh_shape_t.hpp"

template<typename T>
const wmesh_mat_t<T>* wmesh_calculate_shape_restrict(wmesh_int_t 					element_,
						     wmesh_int_t 					ifacet_,
						     wmesh_int_t 					signed_rotation_,
						     
						     wmesh_int_t 					source_shape_family_,
						     wmesh_int_t 					source_shape_degree_,
						     
						     wmesh_int_t 					target_nodes_family_,
						     wmesh_int_t 					target_nodes_degree_);

#if 0

//
// Definition of the shape cell basis restriction over facets.
//
template<typename T>
struct wmesh_shape_restrict_t
{
  wmesh_shape_t 	m_shape;
  wmesh_int_t 		m_num_facets;
  wmesh_int_t 		m_facets_num_dofs[6];
  wmesh_int_t 		m_facets_num_nodes[6];
  wmesh_int_t 		m_facets[6];
  wmesh_int_t 		m_facet_types[6];  
  const wmesh_mat_t<T>	* m_restrict[6][4];
};


#include <iostream>
template<typename T>
inline wmesh_status_t wmesh_shape_restrict_info(const wmesh_shape_restrict_t<T>&self_)
{
  std::cout << "INFO SHAPE RESTRICT" << std::endl;
  std::cout << " - shape " << self_.m_shape.m_element << std::endl;
  std::cout << " - shape_family  " << self_.m_shape.m_family << std::endl;
  std::cout << " - shape_degree  " << self_.m_shape.m_degree << std::endl;
  std::cout << " - shape_ndofs   " << self_.m_shape.m_ndofs << std::endl;
  std::cout << " - num_facets    " << self_.m_num_facets << std::endl;
  for (wmesh_int_t i=0;i<self_.m_num_facets;++i)
    {
      
      std::cout << "   - facet local idx " << i << std::endl;
      std::cout << "     - facet element        "  << self_.m_facets[i] << std::endl;
      std::cout << "     - facet type           "  << self_.m_facet_types[i] << std::endl;
      std::cout << "     - facet num nodes      "  << self_.m_facets_num_nodes[i] << std::endl;
      std::cout << "     - facet num dofs       "  << self_.m_facets_num_dofs[i] << std::endl;
      
    }
  return WMESH_STATUS_SUCCESS;
}

template<typename T>
wmesh_status_t wmesh_shape_restrict_def(wmesh_shape_restrict_t<T>*__restrict__ 	self_,
					wmesh_int_t 				element_,
					wmesh_int_t 				shape_family_,
					wmesh_int_t 				shape_degree_);
template<typename T>
inline const wmesh_mat_t<T>& wmesh_shape_restrict_get(wmesh_shape_restrict_t<T>*__restrict__ 	self_,
						      wmesh_int_t				ifacet_,
						      wmesh_int_t				signed_rotation_)
{
  if (signed_rotation_ > 0)
    {
      return self_->m_restrict[ self_->m_facet_types[ifacet_] ][signed_rotation_-1];
    }
  else 
    {
      return self_->m_restrict[ self_->m_facet_types[ifacet_] ][-signed_rotation_-1];
    }
}

#endif
