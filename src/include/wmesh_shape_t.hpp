#pragma once

#include "wmesh_shape_info_t.hpp"
#include "wmesh_shape_t.h"


struct wmesh_shape_t
{
private:
  
  wmesh_int_t		m_element;
  wmesh_int_t		m_family;
  wmesh_int_t		m_degree;
  wmesh_int_t		m_ndofs;
  wmesh_int_t		m_nodes_family;
  
public:

  inline   wmesh_int_t		get_element() 		const { return this->m_element;}
  inline   wmesh_int_t		get_family() 		const { return this->m_family;}
  inline   wmesh_int_t		get_degree() 		const { return this->m_degree;}
  inline   wmesh_int_t		get_ndofs() 		const { return this->m_ndofs;}
  inline   wmesh_int_t		get_nodes_family() 	const { return this->m_nodes_family;}  

  wmesh_shape_t(wmesh_int_t 			element_,
		wmesh_int_t 			family_,
		wmesh_int_t 			degree_);

  inline wmesh_shape_t(wmesh_int_t 			element_,
		       const wmesh_shape_info_t& 	shape_info_)
    : wmesh_shape_t(element_,shape_info_.get_family(),shape_info_.get_degree())
  {

  };
};
