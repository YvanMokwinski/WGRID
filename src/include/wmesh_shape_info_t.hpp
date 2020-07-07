#pragma once
#include "wmesh-types.hpp"
#include "wmesh_shape_info_t.h"
#include "wmesh-types.h"

struct wmesh_shape_info_t
{
  wmesh_int_t	m_family{};
  wmesh_int_t 	m_degree{};
  wmesh_shape_info_t()
  {
  };

  wmesh_shape_info_t(wmesh_int_t 			family_,
		     wmesh_int_t 			degree_)
    : m_family(family_),
      m_degree(degree_)
  {
  };

  wmesh_shape_info_t(const wmesh_shape_info_t& that_)
    : m_family(that_.m_family),
      m_degree(that_.m_degree)
  {
  };
  
  inline wmesh_int_t get_family() const { return this->m_family; };
  inline wmesh_int_t get_degree() const { return this->m_degree; };
};


