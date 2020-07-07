#pragma once

#include "wmeshspacedg.h"
#include "wmeshspace_t.hpp"

struct wmeshspacedg_t
{
  wmesh_t * 		m_mesh;
  wmesh_int_t 		m_nodes_family;
  wmesh_int_t 		m_degree;
  wmesh_t * 		m_patterns[4];
  wmesh_int_t 		m_ndofs;
  wmesh_int_t 		m_dofs_ptr[4+1];
  wmesh_int_t 		m_dofs_m[4];
};
