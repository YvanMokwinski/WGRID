#ifndef WMESH_TYPES_HPP
#define WMESH_TYPES_HPP

#include "wmesh-types.h"
extern "C"
{

  struct wmesh_int_sparsemat_t
  {
    wmesh_int_t 	m_size;
    wmesh_int_p  	m_m;
    wmesh_int_p  	m_n;
    wmesh_int_p  	m_ld;
    wmesh_int_p  	m_ptr;
    wmesh_int_p  	m_data;
 };
  
  struct wmesh_int_mat_t
  {
    wmesh_int_t m;
    wmesh_int_t n;
    wmesh_int_t ld;
    wmesh_int_p v;
  };
  
  struct wspace_c2n_t
  {
    wmesh_int_mat_t 	m_c2n;
    wmesh_int_t  	m_cell_type;
  };
  
  struct wspace_sparse_c2n_t
  {
    wmesh_int_sparsemat_t 	m_data;
    wmesh_int_t			m_nzones;
    wspace_c2n_t * 		m_zones;
  };

};

#endif
