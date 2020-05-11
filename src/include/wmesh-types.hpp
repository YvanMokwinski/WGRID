#ifndef WMESH_TYPES_HPP
#define WMESH_TYPES_HPP


#include "wmesh-status.h"
extern "C"
{
  
  struct wmesh_int_mat_t
  {
    wmesh_int_t 	m;
    wmesh_int_t 	n;
    wmesh_int_t 	ld;
    wmesh_int_p 	v;
  };

  struct wmesh_int_sparsemat_t
  {
    wmesh_int_t 	m_size;
    wmesh_int_p  	m_ptr;
    wmesh_int_p  	m_m;
    wmesh_int_p  	m_n;
    wmesh_int_p  	m_data;
    wmesh_int_p  	m_ld;    
 };

#define WMESH_INT_SPARSEMAT_FORWARD(_f)				\
  (_f).m_size,							\
    (_f).m_ptr,							\
    (_f).m_m,							\
    (_f).m_n,							\
    (_f).m_data,							\
    (_f).m_ld
  
#define WMESH_INT_SPARSEMAT_CONST_PARAMS				\
  wmesh_int_t 		size_,						\
    const_wmesh_int_p  	ptr_,						\
    const_wmesh_int_p  	m_,						\
    const_wmesh_int_p  	n_,						\
    const_wmesh_int_p  	data_,						\
    const_wmesh_int_p  	ld_
  
#define WMESH_INT_SPARSEMAT_PARAMS		\
  wmesh_int_t 		size_,			\
    const_wmesh_int_p  	ptr_,			\
    const_wmesh_int_p  	m_,			\
    const_wmesh_int_p  	n_,			\
    wmesh_int_p  	data_,			\
    const_wmesh_int_p  	ld_


  
  wmesh_status_t wmesh_int_sparsemat_init(wmesh_int_t 		size_,
					  wmesh_int_p		ptr_,
					  const_wmesh_int_p	m_,
					  const_wmesh_int_p	n_,
					  wmesh_int_p		ld_);

  wmesh_status_t wmesh_int_sparsemat_new(wmesh_int_sparsemat_t* self_,
					 WMESH_INT_SPARSEMAT_PARAMS );
  
  wmesh_status_t wmesh_int_sparsemat_get(wmesh_int_sparsemat_t* self_,
					 wmesh_int_t 		idx_,
					 wmesh_int_mat_t* 	wint_mat_);
  
  wmesh_status_t wmesh_int_sparsemat_free(wmesh_int_sparsemat_t*self_);    
  wmesh_status_t wmesh_int_sparsemat_info(const wmesh_int_sparsemat_t* self_,
					  FILE * out_);
  
  wmesh_status_t wmesh_int_sparsemat_fprintf(const wmesh_int_sparsemat_t* 		self_,
					     FILE * out_);
  
  wmesh_status_t wmesh_int_mat_def	(wmesh_int_mat_t * 	self_,
					 wmesh_int_t   		m,
					 wmesh_int_t 		n,
					 wmesh_int_p 		v,
					 wmesh_int_t 		ld);
  
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
