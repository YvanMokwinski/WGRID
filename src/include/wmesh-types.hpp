#ifndef WMESH_TYPES_HPP
#define WMESH_TYPES_HPP

#include <ostream>
#include "wmesh-status.h"

template<typename T>
struct wmesh_mat_t
{
  wmesh_int_t 		m;
  wmesh_int_t 		n;
  wmesh_int_t 		ld;
  T * __restrict__ 	v;

  
  static void define(wmesh_mat_t*self_,wmesh_int_t m_,wmesh_int_t n_,T*__restrict__ v_,wmesh_int_t ld_)
  {
    self_->m = m_;
    self_->n = n_;
    self_->v = v_;
    self_->ld = ld_;
  };

};



template<typename T>
struct wmesh_cubature_t
{
  wmesh_int_t 		m_element;
  wmesh_int_t 		m_family;
  wmesh_int_t 		m_degree;
  wmesh_int_t 		m_c_storage;
  wmesh_mat_t<T> 	m_c;
  wmesh_mat_t<T> 	m_w;
};


template<typename T>
wmesh_status_t wmesh_cubature_def(wmesh_cubature_t<T>*__restrict__ self_,
				  wmesh_int_t 		element_,
				  wmesh_int_t 		family_,
				  wmesh_int_t 		degree_);










extern "C"
{
struct wmesh_shape_t
{
  wmesh_int_t		m_element;
  wmesh_int_t		m_family;
  wmesh_int_t 		m_degree;
  wmesh_int_t 		m_ndofs;
};

wmesh_status_t wmesh_shape_def(wmesh_shape_t*__restrict__ self_,
			       wmesh_int_t 		element_,
			       wmesh_int_t 		family_,
			       wmesh_int_t 		degree_);
};



template<typename T>
struct wmesh_shape_eval_t
{
  wmesh_shape_t 	m_shape;


  wmesh_int_t 		m_f_storage;
  wmesh_mat_t<T> 	m_f;
  wmesh_int_t 		m_diff_storage;
  wmesh_mat_t<T> 	m_diff[3];

  wmesh_int_t 		m_wf_storage;
  wmesh_mat_t<T> 	m_wf;
  wmesh_int_t 		m_wdiff_storage;
  wmesh_mat_t<T> 	m_wdiff[3];

};



template<typename T>
wmesh_status_t wmesh_shape_eval_def(wmesh_shape_eval_t<T>*__restrict__ 	self_,
				    wmesh_int_t 			element_,				
				    wmesh_int_t 			shape_family_,
				    wmesh_int_t 			shape_degree_,				
				    wmesh_int_t 			nodes_storage_,
				    const wmesh_mat_t<T> * 		nodes_,
				    const wmesh_mat_t<T> * 		weights_ = nullptr);



extern "C"
{
  
  struct wmesh_int_mat_t : wmesh_mat_t<wmesh_int_t>{};
#if 0
  {
    wmesh_int_t 	m;
    wmesh_int_t 	n;
    wmesh_int_t 	ld;
    wmesh_int_p 	v;
  };
#endif
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
