#pragma once

#include "bms.h"

template<typename T>
wmesh_status_t bms_element_geometry(wmesh_int_t 	element_,
				    T*__restrict__ 	c_);
template<typename T>
wmesh_status_t bms_ordering_vertices(wmesh_int_t 	element_,
				     T*__restrict__ 	c_);


template<typename T>
wmesh_status_t bms_ordering_linear_shape(wmesh_int_t 			element_,
					 wmesh_int_t 			c_storage_,
					 wmesh_int_t 			c_m_,
					 wmesh_int_t 			c_n_,
					 const T * __restrict__ 	c_v_,
					 wmesh_int_t 			c_ld_,
					 wmesh_int_t 			ev_storage_,
					 wmesh_int_t 			ev_m_,
					 wmesh_int_t 			ev_n_,
					  T * __restrict__ 	ev_v_,
					 wmesh_int_t 			ev_ld_);


template<typename T>
wmesh_status_t bms_transform(wmesh_int_t 		element_,
			     wmesh_int_t 		r_n_,					 
			     const T*__restrict__	r_,
			     wmesh_int_t		r_inc_,
			     
			     wmesh_int_t		c_storage_,
			     wmesh_int_t		c_m_,
			     wmesh_int_t		c_n_,
			     T * __restrict__		c_v_,					 
			     wmesh_int_t		c_ld_,
			     
			     const_wmesh_int_p 		p_,
			     wmesh_int_t		p_inc_);


template<typename T>
wmesh_status_t bms_sparse_add(wmesh_int_t 		idofs_n_,
			      const_wmesh_int_p 	idofs_,
			      wmesh_int_t 		idofs_inc_,

			      wmesh_int_t 		jdofs_n_,
			      const_wmesh_int_p 	jdofs_,
			      wmesh_int_t 		jdofs_inc_,

			      const T * __restrict__ 	lmat_,
			      wmesh_int_t 		lmat_ld_,
			      
			      wmesh_int_t 		csr_size_,
			      const_wmesh_int_p 	csr_ptr_,
			      const_wmesh_int_p 	csr_ind_,
			      T * __restrict__		csr_val_);
