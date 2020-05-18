#pragma once

template<typename T>
wmesh_status_t bms_transform(wmesh_int_t 		element_,
			     wmesh_int_t 		r_n_,					 
			     const T*__restrict__	r_,
			     wmesh_int_t		r_inc_,
			     
			     wmesh_int_t		b_storage_,
			     wmesh_int_t		b_m_,
			     wmesh_int_t		b_n_,
			     const_wmesh_int_p		b_v_,
			     wmesh_int_t		b_ld_,
			     
			     wmesh_int_t		c_storage_,
			     wmesh_int_t		c_m_,
			     wmesh_int_t		c_n_,
			     T * __restrict__		c_v_,					 
			     wmesh_int_t		c_ld_,
			     
			     wmesh_int_t 		work_n_,
			     wmesh_int_p		work_);