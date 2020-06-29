#pragma once


#include "wmesh_shape_t.hpp"

template<typename T>
struct wmesh_shape_eval_t
{
  wmesh_shape_t 	m_shape;

  wmesh_int_t 		m_f_storage;
  wmesh_mat_t<T> 	m_f;

  wmesh_int_t 		m_nabla_storage;
  wmesh_mat_t<T> 	m_nabla;

  wmesh_int_t 		m_diff_storage;
  wmesh_mat_t<T> 	m_diff[3];

};

template<typename T>
wmesh_status_t wmesh_shape_eval_def(wmesh_shape_eval_t<T>*__restrict__ 	self_,
				    wmesh_int_t 			element_,				
				    wmesh_int_t 			shape_family_,
				    wmesh_int_t 			shape_degree_,				
				    wmesh_int_t 			nodes_storage_,
				    const wmesh_mat_t<T> * 		nodes_,
				    const wmesh_mat_t<T> * 		weights_ = nullptr);

