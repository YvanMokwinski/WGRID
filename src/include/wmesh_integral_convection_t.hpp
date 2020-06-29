#pragma once

#include "wmesh_shape_t.hpp"
#include "wmesh_cubature_t.hpp"
#include "wmesh_shape_eval_t.hpp"
#include "wmesh_cubature_info_t.hpp"

template<typename T>
struct wmesh_integral_convection_t
{
  wmesh_int_t 			m_element;
  wmesh_int_t 			m_topodim;
  wmesh_shape_t 		m_shape_element;
  wmesh_shape_t 		m_shape_velocity;
  wmesh_shape_t 		m_shape_trial;
  wmesh_shape_t 		m_shape_test;
  
  const wmesh_cubature_t<T> * 	m_cubature;
  const wmesh_shape_eval_t<T> * m_shape_eval_element;
  const wmesh_shape_eval_t<T> * m_shape_eval_velocity;
  const wmesh_shape_eval_t<T> * m_shape_eval_trial;
  const wmesh_shape_eval_t<T> * m_shape_eval_test;
  
  wmesh_int_t 			m_eval_trial_nabla_storage;
  const wmesh_mat_t<T> * 	m_eval_trial_nabla;
  
  wmesh_int_t 			m_eval_velocity_storage;
  const wmesh_mat_t<T> * 	m_eval_velocity;
  wmesh_int_t 			m_eval_test_storage;
  const wmesh_mat_t<T> * 	m_eval_test;
  
  wmesh_int_t 			m_build_storage;
  wmesh_mat_t<T> 		m_build;
};

template<typename T>
wmesh_status_t wmesh_integral_convection_def(wmesh_integral_convection_t<T>*__restrict__ self_,
					     wmesh_int_t 			element_,
					     const wmesh_cubature_info_t& 		cubature_info_,
					     const wmesh_shape_info_t& 		shape_info_element_,
					     const wmesh_shape_info_t& 		shape_info_velocity_,
					     const wmesh_shape_info_t& 		shape_info_trial_,
					     const wmesh_shape_info_t& 		shape_info_test_);

template<typename T>
struct
wmesh_integral_convection_data_t
{

  wmesh_int_t 				m_q_velocity_storage;
  wmesh_mat_t<T> 			m_q_velocity;
  wmesh_int_t 				m_q_jacobians_storage;
  wmesh_mat_t<T> 			m_q_jacobians;
  wmesh_mat_t<T> 			m_q_jacobians_det;
};


template<typename T>
wmesh_status_t wmesh_integral_convection_data_def(wmesh_integral_convection_data_t<T>*__restrict__  self_,
						  const wmesh_integral_convection_t<T>& 	parent_);

template<typename T>
wmesh_status_t wmesh_integral_convection_eval(const wmesh_integral_convection_t<T>& 	self_,
					      wmesh_integral_convection_data_t<T>&  	self_data_,
					      wmesh_int_t 				dofs_element_storage_,
					      const wmesh_mat_t<T>& 			dofs_element_,
					      wmesh_int_t 				dofs_velocity_storage_,
					      const wmesh_mat_t<T>& 			dofs_velocity_,
					      const T 					alpha_,
					      wmesh_mat_t<T>& 				local_matrix_);
