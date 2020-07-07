#pragma once

#include "wmesh_shape_t.hpp"
#include "wmesh_cubature_t.hpp"
#include "wmesh_shape_eval_t.hpp"
#include "wmesh_cubature_info_t.hpp"
#include "wmesh_shape_info_t.hpp"
  
template<typename T>
struct wmesh_template_integral_convection_t
{
  struct data_t
  {    
    wmesh_int_t		m_dofs_element_storage;
    wmesh_mat_t<T> 	m_dofs_element;

    wmesh_int_t		m_dofs_velocity_storage;
    wmesh_mat_t<T> 	m_dofs_velocity;
    
    wmesh_int_t		m_q_velocity_storage;
    wmesh_mat_t<T> 	m_q_velocity;

    wmesh_int_t 	m_q_jacobians_storage;
    wmesh_mat_t<T>	m_q_jacobians;

    wmesh_mat_t<T>	m_q_jacobians_det;
    wmesh_mat_t<T> 	m_build_rhs;
    wmesh_mat_t<T> 	m_local_matrix;
    T*__restrict__	m_local_rhs;
    
    data_t(const wmesh_template_integral_convection_t<T>&  parent_);
    
  };

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

  wmesh_int_t 			m_build_storage;
  wmesh_mat_t<T> 		m_build;
  
  wmesh_template_integral_convection_t(wmesh_int_t 				element_,
				       const wmesh_cubature_info_t& 		cubature_info_,
				       const wmesh_shape_info_t& 		shape_info_element_,
				       const wmesh_shape_info_t& 		shape_info_velocity_,
				       const wmesh_shape_info_t& 		shape_info_trial_,
				       const wmesh_shape_info_t& 		shape_info_test_); 
  
  wmesh_status_t eval(data_t&  					data_,
		      const T 					alpha_,
		      wmesh_mat_t<T>& 				local_matrix_) const ;
};
