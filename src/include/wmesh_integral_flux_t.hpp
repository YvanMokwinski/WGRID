#pragma once

#include "wmesh_shape_t.hpp"
#include "wmesh_cubature_t.hpp"
#include "wmesh_shape_eval_t.hpp"
#include "wmesh_cubature_info_t.hpp"
#include "wmesh_shape_info_t.hpp"
  
template<typename T>
struct wmesh_template_integral_flux_t
{
  struct data_t
  {    
    wmesh_int_t		m_dofs_element_storage;
    wmesh_mat_t<T> 	m_dofs_element;

    wmesh_int_t		m_dofs_velocity_storage;
    wmesh_mat_t<T> 	m_dofs_velocity;
    
    wmesh_int_t		m_q_velocity_storage;
    wmesh_mat_t<T> 	m_q_velocity;

    wmesh_int_t 	m_q_element_diff_storage;
    wmesh_mat_t<T>	m_q_element_diff[2];

    wmesh_mat_t<T> 	m_build_rhs;
    wmesh_mat_t<T> 	m_local_matrix;
    wmesh_mat_t<T> 	m_local_rhs;
    wmesh_mat_t<T>	m_q_a;
    wmesh_mat_t<T>	m_q_coo;

    
#if 0
    wmesh_int_t		m_facet;
    wmesh_int_t		m_facet_type;
    
    wmesh_int_t 	m_ielement_type;
    wmesh_int_t		m_ifacet_idx;
    wmesh_int_t		m_ifacet_rotation;

    wmesh_int_t 	m_jelement_type;
    wmesh_int_t		m_jfacet_idx;
    wmesh_int_t		m_jfacet_rotation;
#endif
    
    data_t(const wmesh_template_integral_flux_t<T>&  parent_);    
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
    
  wmesh_int_t 			m_build_residual_storage;
  wmesh_mat_t<T> 		m_build_residual;

  wmesh_template_integral_flux_t(wmesh_int_t 				element_,
				 const wmesh_cubature_info_t& 		cubature_info_,
				 const wmesh_shape_info_t& 		shape_info_element_,
				 const wmesh_shape_info_t& 		shape_info_velocity_,
				 const wmesh_shape_info_t& 		shape_info_trial_,
				 const wmesh_shape_info_t& 		shape_info_test_); 
  
  wmesh_status_t eval(data_t&  					data_,
		      const T 					alpha_,
		      wmesh_mat_t<T>& 				local_matrix_) const ;

  

  wmesh_status_t eval_residual(wmesh_template_integral_flux_t<T>::data_t&  		data_,
								  const T 						alpha_,
								  wmesh_mat_t<T>& 					local_rhs_,
								  T * 							q_a_,
								  wmesh_int_t						q_a_inc_) const;

};
