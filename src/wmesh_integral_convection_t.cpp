
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh.hpp"
#include <chrono>
#include <iostream>
#include "bms.h"
#include "wmesh-utils.hpp"
#include "wmesh-blas.h"
#include "bms_templates.hpp"
#include "wmesh_shape_eval_t.hpp"

#include "wmesh_shape_t.hpp"
#include "wmesh_cubature_t.hpp"
#include "wmesh_shape_eval_t.hpp"

#include "wmesh_cubature_factory_t.hpp"
#include "wmesh_shape_eval_factory_t.hpp"
#include "wmesh_integral_convection_t.hpp"
template<typename T>
wmesh_status_t wmesh_integral_convection_def(wmesh_integral_convection_t<T>*__restrict__ 	self_,
					     wmesh_int_t 					element_,
					     const wmesh_cubature_info_t& 			cubature_info_,					     
					     const wmesh_shape_info_t& 				shape_info_element_,
					     const wmesh_shape_info_t& 				shape_info_velocity_,
					     const wmesh_shape_info_t& 				shape_info_trial_,
					     const wmesh_shape_info_t& 				shape_info_test_)
{

  
  wmesh_status_t status;
  self_->m_element = element_;

  status = wmesh_shape_def(&self_->m_shape_element,
			   element_,
			   shape_info_element_.m_family,
			   shape_info_element_.m_degree);
  WMESH_STATUS_CHECK(status);

  status = wmesh_shape_def(&self_->m_shape_velocity,
			   element_,
			   shape_info_velocity_.m_family,
			   shape_info_velocity_.m_degree);
  WMESH_STATUS_CHECK(status);

  status = wmesh_shape_def(&self_->m_shape_trial,
			   element_,
			   shape_info_trial_.m_family,
			   shape_info_trial_.m_degree);
  WMESH_STATUS_CHECK(status);

  status = wmesh_shape_def(&self_->m_shape_test,
			   element_,
			   shape_info_test_.m_family,
			   shape_info_test_.m_degree);
  WMESH_STATUS_CHECK(status);

  wmesh_int_t topodim;
  status = bms_element2topodim(element_,
			       &topodim);  
  WMESH_STATUS_CHECK(status);
  self_->m_topodim = topodim;

  
  self_->m_cubature 			= wmesh_cubature_factory_t<T>::cubature_instance	(element_,
												 cubature_info_.m_family,
												 cubature_info_.m_degree);
  
  self_->m_shape_eval_velocity  	= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(self_->m_shape_velocity,
												 self_->m_cubature);
  
  self_->m_shape_eval_element  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(self_->m_shape_element,
												 self_->m_cubature);
  
  self_->m_shape_eval_trial  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(self_->m_shape_trial,
												 self_->m_cubature);
  
  self_->m_shape_eval_test  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(self_->m_shape_test,
												 self_->m_cubature);
  
  self_->m_eval_velocity_storage 	= self_->m_shape_eval_velocity->m_f_storage;
  self_->m_eval_velocity 		= &self_->m_shape_eval_velocity->m_f;

  self_->m_eval_test_storage 		= self_->m_shape_eval_test->m_f_storage;
  self_->m_eval_test 			= &self_->m_shape_eval_test->m_f;

  self_->m_eval_trial_nabla_storage 	= self_->m_shape_eval_trial->m_nabla_storage;
  self_->m_eval_trial_nabla 		= &self_->m_shape_eval_trial->m_nabla;
  
  self_->m_build_storage 		= WMESH_STORAGE_INTERLEAVE;
  
  const wmesh_int_t q_n 	= self_->m_cubature->m_w.n;
  const wmesh_int_t test_ndofs 	= self_->m_shape_test.m_ndofs;
  const wmesh_int_t trial_ndofs = self_->m_shape_trial.m_ndofs;

  wmesh_mat_t<T>::alloc(&self_->m_build,
			test_ndofs * trial_ndofs,
			q_n * topodim);

  const wmesh_int_t n1=static_cast<wmesh_int_t>(1);
  for (wmesh_int_t idim=0;idim<topodim;++idim)
    {
      for (wmesh_int_t k=0;k<q_n;++k)
	{
	  const T alpha = self_->m_cubature->m_w.v[self_->m_cubature->m_w.ld*k];
	  xger(&test_ndofs,
	       &trial_ndofs,
	       &alpha,
	       self_->m_eval_test->v + self_->m_eval_test->ld * k,
	       &n1,
	       self_->m_shape_eval_trial->m_diff[idim].v + self_->m_shape_eval_trial->m_diff[idim].ld * k,
	       &n1,
	       self_->m_build.v + self_->m_build.ld * (q_n * idim + k),
	       &test_ndofs);
	}
    }
  
  return WMESH_STATUS_SUCCESS;
}


template
wmesh_status_t wmesh_integral_convection_def<float>(wmesh_integral_convection_t<float>*__restrict__ 	self_,
						    wmesh_int_t 					element_,
						    const wmesh_cubature_info_t& 			cubature_info_,					     
						    const wmesh_shape_info_t& 				shape_info_element_,
						    const wmesh_shape_info_t& 				shape_info_velocity_,
						    const wmesh_shape_info_t& 				shape_info_trial_,
						    const wmesh_shape_info_t& 				shape_info_test_);

template
wmesh_status_t wmesh_integral_convection_def<double>(wmesh_integral_convection_t<double>*__restrict__ 	self_,
						    wmesh_int_t 					element_,
						    const wmesh_cubature_info_t& 			cubature_info_,					     
						    const wmesh_shape_info_t& 				shape_info_element_,
						    const wmesh_shape_info_t& 				shape_info_velocity_,
						    const wmesh_shape_info_t& 				shape_info_trial_,
						    const wmesh_shape_info_t& 				shape_info_test_);


template<typename T>
wmesh_status_t wmesh_integral_convection_data_def(wmesh_integral_convection_data_t<T>*__restrict__  	self_,
						  const wmesh_integral_convection_t<T>& 		parent_)
{
  memset(self_,0,sizeof(wmesh_integral_convection_data_t<T>));
  self_->m_q_velocity_storage 		= WMESH_STORAGE_BLOCK;
  self_->m_q_jacobians_storage 		= WMESH_STORAGE_INTERLEAVE;
  const wmesh_int_t q_n 		= parent_.m_cubature->m_w.n;
  const wmesh_int_t ndofs_velocity 	= parent_.m_shape_velocity.m_ndofs;
  const wmesh_int_t topodim 		= parent_.m_topodim;
  wmesh_mat_t<T>::alloc(&self_->m_q_velocity,q_n,ndofs_velocity);
  wmesh_mat_t<T>::alloc(&self_->m_q_jacobians,topodim*topodim,q_n);
  wmesh_mat_t<T>::alloc(&self_->m_q_jacobians_det,1,q_n);
  return WMESH_STATUS_SUCCESS;
};

template
wmesh_status_t wmesh_integral_convection_data_def<float>(wmesh_integral_convection_data_t<float>*__restrict__  	self_,
							 const wmesh_integral_convection_t<float>& 		parent_);
template
wmesh_status_t wmesh_integral_convection_data_def<double>(wmesh_integral_convection_data_t<double>*__restrict__  	self_,
							 const wmesh_integral_convection_t<double>& 		parent_);

template<typename T>
wmesh_status_t wmesh_integral_convection_eval(const wmesh_integral_convection_t<T>& 	self_,
					      wmesh_integral_convection_data_t<T>&  	self_data_,
					      wmesh_int_t 				dofs_element_storage_,
					      const wmesh_mat_t<T>& 			dofs_element_,
					      wmesh_int_t 				dofs_velocity_storage_,
					      const wmesh_mat_t<T>& 			dofs_velocity_,
					      const T 					alpha_,
					      wmesh_mat_t<T>& 				local_matrix_)
{
  
  //
  // Eval velocity. We want it the result to be block here.
  //
  const wmesh_int_t n1=static_cast<wmesh_int_t>(1);
  const T r1=static_cast<T>(1);
  const T r0 = static_cast<T>(0);
  const wmesh_int_t q_n = self_.m_cubature->m_w.n;
  const wmesh_int_t topodim = self_.m_topodim;
  wmesh_status_t status;
  
  wmesh_mat_gemm(self_.m_shape_eval_velocity->m_f_storage,
		 dofs_velocity_storage_,
		 self_data_.m_q_velocity_storage,
		 
		 static_cast<T>(1),
		 self_.m_shape_eval_velocity->m_f,
		 dofs_velocity_,
		 static_cast<T>(0),
		 self_data_.m_q_velocity);

  
  //
  // Compute the jacobian.
  //
  status =  bms_element_jacobians(self_.m_element,
				  dofs_element_storage_,
				  dofs_element_,
				  self_.m_shape_eval_element->m_diff_storage,
				  &self_.m_shape_eval_element->m_diff[0],
				  self_data_.m_q_jacobians,
				  self_data_.m_q_jacobians_det);
  WMESH_STATUS_CHECK(status);
  
  //
  // Transform the velocity.
  //
  T tmp[3];
  T tmpu[3];
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      const T * J    = self_data_.m_q_jacobians.v + self_data_.m_q_jacobians.ld * k;
      const T * detJ = self_data_.m_q_jacobians_det.v + self_data_.m_q_jacobians_det.ld * k;
      if (self_data_.m_q_velocity_storage==WMESH_STORAGE_BLOCK)
	{
	  for (wmesh_int_t i=0;i<topodim;++i)
	    {
	      tmpu[i] = *(self_data_.m_q_velocity.v + self_data_.m_q_velocity.ld*i+k);
	    }
	}
      else
	{
	  for (wmesh_int_t i=0;i<topodim;++i)
	    {
	      tmpu[i] = *(self_data_.m_q_velocity.v + self_data_.m_q_velocity.ld*k+i);
	    }
	}
      
      xgemv("N",
	    &topodim,
	    &topodim,
	    detJ,
	    J,
	    &topodim,
	    tmpu,
	    &n1,
	    &r0,
	    tmp,
	    &n1);
      
      if (self_data_.m_q_velocity_storage==WMESH_STORAGE_BLOCK)
	{
	  for (wmesh_int_t i=0;i<topodim;++i)
	    {
	      *(self_data_.m_q_velocity.v + self_data_.m_q_velocity.ld*i+k) = tmp[i];
	    }
	}
      else
	{
	  for (wmesh_int_t i=0;i<topodim;++i)
	    {
	      *(self_data_.m_q_velocity.v + self_data_.m_q_velocity.ld*k+i) = tmp[i];
	    }

	}
    }

  
  //
  // Build
  //
  if (local_matrix_.ld == local_matrix_.m)
    {
      xgemv("N",
	    &self_.m_build.m,
	    &self_.m_build.n,
	    &r1,
	    self_.m_build.v,
	    &self_.m_build.ld,
	    self_data_.m_q_velocity.v,
	    &n1,
	    &alpha_,
	    local_matrix_.v,
	    &n1);
    }
  else
    {
      std::cerr << "need a correction." << std::endl;
      exit(1);
    }
  
  return WMESH_STATUS_SUCCESS;
}

template
wmesh_status_t wmesh_integral_convection_eval<float>(const wmesh_integral_convection_t<float>& 	self_,
						     wmesh_integral_convection_data_t<float>&  	self_data_,
						     wmesh_int_t 				dofs_element_storage_,
						     const wmesh_mat_t<float>& 			dofs_element_,
						     wmesh_int_t 				dofs_velocity_storage_,
						     const wmesh_mat_t<float>& 			dofs_velocity_,
						     const float				alpha_,
						     wmesh_mat_t<float>& 			local_matrix_);

template
wmesh_status_t wmesh_integral_convection_eval<double>(const wmesh_integral_convection_t<double>& 	self_,
						      wmesh_integral_convection_data_t<double>&  	self_data_,
						      wmesh_int_t 				dofs_element_storage_,
						      const wmesh_mat_t<double>& 			dofs_element_,
						      wmesh_int_t 				dofs_velocity_storage_,
						      const wmesh_mat_t<double>& 			dofs_velocity_,
						      const double				alpha_,
						     wmesh_mat_t<double>& 			local_matrix_);
