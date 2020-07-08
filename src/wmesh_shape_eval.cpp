
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_t.hpp"
#include <chrono>
#include <iostream>
#include "bms.h"
#include "wmesh-utils.hpp"
#include "wmesh-blas.h"
#include "bms_templates.hpp"
#include "wmesh_integral_convection_t.hpp"


template<typename T>
static wmesh_status_t wmesh_shape_eval_def_init(wmesh_shape_eval_t<T>*__restrict__ 	self_,
						wmesh_int_t 				nodes_storage_,
						const wmesh_mat_t<T> * 		nodes_)
{
  WMESH_CHECK_POINTER(self_);
  WMESH_CHECK_POINTER(nodes_);

  
  wmesh_int_t iw_n,rw_n;
  wmesh_status_t status;

  const wmesh_int_t element = self_->m_shape.get_element();
  const wmesh_int_t family  = self_->m_shape.get_family();
  const wmesh_int_t degree  = self_->m_shape.get_degree();

  wmesh_int_t topodim;
  status = bms_element2topodim(element,
			       &topodim);  
  WMESH_STATUS_CHECK(status);

  status = bms_shape_buffer_size(element,
				 self_->m_shape.get_family(),
				 self_->m_shape.get_degree(),
				 &iw_n,
				 &rw_n);
  WMESH_STATUS_CHECK(status);

  iw_n = 20000;
  rw_n = 20000;
  wmesh_int_p iw = (iw_n > 0) ? (wmesh_int_p)malloc(sizeof(wmesh_int_t)*iw_n) : nullptr;
  if (iw_n > 0 && !iw)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
    }
  
  T * rw = (rw_n > 0) ? (T*)malloc(sizeof(T)*rw_n) : nullptr;
  if (rw_n > 0 && !rw)
    {
      if (iw)
	{
	  free(iw);	  
	}
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
    }
  
  wmesh_int_t diff[3] = {0,0,0};  
  status = bms_template_shape(element,
			      family,
			      degree,
			    
			      diff,
			    
			      nodes_storage_,			      
			      nodes_->m,
			      nodes_->n,
			      nodes_->v,
			      nodes_->ld,
			      
			      WMESH_STORAGE_INTERLEAVE,			    
			      self_->m_f.m,
			      self_->m_f.n,
			      self_->m_f.v,
			      self_->m_f.ld,
			    
			      iw_n,
			      iw,
			      rw_n,
			      rw);
  if (WMESH_STATUS_SUCCESS != status)
    {
      if (rw)
	{
	  free(rw);
	}
      if (iw)
	{
	  free(iw);
	}
      WMESH_STATUS_CHECK(status);
    }
	
  for (wmesh_int_t j=0;j<topodim;++j)
    {
      diff[j] = 1;
      status = bms_template_shape(element,
				  family,
				  degree,
			     
				  diff,
			     
				  nodes_storage_,
				  nodes_->m,
				  nodes_->n,
				  nodes_->v,
				  nodes_->ld,
				  
				  WMESH_STORAGE_INTERLEAVE,			    
				  self_->m_diff[j].m,
				  self_->m_diff[j].n,
				  self_->m_diff[j].v,
				  self_->m_diff[j].ld,
			      
				  iw_n,
				  iw,
				  rw_n,
				  rw);
      if (WMESH_STATUS_SUCCESS != status)
	{
	  if (rw)
	    {
	      free(rw);
	    }
	  if (iw)
	    {
	      free(iw);
	    }
	  WMESH_STATUS_CHECK(status);
	}
      diff[j] = 0;
    }
	
  if (rw)
    {
      free(rw);
    }
  if (iw)
    {
      free(iw);
    }

  
  return WMESH_STATUS_SUCCESS;  
}


template<typename T>
wmesh_shape_eval_t<T>::wmesh_shape_eval_t(wmesh_int_t 				element_,				
					  wmesh_int_t 				shape_family_,
					  wmesh_int_t 				shape_degree_,				
					  wmesh_int_t 				nodes_storage_,
					  const wmesh_mat_t<T> * 		nodes_,
					  const wmesh_mat_t<T> * 		weights_)
  : m_shape(element_,shape_family_,shape_degree_)
{
  

  wmesh_status_t status;
  const wmesh_int_t
    element 		= element_,
    degree		= shape_degree_,
    num_nodes 		= (nodes_storage_ == WMESH_STORAGE_INTERLEAVE) ? nodes_->n : nodes_->m;
  
  wmesh_int_t
    topodim,
    num_dofs_per_element;
  
  status = bms_element2topodim(element,&topodim);
  WMESH_STATUS_CHECK_EXIT(status);

  status = bms_ndofs(element,
		     degree,
		     &num_dofs_per_element);
  WMESH_STATUS_CHECK_EXIT(status);
  
  //
  //
  //
  {
    T*__restrict__ f_ptr = (T*)malloc(sizeof(T)*num_dofs_per_element*num_nodes);
  
    if (!f_ptr)
      {
	WMESH_STATUS_CHECK_EXIT(WMESH_STATUS_ERROR_MEMORY);
      }
    
    wmesh_mat_t<T>::alloc(&this->m_nabla,num_dofs_per_element, num_nodes * topodim);

    T*__restrict__ diff_ptr[3];
    for (wmesh_int_t j=0;j<topodim;++j)
      {
	diff_ptr[j] = this->m_nabla.v + this->m_nabla.ld * ( num_nodes * j );
      }

    this->m_nabla_storage = WMESH_STORAGE_INTERLEAVE;
    this->m_diff_storage = WMESH_STORAGE_INTERLEAVE;
    if (!f_ptr)
      {
	free(f_ptr);
	f_ptr = nullptr;
      }

    wmesh_mat_t<T>::define(&this->m_f,
			   num_dofs_per_element,
			   num_nodes,
			   f_ptr,
			   num_dofs_per_element);
  
    for (wmesh_int_t j=0;j<topodim;++j)
      {
	wmesh_mat_t<T>::define(&this->m_diff[j],
			       num_dofs_per_element,
			       num_nodes,
			       diff_ptr[j],
			       num_dofs_per_element);
      }
  }
  

  //
  //
  //

  status = wmesh_shape_eval_def_init(this,
				     nodes_storage_,
				     nodes_);
  WMESH_STATUS_CHECK_EXIT(status);
};

template struct wmesh_shape_eval_t<float>;
template struct wmesh_shape_eval_t<double>;

#if 0
template<typename T>
wmesh_status_t wmesh_shape_eval_def(wmesh_shape_eval_t<T>*__restrict__ 	self_,
				    wmesh_int_t 			element_,				
				    wmesh_int_t 			shape_family_,
				    wmesh_int_t 			shape_degree_,				
				    wmesh_int_t 			nodes_storage_,
				    const wmesh_mat_t<T> * 		nodes_,
				    const wmesh_mat_t<T> * 		weights_);


template
wmesh_status_t wmesh_shape_eval_def<float>(wmesh_shape_eval_t<float>*__restrict__ 	self_,
					   wmesh_int_t 			element_,				
					   wmesh_int_t 			shape_family_,
					   wmesh_int_t 			shape_degree_,				
					   wmesh_int_t 			nodes_storage_,
					   const wmesh_mat_t<float> * 		nodes_,
					   const wmesh_mat_t<float> * 		weights_ = nullptr);

template
wmesh_status_t wmesh_shape_eval_def<double>(wmesh_shape_eval_t<double>*__restrict__ 	self_,
					   wmesh_int_t 			element_,				
					   wmesh_int_t 			shape_family_,
					   wmesh_int_t 			shape_degree_,				
					   wmesh_int_t 			nodes_storage_,
					   const wmesh_mat_t<double> * 	nodes_,
					   const wmesh_mat_t<double> * 	weights_ = nullptr);
#endif
