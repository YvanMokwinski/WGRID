#include "wmesh-types.hpp"
#include "bms.h"
#include <string.h>

extern "C"
{

  wmesh_status_t wmesh_shape_info_def(wmesh_shape_info_t*__restrict__ 	self_,
				      wmesh_int_t 			family_,
				      wmesh_int_t 			degree_)
  {
    WMESH_CHECK_POINTER(self_);
    memset(self_,0,sizeof(wmesh_shape_info_t));
    self_->m_family 	= family_;
    self_->m_degree 	= degree_;
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t wmesh_shape_def(wmesh_shape_t*__restrict__ self_,
				 wmesh_int_t 		element_,
				 wmesh_int_t 		family_,
				 wmesh_int_t 		degree_)
  {
    WMESH_CHECK_POINTER(self_);
    memset(self_,0,sizeof(wmesh_shape_t));
    self_->m_element 	= element_;
    self_->m_family 	= family_;
    self_->m_degree 	= degree_;

    wmesh_status_t status = bms_ndofs(element_,
				      degree_,
				      &self_->m_ndofs);
    WMESH_STATUS_CHECK(status);    
    switch(family_)
      {
      case WMESH_SHAPE_FAMILY_LAGRANGE:
	{
	  self_->m_nodes_family = WMESH_NODES_FAMILY_LAGRANGE;
	  break;
	}
      case WMESH_SHAPE_FAMILY_LEGENDRE:
	{
	  self_->m_nodes_family = WMESH_NODES_FAMILY_GAUSSLOBATTO;
	  break;
	}
      case WMESH_SHAPE_FAMILY_ORTHOGONAL:
	{
	  self_->m_nodes_family = WMESH_NODES_FAMILY_GAUSSLOBATTO;
	  break;
	}
      }
    
    return WMESH_STATUS_SUCCESS;
  }
  
};


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



template<typename T>
wmesh_status_t wmesh_shape_calculate_eval(const wmesh_shape_t& 		shape_,
					  wmesh_int_t 			nodes_storage_,
					  const wmesh_mat_t<T>& 	nodes_,
					  wmesh_mat_t<T>& 		eval_)
{
  
  wmesh_status_t status;
  
  const wmesh_int_t
    num_nodes 		= (nodes_storage_ == WMESH_STORAGE_INTERLEAVE) ? nodes_.n : nodes_.m;

#if 0
  wmesh_mat_t<T>::alloc(&eval_,
			shape_.m_ndofs,
			num_nodes);
#endif
  
  wmesh_int_t iw_n,rw_n;  
  status = bms_shape_buffer_size(shape_.m_element,
				 shape_.m_family,
				 shape_.m_degree,
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
  status = bms_template_shape(shape_.m_element,
			      shape_.m_family,
			      shape_.m_degree,
			    
			      diff,
			    
			      nodes_storage_,			      
			      WMESH_MAT_FORWARD(nodes_),
			      
			      WMESH_STORAGE_INTERLEAVE,			    
			      WMESH_MAT_FORWARD(eval_),
			      
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


template
wmesh_status_t wmesh_shape_calculate_eval<float>(const wmesh_shape_t& 		shape_,
					  wmesh_int_t 			nodes_storage_,
					  const wmesh_mat_t<float>& 	nodes_,
					  wmesh_mat_t<float>& 		eval_);
template
wmesh_status_t wmesh_shape_calculate_eval<double>(const wmesh_shape_t& 		shape_,
					  wmesh_int_t 			nodes_storage_,
					  const wmesh_mat_t<double>& 	nodes_,
					  wmesh_mat_t<double>& 		eval_);
