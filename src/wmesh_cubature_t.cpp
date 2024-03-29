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
#include "wmesh_cubature_t.hpp"
#include "wmesh_cubature_info_t.hpp"

template<typename T>
wmesh_cubature_t<T>::wmesh_cubature_t(wmesh_int_t 		element_,
				   wmesh_int_t 		family_,
				   wmesh_int_t 		degree_)
{
  wmesh_int_t rw_n;
  T * __restrict__ rw;
  

  this->m_element 	= element_;
  this->m_family 	= family_;
  this->m_degree 	= degree_;
  wmesh_status_t status;
  
  wmesh_int_t topodim;
  status = bms_element2topodim(element_,
			       &topodim);

  wmesh_int_t num_nodes1d;
  status = bms_cubature_num_nodes(WMESH_ELEMENT_EDGE,
				  family_,
				  degree_,
				  &num_nodes1d);
  WMESH_STATUS_CHECK_EXIT(status);

  wmesh_int_t num_nodes;    
  status = bms_cubature_num_nodes(element_,
				  family_,
				  degree_,
				  &num_nodes);
  WMESH_STATUS_CHECK_EXIT(status);

  this->m_c_storage = WMESH_STORAGE_INTERLEAVE;	      
  
  wmesh_mat_t<T>::define(&this->m_c,
			 topodim,
			 num_nodes,
			 (T*)malloc(sizeof(T) * topodim * num_nodes),
			 topodim);
  
  wmesh_mat_t<T>::define(&this->m_w,
			 1,
			 num_nodes,
			 (T*)malloc(sizeof(T) * 1 * num_nodes),
			 1);
  
  status = bms_cubature_buffer_size(element_,
				    family_,
				    num_nodes1d,			
				    &rw_n);
  WMESH_STATUS_CHECK_EXIT(status);
  
  rw = (rw_n > 0) ? (T*)malloc(sizeof(T)*rw_n) : nullptr;
  if (rw_n > 0 && !rw)
    {
      WMESH_STATUS_CHECK_EXIT(WMESH_STATUS_ERROR_MEMORY);
    }

  status = bms_template_cubature(element_,
				 family_,
				 num_nodes1d,
				 
				 this->m_c_storage,
				 this->m_c.m,
				 this->m_c.n,
				 this->m_c.v,
				 this->m_c.ld,
				 
				 this->m_w.n,
				 this->m_w.v,
				 this->m_w.ld,
				 
				 rw_n,
				 rw);
  
  WMESH_STATUS_CHECK_EXIT(status);
  if (rw)
    {
      free(rw);
    }
  rw = nullptr;  

};

template
struct wmesh_cubature_t<double>;

template
struct wmesh_cubature_t<float>;

#if 0
template<typename T>
wmesh_status_t wmesh_cubature_def(wmesh_cubature_t<T>**__restrict__ 	self__,
				  wmesh_int_t 				element_,
				  wmesh_int_t 				family_,
				  wmesh_int_t 				degree_)
{
  self__[0] = new wmesh_cubature<T>(element_,
				    family_,
				    degree_)
  return WMESH_STATUS_SUCCESS;
}

template
wmesh_status_t wmesh_cubature_def<float>(wmesh_cubature_t<float>*__restrict__ 	self_,
					 wmesh_int_t 				element_,
					 wmesh_int_t 				family_,
					 wmesh_int_t 				degree_);

template
wmesh_status_t wmesh_cubature_def<double>(wmesh_cubature_t<double>*__restrict__ 	self_,
					  wmesh_int_t 					element_,
					  wmesh_int_t 					family_,
					  wmesh_int_t 					degree_);
#endif
