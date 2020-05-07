
#include "wmesh-types.hpp"
#include <stdlib.h>
#include <string.h>
#include "wmesh-status.h"
extern "C"
{
  

  wmesh_status_t wmesh_int_sparsemat_fprintf(const wmesh_int_sparsemat_t* 		self_,
					  FILE * out_)
  {
    if (!self_ || !out_)
      {
	return WMESH_STATUS_INVALID_POINTER;
      }
    
    for (wmesh_int_t k=0;k<self_->m_size;++k)
      {
	if (self_->m_n[k]>0)
	  {
	    for (wmesh_int_t j=0;j<self_->m_n[k];++j)
	      {
		for (wmesh_int_t i=0;i<self_->m_m[k];++i)
		  {
		    fprintf(out_," " WMESH_INT_FORMAT,self_->m_data[self_->m_ptr[k] + (self_->m_ld[k]*j+i) ]);
		  }
		fprintf(out_,"\n");
	      }
	  }
      }
    
    return WMESH_STATUS_SUCCESS;
  };

  wmesh_status_t wmesh_int_sparsemat_get(wmesh_int_sparsemat_t* self_,
					 wmesh_int_t idx_,
					 wmesh_int_mat_t* wint_mat_)
  {
    wmesh_int_mat_def(wint_mat_,
		      self_->m_m[idx_],
		      self_->m_n[idx_],
		      &self_->m_data[self_->m_ptr[idx_]],
		 self_->m_ld[idx_]);
    return WMESH_STATUS_SUCCESS;
  }
  
  wmesh_status_t wmesh_int_sparsemat_new(wmesh_int_sparsemat_t*self_,
					 WMESH_INT_SPARSEMAT_PARAMS)
  {
    self_->m_size = size_;
    self_->m_m = (wmesh_int_t*)malloc(sizeof(wmesh_int_t)*self_->m_size);
    self_->m_n = (wmesh_int_t*)malloc(sizeof(wmesh_int_t)*self_->m_size);
    self_->m_ld = (wmesh_int_t*)malloc(sizeof(wmesh_int_t)*self_->m_size);
    self_->m_ptr = (wmesh_int_t*)malloc(sizeof(wmesh_int_t)*(self_->m_size+1));
    self_->m_data = data_;
    for (wmesh_int_t i=0;i<self_->m_size;++i)
      {
	self_->m_m[i] = m_[i];
      }
    
    for (wmesh_int_t i=0;i<self_->m_size;++i)
      {
	self_->m_n[i] = n_[i];
      }
    
    for (wmesh_int_t i=0;i<self_->m_size;++i)
      {
	self_->m_ld[i] = ld_[i];
      }

    for (wmesh_int_t i=0;i<=self_->m_size;++i)
      {
	self_->m_ptr[i] = ptr_[i];
      }
    
    return WMESH_STATUS_SUCCESS;
  }
  
  wmesh_status_t wmesh_int_sparsemat_free(wmesh_int_sparsemat_t*self_)
  {
    if (nullptr != self_)
      {
	if (nullptr != self_->m_ptr)
	  {
	    free(self_->m_ptr);
	  }
	if (nullptr != self_->m_ld)
	  {
	    free(self_->m_ld);
	  }
	if (nullptr != self_->m_n)
	  {
	    free(self_->m_m);
	  }
	if (nullptr != self_->m_m)
	  {
	    free(self_->m_m);
	  }
	memset(self_,0,sizeof(wmesh_int_sparsemat_t));
      }	  
    return WMESH_STATUS_SUCCESS;
  }
  

}
