#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh.hpp"

#include <chrono>
#include <iostream>
#include "wmesh.hpp"
#include "GenericEncoding.hpp"
#include "wmesh_utils.hpp"
#include "wmesh-utils.hpp"

using namespace std::chrono;

wmesh_status_t wmesh_sfc_fill333(wmesh_int_t 			c2n_m_,
				 wmesh_int_t 			c2n_n_,
				 const_wmesh_int_p		c2n_v_,
				 wmesh_int_t 			c2n_ld_,
				 
				 wmesh_int_t 			coo_m_,
				 wmesh_int_t 			coo_n_,
				 const double *__restrict__ 	coo_v_,
				 wmesh_int_t 			coo_ld_,
				 
				 unsigned long long int*		sfc_v_,
				 wmesh_int_t 			sfc_ld_,
				 
			      const double *__restrict__ 	box_,
			      size_t 				work_size_,
			      void * 				work_)
{
  
  if (work_size_ < sizeof(double)*coo_m_ )
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }
  
  double * center = (double*)work_;
  for (wmesh_int_t j=0;j<c2n_n_;++j)
    {
      for (wmesh_int_t k=0;k<coo_m_;++k)
	{		  
	  center[k] = 0.0;
	}
      
      for (wmesh_int_t i=0;i<c2n_m_;++i)
	{
	  wmesh_int_t v = c2n_v_[j*c2n_ld_ + i];	  
	  for (wmesh_int_t k=0;k<coo_m_;++k)
	    {		  
	      center[k] += coo_v_[coo_ld_*(v-1)+k];
	    }
	}
      
      for (wmesh_int_t k=0;k<coo_m_;++k)
	{		  
	  center[k] /= c2n_m_;
	}
      
      sfc_v_[sfc_ld_ * j] = hilbert_coordinate(center,
					       box_,
					       31);
    }
  return WMESH_STATUS_SUCCESS;
}

template <typename T>
wmesh_status_t wmesh_sfc_fill(wmesh_int_t 			c2n_size_,
			      const_wmesh_int_p			c2n_ptr_,
			      const_wmesh_int_p			c2n_m_,
			      const_wmesh_int_p			c2n_n_,
			      const_wmesh_int_p			c2n_v_,
			      const_wmesh_int_p			c2n_ld_,
			      
			      wmesh_int_t 			coo_m_,
			      wmesh_int_t 			coo_n_,
			      const T *__restrict__ 		coo_v_,
			      wmesh_int_t 			coo_ld_,
			      
			      unsigned long long int*		sfc_v_,
			      wmesh_int_t 			sfc_ld_,			      
			      const T *__restrict__ 		box_)
{
  T center[3];
  unsigned long long int idx = 0;
  for (wmesh_int_t l=0;l<c2n_size_;++l)
    {
      for (wmesh_int_t j=0;j<c2n_n_[l];++j)
	{	  
	  for (wmesh_int_t k=0;k<coo_m_;++k)
	    {		  
	      center[k] = 0.0;
	    }
      
	  for (wmesh_int_t i=0;i<c2n_m_[l];++i)
	    {
	      wmesh_int_t v = c2n_v_[c2n_ptr_[l] + j*c2n_ld_[l] + i] - 1;	  
	      for (wmesh_int_t k=0;k<coo_m_;++k)
		{		  
		  center[k] += coo_v_[coo_ld_* v +k];
		}
	    }
      
	  for (wmesh_int_t k=0;k<coo_m_;++k)
	    {		  
	      center[k] /= c2n_m_[l];
	    }
      
	  sfc_v_[ idx * sfc_ld_ + 0 ] = hilbert_coordinate(center,
							   box_,
							   31);
	  ++idx;
	}
    }
  return WMESH_STATUS_SUCCESS;
}


extern "C"
{
  wmesh_status_t bms_partitioning_sfc(wmesh_int_t 				nparts_,
				      wmesh_int_t 				num_cells_,
				      wmesh_int_p 				p_v_,				      
				      wmesh_int_t  				p_ld_,				      

				      wmesh_int_t 				c2n_size_,
				      const_wmesh_int_p 			c2n_ptr_,
				      const_wmesh_int_p 			c2n_m_,
				      const_wmesh_int_p 			c2n_n_,
				      const_wmesh_int_p 			c2n_v_,
				      const_wmesh_int_p 			c2n_ld_,
				      
				      wmesh_int_t 				coo_m_,
				      wmesh_int_t 				coo_n_,
				      const double *__restrict__ 		coo_v_,
				      wmesh_int_t 				coo_ld_,
				      
				      size_t*__restrict__                       work_size_,
				      void *__restrict__                        work_)
  {
    if (num_cells_ == 0)
      {
	return WMESH_STATUS_SUCCESS;
      }
    
    WMESH_CHECK_POINTER(p_v_);
    WMESH_CHECK_POINTER(c2n_ptr_);
    WMESH_CHECK_POINTER(c2n_m_);
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(c2n_v_);
    WMESH_CHECK_POINTER(c2n_ld_);
    WMESH_CHECK_POINTER(coo_v_);
    WMESH_CHECK_POINTER(work_size_);
    
    if (work_size_[0] < 0)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
    if (!work_ && work_size_[0] == 0)
      {
	work_size_[0] = sizeof(unsigned long long int)*2*num_cells_;
	return WMESH_STATUS_SUCCESS;
      }
    if (!work_ && work_size_[0] > 0)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
      }
    if (work_size_[0] < sizeof(unsigned long long int)*2*num_cells_)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
      }

    
    unsigned long long int * __restrict__ p = (unsigned long long int * __restrict__)work_;
    for (wmesh_int_t i=0;i<num_cells_;++i)
      {
	p[2*i+1] = i;
      }

    wmesh_status_t status;
    double box[6];
    status = wmesh_raw_bounding_box(coo_m_,
				    coo_n_,
				    coo_v_,
				    coo_ld_,
				    box);

    WMESH_STATUS_CHECK(status);
    
    status = wmesh_sfc_fill(c2n_size_,
			    c2n_ptr_,
			    c2n_m_,
			    c2n_n_,
			    c2n_v_,
			    c2n_ld_,
			    
			    coo_m_,
			    coo_n_,
			    coo_v_,
			    coo_ld_,
			    
			    p,
			    2,			      
			    box);

    WMESH_STATUS_CHECK(status);

    //
    // Sort.
    //
    qsort(p,
	  num_cells_,
	  sizeof(unsigned long long int)*2,
	  wmesh_qsort_increasing_predicate<unsigned long long int>);
    
    {
      wmesh_int_t idx = 0;
      for (wmesh_int_t i=0;i<nparts_;++i)
	{
	  for (wmesh_int_t j=0;j<num_cells_ / nparts_;++j)
	    {
	      p[2*idx+0] = i;
	      ++idx;
	    }
	}
      for (wmesh_int_t j=0;j<num_cells_ % nparts_;++j)
	{
	  p[2*idx+0] = nparts_-1;
	  ++idx;
	}
    }
    
    for (wmesh_int_t i=0;i<num_cells_;++i)
      {
	p_v_[p_ld_ * p[2*i+1] + 0] = p[2*i+0];
      }
    
    return WMESH_STATUS_SUCCESS;
  }
  
};
