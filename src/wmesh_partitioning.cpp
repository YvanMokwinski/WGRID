#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"

#include <chrono>
#include <iostream>
#include "wmesh.hpp"
#include "GenericEncoding.hpp"
#include "wmesh_utils.hpp"

using namespace std::chrono;
static int hilbert_sort_predicate(const void * a_,
			   const void * b_)
{
  const unsigned long long int * a = (const unsigned long long int*)a_;
  const unsigned long long int * b = (const unsigned long long int*)b_;
  if (*a < *b)
    {
      return -1;
    }
  else if (*a > *b)    
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

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


#if 0
wmesh_status_t wmesh_partitioning_sfc(wmesh_int_t 			nparts_,
				      wmesh_int_t 			num_cells_,
				      wmesh_int_t * 			part_id_,
				      
				      wmesh_int_t 			c2n_size_,
				      const_wmesh_int_p 		c2n_ptr_,
				      const_wmesh_int_p 		c2n_m_,
				      const_wmesh_int_p 		c2n_n_,
				      const_wmesh_int_p 		c2n_v_,
				      const_wmesh_int_p 		c2n_ld_,
				      
				      wmesh_int_t 			coo_m_,
				      wmesh_int_t 			coo_n_,
				      const double *__restrict__ 	coo_v_,
				      wmesh_int_t 			coo_ld_,
				      
				      size_t 				work_size_,
				      void *__restrict__ 		work_)
{
  
  double box[6];
  double center[3];
  
  
  wmesh_status_t status;
  if (work_size_ < 2 * num_cells_ * sizeof(unsigned long long int))
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }
    
  status = wmesh_raw_bounding_box(coo_m_,
				  coo_n_,
				  coo_v_,
				  coo_ld_,
				  box,
				  work_size_ / sizeof(double),
				  (double*)work_);
  
  unsigned long long int*__restrict__	sfc_v	= (unsigned long long int*__restrict__)work_;  
  for (wmesh_int_t l=0;l<c2n_size_;++l)
    {
      status = wmesh_sfc_fill(c2n_m_[l],
			      c2n_n_[l],
			      c2n_v_ + c2n_ptr_[l],
			      c2n_ld_[l],
			      
			      coo_m_,
			      coo_n_,
			      coo_v_,
			      coo_ld_,
			      
			      sfc_v + 2 * (c2n_ptr_[l] / c2n_ld_[l]),
			      2,
			      
			      box,
			      work_size_,
			      work_);
      
      WMESH_STATUS_CHECK(status);
    }
  
  //
  // Fill ids.
  //
  {
    unsigned long long int idx = 0;
    for (wmesh_int_t l=0;l<c2n_size_;++l)
      for (int j=0;j<c2n_n_[l];++j)
	{
	  sfc_v[2 * idx + 1] = GenericEncoding<unsigned long long int,2>::Encod(idx+1, (unsigned long long int )l);
	  ++idx;
	}
  }

  //
  // Sort.
  //
  for (int i=0;i<self_->m_num_cells;++i)
    {
      std::cout << sfc_v[2*i+0] <<  " " << sfc_v[2*i+1] << std::endl;
    }
  
  qsort(sfc_v,
	num_cells_,
	sizeof(unsigned long long int)*2,
	hilbert_sort_predicate);

  //
  // 
  //
  {
    unsigned long long int idx = 0;
    unsigned long long int shift = 0;
    unsigned long long int N = num_cells_ / nparts_ + num_cells_ % nparts_ ;

    for (unsigned long long int i=0;i<N;++i)
      {
	sfc_v[2*i+0] = 0;
	sfc_v[2*i+1] = GenericEncoding<unsigned long long int,2>::Up(sfc_v[2*i+1])-1;	    
      }    
    for (wmesh_int_t ipart=1;ipart<nparts_;++ipart)
      {
	for (unsigned long long int i=N;i<N + num_cells_ / nparts_;++i)
	  {
	    sfc_v[2*i+0] = ipart;	    
	    sfc_v[2*i+1] = GenericEncoding<unsigned long long int,2>::Up(sfc_v[2*i+1])-1;	    
	  }
	N += num_cells_ / nparts_;
      }
    

    for (wmesh_int_t i=0;i<num_cells_;++i)
      {
	part_id_[sfc_v[2*i+1]] = sfc_v[2*i+0];
      }

  }
  return WMESH_STATUS_SUCCESS;
}
#endif


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
    
    WMESH_POINTER_CHECK(p_v_);
    WMESH_POINTER_CHECK(c2n_ptr_);
    WMESH_POINTER_CHECK(c2n_m_);
    WMESH_POINTER_CHECK(c2n_n_);
    WMESH_POINTER_CHECK(c2n_v_);
    WMESH_POINTER_CHECK(c2n_ld_);
    WMESH_POINTER_CHECK(coo_v_);
    WMESH_POINTER_CHECK(work_size_);
    
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
	  hilbert_sort_predicate);
    
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

  
  wmesh_status_t wmesh_partitioning(wmesh_t*	self_,
				    wmesh_int_t nparts_)
  {

    size_t  work_size = 0;
    void *  work = nullptr;
    wmesh_status_t status;
    status = bms_partitioning_sfc(nparts_,
				  self_->m_num_cells,
				  self_->m_c_c.m_data,
				  1,
				  self_->m_c2n.m_size,
				  self_->m_c2n.m_ptr,
				  self_->m_c2n.m_m,
				  self_->m_c2n.m_n,
				  self_->m_c2n.m_data,
				  self_->m_c2n.m_ld,
				  
				  3,
				  self_->m_num_nodes,
				  self_->m_coo,
				  3,
				  
				  &work_size,
				  work);
    WMESH_STATUS_CHECK(status);

    work = (void*)malloc(work_size);
    status = bms_partitioning_sfc(nparts_,
				  self_->m_num_cells,
				  self_->m_c_c.m_data,
				  1,
				  self_->m_c2n.m_size,
				  self_->m_c2n.m_ptr,
				  self_->m_c2n.m_m,
				  self_->m_c2n.m_n,
				  self_->m_c2n.m_data,
				  self_->m_c2n.m_ld,
				  
				  3,
				  self_->m_num_nodes,
				  self_->m_coo,
				  3,
				  
				  &work_size,
				  work);
    free(work);
    work = nullptr;
    if (WMESH_STATUS_SUCCESS != status)
      {
	WMESH_STATUS_CHECK(status);
      }

#if 0
    double box[6];
    std::cout << "num_cells " << self_->m_num_cells << std::endl;

    
    status = wmesh_raw_bounding_box(3,
				    self_->m_num_nodes,
				    self_->m_coo,
				    3,
				    box);    
    WMESH_STATUS_CHECK(status);
    
    status = wmesh_sfc_fill(self_->m_c2n.m_size,
			    self_->m_c2n.m_ptr,
			    self_->m_c2n.m_m,
			    self_->m_c2n.m_n,
			    self_->m_c2n.m_data,
			    self_->m_c2n.m_ld,
			    
			    3,
			    self_->m_num_nodes,
			    self_->m_coo,
			    3,
			    
			    p,
			    2,			      
			    box);

    WMESH_STATUS_CHECK(status);
    
    {
      wmesh_int_t idx = 0;
      for (wmesh_int_t i=0;i<nparts_;++i)
	{
	  for (wmesh_int_t j=0;j<self_->m_num_cells / nparts_;++j)
	    {
	      p[2*idx+0] = i;
	      ++idx;
	    }
	}
      for (wmesh_int_t j=0;j<self_->m_num_cells % nparts_;++j)
	{
	  p[2*idx+0] = nparts_-1;
	  ++idx;
	}
    }
    
    for (wmesh_int_t i=0;i<self_->m_num_cells;++i)
      {
	self_->m_c_c.m_data[p[2*i+1]] = p[2*i+0];
      }
#endif    
    return WMESH_STATUS_SUCCESS;
  }
  
};
