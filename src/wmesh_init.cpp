
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "bms.h"

#include <chrono>
#include <iostream>

#define WINT_SPARSE_MAT_PARAMS(_f)				\
    _f.m_ptr,							\
      _f.m_m,							\
      _f.m_n,							\
      _f.m_data,						\
      _f.m_ld

using namespace std::chrono;
extern "C"
{
  
  wmesh_status_t wmesh_init_bf2n(wmesh_int_t 			num_triangles_,
				 wmesh_int_t 			num_quadrilaterals_,
				 wmesh_int_sparsemat_t*		bf2n_)
  {
    wmesh_int_t bf2n_m[2] {3,4};
    wmesh_int_t bf2n_ld[2] {3,4};
    wmesh_int_t bf2n_ptr[2+1];
    wmesh_int_t bf2n_n[2];

    bf2n_n[0] = num_triangles_;
    bf2n_n[1] = num_quadrilaterals_;
    
    bf2n_ptr[0] = 0;
    bf2n_ptr[1] = bf2n_ptr[0] + bf2n_n[0] * bf2n_ld[0];
    bf2n_ptr[2] = bf2n_ptr[1] + bf2n_n[1] * bf2n_ld[1];
    wmesh_int_t * __restrict__ bf2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*bf2n_ptr[2]);
    if (!bf2n_v)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    return wmesh_int_sparsemat_new(bf2n_,
				   2,				   
				   bf2n_ptr,
				   bf2n_m,
				   bf2n_n,
				   bf2n_v,
				   bf2n_ld);
  }

  
   wmesh_status_t wmesh_init_c2f(const wmesh_int_sparsemat_t*	c2n_,
				       wmesh_int_sparsemat_t*		c2f_)
  {
    wmesh_int_t c2f_m[4] {4,5,5,6};
    wmesh_int_t c2f_ld[4] {4,5,5,6};
    wmesh_int_t c2f_ptr[4+1];
    wmesh_int_t c2f_n[4];      
    for (wmesh_int_t i=0;i<4;++i) c2f_n[i] = c2n_->m_n[i];
    c2f_ptr[0] = 0;
    c2f_ptr[1] = c2f_ptr[0] + c2f_n[0] * c2f_ld[0];
    c2f_ptr[2] = c2f_ptr[1] + c2f_n[1] * c2f_ld[1];
    c2f_ptr[3] = c2f_ptr[2] + c2f_n[2] * c2f_ld[2];
    c2f_ptr[4] = c2f_ptr[3] + c2f_n[3] * c2f_ld[3];
    wmesh_int_t * __restrict__ c2f_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*c2f_ptr[4]);
    if (!c2f_v)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }    
    return wmesh_int_sparsemat_new(c2f_,
				   4,
				   c2f_ptr,
				   c2f_m,
				   c2f_n,
				   c2f_v,
				   c2f_ld);
  }
  


  wmesh_status_t wmesh_init_c2f_q(wmesh_int_sparsemat_t*	c2f_,
					 wmesh_int_sparsemat_t*	c2f_q_)
  {
    wmesh_int_t c2f_q_m[4]  {0,1,3,6};
    wmesh_int_t c2f_q_ptr[4+1];        
    c2f_q_ptr[0] = c2f_->m_ptr[0] + 4; 
    c2f_q_ptr[1] = c2f_->m_ptr[1] + 4; 
    c2f_q_ptr[2] = c2f_->m_ptr[2] + 2; 
    c2f_q_ptr[3] = c2f_->m_ptr[3] + 0;
    c2f_q_ptr[4] = c2f_->m_ptr[4];
    return wmesh_int_sparsemat_new(c2f_q_,
				   4,
				   c2f_q_ptr,
				   c2f_q_m,
				   c2f_->m_n,
				   c2f_->m_data,
				   c2f_->m_ld);
  }

   wmesh_status_t wmesh_init_c2f_t(wmesh_int_sparsemat_t*	c2f_,
					 wmesh_int_sparsemat_t*	c2f_t_)
  {
    wmesh_int_t c2f_t_m[4]  {4,4,2,0};
    return wmesh_int_sparsemat_new(c2f_t_,
				   4,
				   c2f_->m_ptr,
				   c2f_t_m,
				   c2f_->m_n,
				   c2f_->m_data,
				   c2f_->m_ld);
  }

  

  wmesh_status_t wmesh_init_c2e(const wmesh_int_sparsemat_t*	c2n_,
				wmesh_int_sparsemat_t*		c2e_,
				wmesh_int_t 			topology_dimension_)
  {
    if (topology_dimension_==3)
      {
	wmesh_int_t c2e_m[4] {6,8,9,12};
	wmesh_int_t c2e_ptr[4+1];
	wmesh_int_t c2e_n[4];      
	for (wmesh_int_t i=0;i<4;++i) c2e_n[i] = c2n_->m_n[i];
	c2e_ptr[0] = 0;
	c2e_ptr[1] = c2e_ptr[0] + c2e_n[0] * c2e_m[0];
	c2e_ptr[2] = c2e_ptr[1] + c2e_n[1] * c2e_m[1];
	c2e_ptr[3] = c2e_ptr[2] + c2e_n[2] * c2e_m[2];
	c2e_ptr[4] = c2e_ptr[3] + c2e_n[3] * c2e_m[3];
	wmesh_int_t * __restrict__ c2e_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*c2e_ptr[4]);
	if (!c2e_v)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }
	return wmesh_int_sparsemat_new(c2e_,
				       4,
				       c2e_ptr,
				       c2e_m,
				       c2e_n,
				       c2e_v,
				       c2e_m);
      }
    else if (topology_dimension_ == 2)
      {
	wmesh_int_t c2e_m[2] {3,4};
	wmesh_int_t c2e_ptr[2+1];
	wmesh_int_t c2e_n[2];      

	for (wmesh_int_t i=0;i<2;++i) c2e_n[i] = c2n_->m_n[i];

	c2e_ptr[0] = 0;
	c2e_ptr[1] = c2e_ptr[0] + c2e_n[0] * c2e_m[0];
	c2e_ptr[2] = c2e_ptr[1] + c2e_n[1] * c2e_m[1];

	wmesh_int_t * __restrict__ c2e_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*c2e_ptr[2]);
	if (!c2e_v)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }
	return wmesh_int_sparsemat_new(c2e_,
				       2,
				       c2e_ptr,
				       c2e_m,
				       c2e_n,
				       c2e_v,
				       c2e_m);
      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
      }
  };


};
