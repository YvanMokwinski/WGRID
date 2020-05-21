#include <stdlib.h>
#include <string.h>

#include "wmesh-status.h"
#include "wmesh.hpp"
#include "bms.h"

extern "C"
{

  wmesh_status_t wmesh_build_s_e2n(wmesh_int_sparsemat_t* 	self_,
				   wmesh_int_t 			dim_)
  {

    wmesh_int_t s_e2n_size;
    wmesh_int_t s_e2n_ptr[4+1];
    wmesh_int_t s_e2n_m[4];
    wmesh_int_t s_e2n_n[4];
    wmesh_int_t v[128];
    wmesh_int_t s_e2n_ld[4];

    wmesh_status_t status = bms_s_e2n(dim_,
				      &s_e2n_size,
				      s_e2n_ptr,
				      s_e2n_m,
				      s_e2n_n,
				      v,
				      s_e2n_ld);
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_t * __restrict__ s_e2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*s_e2n_ptr[s_e2n_size]);
    memcpy(s_e2n_v,v,sizeof(wmesh_int_t)*s_e2n_ptr[s_e2n_size]);
    return wmesh_int_sparsemat_new(self_,
				   s_e2n_size,
				   s_e2n_ptr,
				   s_e2n_m,
				   s_e2n_n,
				   s_e2n_v,
				   s_e2n_ld);	
  }

  wmesh_status_t wmesh_build_s_q2n(wmesh_int_sparsemat_t* 	self_)
  {


    wmesh_int_t s_q2n_size;
    wmesh_int_t s_q2n_ptr[4+1];
    wmesh_int_t s_q2n_m[4];
    wmesh_int_t s_q2n_n[4];
    wmesh_int_t v[128];
    wmesh_int_t s_q2n_ld[4];

    wmesh_status_t status = bms_s_q2n(&s_q2n_size,
				      s_q2n_ptr,
				      s_q2n_m,
				      s_q2n_n,
				      v,
				      s_q2n_ld);
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_t * __restrict__ s_q2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*s_q2n_ptr[s_q2n_size]);
    memcpy(s_q2n_v,v,sizeof(wmesh_int_t)*s_q2n_ptr[s_q2n_size]);
    return wmesh_int_sparsemat_new(self_,
				   s_q2n_size,
				   s_q2n_ptr,
				   s_q2n_m,
				   s_q2n_n,
				   s_q2n_v,
				   s_q2n_ld);	

    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
  }
  

  wmesh_status_t wmesh_build_s_t2n(wmesh_int_sparsemat_t* self_)
  {

    wmesh_int_t s_t2n_size;
    wmesh_int_t s_t2n_ptr[4+1];
    wmesh_int_t s_t2n_m[4];
    wmesh_int_t s_t2n_n[4];
    wmesh_int_t v[128];
    wmesh_int_t s_t2n_ld[4];
    wmesh_status_t status = bms_s_t2n(&s_t2n_size,
				      s_t2n_ptr,
				      s_t2n_m,
				      s_t2n_n,
				      v,
				      s_t2n_ld);
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_t * __restrict__ s_t2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*s_t2n_ptr[s_t2n_size]);
    memcpy(s_t2n_v,v,sizeof(wmesh_int_t)*s_t2n_ptr[s_t2n_size]);
    return wmesh_int_sparsemat_new(self_,
				   s_t2n_size,
				   s_t2n_ptr,
				   s_t2n_m,
				   s_t2n_n,
				   s_t2n_v,
				   s_t2n_ld);	

  }

}
