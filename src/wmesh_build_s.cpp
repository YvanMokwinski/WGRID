
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"


extern "C"
{

  wmesh_status_t wmesh_build_s_e2n(wmesh_int_sparsemat_t* 	self_,
				   wmesh_int_t 			dim_)
  {
#ifndef NDEBUG
    WMESH_CHECK_POINTER(self_);
#endif    
    if (dim_==3)
      {
	wmesh_int_t s_e2n_m[4] {2,2,2,2};
	wmesh_int_t s_e2n_ld[4] {2,2,2,2};
	wmesh_int_t s_e2n_n[4]{6,8,9,12};
	wmesh_int_t s_e2n_ptr[4+1]{0,12,28,46,70};
	wmesh_int_t v[70]{
	  // tet
	  1,2,	  
	    2,0,
	    0,1,
	    2,3,
	    0,3,
	    1,3,
	    // pyr
	    0,1,
	    1,2,
	    2,3,
	    3,0,
	    0,4,
	    1,4,
	    2,4,
	    3,4,
	    //wedge
	    0,1,
	    1,2,
	    2,0,
	    3,4,
	    4,5,
	    5,3,
	    0,3,
	    1,4,
	    2,5,
	    // hex
	    0,1,
	    1,2,
	    2,3,
	    3,0,
	    4,5,
	    5,6,
	    6,7,
	    7,4,
	    0,4,
	    1,5,
	    2,6,
	    3,7
	    };
      
	wmesh_int_t * __restrict__ s_e2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*s_e2n_ptr[4]);
	memcpy(s_e2n_v,v,sizeof(wmesh_int_t)*s_e2n_ptr[4]);
	return wmesh_int_sparsemat_new(self_,
				       4,
				       s_e2n_ptr,
				       s_e2n_m,
				       s_e2n_n,
				       s_e2n_v,
				       s_e2n_ld);	
      }
    else if (dim_==2)
      {
	wmesh_int_t s_e2n_m[2] 		{2,2};
	wmesh_int_t s_e2n_ld[2] 	{2,2};
	wmesh_int_t s_e2n_n[2]		{3,4};
	wmesh_int_t s_e2n_ptr[2+1]	{0,6,14};
	wmesh_int_t v[14]{ 0,1, 1,2, 2,0, 0,1, 1,2, 2,3, 3,0};
	wmesh_int_t * __restrict__ s_e2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*s_e2n_ptr[2]);
	memcpy(s_e2n_v,v,sizeof(wmesh_int_t)*s_e2n_ptr[2]);

	return wmesh_int_sparsemat_new(self_,
				       2,
				       s_e2n_ptr,
				       s_e2n_m,
				       s_e2n_n,
				       s_e2n_v,
				       s_e2n_ld);
	
      }
    else if (dim_==1)
      {
	wmesh_int_t s_e2n_m[1] 		{2};
	wmesh_int_t s_e2n_ld[1] 	{2};
	wmesh_int_t s_e2n_n[1]		{1};
	wmesh_int_t s_e2n_ptr[1+1]	{0, 2};
	wmesh_int_t v[2]{ 0, 1};
	wmesh_int_t * __restrict__ s_e2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*s_e2n_ptr[1]);
	memcpy(s_e2n_v,v,sizeof(wmesh_int_t)*s_e2n_ptr[1]);

	return wmesh_int_sparsemat_new(self_,
				       1,
				       s_e2n_ptr,
				       s_e2n_m,
				       s_e2n_n,
				       s_e2n_v,
				       s_e2n_ld);
	
	return WMESH_STATUS_SUCCESS;
      }

    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
  };

  wmesh_status_t wmesh_build_s_q2n(wmesh_int_sparsemat_t* 	self_,
				   wmesh_int_t 			dim_)
  {
#ifndef NDEBUG
    WMESH_CHECK_POINTER(self_);
#endif    
    if (3 == dim_)
      {
	wmesh_int_t s_q2n_m[4]  {4,4,4,4};
	wmesh_int_t s_q2n_ld[4] {4,4,4,4};
	wmesh_int_t s_q2n_n[4]  {0,1,3,6};
	wmesh_int_t s_q2n_ptr[4+1]{0,0,4,16,40};
	wmesh_int_t v[40]{
	  // tet
	  // pyr
	  0,3,2, 1,
	    //wedge
	    0,2,5,3,
	    4,5,2,1,
	    0,3,4,1,
	    // hex
	    0,3,2,1,
	    4,5,6,7,	  
	    0,1,5,4,
	    1,2,6,5,
	    2,3,7,6,
	    3,0,4,7
	    };
	
	wmesh_int_t * __restrict__ s_q2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*s_q2n_ptr[4]);
	memcpy(s_q2n_v,v,sizeof(wmesh_int_t)*s_q2n_ptr[4]);
	return wmesh_int_sparsemat_new(self_,
				       4,
				       s_q2n_ptr,
				       s_q2n_m,
				       s_q2n_n,
				       s_q2n_v,				       
				       s_q2n_ld);
      }
    else if (2==dim_)
      {
	wmesh_int_t s_e2n_m[1] 		{4};
	wmesh_int_t s_e2n_ld[1] 	{4};
	wmesh_int_t s_e2n_n[1]		{1};
	wmesh_int_t s_e2n_ptr[1+1]	{0, 4};
	wmesh_int_t v[4]{ 0, 1, 2, 3};
	wmesh_int_t * __restrict__ s_e2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*s_e2n_ptr[1]);
	memcpy(v,s_e2n_v,sizeof(wmesh_int_t)*s_e2n_ptr[1]);

	return wmesh_int_sparsemat_new(self_,
				       1,
				       s_e2n_ptr,
				       s_e2n_m,
				       s_e2n_n,
				       s_e2n_v,
				       s_e2n_ld);
	
	return WMESH_STATUS_SUCCESS;
      }
    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
  }
  

  wmesh_status_t wmesh_build_s_t2n(wmesh_int_sparsemat_t* 	self_,
				   wmesh_int_t dim_)
  {
#ifndef NDEBUG
    WMESH_CHECK_POINTER(self_);
#endif    
    if (3==dim_)
      {
	wmesh_int_t s_t2n_m[4]  {3,3,3,3};
	wmesh_int_t s_t2n_ld[4] {3,3,3,3};
	wmesh_int_t s_t2n_n[4]  {4,4,2,0};
	wmesh_int_t s_t2n_ptr[4+1]{0,12,24,30,30};
	wmesh_int_t v[30]{
	  // tet
	  1,2,3,
	    2,0,3,
	    0,1,3,
	    0,2,1,
	    // pyr
	    0,1,4, 
	    1,2,4, 
	    2,3,4, 
	    3,0,4, 
	    //wedge
	    0,1,2,
	    3,5,4							     
	    // hex
	    };
    
	wmesh_int_p s_t2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*s_t2n_ptr[4]);
	memcpy(s_t2n_v,v,sizeof(wmesh_int_t)*s_t2n_ptr[4]);
	return wmesh_int_sparsemat_new(self_,
				       4,
				       s_t2n_ptr,
				       s_t2n_m,
				       s_t2n_n,
				       s_t2n_v,
				       s_t2n_ld);
      }
    else if (2 == dim_)
      {
	wmesh_int_t s_e2n_m[1] 		{4};
	wmesh_int_t s_e2n_ld[1] 	{4};
	wmesh_int_t s_e2n_n[1]		{1};
	wmesh_int_t s_e2n_ptr[1+1]	{0, 4};
	wmesh_int_t v[4]{ 0, 1, 2, 3};
	wmesh_int_p s_e2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*s_e2n_ptr[1]);
	memcpy(v,s_e2n_v,sizeof(wmesh_int_t)*s_e2n_ptr[1]);

	return wmesh_int_sparsemat_new(self_,
				       1,
				       s_e2n_ptr,
				       s_e2n_m,
				       s_e2n_n,
				       s_e2n_v,
				       s_e2n_ld);
	
	return WMESH_STATUS_SUCCESS;
      }
    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);

  };
  
};
