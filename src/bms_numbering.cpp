#include "wmesh.h"
#include "wmesh-status.h"
#include "bms.h"
#include <string.h>
extern "C"
{




  
#define bms_s_init						\
  memcpy(s_size_,&s_size,sizeof(wmesh_int_t));			\
  memcpy(s_ptr_,s_ptr,sizeof(wmesh_int_t)*(s_size+1));		\
  memcpy(s_m_,s_m,sizeof(wmesh_int_t)*s_size);			\
  memcpy(s_n_,s_n,sizeof(wmesh_int_t)*s_size);			\
  memcpy(s_v_,s_v,sizeof(wmesh_int_t)*s_ptr[s_size]);		\
  memcpy(s_ld_,s_ld,sizeof(wmesh_int_t)*s_size)
  
  //
  // s_v_ being 128 is enough.
  //
  wmesh_status_t bms_s_e2n(wmesh_int_t 	dim_,
			   wmesh_int_p	s_size_,
			   wmesh_int_p	s_ptr_,
			   wmesh_int_p	s_m_,
			   wmesh_int_p	s_n_,
			   wmesh_int_p	s_v_,
			   wmesh_int_p	s_ld_)
  {

    WMESH_CHECK_POINTER(s_size_);
    WMESH_CHECK_POINTER(s_ptr_);
    WMESH_CHECK_POINTER(s_m_);
    WMESH_CHECK_POINTER(s_n_);
    WMESH_CHECK_POINTER(s_v_);
    WMESH_CHECK_POINTER(s_ld_);
    if (dim_==3)
      {
	wmesh_int_t s_size = 4;
	wmesh_int_t s_m[4] {2,2,2,2};
	wmesh_int_t s_n[4]{6,8,9,12};
	wmesh_int_t s_ptr[4+1]{0,12,28,46,70};
	wmesh_int_t s_ld[4] {2,2,2,2};
	wmesh_int_t s_v[70]{
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

	bms_s_init;      
	return WMESH_STATUS_SUCCESS;
      }
    else if (dim_==2)
      {
	wmesh_int_t s_size = 2;
	wmesh_int_t s_m[2] 		{2,2};
	wmesh_int_t s_ld[2] 	{2,2};
	wmesh_int_t s_n[2]		{3,4};
	wmesh_int_t s_ptr[2+1]	{0,6,14};
	wmesh_int_t s_v[14]{ 0,1, 1,2, 2,0, 0,1, 1,2, 2,3, 3,0};
	bms_s_init;
	return WMESH_STATUS_SUCCESS;
      }
    else if (dim_==1)
      {
	wmesh_int_t s_size = 1;
	wmesh_int_t s_m[1] 		{2};
	wmesh_int_t s_ld[1] 	{2};
	wmesh_int_t s_n[1]		{1};
	wmesh_int_t s_ptr[1+1]	{0, 2};
	wmesh_int_t s_v[2]{ 0, 1};
	bms_s_init;      
	return WMESH_STATUS_SUCCESS;
      }
    
    WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
  };


  wmesh_status_t bms_s_q2n(wmesh_int_p	s_size_,
			   wmesh_int_p	s_ptr_,
			   wmesh_int_p	s_m_,
			   wmesh_int_p	s_n_,
			   wmesh_int_p	s_v_,
			   wmesh_int_p	s_ld_)
  {

    WMESH_CHECK_POINTER(s_size_);
    WMESH_CHECK_POINTER(s_ptr_);
    WMESH_CHECK_POINTER(s_m_);
    WMESH_CHECK_POINTER(s_n_);
    WMESH_CHECK_POINTER(s_v_);
    WMESH_CHECK_POINTER(s_ld_);

    wmesh_int_t s_size = 4;
    wmesh_int_t s_m[4]  {4,4,4,4};
    wmesh_int_t s_ld[4] {4,4,4,4};
    wmesh_int_t s_n[4]  {0,1,3,6};
    wmesh_int_t s_ptr[4+1]{0,0,4,16,40};
    wmesh_int_t s_v[40]{
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
    
    bms_s_init;
    return WMESH_STATUS_SUCCESS;

  };


  wmesh_status_t bms_s_t2n(wmesh_int_p	s_size_,
			   wmesh_int_p	s_ptr_,
			   wmesh_int_p	s_m_,
			   wmesh_int_p	s_n_,
			   wmesh_int_p	s_v_,
			   wmesh_int_p	s_ld_)
  {
    WMESH_CHECK_POINTER(s_size_);
    WMESH_CHECK_POINTER(s_ptr_);
    WMESH_CHECK_POINTER(s_m_);
    WMESH_CHECK_POINTER(s_n_);
    WMESH_CHECK_POINTER(s_v_);
    WMESH_CHECK_POINTER(s_ld_);
    
    wmesh_int_t s_size = 4;
    wmesh_int_t s_m[4]  {3,3,3,3};
    wmesh_int_t s_ld[4] {3,3,3,3};
    wmesh_int_t s_n[4]  {4,4,2,0};
    wmesh_int_t s_ptr[4+1]{0,12,24,30,30};
    wmesh_int_t s_v[30]{
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
    
    bms_s_init;
    return WMESH_STATUS_SUCCESS;
  }
  
#undef bms_s_init
};
