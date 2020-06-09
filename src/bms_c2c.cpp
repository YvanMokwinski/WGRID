#include <limits>
#include <iostream>

#include "bms.h"

#include "GenericEncoding.hpp"
extern "C"
{

  wmesh_status_t bms_c2c_buffer_size(wmesh_int_t		c2n_size_,
						const_wmesh_int_p 	c2n_n_,
						wmesh_int_p 		work_n_)
  {
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(work_n_);
    wmesh_status_t status;
    wmesh_int_t n, work_n = 0;
    status = bms_c2c_t_buffer_size(c2n_size_,
					      c2n_n_,
					      &n);
    WMESH_STATUS_CHECK(status);    
    work_n = (work_n < n) ? n : work_n;

    status = bms_c2c_q_buffer_size(c2n_size_,
					      c2n_n_,
					      &n);
    WMESH_STATUS_CHECK(status);
    
    work_n = (work_n < n) ? n : work_n;

    work_n_[0] = work_n;
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t bms_c2c(wmesh_int_t		c2n_size_,
			 const_wmesh_int_p 	c2n_ptr_,
			 const_wmesh_int_p 	c2n_m_,
			 const_wmesh_int_p 	c2n_n_,
			 const_wmesh_int_p 	c2n_v_,
			 const_wmesh_int_p 	c2n_ld_,
				     
			 wmesh_int_t		c2c_t_size_,
			 const_wmesh_int_p 	c2c_t_ptr_,
			 const_wmesh_int_p 	c2c_t_m_,
			 const_wmesh_int_p 	c2c_t_n_,
			 wmesh_int_p		c2c_t_v_,
			 const_wmesh_int_p 	c2c_t_ld_,
				     
			 wmesh_int_t		c2c_q_size_,
			 const_wmesh_int_p 	c2c_q_ptr_,
			 const_wmesh_int_p 	c2c_q_m_,
			 const_wmesh_int_p 	c2c_q_n_,
			 wmesh_int_p		c2c_q_v_,
			 const_wmesh_int_p 	c2c_q_ld_,
				     
			 wmesh_int_t 		s_t2n_size_,
			 const_wmesh_int_p 	s_t2n_ptr_,
			 const_wmesh_int_p 	s_t2n_m_,
			 const_wmesh_int_p 	s_t2n_n_,
			 const_wmesh_int_p 	s_t2n_v_,
			 const_wmesh_int_p 	s_t2n_ld_,
			 
			 wmesh_int_t 		s_q2n_size_,
			 const_wmesh_int_p 	s_q2n_ptr_,
			 const_wmesh_int_p 	s_q2n_m_,
			 const_wmesh_int_p 	s_q2n_n_,
			 const_wmesh_int_p 	s_q2n_v_,
			 const_wmesh_int_p 	s_q2n_ld_,
			 
			 wmesh_int_t 		work_n_,
			 wmesh_int_p 		work_,
			 wmesh_int_p 		num_faces_,
			 wmesh_int_p 		num_bfaces_)
  {

    wmesh_status_t 	status;
    wmesh_int_t 	required_work_n;

    status =  bms_c2c_buffer_size(c2n_size_,
				  c2n_n_,
				  &required_work_n);
    WMESH_STATUS_CHECK(status);
    if (work_n_ < required_work_n)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
      }

    //
    // Enumerates triangles faces.
    //
    status = bms_c2c_t(c2n_size_,
		       c2n_ptr_,
		       c2n_m_,
		       c2n_n_,
		       c2n_v_,
		       c2n_ld_,

		       c2c_t_size_,
		       c2c_t_ptr_,
		       c2c_t_m_,
		       c2c_t_n_,
		       c2c_t_v_,
		       c2c_t_ld_,
		       
		       s_t2n_size_,
		       s_t2n_ptr_,
		       s_t2n_m_,
		       s_t2n_n_,
		       s_t2n_v_,
		       s_t2n_ld_,
		       
		       &num_bfaces_[0],
		       work_n_,
		       work_);    
    WMESH_STATUS_CHECK(status);
    num_faces_[0] = ((c2n_n_[0]*4+c2n_n_[1]*4+c2n_n_[2]*2) + num_bfaces_[0])/2;
    
    //
    // Enumerates quadrilaterals faces.
    //    
    status =  bms_c2c_q(c2n_size_,
			c2n_ptr_,
			c2n_m_,
			c2n_n_,
			c2n_v_,
			c2n_ld_,
				    
			c2c_q_size_,
			c2c_q_ptr_,
			c2c_q_m_,
			c2c_q_n_,
			c2c_q_v_,
			c2c_q_ld_,
				    
			s_q2n_size_,
			s_q2n_ptr_,
			s_q2n_m_,
			s_q2n_n_,
			s_q2n_v_,
			s_q2n_ld_,
				    
			&num_bfaces_[1],
			work_n_,
			work_);    
    WMESH_STATUS_CHECK(status);
    num_faces_[1] = ((c2n_n_[1]*1+c2n_n_[2]*3+c2n_n_[3]*6) + num_bfaces_[1])/2;    
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t bms_c2c_cindex(wmesh_int_t 	c_,
				wmesh_int_p 	cindex_)
  {
    cindex_[0] = GenericEncoding<wmesh_int_t,2>::Up(c_) - 1;
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t bms_c2c_ctype(wmesh_int_t 	c_,
			       wmesh_int_p 	ctype_)
  {
    ctype_[0] = GenericEncoding<wmesh_int_t,2>::Low(c_);
    return WMESH_STATUS_SUCCESS;
  }

}
