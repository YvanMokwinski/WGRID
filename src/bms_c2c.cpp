
#include "bms.h"

extern "C"
{

  wmesh_status_t wbms_c2c_calculate_buffer_size(wmesh_int_t		c2n_size_,
						const_wmesh_int_p 	c2n_n_,
						wmesh_int_p 		work_n_)
  {
    WMESH_POINTER_CHECK(c2n_n_);
    WMESH_POINTER_CHECK(work_n_);
    wmesh_status_t status;
    wmesh_int_t n, work_n = 0;
    status = wbms_c2c_t_calculate_buffer_size(c2n_size_,
					      c2n_n_,
					      &n);
    WMESH_STATUS_CHECK(status);    
    work_n = (work_n < n) ? n : work_n;

    status = wbms_c2c_q_calculate_buffer_size(c2n_size_,
					      c2n_n_,
					      &n);
    WMESH_STATUS_CHECK(status);
    
    work_n = (work_n < n) ? n : work_n;

    work_n_[0] = work_n;
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t wbms_c2c_calculate(wmesh_int_t		c2n_size_,
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
				     
				    wmesh_int_t 	s_t2n_size_,
				    const_wmesh_int_p 	s_t2n_ptr_,
				    const_wmesh_int_p 	s_t2n_m_,
				    const_wmesh_int_p 	s_t2n_n_,
				    const_wmesh_int_p 	s_t2n_v_,
				    const_wmesh_int_p 	s_t2n_ld_,
			 
				    wmesh_int_t 	s_q2n_size_,
				    const_wmesh_int_p 	s_q2n_ptr_,
				    const_wmesh_int_p 	s_q2n_m_,
				    const_wmesh_int_p 	s_q2n_n_,
				    const_wmesh_int_p 	s_q2n_v_,
				    const_wmesh_int_p 	s_q2n_ld_,
			 
				    wmesh_int_t 	work_n_,
				    wmesh_int_p 	work_,
				    wmesh_int_p 	num_faces_,
				    wmesh_int_p 	num_bfaces_)
  {

    wmesh_status_t 	status;
    wmesh_int_t 	required_work_n;

    status =  wbms_c2c_calculate_buffer_size(c2n_size_,
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
    status = wbms_c2c_t_calculate(c2n_size_,
				  c2n_ptr_,
				  c2n_m_,
				  c2n_n_,
				  c2n_v_,
				  c2n_ld_,
				  
				  c2c_t_ptr_,
				  c2c_t_m_,
				  c2c_t_n_,
				  c2c_t_v_,
				  c2c_t_ld_,
				  
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
    status =  wbms_c2c_q_calculate(c2n_size_,
				   c2n_ptr_,
				   c2n_m_,
				   c2n_n_,
				   c2n_v_,
				   c2n_ld_,
				    
				   c2c_q_ptr_,
				   c2c_q_m_,
				   c2c_q_n_,
				   c2c_q_v_,
				   c2c_q_ld_,
				    
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
  
}
#if 0
  wmesh_status_t bms_c2c2d_buffer_size(wmesh_int_t		c2n_size_,
				       const_wmesh_int_p 	c2n_n_,
				       wmesh_int_p 		work_n_)
  {
    WMESH_POINTER_CHECK(c2n_n_);
    WMESH_POINTER_CHECK(work_n_);
    wmesh_status_t status;
    wmesh_int_t n, work_n = 0;
    
    status = wmesh_c2c_e_calculate_buffer_size(c2n_size_,
					       c2n_n_,
					       &n);
    WMESH_STATUS_CHECK(status);
    
    work_n = (work_n < n) ? n : work_n;
    
    work_n_[0] = work_n;
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t bms_c2c2d(wmesh_int_t		c2n_size_,
			   const_wmesh_int_p 	c2n_ptr_,
			   const_wmesh_int_p 	c2n_m_,
			   const_wmesh_int_p 	c2n_n_,
			   const_wmesh_int_p 	c2n_v_,
			   const_wmesh_int_p 	c2n_ld_,
			   
			   wmesh_int_t		c2c_size_,
			   const_wmesh_int_p 	c2c_ptr_,
			   const_wmesh_int_p 	c2c_m_,
			   const_wmesh_int_p 	c2c_n_,
			   wmesh_int_p		c2c_v_,
			   const_wmesh_int_p 	c2c_ld_,
			   
			 
			   wmesh_int_t 		s_e2n_size_,
			   const_wmesh_int_p 	s_e2n_ptr_,
			   const_wmesh_int_p 	s_e2n_m_,
			   const_wmesh_int_p 	s_e2n_n_,
			   const_wmesh_int_p 	s_e2n_v_,
			   const_wmesh_int_p 	s_e2n_ld_,
			 
			   wmesh_int_t 		work_n_,
			   wmesh_int_p 		work_,
			   
			   wmesh_int_t 		num_edges_[2],
			   wmesh_int_t 		num_bedges_[2])
  {
    wmesh_status_t 	status;
    wmesh_int_t 	required_work_n;

    status =  bms_c2c2d_buffer_size(c2n_size_,
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
    status = wmesh_c2c_e_calculate(c2n_size_,
				   c2n_ptr_,
				   c2n_m_,
				   c2n_n_,
				   c2n_v_,
				   c2n_ld_,
				   
				   c2c_e_ptr,
				   c2c_e_m,
				   c2c_e_n,
				   c2c_e_v,
				   c2c_e_ld,
				   
				   s_e2n_ptr_,
				   s_e2n_m_,
				   s_e2n_n_,
				   s_e2n_v_,
				   s_e2n_ld_,
				   
				   &num_bedges_[0],
				   work_n_,
				   work_);    
    WMESH_STATUS_CHECK(status);
    num_edges_[0] = ((c2n_n_[0]*3+c2n_n_[1]*4) + num_bfaces_[0])/2;
    return WMESH_STATUS_SUCCESS;
  };    
#endif
