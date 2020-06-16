#include <limits>
#include <iostream>
#include <array>

#include "bms_c2c_x_template.hpp"

extern "C"
{
  wmesh_status_t  bms_c2c_q_buffer_size(wmesh_int_t		c2n_size_,
					const_wmesh_int_p	c2n_n_,
					wmesh_int_p 		work_n_)
  {
    return bms_c2c_x_buffer_size<WMESH_ELEMENT_QUADRILATERAL>(c2n_size_,
							 c2n_n_,
							 work_n_);
  }
  
  wmesh_status_t  bms_c2c_q		(wmesh_int_t		c2n_size_,
					 const_wmesh_int_p 	c2n_ptr_,
					 const_wmesh_int_p 	c2n_m_,
					 const_wmesh_int_p 	c2n_n_,
					 const_wmesh_int_p 	c2n_v_,
					 const_wmesh_int_p 	c2n_ld_,
					 
					 wmesh_int_t		c2c_q_size_,
					 const_wmesh_int_p 	c2c_q_ptr_,
					 const_wmesh_int_p 	c2c_q_m_,
					 const_wmesh_int_p 	c2c_q_n_,
					 wmesh_int_p 		c2c_q_v_,
					 const_wmesh_int_p 	c2c_q_ld_,
					 
					 wmesh_int_t		s_q2n_size_,
					 const_wmesh_int_p 	s_q2n_ptr_,
					 const_wmesh_int_p	s_q2n_m_,
					 const_wmesh_int_p	s_q2n_n_,
					 const_wmesh_int_p 	s_q2n_v_,
					 const_wmesh_int_p	s_q2n_ld_,
					 
					 wmesh_int_p		num_bfacets_,
					 wmesh_int_t		work_n_,
					 wmesh_int_p		work_)
  {
    return bms_c2c_x<WMESH_ELEMENT_QUADRILATERAL>(c2n_size_,
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
						  true,
						  num_bfacets_,
						  work_n_,
						  work_);
#if 0
    WMESH_CHECK_POINTER(c2n_n_);    
    WMESH_CHECK_POINTER(c2n_m_);
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(c2n_v_);
    WMESH_CHECK_POINTER(c2n_ld_);
  
    WMESH_CHECK_POINTER(c2c_q_m_);
    WMESH_CHECK_POINTER(c2c_q_n_);
    WMESH_CHECK_POINTER(c2c_q_v_);
    WMESH_CHECK_POINTER(c2c_q_ld_);

    WMESH_CHECK_POINTER(s_q2n_m_);
    WMESH_CHECK_POINTER(s_q2n_n_);
    WMESH_CHECK_POINTER(s_q2n_v_);
    WMESH_CHECK_POINTER(s_q2n_ld_);

    WMESH_CHECK_POINTER(num_bfacets_);
    WMESH_CHECK_POINTER(work_);

    bool match_mode = true;
    wmesh_status_t status;
    wmesh_int_t required_work_n;
    status = bms_c2c_x_buffer_size<WMESH_ELEMENT_QUADRILATERAL>(c2n_size_,
								c2n_n_,
								&required_work_n);
    
    if (work_n_ < required_work_n)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
      }
    work_[0] = work_n_ - 1;

    
    //    std::cout << "triangle index before " << triangle_idx_[0] << std::endl;
    static constexpr const wmesh_int_t s_default_hash = - std::numeric_limits<wmesh_int_t>::max();
    const wmesh_int_t hash_size = work_[0];
    wmesh_int_t * hash_link = work_+1;
    for (wmesh_int_t i=0;i<hash_size;++i)
      {
	hash_link[i] = s_default_hash;
      }

    for (wmesh_int_t i=c2n_size_-1;i>=0;--i)
      {
	if (c2n_n_[i]>0)
	  {
	    status = bms_c2c_x_hash<WMESH_ELEMENT_QUADRILATERAL>(i,
								 c2n_m_[i],
								 c2n_n_[i],
								 c2n_v_ + c2n_ptr_[i],
								 c2n_ld_[i],
								 c2c_q_m_[i],
								 c2c_q_n_[i],
								 c2c_q_v_ + c2c_q_ptr_[i],
								 c2c_q_ld_[i],
								 s_q2n_m_[i],
								 s_q2n_n_[i],
								 s_q2n_v_ + s_q2n_ptr_[i],
								 s_q2n_ld_[i],
								 work_n_,
								 work_);
	    WMESH_STATUS_CHECK(status);
	  }
      }
    
    for (wmesh_int_t i=0;i<c2n_size_;++i)
      {
	status = bms_c2c_x_search<WMESH_ELEMENT_QUADRILATERAL>(i,
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

							       match_mode,
							       num_bfacets_,
							       work_n_,
							       work_);
	WMESH_STATUS_CHECK(status);
      }

    return WMESH_STATUS_SUCCESS;
#endif
  }

};
