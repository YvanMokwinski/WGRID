
#include <limits>
#include <iostream>
#include <array>
#include "wmesh.hpp"

#include "GenericEncoding.hpp"
#include "bms_c2x_template.hpp"

extern "C"
{
  wmesh_status_t  bms_c2f_q_buffer_size(wmesh_int_t		c2n_size_,
					const_wmesh_int_p	c2n_n_,
					wmesh_int_p 		work_n_)
  {
    return bms_c2x_buffer_size<WMESH_ELEMENT_QUADRILATERAL>(c2n_size_,
							    c2n_n_,
							    work_n_);
  };
  
  wmesh_status_t  bms_c2f_q(wmesh_int_t 		c2n_size_,
						
			    const_wmesh_int_p 	c2n_ptr_,
			    const_wmesh_int_p 	c2n_m_,
			    const_wmesh_int_p 	c2n_n_,
			    const_wmesh_int_p 	c2n_v_,
			    const_wmesh_int_p 	c2n_ld_,
						
			    wmesh_int_t 		c2q_size_,
			    const_wmesh_int_p 	c2q_ptr_,
			    const_wmesh_int_p 	c2q_m_,
			    const_wmesh_int_p 	c2q_n_,
			    wmesh_int_p 		c2q_v_,
			    const_wmesh_int_p 	c2q_ld_,
						
			    wmesh_int_t 		s_q2n_size_,
			    const_wmesh_int_p 	s_q2n_ptr_,
			    const_wmesh_int_p	s_q2n_m_,
			    const_wmesh_int_p	s_q2n_n_,
			    const_wmesh_int_p 	s_q2n_v_,
			    const_wmesh_int_p	s_q2n_ld_,
				       
			    wmesh_int_p		idx_,
			    wmesh_int_t		work_n_,
			    wmesh_int_p		work_)
  {
    
    bool match_mode = true;
    return bms_c2x<WMESH_ELEMENT_QUADRILATERAL>(c2n_size_,
						
						c2n_ptr_,
						c2n_m_,
						c2n_n_,
						c2n_v_,
						c2n_ld_,
						
						c2q_size_,
						c2q_ptr_,
						c2q_m_,
						c2q_n_,
						c2q_v_,
						c2q_ld_,
						
						s_q2n_size_,
						s_q2n_ptr_,
						s_q2n_m_,
						s_q2n_n_,
						s_q2n_v_,
						s_q2n_ld_,
						
						match_mode,
						idx_,
						work_n_,
						work_);
  }
};

