#include <limits>
#include <iostream>
#include <array>
#include "wmesh.hpp"
#include "GenericEncoding.hpp"
#include "bms_c2c_x_template.hpp"

extern "C"
{
  wmesh_status_t  bms_c2c_e_buffer_size(wmesh_int_t		c2n_size_,
					const_wmesh_int_p	c2n_n_,
					wmesh_int_p 		work_n_)
  {
    return bms_c2c_x_buffer_size<WMESH_ELEMENT_EDGE>(c2n_size_,
						     c2n_n_,
						     work_n_);
  }
  
  wmesh_status_t  bms_c2c_e(wmesh_int_t		c2n_size_,
			    const_wmesh_int_p 	c2n_ptr_,
			    const_wmesh_int_p 	c2n_m_,
			    const_wmesh_int_p 	c2n_n_,
			    const_wmesh_int_p 	c2n_v_,
			    const_wmesh_int_p 	c2n_ld_,
			    
			    wmesh_int_t		c2c_e_size_,
			    const_wmesh_int_p 	c2c_e_ptr_,
			    const_wmesh_int_p 	c2c_e_m_,
			    const_wmesh_int_p 	c2c_e_n_,
			    wmesh_int_p 	c2c_e_v_,
			    const_wmesh_int_p 	c2c_e_ld_,
			    
			    wmesh_int_t		s_e2n_size_,
			    const_wmesh_int_p 	s_e2n_ptr_,
			    const_wmesh_int_p	s_e2n_m_,
			    const_wmesh_int_p	s_e2n_n_,
			    const_wmesh_int_p 	s_e2n_v_,
			    const_wmesh_int_p	s_e2n_ld_,
			    
			    wmesh_int_p		num_bfacets_,
			    wmesh_int_t		work_n_,
			    wmesh_int_p		work_)
  {
    //
    // if true, the algorithm won't search for more than one neighbor.
    // Meaning that the surface needs to be manifold.
    //
    bool match_mode = true;
    return bms_c2c_x<WMESH_ELEMENT_EDGE>(c2n_size_,
					 c2n_ptr_,
					 c2n_m_,
					 c2n_n_,
					 c2n_v_,
					 c2n_ld_,
					 
					 c2c_e_size_,
					 c2c_e_ptr_,
					 c2c_e_m_,
					 c2c_e_n_,
					 c2c_e_v_,
					 c2c_e_ld_,
					 
					 s_e2n_size_,
					 s_e2n_ptr_,
					 s_e2n_m_,
					 s_e2n_n_,
					 s_e2n_v_,
					 s_e2n_ld_,
					 
					 match_mode,
					 num_bfacets_,
					 work_n_,
					 work_);
    
  }

};
