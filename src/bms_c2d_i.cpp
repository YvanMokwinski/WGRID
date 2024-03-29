#include <limits>
#include <iostream>
#include <array>
#include "wmesh_t.hpp"

static inline wmesh_status_t bms_c2d_i_calculate(wmesh_int_t 			c2d_i_ptr_,
						 wmesh_int_t 			c2d_i_m_,
						 wmesh_int_t 			c2d_i_n_,
						 wmesh_int_p			c2d_i_v_,
						 wmesh_int_t 			c2d_i_ld_,
						 wmesh_int_t 			dof_idx_)
{
  
  for (wmesh_int_t cell_idx = 0;cell_idx < c2d_i_n_;++cell_idx)
    {	  
      for (wmesh_int_t ldof=0;ldof<c2d_i_m_;++ldof)
	{
	  c2d_i_v_[c2d_i_ptr_ + cell_idx * c2d_i_ld_ + ldof] = dof_idx_ + cell_idx * c2d_i_m_ + ldof;
	}
    }
  
  return WMESH_STATUS_SUCCESS;

  
};

extern "C"
{
  wmesh_status_t  bms_c2d_i(wmesh_int_t 	c2d_i_size_,
			    const_wmesh_int_p 	c2d_i_ptr_,
			    const_wmesh_int_p 	c2d_i_m_,
			    const_wmesh_int_p 	c2d_i_n_,
			    wmesh_int_p 	c2d_i_v_,
			    const_wmesh_int_p 	c2d_i_ld_,
			    wmesh_int_p 	dof_idx_origin_)
  {

    WMESH_CHECK_POINTER(c2d_i_ptr_);    
    WMESH_CHECK_POINTER(c2d_i_m_);    
    WMESH_CHECK_POINTER(c2d_i_n_);
    WMESH_CHECK_POINTER(c2d_i_v_);
    WMESH_CHECK_POINTER(c2d_i_ld_);
    WMESH_CHECK_POINTER(dof_idx_origin_);

    wmesh_int_t dof_idx_origin = dof_idx_origin_[0];
    for (wmesh_int_t i=0;i<c2d_i_size_;++i)
      {	
	wmesh_status_t status = bms_c2d_i_calculate(c2d_i_ptr_[i],
									c2d_i_m_[i],
									c2d_i_n_[i],
									c2d_i_v_,
									c2d_i_ld_[i],
									dof_idx_origin);
	dof_idx_origin += c2d_i_m_[i] * c2d_i_n_[i];
	WMESH_STATUS_CHECK(status);
      }
    
    dof_idx_origin_[0] = dof_idx_origin;
    return WMESH_STATUS_SUCCESS;
  }
};
