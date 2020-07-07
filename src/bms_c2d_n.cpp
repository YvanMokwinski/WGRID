#include <limits>
#include <iostream>
#include "wmesh_t.hpp"

static inline wmesh_status_t bms_c2d_n_calculate(wmesh_int_t 			c2n_ptr_,
						 wmesh_int_t 			c2n_m_,
						 wmesh_int_t 			c2n_n_,
						 const_wmesh_int_p		c2n_v_,
						 wmesh_int_t 			c2n_ld_,
						 
						 
						 wmesh_int_t 			c2d_n_ptr_,
						 wmesh_int_t 			c2d_n_m_,
						 wmesh_int_t 			c2d_n_n_,
						 wmesh_int_p			c2d_n_v_,
						 wmesh_int_t 			c2d_n_ld_,
						 
						 wmesh_int_t 			num_nodes_,
						 wmesh_int_t 			ndofs_per_node_,
						 wmesh_int_t 			dof_idx_)
{
  
  wmesh_int_t
    c2n[8];
  
  for (wmesh_int_t cell_idx = 0;cell_idx < c2n_n_;++cell_idx)
    {
      
      //
      // Extract nodes.
      //
      get_c2n(c2n_m_,
	      c2n_v_ + c2n_ptr_,
	      c2n_ld_,
	      cell_idx,
	      c2n);

      for (wmesh_int_t lnode_idx = 0;lnode_idx<c2n_m_;++lnode_idx)
	{
	  wmesh_int_t node_idx = c2n[lnode_idx]-1;
	  for (wmesh_int_t ldof=0;ldof<ndofs_per_node_;++ldof)
	    {
	      c2d_n_v_[c2d_n_ptr_ + cell_idx * c2d_n_ld_ + lnode_idx * ndofs_per_node_ + ldof] = dof_idx_ + node_idx * ndofs_per_node_ + ldof;
	    }
	}
    }
  
  return WMESH_STATUS_SUCCESS;

};


extern "C"
{
  wmesh_status_t  bms_c2d_n(wmesh_int_t 	c2n_size_,
			    const_wmesh_int_p 	c2n_ptr_,
			    const_wmesh_int_p 	c2n_m_,
			    const_wmesh_int_p 	c2n_n_,
			    const_wmesh_int_p 	c2n_v_,
			    const_wmesh_int_p 	c2n_ld_,
			    
			    wmesh_int_t 	c2d_n_size_,
			    const_wmesh_int_p 	c2d_n_ptr_,
			    const_wmesh_int_p 	c2d_n_m_,
			    const_wmesh_int_p 	c2d_n_n_,
			    wmesh_int_p 	c2d_n_v_,
			    const_wmesh_int_p 	c2d_n_ld_,
			    
			    wmesh_int_t 	num_nodes_,
			    wmesh_int_t 	num_dofs_per_node_,
			    wmesh_int_t 	dof_idx_origin_)
  {
    WMESH_CHECK_POINTER(c2n_ptr_);    
    WMESH_CHECK_POINTER(c2n_m_);    
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(c2n_v_);
    WMESH_CHECK_POINTER(c2n_ld_);
    
    WMESH_CHECK_POINTER(c2d_n_ptr_);    
    WMESH_CHECK_POINTER(c2d_n_m_);    
    WMESH_CHECK_POINTER(c2d_n_n_);
    WMESH_CHECK_POINTER(c2d_n_v_);
    WMESH_CHECK_POINTER(c2d_n_ld_);

    for (wmesh_int_t cell_type=0;cell_type<c2n_size_;++cell_type)
      {	
	wmesh_status_t status = bms_c2d_n_calculate(c2n_ptr_[cell_type],
							     c2n_m_[cell_type],
							     c2n_n_[cell_type],
							     c2n_v_,
							     c2n_ld_[cell_type],
							     
							     c2d_n_ptr_[cell_type],
							     c2d_n_m_[cell_type],
							     c2d_n_n_[cell_type],
							     c2d_n_v_,
							     c2d_n_ld_[cell_type],
							     
							     num_nodes_,
							     num_dofs_per_node_,
							     dof_idx_origin_);
	WMESH_STATUS_CHECK(status);
      }
    return WMESH_STATUS_SUCCESS;
  }
};
