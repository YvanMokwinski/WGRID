#include <limits>
#include <iostream>
#include "wmesh_t.hpp"



static inline wmesh_status_t bms_c2d_e_calculate(wmesh_int_t 			c2n_ptr_,
								  wmesh_int_t 			c2n_m_,
								  wmesh_int_t 			c2n_n_,
								  const_wmesh_int_p		c2n_v_,
								  wmesh_int_t 			c2n_ld_,
								  
								  wmesh_int_t 			c2e_ptr_,
								  wmesh_int_t 			c2e_m_,
								  wmesh_int_t 			c2e_n_,
								  const_wmesh_int_p		c2e_v_,
								  wmesh_int_t 			c2e_ld_,
								  
								  wmesh_int_t 			c2d_e_ptr_,
								  wmesh_int_t 			c2d_e_m_,
								  wmesh_int_t 			c2d_e_n_,
								  wmesh_int_p			c2d_e_v_,
								  wmesh_int_t 			c2d_e_ld_,
								  
								  wmesh_int_t 			s_e2n_ptr_,
								  wmesh_int_t 			s_e2n_m_,
								  wmesh_int_t 			s_e2n_n_,
								  const_wmesh_int_p 		s_e2n_v_,
								  wmesh_int_t 			s_e2n_ld_,

								  wmesh_int_t 			num_edges_,
								  wmesh_int_t 			ndofs_per_edge_,
								  wmesh_int_t 			dof_idx_)
{

  wmesh_int_t
    e2n[2],
    c2e[12],
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

      //
      // Extract edge indexing.
      //
      get_c2e(c2e_m_,
	      c2e_v_ + c2e_ptr_,
	      c2e_ld_,
	      cell_idx,
	      c2e);

      //
      // Loop over the edges.
      //
      for (wmesh_int_t ledge_idx = 0;ledge_idx<c2e_m_;++ledge_idx)
	{
	  wmesh_int_t edge_idx = c2e[ledge_idx];
	  
	  get_e2n(c2n,
		  ledge_idx,
		  e2n,
		  s_e2n_m_,
		  s_e2n_n_,
		  s_e2n_v_ + s_e2n_ptr_,
		  s_e2n_ld_);
	  
	  if (e2n[0] < e2n[1])
	    {
	      for (wmesh_int_t ldof=0;ldof<ndofs_per_edge_;++ldof)
		{
		  c2d_e_v_[c2d_e_ptr_ + cell_idx * c2d_e_ld_ + ledge_idx * ndofs_per_edge_ + ldof] = dof_idx_ + edge_idx * ndofs_per_edge_ + ldof;
		}
	    }
	  else
	    {
	      for (wmesh_int_t ldof=0;ldof<ndofs_per_edge_;++ldof)
		{
		  c2d_e_v_[c2d_e_ptr_ + cell_idx * c2d_e_ld_ + ledge_idx * ndofs_per_edge_ + ldof] = dof_idx_ + edge_idx * ndofs_per_edge_ + ndofs_per_edge_ - 1 - ldof;
		}
	    }
	}
    }
  
  return WMESH_STATUS_SUCCESS;

};


extern "C"
{
  wmesh_status_t  bms_c2d_e(wmesh_int_t 	c2n_size_,
					     const_wmesh_int_p 	c2n_ptr_,
					     const_wmesh_int_p 	c2n_m_,
					     const_wmesh_int_p 	c2n_n_,
					     const_wmesh_int_p 	c2n_v_,
					     const_wmesh_int_p 	c2n_ld_,
					     
					     wmesh_int_t 	c2e_size_,
					     const_wmesh_int_p 	c2e_ptr_,
					     const_wmesh_int_p 	c2e_m_,
					     const_wmesh_int_p 	c2e_n_,
					     const_wmesh_int_p 	c2e_v_,
					     const_wmesh_int_p 	c2e_ld_,
					     
					     wmesh_int_t 	c2d_e_size_,
					     const_wmesh_int_p 	c2d_e_ptr_,
					     const_wmesh_int_p 	c2d_e_m_,
					     const_wmesh_int_p 	c2d_e_n_,
					     wmesh_int_p 	c2d_e_v_,
					     const_wmesh_int_p 	c2d_e_ld_,

					     wmesh_int_t 	s_e2n_size_,
					     const_wmesh_int_p 	s_e2n_ptr_,
					     const_wmesh_int_p	s_e2n_m_,
					     const_wmesh_int_p	s_e2n_n_,
					     const_wmesh_int_p 	s_e2n_v_,
					     const_wmesh_int_p	s_e2n_ld_,
					     
					     wmesh_int_t 	num_edges_,
					     wmesh_int_t 	num_dofs_per_edge_,
					     wmesh_int_t 	dof_idx_origin_)
  {
    WMESH_CHECK_POINTER(c2n_ptr_);    
    WMESH_CHECK_POINTER(c2n_m_);    
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(c2n_v_);
    WMESH_CHECK_POINTER(c2n_ld_);

    WMESH_CHECK_POINTER(c2e_ptr_);    
    WMESH_CHECK_POINTER(c2e_m_);    
    WMESH_CHECK_POINTER(c2e_n_);
    WMESH_CHECK_POINTER(c2e_v_);
    WMESH_CHECK_POINTER(c2e_ld_);

    WMESH_CHECK_POINTER(c2d_e_ptr_);    
    WMESH_CHECK_POINTER(c2d_e_m_);    
    WMESH_CHECK_POINTER(c2d_e_n_);
    WMESH_CHECK_POINTER(c2d_e_v_);
    WMESH_CHECK_POINTER(c2d_e_ld_);

    WMESH_CHECK_POINTER(s_e2n_ptr_);    
    WMESH_CHECK_POINTER(s_e2n_m_);    
    WMESH_CHECK_POINTER(s_e2n_n_);
    WMESH_CHECK_POINTER(s_e2n_v_);
    WMESH_CHECK_POINTER(s_e2n_ld_);
    
    for (wmesh_int_t cell_type=0;cell_type<c2n_size_;++cell_type)
      {	
	wmesh_status_t status = bms_c2d_e_calculate(c2n_ptr_[cell_type],
								     c2n_m_[cell_type],
								     c2n_n_[cell_type],
								     c2n_v_,
								     c2n_ld_[cell_type],
								     
								     c2e_ptr_[cell_type],
								     c2e_m_[cell_type],
								     c2e_n_[cell_type],
								     c2e_v_,
								     c2e_ld_[cell_type],

								     c2d_e_ptr_[cell_type],
								     c2d_e_m_[cell_type],
								     c2d_e_n_[cell_type],
								     c2d_e_v_,
								     c2d_e_ld_[cell_type],

								     s_e2n_ptr_[cell_type],
								     s_e2n_m_[cell_type],
								     s_e2n_n_[cell_type],
								     s_e2n_v_,
								     s_e2n_ld_[cell_type],
								     
								     num_edges_,
								     num_dofs_per_edge_,
								     dof_idx_origin_);
	WMESH_STATUS_CHECK(status);
      }
    return WMESH_STATUS_SUCCESS;
  }
};
