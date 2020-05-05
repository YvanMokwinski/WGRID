#include <limits>
#include <iostream>
#include <array>
#include "wmesh.hpp"

static  inline void get_c2n(wmesh_int_t 		numCellNodes_,
			    const_wmesh_int_p 	cellsToNodes_,
			    wmesh_int_t 			cellsToNodesLd_,
			    wmesh_int_t 			cellIndex_,
			    wmesh_int_p 	cnc_)
{
  for (wmesh_int_t localNodeIndex=0;localNodeIndex<numCellNodes_;++localNodeIndex)
    {
      cnc_[localNodeIndex] = cellsToNodes_[cellIndex_*cellsToNodesLd_+localNodeIndex];      
    }
};

static  inline void get_c2e(wmesh_int_t 	m_,
			    const_wmesh_int_p 	c2e_,
			    wmesh_int_t 	c2e_ld_,
			    wmesh_int_t 	cell_idx_,
			    wmesh_int_p 	cnc_)
{
  for (wmesh_int_t i=0;i<m_;++i)
    {
      cnc_[i] = c2e_[cell_idx_*c2e_ld_+i];      
    }
};

static  inline void get_e2n(const_wmesh_int_p	c2n_,
			    const wmesh_int_t 	localEdgeIndex_,
			    wmesh_int_p		e2n_,
			    wmesh_int_t 	s_e2n_m_,
			    wmesh_int_t 	s_e2n_n_,
			    const_wmesh_int_p 	s_e2n_v_,
			    wmesh_int_t 	s_e2n_ld_)		    
{
  e2n_[0] = c2n_[s_e2n_v_[s_e2n_ld_ * localEdgeIndex_+ 0]];
  e2n_[1] = c2n_[s_e2n_v_[s_e2n_ld_ * localEdgeIndex_+ 1]];
};


static inline wmesh_status_t wmesh_space_indexing_edges_calculate(wmesh_int_t 			c2n_ptr_,
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
  wmesh_status_t  wmesh_space_indexing_edges(wmesh_int_t 	ntypes_,
					     
					     const_wmesh_int_p 	c2n_ptr_,
					     const_wmesh_int_p 	c2n_m_,
					     const_wmesh_int_p 	c2n_n_,
					     const_wmesh_int_p 	c2n_v_,
					     const_wmesh_int_p 	c2n_ld_,
					     
					     const_wmesh_int_p 	c2e_ptr_,
					     const_wmesh_int_p 	c2e_m_,
					     const_wmesh_int_p 	c2e_n_,
					     const_wmesh_int_p 	c2e_v_,
					     const_wmesh_int_p 	c2e_ld_,
					     
					     const_wmesh_int_p 	c2d_e_ptr_,
					     const_wmesh_int_p 	c2d_e_m_,
					     const_wmesh_int_p 	c2d_e_n_,
					     wmesh_int_p 	c2d_e_v_,
					     const_wmesh_int_p 	c2d_e_ld_,

					     const_wmesh_int_p 	s_e2n_ptr_,
					     const_wmesh_int_p	s_e2n_m_,
					     const_wmesh_int_p	s_e2n_n_,
					     const_wmesh_int_p 	s_e2n_v_,
					     const_wmesh_int_p	s_e2n_ld_,

					     wmesh_int_t 	num_edges_,
					     wmesh_int_t 	num_dofs_per_edge_,
					     wmesh_int_t 	dof_idx_origin_)
  {
    WMESH_POINTER_CHECK(c2n_ptr_);    
    WMESH_POINTER_CHECK(c2n_m_);    
    WMESH_POINTER_CHECK(c2n_n_);
    WMESH_POINTER_CHECK(c2n_v_);
    WMESH_POINTER_CHECK(c2n_ld_);

    WMESH_POINTER_CHECK(c2e_ptr_);    
    WMESH_POINTER_CHECK(c2e_m_);    
    WMESH_POINTER_CHECK(c2e_n_);
    WMESH_POINTER_CHECK(c2e_v_);
    WMESH_POINTER_CHECK(c2e_ld_);

    WMESH_POINTER_CHECK(c2d_e_ptr_);    
    WMESH_POINTER_CHECK(c2d_e_m_);    
    WMESH_POINTER_CHECK(c2d_e_n_);
    WMESH_POINTER_CHECK(c2d_e_v_);
    WMESH_POINTER_CHECK(c2d_e_ld_);

    WMESH_POINTER_CHECK(s_e2n_ptr_);    
    WMESH_POINTER_CHECK(s_e2n_m_);    
    WMESH_POINTER_CHECK(s_e2n_n_);
    WMESH_POINTER_CHECK(s_e2n_v_);
    WMESH_POINTER_CHECK(s_e2n_ld_);
    
    for (wmesh_int_t cell_type=0;cell_type<ntypes_;++cell_type)
      {	
	wmesh_status_t status = wmesh_space_indexing_edges_calculate(c2n_ptr_[cell_type],
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
