#include <limits>
#include <iostream>
#include <array>
#include "wmesh.hpp"

#include "GenericEncoding.hpp"

//!
//!
//!
//static wmesh_int_t count=0;
struct wmesh_compare_edge
{

public: 
  static inline bool are_same(const wmesh_int_t thisEdgeToNodes_[],
			      const wmesh_int_t thatEdgeToNodes_[]) noexcept
  {
    //    ++count;
    return ( ( (thisEdgeToNodes_[0] == thatEdgeToNodes_[0]) && (thisEdgeToNodes_[1] == thatEdgeToNodes_[1]) )
	     ||
	     ( (thisEdgeToNodes_[1] == thatEdgeToNodes_[0]) && (thisEdgeToNodes_[0] == thatEdgeToNodes_[1]) ) );    
  };
  
public: 
  static inline wmesh_int_t hash(const wmesh_int_t thisEdgeToNodes_[],
				 const wmesh_int_t bound_) noexcept
  {
    wmesh_int_t a = (thisEdgeToNodes_[0] < thisEdgeToNodes_[1]) ? thisEdgeToNodes_[0] : thisEdgeToNodes_[1];
    wmesh_int_t b = (thisEdgeToNodes_[0] < thisEdgeToNodes_[1]) ? thisEdgeToNodes_[1] : thisEdgeToNodes_[0];
    //    const wmesh_int_t h = ( (thisEdgeToNodes_[0] < thisEdgeToNodes_[1]) ? thisEdgeToNodes_[0] : thisEdgeToNodes_[1] ) % bound_;
    const wmesh_int_t h = (53 * a + 79*b) % bound_;
    // const wmesh_int_t h = (thisEdgeToNodes_[0] + thisEdgeToNodes_[1]) % bound_;
    return (h<0)
      ? -h
      : h;
  };  
  
};

static  inline void half_edge_decod(wmesh_int_t halfEdgeIndex_,
				    wmesh_int_t num_edges_per_cell_,
				    wmesh_int_t ndofs_per_edge_,
				    wmesh_int_p cell_idx_,
				    wmesh_int_p ledge_idx_)
{
  *cell_idx_  = halfEdgeIndex_ / num_edges_per_cell_;
  *ledge_idx_ = (halfEdgeIndex_ % num_edges_per_cell_) / ndofs_per_edge_;
};




wmesh_status_t wmesh_indexing_space_edge_hashing(wmesh_int_t 		ntypes_,

						 const_wmesh_int_p	c2n_ptr_,
						 const_wmesh_int_p	c2n_m_,
						 const_wmesh_int_p	c2n_n_,
						 const_wmesh_int_p 	c2n_v_,
						 const_wmesh_int_p	c2n_ld_,
						 
						 const_wmesh_int_p	c2d_e_ptr_,
						 const_wmesh_int_p	c2d_e_m_,
						 const_wmesh_int_p	c2d_e_n_,
						 wmesh_int_p 		c2d_e_v_,
						 const_wmesh_int_p 	c2d_e_ld_,
						 
						 const_wmesh_int_p 	s_e2n_ptr_,
						 const_wmesh_int_p 	s_e2n_m_,
						 const_wmesh_int_p 	s_e2n_n_,
						 const_wmesh_int_p 	s_e2n_v_,
						 const_wmesh_int_p 	s_e2n_ld_,

						 wmesh_int_t 		ndofs_per_edge_,
						 wmesh_int_t		work_n_,
						 wmesh_int_p		work_)
{  
  static constexpr const wmesh_int_t s_default_hash = - std::numeric_limits<wmesh_int_t>::max();
  const wmesh_int_t hash_size = work_[0];
  wmesh_int_t * hash_link = work_+1;
  for (wmesh_int_t i=0;i<hash_size;++i)
    {
      hash_link[i] = s_default_hash;
    }

  wmesh_int_t edgeToNodes[2];	  
  wmesh_int_t cellToNodes[8];  
  for (wmesh_int_t itype=ntypes_-1;itype>=0;--itype)
    {
      wmesh_int_t c2n_ptr = c2n_ptr_[itype];
      wmesh_int_t c2n_n = c2n_n_[itype];
      wmesh_int_t c2n_m = c2n_m_[itype];
      wmesh_int_t c2n_ld = c2n_ld_[itype];

      wmesh_int_t nedges = c2n_m / ndofs_per_edge_;
      
      wmesh_int_t c2d_e_ptr = c2d_e_ptr_[itype];
      //      wmesh_int_t c2d_e_n = c2d_e_n_[itype];
      //      wmesh_int_t c2d_e_m = c2d_e_m_[itype];
      wmesh_int_t c2d_e_ld = c2d_e_ld_[itype];
      if (c2n_n>0)
	{
	  for (wmesh_int_t cell_idx=c2n_n-1;cell_idx>=0;--cell_idx)
	    {
	      
	      //
	      // Copy connectivity.
	      //
	      get_c2n(c2n_m,
		      c2n_v_ + c2n_ptr,
		      c2n_ld,
		      cell_idx,
		      cellToNodes);
	      
	      wmesh_int_t at = c2d_e_ld * cell_idx;
	      for (wmesh_int_t ledge_idx=nedges-1;ledge_idx>=0;--ledge_idx)
		{ 
		  //
		  // Extract the edgeToNodes.
		  //
		  get_e2n(cellToNodes,
			  ledge_idx,
			  edgeToNodes,
			  s_e2n_m_[itype],
			  s_e2n_n_[itype],
			  s_e2n_v_ + s_e2n_ptr_[itype],
			  s_e2n_ld_[itype]);
	  
		  //
		  // Compute hash value.
		  //
		  const wmesh_int_t hashValue = wmesh_compare_edge::hash(edgeToNodes,
									 hash_size);
		  
		  //
		  // Assign the last half edge index with the same hash value.
		  //
		  c2d_e_v_[c2d_e_ptr + at + ledge_idx * ndofs_per_edge_ + 0] = hash_link[hashValue];
		  
		  //
		  // Update the last half edge index with respect to the hash value.
		  //
		  hash_link[hashValue] = - GenericEncoding<wmesh_int_t,2>::Encod(at + ledge_idx * ndofs_per_edge_ + 0 +  1, itype);
		}	
	    }
	}
    }
  return WMESH_STATUS_SUCCESS;
};


wmesh_status_t wmesh_indexing_space_edges_calculate(wmesh_int_t				ntypes_,					  
						    
						    const_wmesh_int_p 			c2n_ptr_,
						    const_wmesh_int_p 			c2n_m_,
						    const_wmesh_int_p 			c2n_n_,
						    const_wmesh_int_p 			c2n_v_,
						    const_wmesh_int_p 			c2n_ld_,					  
						    
						    const_wmesh_int_p 			c2d_e_ptr_,
						    const_wmesh_int_p 			c2d_e_m_,
						    const_wmesh_int_p 			c2d_e_n_,
						    wmesh_int_p				c2d_e_v_,
						    const_wmesh_int_p 			c2d_e_ld_,					  
						    
						    const_wmesh_int_p 			s_e2n_ptr_,
						    const_wmesh_int_p 			s_e2n_m_,
						    const_wmesh_int_p 			s_e2n_n_,
						    const_wmesh_int_p 			s_e2n_v_,
						    const_wmesh_int_p 			s_e2n_ld_,
						    
						    wmesh_int_p				edge_idx_,
						    wmesh_int_p				dof_idx_,
						    wmesh_int_t				ndofs_per_edge_,
						    wmesh_int_t				work_n_,
						    wmesh_int_p				work_)
{

  WMESH_POINTER_CHECK(c2n_ptr_);
  WMESH_POINTER_CHECK(c2n_m_);
  WMESH_POINTER_CHECK(c2n_n_);
  WMESH_POINTER_CHECK(c2n_v_);
  WMESH_POINTER_CHECK(c2n_ld_);
  
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
  
  WMESH_POINTER_CHECK(dof_idx_);
  WMESH_POINTER_CHECK(edge_idx_);
  WMESH_POINTER_CHECK(work_);
  
  wmesh_int_t
    e2n[2],
    c2n[8],
    dof_idx = dof_idx_[0],
    edge_idx = edge_idx_[0];

  wmesh_int_t tested_cell_idx;
  wmesh_int_t tested_ledge_idx;
  wmesh_int_t tested_e2n[2];
  wmesh_int_t tested_c2n[8];

  //
  // For each types.
  //


  for (wmesh_int_t itype=0;itype < ntypes_;++itype)
    {
      static constexpr const wmesh_int_t default_hash = - std::numeric_limits<wmesh_int_t>::max();
      for (wmesh_int_t cellIndex = 0;cellIndex < c2n_n_[itype];++cellIndex)
	{
	  //
	  // Extract nodes.
	  //
	  get_c2n(c2n_m_[itype],
		  c2n_v_ + c2n_ptr_[itype],
		  c2n_ld_[itype],
		  cellIndex,
		  c2n);
	  //
	  // Reverse
	  //
	  for (int ledge_idx = 0;ledge_idx<c2d_e_m_[itype];++ledge_idx)
	    {
	      //	  const wmesh_int_t start_half = -(c2d_e_ld_*cellIndex + 1*ledge_idx + 0);
	      const wmesh_int_t start_half = - GenericEncoding<wmesh_int_t,2>::Encod(c2d_e_ld_[itype] * cellIndex + ledge_idx * ndofs_per_edge_ + 0 + 1, itype);
	      wmesh_int_t next_half = c2d_e_v_[c2d_e_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-start_half)] + GenericEncoding<wmesh_int_t,2>::Up(-start_half)-1];
	      if (next_half < 0)
		{
		  get_e2n(c2n,
			  ledge_idx,
			  e2n,
			  s_e2n_m_[itype],
			  s_e2n_n_[itype],
			  s_e2n_v_ + s_e2n_ptr_[itype],
			  s_e2n_ld_[itype]);

		  //   std::cout << "starting with " << e2n[0] << " " << e2n[1] << std::endl;
		  ++edge_idx;
	      
		  //
		  // Assign the interior dofs to the current new edge
		  //
		  if (e2n[0] < e2n[1])
		    {
		      for (wmesh_int_t k=0;k<ndofs_per_edge_;++k)
			{
			  c2d_e_v_[c2d_e_ptr_[itype] + c2d_e_ld_[itype] * cellIndex  + ndofs_per_edge_ * ledge_idx + k]
			    = dof_idx + ((edge_idx-1) * ndofs_per_edge_ + k);
			}
		    }
		  else
		    {
		      for (wmesh_int_t k=0;k<ndofs_per_edge_;++k)
			{
			  c2d_e_v_[c2d_e_ptr_[itype] + c2d_e_ld_[itype] * cellIndex  + ndofs_per_edge_ * ledge_idx + ndofs_per_edge_-1-k]
			    = dof_idx + ((edge_idx-1) * ndofs_per_edge_ + k);
			}
		    }
		  
		  wmesh_int_t last_half_non_matched = start_half;
		  while (next_half < 0 && next_half != default_hash)
		    {		    
		      wmesh_int_t tested_cell_type = GenericEncoding<wmesh_int_t,2>::Low(-next_half);
		      wmesh_int_t tested_at = GenericEncoding<wmesh_int_t,2>::Up(-next_half)-1;		  

		      half_edge_decod(tested_at,
				      c2d_e_ld_[tested_cell_type],
				      ndofs_per_edge_,
				      &tested_cell_idx,
				      &tested_ledge_idx);
		  
		      //
		      // We need to extract the next half edge before we overwrite it.
		      //
		      const wmesh_int_t next_halfBackup = next_half;
		      next_half = c2d_e_v_[c2d_e_ptr_[tested_cell_type] + tested_at];

		      get_c2n(c2n_m_[tested_cell_type],
			      c2n_v_+c2n_ptr_[tested_cell_type],
			      c2n_ld_[tested_cell_type],
			      tested_cell_idx,
			      tested_c2n);
		  
		      get_e2n(tested_c2n,
			      tested_ledge_idx,
			      tested_e2n,
			      s_e2n_m_[tested_cell_type],
			      s_e2n_n_[tested_cell_type],
			      s_e2n_v_ + s_e2n_ptr_[tested_cell_type],
			      s_e2n_ld_[tested_cell_type]);

		      if (wmesh_compare_edge::are_same(e2n,
						       tested_e2n))
			{
			  //
			  // We need to rebuild the linked list.
			  //
			  if (last_half_non_matched != start_half)
			    {
			      const wmesh_int_t shift = c2d_e_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-last_half_non_matched)];
			      c2d_e_v_[ shift + ( GenericEncoding<wmesh_int_t,2>::Up(-last_half_non_matched) - 1) ] = next_half;
			    }
		      
			  //  const wmesh_int_t at = c2d_e_ld_[tested_cell_type] * tested_cell_idx + tested_ledge_idx;
			  //  c2d_e_v_[c2d_e_ptr_[tested_cell_type] + at] = edge_idx-1;
			  if (tested_e2n[0] < tested_e2n[1])
			    {
			      for (wmesh_int_t k=0;k<ndofs_per_edge_;++k)
				{
				  c2d_e_v_[c2d_e_ptr_[tested_cell_type] + c2d_e_ld_[tested_cell_type] * tested_cell_idx  + ndofs_per_edge_ * tested_ledge_idx + k]
				    = dof_idx + ((edge_idx-1) * ndofs_per_edge_ + k);
				}
			    }
			  else
			    {
			      for (wmesh_int_t k=0;k<ndofs_per_edge_;++k)
				{
				  c2d_e_v_[c2d_e_ptr_[tested_cell_type] + c2d_e_ld_[tested_cell_type] * tested_cell_idx  + ndofs_per_edge_ * tested_ledge_idx + ndofs_per_edge_ - 1 - k]
				    = dof_idx + ((edge_idx-1) * ndofs_per_edge_ + k);
				}
			    }
			  
			}
		      else
			{
			  //
			  // this is NOT the same edge.
			  //
			  last_half_non_matched = next_halfBackup;
			}
		    }		
		}
	    }
	}
    }
  
  *dof_idx_  += (edge_idx-edge_idx_[0]) * ndofs_per_edge_;
  *edge_idx_  = edge_idx;  
  return WMESH_STATUS_SUCCESS;
}


extern "C"
{
  wmesh_status_t  wmesh_indexing_space(wmesh_int_t 	ntypes_,
				       
				       const_wmesh_int_p c2n_ptr_,
				       const_wmesh_int_p c2n_m_,
				       const_wmesh_int_p c2n_n_,
				       const_wmesh_int_p c2n_v_,
				       const_wmesh_int_p c2n_ld_,
				       
				       const_wmesh_int_p c2d_e_ptr_,
				       const_wmesh_int_p c2d_e_m_,
				       const_wmesh_int_p c2d_e_n_,
				       wmesh_int_p 	 c2d_e_v_,
				       const_wmesh_int_p c2d_e_ld_,

				       const_wmesh_int_p c2d_t_ptr_,
				       const_wmesh_int_p c2d_t_m_,
				       const_wmesh_int_p c2d_t_n_,
				       wmesh_int_p 	 c2d_t_v_,
				       const_wmesh_int_p c2d_t_ld_,

				       const_wmesh_int_p c2d_q_ptr_,
				       const_wmesh_int_p c2d_q_m_,
				       const_wmesh_int_p c2d_q_n_,
				       wmesh_int_p 	 c2d_q_v_,
				       const_wmesh_int_p c2d_q_ld_,

				       const_wmesh_int_p s_e2n_ptr_,
				       const_wmesh_int_p s_e2n_m_,
				       const_wmesh_int_p s_e2n_n_,
				       const_wmesh_int_p s_e2n_v_,
				       const_wmesh_int_p s_e2n_ld_,

				       const_wmesh_int_p s_t2n_ptr_,
				       const_wmesh_int_p s_t2n_m_,
				       const_wmesh_int_p s_t2n_n_,
				       const_wmesh_int_p s_t2n_v_,
				       const_wmesh_int_p s_t2n_ld_,

				       const_wmesh_int_p s_q2n_ptr_,
				       const_wmesh_int_p s_q2n_m_,
				       const_wmesh_int_p s_q2n_n_,
				       const_wmesh_int_p s_q2n_v_,
				       const_wmesh_int_p s_q2n_ld_,
				       
				       wmesh_int_p	dof_idx_,
				       wmesh_int_t 	ndofs_per_edge_,
				       wmesh_int_p	work_n_,
				       wmesh_int_p	work_)
  {
    WMESH_POINTER_CHECK(c2n_n_);    
    wmesh_int_t mx_num_cells = 0;
    for (wmesh_int_t i=0;i<ntypes_;++i)
      {
	//
	// higher it is, better it is.
	//
	mx_num_cells = std::max(mx_num_cells, c2n_n_[i]);
      }
    
    WMESH_POINTER_CHECK(work_n_);
    if (work_ == nullptr)
      {
	work_n_[0] = 1 + std::max(mx_num_cells,(wmesh_int_t)128);
	return WMESH_STATUS_SUCCESS;
      }
    else if (work_n_[0] < 1 + std::max(mx_num_cells,(wmesh_int_t)128))
      {
	return WMESH_STATUS_ERROR_WORKSPACE; 
      }    
    work_[0] = std::max(mx_num_cells,(wmesh_int_t)128);

    WMESH_POINTER_CHECK(dof_idx_);
    
    //    std::cout << "edge index before " << edge_idx_[0] << std::endl;
    wmesh_int_t status;

    wmesh_int_t num_edges = 0;
    
    //
    // Indexing dofs on nodes.
    //



    
    //
    // Indexing dofs on edges.
    //
    {
      static constexpr const wmesh_int_t s_default_hash = - std::numeric_limits<wmesh_int_t>::max();  
      const wmesh_int_t hash_size = work_[0];
      wmesh_int_t * hash_link = work_+1;
      for (wmesh_int_t i=0;i<hash_size;++i)
	{
	  hash_link[i] = s_default_hash;
	}
      
      status = wmesh_indexing_space_edge_hashing(ntypes_,
						 c2n_ptr_,
						 c2n_m_,
						 c2n_n_,
						 c2n_v_,
						 c2n_ld_,

						 c2d_e_ptr_,
						 c2d_e_m_,
						 c2d_e_n_,
						 c2d_e_v_,
						 c2d_e_ld_,

						 s_e2n_ptr_,
						 s_e2n_m_,
						 s_e2n_n_,
						 s_e2n_v_,
						 s_e2n_ld_,

						 ndofs_per_edge_,
						 work_n_[0],
						 work_);
      
      WMESH_STATUS_CHECK(status);

      status = wmesh_indexing_space_edges_calculate(ntypes_,
					  
						    c2n_ptr_,
						    c2n_m_,
						    c2n_n_,
						    c2n_v_,
						    c2n_ld_,
					  
						    c2d_e_ptr_,
						    c2d_e_m_,
						    c2d_e_n_,
						    c2d_e_v_,
						    c2d_e_ld_,
						    
						    s_e2n_ptr_,
						    s_e2n_m_,
						    s_e2n_n_,
						    s_e2n_v_,
						    s_e2n_ld_,
						    
						    &num_edges,
						    dof_idx_,
						    ndofs_per_edge_,
						    work_n_[0],
						    work_);
      WMESH_STATUS_CHECK(status);
    }



    
    //
    // Indexing dofs on triangles.
    //

    //
    // Indexing dofs on quadrilaterals.
    //

    //
    // Indexing dofs on interior.
    //

    return WMESH_STATUS_SUCCESS;    
  }

};
