#if 1
#include <limits>
#include <iostream>
#include <array>
#include "wmesh.hpp"

#include "GenericEncoding.hpp"
#include "bms_c2x_template.hpp"

extern "C"
{
  wmesh_status_t  wmesh_indexing_edges_buffer_size(wmesh_int_t		c2n_size_,
						   const_wmesh_int_p	c2n_n_,
						   wmesh_int_p 		work_n_)
  {
    return bms_c2x_buffer_size<WMESH_ELEMENT_EDGE>(c2n_size_,
						   c2n_n_,
						   work_n_);
  };
  
  wmesh_status_t  wmesh_indexing_edges(wmesh_int_t 		c2n_size_,
				       
				       const_wmesh_int_p 	c2n_ptr_,
				       const_wmesh_int_p 	c2n_m_,
				       const_wmesh_int_p 	c2n_n_,
				       const_wmesh_int_p 	c2n_v_,
				       const_wmesh_int_p 	c2n_ld_,
				       
				       wmesh_int_t 		c2e_size_,
				       const_wmesh_int_p 	c2e_ptr_,
				       const_wmesh_int_p 	c2e_m_,
				       const_wmesh_int_p 	c2e_n_,
				       wmesh_int_p 		c2e_v_,
				       const_wmesh_int_p 	c2e_ld_,
				       
				       wmesh_int_t 		s_e2n_size_,
				       const_wmesh_int_p 	s_e2n_ptr_,
				       const_wmesh_int_p	s_e2n_m_,
				       const_wmesh_int_p	s_e2n_n_,
				       const_wmesh_int_p 	s_e2n_v_,
				       const_wmesh_int_p	s_e2n_ld_,
				       
				       wmesh_int_p		edge_idx_,
				       wmesh_int_t		work_n_,
				       wmesh_int_p		work_)
  {
    
    bool match_mode = false;
    return bms_c2x<WMESH_ELEMENT_EDGE>(c2n_size_,
				       
				       c2n_ptr_,
				       c2n_m_,
				       c2n_n_,
				       c2n_v_,
				       c2n_ld_,
				       
				       c2e_size_,
				       c2e_ptr_,
				       c2e_m_,
				       c2e_n_,
				       c2e_v_,
				       c2e_ld_,
				       
				       s_e2n_size_,
				       s_e2n_ptr_,
				       s_e2n_m_,
				       s_e2n_n_,
				       s_e2n_v_,
				       s_e2n_ld_,

				       match_mode,
				       edge_idx_,
				       work_n_,
				       work_);
  }
};

#else
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
				    wmesh_int_p cell_idx_,
				    wmesh_int_p ledge_idx_)
{
  *cell_idx_  = halfEdgeIndex_ / num_edges_per_cell_;
  *ledge_idx_ = halfEdgeIndex_ % num_edges_per_cell_;
};





wmesh_status_t wmesh_indexing_edges_hash(wmesh_int_t 			cell_type_,
						wmesh_int_t 			c2n_m_,
						wmesh_int_t 			c2n_n_,
						const_wmesh_int_p 	c2n_v_,
						wmesh_int_t 			c2n_ld_,
						wmesh_int_t 			c2e_m_,
						wmesh_int_t 			c2e_n_,
						wmesh_int_t* __restrict__ 	c2e_v_,
						wmesh_int_t 			c2e_ld_,				     
						wmesh_int_t 			s_e2n_m_,
						wmesh_int_t 			s_e2n_n_,
						const_wmesh_int_p s_e2n_v_,
						wmesh_int_t 			s_e2n_ld_,
						wmesh_int_t			work_n_,
						wmesh_int_t*			work_)
{  
  //  static constexpr const wmesh_int_t s_default_hash = - std::numeric_limits<wmesh_int_t>::max();
  const wmesh_int_t hash_size = work_[0];
  wmesh_int_t * hash_link = work_+1;

  wmesh_int_t edgeToNodes[2];	  
  wmesh_int_t cellToNodes[8];  
  for (wmesh_int_t cell_idx=c2n_n_-1;cell_idx>=0;--cell_idx)
    {

      //
      // Copy connectivity.
      //
      get_c2n(c2n_m_,
	      c2n_v_,
	      c2n_ld_,
	      cell_idx,
	      cellToNodes);
      
      wmesh_int_t at = c2e_ld_ * cell_idx;
      for (wmesh_int_t ledge_idx=c2e_m_-1;ledge_idx>=0;--ledge_idx)
	{ 
	  //
	  // Extract the edgeToNodes.
	  //
	  get_e2n(cellToNodes,
		  ledge_idx,
		  edgeToNodes,
		  s_e2n_m_,
		  s_e2n_n_,
		  s_e2n_v_,
		  s_e2n_ld_);
	  
	  //
	  // Compute hash value.
	  //
	  const wmesh_int_t hashValue = wmesh_compare_edge::hash(edgeToNodes,
								 hash_size);
	  
	  //
	  // Assign the last half edge index with the same hash value.
	  //
	  c2e_v_[at + ledge_idx] = hash_link[hashValue];
	  
	  //
	  // Update the last half edge index with respect to the hash value.
	  //
	  hash_link[hashValue] = - GenericEncoding<wmesh_int_t,2>::Encod(at + ledge_idx + 1, cell_type_);
	}	
    }
#if 0
  static constexpr const wmesh_int_t default_hash = - std::numeric_limits<wmesh_int_t>::max();
  wmesh_int_t j=0;
  for (wmesh_int_t i=0;i<hash_size;++i)
    if (hash_link[i]==default_hash)
      {
	++j;
      }
  std::cout << "  j " <<j<< std::endl;
#endif
  return WMESH_STATUS_SUCCESS;
};


wmesh_status_t wmesh_indexing_edges_calculate(wmesh_int_t 				cell_type_,
					      
					      const_wmesh_int_p 			c2n_ptr_,
					      const_wmesh_int_p 			c2n_m_,
					      const_wmesh_int_p 			c2n_n_,
					      const_wmesh_int_p 			c2n_v_,
					      const_wmesh_int_p 			c2n_ld_,
					      
					      const_wmesh_int_p 			c2e_ptr_,
					      const_wmesh_int_p 			c2e_m_,
					      const_wmesh_int_p 			c2e_n_,
					      wmesh_int_p 				c2e_v_,
					      const_wmesh_int_p 			c2e_ld_,
					      
					      const_wmesh_int_p 			s_e2n_ptr_,
					      const_wmesh_int_p 			s_e2n_m_,
					      const_wmesh_int_p 			s_e2n_n_,
					      const_wmesh_int_p 			s_e2n_v_,
					      const_wmesh_int_p 			s_e2n_ld_,
					      
					      wmesh_int_p 				edge_idx_,
					      wmesh_int_t				work_n_,
					      wmesh_int_p				work_)
{

  WMESH_CHECK_POINTER(c2n_m_);
  WMESH_CHECK_POINTER(c2n_n_);
  WMESH_CHECK_POINTER(c2n_v_);
  WMESH_CHECK_POINTER(c2n_ld_);
  
  WMESH_CHECK_POINTER(c2e_m_);
  WMESH_CHECK_POINTER(c2e_n_);
  WMESH_CHECK_POINTER(c2e_v_);
  WMESH_CHECK_POINTER(c2e_ld_);
  
  WMESH_CHECK_POINTER(s_e2n_m_);
  WMESH_CHECK_POINTER(s_e2n_n_);
  WMESH_CHECK_POINTER(s_e2n_v_);
  WMESH_CHECK_POINTER(s_e2n_ld_);
  
  WMESH_CHECK_POINTER(edge_idx_);
  WMESH_CHECK_POINTER(work_);
  
  wmesh_int_t
    e2n[2],
    c2n[8],
    edge_idx = edge_idx_[0];

  wmesh_int_t tested_cell_idx;
  wmesh_int_t tested_ledge_idx;
  wmesh_int_t tested_e2n[2];
  wmesh_int_t tested_c2n[8];
  
  static constexpr const wmesh_int_t default_hash = - std::numeric_limits<wmesh_int_t>::max();
  for (wmesh_int_t cellIndex = 0;cellIndex < c2n_n_[cell_type_];++cellIndex)
    {
      //
      // Extract nodes.
      //
      get_c2n(c2n_m_[cell_type_],
	      c2n_v_ + c2n_ptr_[cell_type_],
	      c2n_ld_[cell_type_],
	      cellIndex,
	      c2n);
      
      //
      // Reverse
      //
      for (int ledge_idx = 0;ledge_idx<c2e_m_[cell_type_];++ledge_idx)
	{
	  //	  const wmesh_int_t start_half = -(c2e_ld_*cellIndex + 1*ledge_idx + 0);
	  const wmesh_int_t start_half = - GenericEncoding<wmesh_int_t,2>::Encod(c2e_ld_[cell_type_] * cellIndex + ledge_idx + 1, cell_type_);
	  wmesh_int_t next_half = c2e_v_[c2e_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-start_half)] + GenericEncoding<wmesh_int_t,2>::Up(-start_half)-1];
	  if (next_half < 0)
	    {
	      get_e2n(c2n,
		      ledge_idx,
		      e2n,
		      s_e2n_m_[cell_type_],
		      s_e2n_n_[cell_type_],
		      s_e2n_v_ + s_e2n_ptr_[cell_type_],
		      s_e2n_ld_[cell_type_]);

	      //   std::cout << "starting with " << e2n[0] << " " << e2n[1] << std::endl;
	      ++edge_idx;
	      
	      //
	      // Assign the interior dofs to the current new edge
	      //
	      c2e_v_[c2e_ptr_[cell_type_] + c2e_ld_[cell_type_] * cellIndex  + ledge_idx] = edge_idx-1;
	      
	      wmesh_int_t last_half_non_matched = start_half;
	      while (next_half < 0 && next_half != default_hash)
		{		    
		  wmesh_int_t tested_cell_type = GenericEncoding<wmesh_int_t,2>::Low(-next_half);
		  wmesh_int_t tested_idx = GenericEncoding<wmesh_int_t,2>::Up(-next_half)-1;		  

		  half_edge_decod(tested_idx,
				  c2e_ld_[tested_cell_type],
				  &tested_cell_idx,
				  &tested_ledge_idx);
		  
		  //
		  // We need to extract the next half edge before we overwrite it.
		  //
		  const wmesh_int_t next_halfBackup = next_half;
		  next_half = c2e_v_[c2e_ptr_[tested_cell_type] + tested_idx];

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
			  const wmesh_int_t shift = c2e_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-last_half_non_matched)];
			  c2e_v_[ shift + ( GenericEncoding<wmesh_int_t,2>::Up(-last_half_non_matched) - 1) ] = next_half;
			}
		      
		      const wmesh_int_t at = c2e_ld_[tested_cell_type] * tested_cell_idx + tested_ledge_idx;
		      c2e_v_[c2e_ptr_[tested_cell_type] + at] = edge_idx-1;
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

  *edge_idx_ = edge_idx;  
  return WMESH_STATUS_SUCCESS;
};


extern "C"
{
  wmesh_status_t  wmesh_indexing_edges_buffer_size(wmesh_int_t		c2n_size_,
						   const_wmesh_int_p	c2n_n_,
						   wmesh_int_p 		work_n_)
  {
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(work_n_);
    wmesh_int_t mx_num_cells = 0;
    for (wmesh_int_t i=0;i<c2n_size_;++i)
      {
	//
	// higher it is, better it is.
	//
	mx_num_cells = std::max(mx_num_cells, c2n_n_[i] );
      }
    work_n_[0] = 1 + std::max(mx_num_cells,(wmesh_int_t)128);
    return WMESH_STATUS_SUCCESS;
  };
  
  
  wmesh_status_t  wmesh_indexing_edges(wmesh_int_t 		c2n_size_,
				       
				       const_wmesh_int_p 	c2n_ptr_,
				       const_wmesh_int_p 	c2n_m_,
				       const_wmesh_int_p 	c2n_n_,
				       const_wmesh_int_p 	c2n_v_,
				       const_wmesh_int_p 	c2n_ld_,
				       
				       wmesh_int_t 		c2e_size_,
				       const_wmesh_int_p 	c2e_ptr_,
				       const_wmesh_int_p 	c2e_m_,
				       const_wmesh_int_p 	c2e_n_,
				       wmesh_int_p 		c2e_v_,
				       const_wmesh_int_p 	c2e_ld_,
				       
				       wmesh_int_t 		s_e2n_size_,
				       const_wmesh_int_p 	s_e2n_ptr_,
				       const_wmesh_int_p	s_e2n_m_,
				       const_wmesh_int_p	s_e2n_n_,
				       const_wmesh_int_p 	s_e2n_v_,
				       const_wmesh_int_p	s_e2n_ld_,
				       
				       wmesh_int_p		edge_idx_,
				       wmesh_int_t		work_n_,
				       wmesh_int_p		work_)
  {

    WMESH_CHECK_POINTER(c2n_m_);
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(c2n_v_);
    WMESH_CHECK_POINTER(c2n_ld_);
  
    WMESH_CHECK_POINTER(c2e_m_);
    WMESH_CHECK_POINTER(c2e_n_);
    WMESH_CHECK_POINTER(c2e_v_);
    WMESH_CHECK_POINTER(c2e_ld_);

    WMESH_CHECK_POINTER(s_e2n_m_);
    WMESH_CHECK_POINTER(s_e2n_n_);
    WMESH_CHECK_POINTER(s_e2n_v_);
    WMESH_CHECK_POINTER(s_e2n_ld_);
    WMESH_CHECK_POINTER(edge_idx_);
    WMESH_CHECK_POINTER(work_);

    wmesh_status_t status;
    wmesh_int_t required_work_n;
    status = wmesh_indexing_edges_buffer_size(c2n_size_,
					      c2n_n_,
					      &required_work_n);
    WMESH_STATUS_CHECK(status);
    if (work_n_ < required_work_n)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
      }
    work_[0] = work_n_ - 1;
    
    //    std::cout << "edge index before " << edge_idx_[0] << std::endl;
    static constexpr const wmesh_int_t s_default_hash = - std::numeric_limits<wmesh_int_t>::max();

    const wmesh_int_t hash_size = work_[0];
    wmesh_int_t * hash_link = work_+1;
    for (wmesh_int_t i=0;i<hash_size;++i)
      {
	hash_link[i] = s_default_hash;
      }
    
    //
    //
    //
    //    std::cout << "edge hashing ..." << std::endl;
    for (wmesh_int_t i=4-1;i>=0;--i)
      {
	if (c2n_n_[i]>0)
	  {
	    status = wmesh_indexing_edges_hash(i,
					       c2n_m_[i],
					       c2n_n_[i],
					       c2n_v_ + c2n_ptr_[i],
					       c2n_ld_[i],

					       c2e_m_[i],
					       c2e_n_[i],
					       c2e_v_ + c2e_ptr_[i],
					       c2e_ld_[i],
					       s_e2n_m_[i],
					       s_e2n_n_[i],
					       s_e2n_v_ + s_e2n_ptr_[i],
					       s_e2n_ld_[i],
					       work_n_,
					       work_);	    
	    WMESH_STATUS_CHECK(status);
	  }
      }
    
    
#if 0
    for (int k=3;k>=0;--k)
      {
	if (c2e_n_[k]>0)
	  {
	    for (int j=c2e_n_[k]-1;j>=0;--j)
	      {
		for (int i=c2e_m_[k]-1;i>=0;--i)
		  {
		    std::cout << "edge ("<< k << "," << j << "," << i << ")" <<  std::endl;
		    auto s = c2e_v_[c2e_ptr_[k]+j*c2e_ld_[k]+i];		    
		    static constexpr const wmesh_int_t d = - std::numeric_limits<wmesh_int_t>::max();
		    while (s < 0 && s != d)
		      {
			auto ltype = GenericEncoding<wmesh_int_t,2>::Low(-s);
			auto lidx = (GenericEncoding<wmesh_int_t,2>::Up(-s)-1)/c2e_ld_[ltype];
			auto llidx = (GenericEncoding<wmesh_int_t,2>::Up(-s)-1)%c2e_ld_[ltype];
			std::cout << "          linked edge ("<< ltype << "," << lidx << "," << llidx << ")" <<  std::endl;
			s = c2e_v_[c2e_ptr_[ltype]+lidx*c2e_ld_[ltype]+llidx];		    
		      }		    
		  }
	      }
	  }
      }
#endif


    
#if 0
    //    for (int i=0;i<c2e_ptr_[4];++i) std::cout << " " << c2e_v_[i]<<std::endl;
    for (int i=c2e_ptr_[4]-1;i>=0;--i)
      {
	static constexpr const wmesh_int_t s = - std::numeric_limits<wmesh_int_t>::max();
	wmesh_int_t h = c2e_v_[i];
	if (h!=s)
	  {
	    wmesh_int_t cell_type = GenericEncoding<wmesh_int_t,2>::Low(-h-1);
	    wmesh_int_t idx  = GenericEncoding<wmesh_int_t,2>::Up(-h-1);
	    std::cout << " " << cell_type << " " << idx % c2e_ld_[cell_type] << std::endl;
	  }
	else
	  {
	    std::cout << " " << h << std::endl;
	    
	  }
      }
#endif    
    //    std::cout << "edge hashing done." << std::endl;
    //    exit(1);
#ifndef NDEBUG
    wmesh_int_t edge_idx_bak = edge_idx_[0];
#endif
    for (int i=0;i<4;++i)
      {
	status = wmesh_indexing_edges_calculate(i,
						c2n_ptr_,
						c2n_m_,
						c2n_n_,
						c2n_v_,
						c2n_ld_,
							   
						c2e_ptr_,
						c2e_m_,
						c2e_n_,
						c2e_v_,
						c2e_ld_,
							   
						s_e2n_ptr_,
						s_e2n_m_,
						s_e2n_n_,
						s_e2n_v_,
						s_e2n_ld_,
							   
						edge_idx_,
						work_n_,
						work_);
	WMESH_STATUS_CHECK(status);
      }

#ifndef NDEBUG
    std::cerr << "(FILE="
	      << __FILE__
	      << ",Line="
	      << __LINE__
	      << ") num edges "
	      << edge_idx_[0] - edge_idx_bak
	      << std::endl;
#endif
    return WMESH_STATUS_SUCCESS;
  }
};
#endif
