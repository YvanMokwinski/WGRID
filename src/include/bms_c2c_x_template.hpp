#pragma once
#include <limits>
#include <iostream>
#include <array>
#include "wmesh.hpp"

#include "GenericEncoding.hpp"
#include "bms_compare.hpp"


template<wmesh_int_t XTYPE>
wmesh_status_t  bms_c2c_x_buffer_size(wmesh_int_t		c2n_size_,
				      const_wmesh_int_p 	c2n_n_,						
				      wmesh_int_t*__restrict__ 	size_)
{
  WMESH_CHECK_POINTER(c2n_n_);    
  WMESH_CHECK_POINTER(size_);    
  wmesh_int_t mx_num_cells = 0;
  for (wmesh_int_t i=0;i<c2n_size_;++i)
    {
      //
      // higher it is, better it is.
      //
      mx_num_cells = std::max(mx_num_cells, c2n_n_[i] );
    }
    
  size_[0] = (1 + std::max(mx_num_cells,(wmesh_int_t)128));
  return WMESH_STATUS_SUCCESS;
}

static inline void half_decod(wmesh_int_t half_idx_,
			      wmesh_int_t num_entities_per_cell_,
			      wmesh_int_p cell_idx_,
			      wmesh_int_p local_idx_)
{
  *cell_idx_  = half_idx_ / num_entities_per_cell_;
  *local_idx_ = half_idx_ % num_entities_per_cell_;
};

template<wmesh_int_t XTYPE>
wmesh_status_t bms_c2c_x_hash(wmesh_int_t 		cell_type_,
			      wmesh_int_t 		c2n_m_,
			      wmesh_int_t 		c2n_n_,
			      const_wmesh_int_p 	c2n_v_,
			      wmesh_int_t 		c2n_ld_,
			      wmesh_int_t 		c2c_x_m_,
			      wmesh_int_t 		c2c_x_n_,
			      wmesh_int_p 		c2c_x_v_,
			      wmesh_int_t 		c2c_x_ld_,				     
			      wmesh_int_t 		s_x2n_m_,
			      wmesh_int_t 		s_x2n_n_,
			      const_wmesh_int_p 		s_x2n_v_,
			      wmesh_int_t 		s_x2n_ld_,
			      wmesh_int_t		work_n_,
			      wmesh_int_p		work_)
{  
  const wmesh_int_t 	hash_size = work_[0];
  wmesh_int_t * 	hash_link = work_+1;
  wmesh_int_t 		hash_dna[4];   
  wmesh_int_t 		x2n[4];	  
  wmesh_int_t 		c2n[8];  
  for (wmesh_int_t cell_idx=c2n_n_-1;cell_idx>=0;--cell_idx)
    {
       
      //
      // Copy connectivity.
      //
      get_c2n(c2n_m_,
	      c2n_v_,
	      c2n_ld_,
	      cell_idx,
	      c2n);
       
      wmesh_int_t at = c2c_x_ld_ * cell_idx;
      for (wmesh_int_t t_lidx=c2c_x_m_-1;t_lidx>=0;--t_lidx)
	{ 
	  //
	  // Extract.
	  //
	  get_x2n(c2n,
		  t_lidx,
		  x2n,
		  s_x2n_m_,
		  s_x2n_n_,
		  s_x2n_v_,
		  s_x2n_ld_);
	   
	  //
	  // Compute hash value.
	  //
	  
	  bms_compare<XTYPE>::hash_dna(x2n,
				       hash_dna);
	   
	   
	  const wmesh_int_t hashValue = bms_compare<XTYPE>::hash(hash_dna,
								 hash_size);
	   
	  //
	  // Assign the last half quadrilateral index with the same hash value.
	  //
	  c2c_x_v_[at + t_lidx] = hash_link[hashValue];
	  
	  //
	  // Update the last half quadrilateral index with respect to the hash value.
	  //
	  hash_link[hashValue] = - GenericEncoding<wmesh_int_t,2>::Encod(at + t_lidx + 1, cell_type_);
	}	
    }
   
  return WMESH_STATUS_SUCCESS;
};


template<wmesh_int_t XTYPE>
wmesh_status_t bms_c2c_x_search(wmesh_int_t 				cell_type_,
				   
				const_wmesh_int_p 			c2n_ptr_,
				const_wmesh_int_p 			c2n_m_,
				const_wmesh_int_p 			c2n_n_,
				const_wmesh_int_p 			c2n_v_,
				const_wmesh_int_p 			c2n_ld_,
				   
				const_wmesh_int_p 			c2c_x_ptr_,
				const_wmesh_int_p 			c2c_x_m_,
				const_wmesh_int_p 			c2c_x_n_,
				wmesh_int_p 				c2c_x_v_,
				const_wmesh_int_p 			c2c_x_ld_,
				   
				const_wmesh_int_p 			s_x2n_ptr_,
				const_wmesh_int_p 			s_x2n_m_,
				const_wmesh_int_p 			s_x2n_n_,
				const_wmesh_int_p 			s_x2n_v_,
				const_wmesh_int_p 			s_x2n_ld_,
				   
				bool 					match_mode_,
				wmesh_int_p				num_bfacets_,
				wmesh_int_t				work_n_,
				wmesh_int_p				work_)
{
  
  WMESH_CHECK_POINTER(c2n_m_);
  WMESH_CHECK_POINTER(c2n_n_);
  WMESH_CHECK_POINTER(c2n_v_);
  WMESH_CHECK_POINTER(c2n_ld_);
  
  WMESH_CHECK_POINTER(c2c_x_m_);
  WMESH_CHECK_POINTER(c2c_x_n_);
  WMESH_CHECK_POINTER(c2c_x_v_);
  WMESH_CHECK_POINTER(c2c_x_ld_);
  
  WMESH_CHECK_POINTER(s_x2n_m_);
  WMESH_CHECK_POINTER(s_x2n_n_);
  WMESH_CHECK_POINTER(s_x2n_v_);
  WMESH_CHECK_POINTER(s_x2n_ld_);
  
  WMESH_CHECK_POINTER(num_bfacets_);
  WMESH_CHECK_POINTER(work_);
  
  wmesh_int_t
    x2n[4],
    c2n[8],
    num_bfacets = num_bfacets_[0];

  wmesh_int_t tested_cell_idx;
  wmesh_int_t tested_t_lidx;
  wmesh_int_t tested_x2n[4];
  wmesh_int_t tested_c2n[8];
  wmesh_int_t x2n_dna[4];
  wmesh_int_t tested_x2n_dna[4];
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
      for (int t_lidx = 0;t_lidx<c2c_x_m_[cell_type_];++t_lidx)
	{
	  //	  const wmesh_int_t start_half = -(c2c_x_ld_*cellIndex + 1*t_lidx + 0);
	  const wmesh_int_t start_half = - GenericEncoding<wmesh_int_t,2>::Encod(c2c_x_ld_[cell_type_] * cellIndex + t_lidx + 1, cell_type_);
	  wmesh_int_t next_half = c2c_x_v_[c2c_x_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-start_half)] + GenericEncoding<wmesh_int_t,2>::Up(-start_half)-1];
	  if (next_half < 0)
	    {
	      get_x2n(c2n,
		      t_lidx,
		      x2n,
		      s_x2n_m_[cell_type_],
		      s_x2n_n_[cell_type_],
		      s_x2n_v_ + s_x2n_ptr_[cell_type_],
		      s_x2n_ld_[cell_type_]);

	      bms_compare<XTYPE>::hash_dna(x2n,
							    x2n_dna);

	      bool matched = false;
	      //   std::cout << "starting with " << x2n[0] << " " << x2n[1] << std::endl;
	      //
	      // Assign the interior dofs to the current new triangle
	      //
	      c2c_x_v_[c2c_x_ptr_[cell_type_] + c2c_x_ld_[cell_type_] * cellIndex  + t_lidx] = 0;	      
	      wmesh_int_t last_half_non_matched = start_half;
	      while (next_half < 0 && next_half != default_hash)
		{		    
		  wmesh_int_t tested_cell_type = GenericEncoding<wmesh_int_t,2>::Low(-next_half);
		  wmesh_int_t tested_idx = GenericEncoding<wmesh_int_t,2>::Up(-next_half)-1;		  
		  half_decod(tested_idx,
			     c2c_x_ld_[tested_cell_type],
			     &tested_cell_idx,
			     &tested_t_lidx);
		  //
		  // We need to extract the next half triangle before we overwrite it.
		  //
		  const wmesh_int_t next_halfBackup = next_half;
		  next_half = c2c_x_v_[c2c_x_ptr_[tested_cell_type] + tested_idx];

		  get_c2n(c2n_m_[tested_cell_type],
			  c2n_v_+c2n_ptr_[tested_cell_type],
			  c2n_ld_[tested_cell_type],
			  tested_cell_idx,
			  tested_c2n);
		  
		  get_x2n(tested_c2n,
			  tested_t_lidx,
			  tested_x2n,
			  s_x2n_m_[tested_cell_type],
			  s_x2n_n_[tested_cell_type],
			  s_x2n_v_ + s_x2n_ptr_[tested_cell_type],
			  s_x2n_ld_[tested_cell_type]);
		  

		  bms_compare<XTYPE>::hash_dna			(tested_x2n,
								 tested_x2n_dna);
	      
	      bool are_same = bms_compare<XTYPE>::comp_dna	(x2n_dna,
								 tested_x2n_dna);
	      
	      if (are_same)
		{
		  matched = true;
		  //
		  // We need to rebuild the linked list.
		  //
		  if (last_half_non_matched != start_half)
		    {
		      const wmesh_int_t shift = c2c_x_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-last_half_non_matched)];
		      c2c_x_v_[ shift + ( GenericEncoding<wmesh_int_t,2>::Up(-last_half_non_matched) - 1) ] = next_half;
		    }

		      
		  const wmesh_int_t at = c2c_x_ld_[tested_cell_type] * tested_cell_idx + tested_t_lidx;
		  c2c_x_v_[c2c_x_ptr_[tested_cell_type] + at] = GenericEncoding<wmesh_int_t,2>::Encod(cellIndex + 1,cell_type_);
		  c2c_x_v_[c2c_x_ptr_[cell_type_] + c2c_x_ld_[cell_type_] * cellIndex  + t_lidx] = GenericEncoding<wmesh_int_t,2>::Encod(tested_cell_idx + 1,tested_cell_type);
		  break;
		}
	      else
		{
		      
		  //
		  // this is NOT the same triangle.
		  //
		  last_half_non_matched = next_halfBackup;
		}
		}
	      
	      if (false == matched)
		{
		  ++num_bfacets;
		} 
	    }
	}
    }

  *num_bfacets_ = num_bfacets;  
  return WMESH_STATUS_SUCCESS;
};


template<wmesh_int_t XTYPE>
wmesh_status_t  bms_c2c_x(wmesh_int_t		c2n_size_,
			  const_wmesh_int_p 	c2n_ptr_,
			  const_wmesh_int_p 	c2n_m_,
			  const_wmesh_int_p 	c2n_n_,
			  const_wmesh_int_p 	c2n_v_,
			  const_wmesh_int_p 	c2n_ld_,
			  
			  wmesh_int_t		c2c_x_size_,
			  const_wmesh_int_p 	c2c_x_ptr_,
			  const_wmesh_int_p 	c2c_x_m_,
			  const_wmesh_int_p 	c2c_x_n_,
			  wmesh_int_p		c2c_x_v_,
			  const_wmesh_int_p 	c2c_x_ld_,
			  
			  wmesh_int_t		s_x2n_size_,
			  const_wmesh_int_p 	s_x2n_ptr_,
			  const_wmesh_int_p	s_x2n_m_,
			  const_wmesh_int_p	s_x2n_n_,
			  const_wmesh_int_p 	s_x2n_v_,
			  const_wmesh_int_p	s_x2n_ld_,
			  
			  bool			match_mode_,
			  wmesh_int_p		idx_,
			  wmesh_int_t		work_n_,
			  wmesh_int_p		work_)
{
  
  WMESH_CHECK_POINTER(c2n_ptr_);
  WMESH_CHECK_POINTER(c2n_m_);
  WMESH_CHECK_POINTER(c2n_n_);
  WMESH_CHECK_POINTER(c2n_v_);
  WMESH_CHECK_POINTER(c2n_ld_);
  
  WMESH_CHECK_POINTER(c2c_x_ptr_);
  WMESH_CHECK_POINTER(c2c_x_m_);
  WMESH_CHECK_POINTER(c2c_x_n_);
  WMESH_CHECK_POINTER(c2c_x_v_);
  WMESH_CHECK_POINTER(c2c_x_ld_);
    
  WMESH_CHECK_POINTER(s_x2n_ptr_);
  WMESH_CHECK_POINTER(s_x2n_m_);
  WMESH_CHECK_POINTER(s_x2n_n_);
  WMESH_CHECK_POINTER(s_x2n_v_);
  WMESH_CHECK_POINTER(s_x2n_ld_);

  WMESH_CHECK_POINTER(idx_);
  WMESH_CHECK_POINTER(work_);
    
  wmesh_int_t required_work_n;
  wmesh_status_t status;

    
  status = bms_c2c_x_buffer_size<XTYPE>(c2n_size_,
					c2n_n_,						
					&required_work_n);
    
  WMESH_STATUS_CHECK(status);
  if (required_work_n < work_n_)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
    }
    
  work_[0] = work_n_ - 1;    

  //
  // Initialize the hash links.
  //
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
	  status = bms_c2c_x_hash<XTYPE>(i,
					 c2n_m_[i],
					 c2n_n_[i],
					 c2n_v_ + c2n_ptr_[i],
					 c2n_ld_[i],
					 c2c_x_m_[i],
					 c2c_x_n_[i],
					 c2c_x_v_ + c2c_x_ptr_[i],
					 c2c_x_ld_[i],
					 s_x2n_m_[i],
					 s_x2n_n_[i],
					 s_x2n_v_ + s_x2n_ptr_[i],
					 s_x2n_ld_[i],
					 work_n_,
					 work_);   
	  WMESH_STATUS_CHECK(status);
	}
    }
    
#ifndef NDEBUG
  wmesh_int_t idx_bak = idx_[0];
#endif
    
  for (int cell_type=0;cell_type<c2n_size_;++cell_type)
    {
      status = bms_c2c_x_search<XTYPE>(cell_type,
					  
				       c2n_ptr_,
				       c2n_m_,
				       c2n_n_,
				       c2n_v_,
				       c2n_ld_,
					  
				       c2c_x_ptr_,
				       c2c_x_m_,
				       c2c_x_n_,
				       c2c_x_v_,
				       c2c_x_ld_,
					  
				       s_x2n_ptr_,
				       s_x2n_m_,
				       s_x2n_n_,
				       s_x2n_v_,
				       s_x2n_ld_,

				       match_mode_,
				       idx_,
				       work_n_,
				       work_);
      WMESH_STATUS_CHECK(status);
    }
    
#ifndef NDEBUG
  std::cerr << "(FILE="
	    << __FILE__
	    << ",Line="
	    << __LINE__
	    << ") num x "
	    << idx_[0] - idx_bak
	    << std::endl;
#endif

  return WMESH_STATUS_SUCCESS;
};


