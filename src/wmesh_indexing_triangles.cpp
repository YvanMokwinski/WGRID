#if 1
#include <limits>
#include <iostream>
#include <array>
#include "wmesh.hpp"

#include "GenericEncoding.hpp"
#include "bms_c2x_template.hpp"

extern "C"
{
  wmesh_status_t  wmesh_indexing_triangles_buffer_size(wmesh_int_t		c2n_size_,
						   const_wmesh_int_p	c2n_n_,
						   wmesh_int_p 		work_n_)
  {
    return bms_c2x_buffer_size<WMESH_ELEMENT_TRIANGLE>(c2n_size_,
						   c2n_n_,
						   work_n_);
  };
  
  wmesh_status_t  wmesh_indexing_triangles(wmesh_int_t 		c2n_size_,
					   
					   const_wmesh_int_p 	c2n_ptr_,
					   const_wmesh_int_p 	c2n_m_,
					   const_wmesh_int_p 	c2n_n_,
					   const_wmesh_int_p 	c2n_v_,
					   const_wmesh_int_p 	c2n_ld_,
					   
					   wmesh_int_t 		c2t_size_,
					   const_wmesh_int_p 	c2t_ptr_,
					   const_wmesh_int_p 	c2t_m_,
					   const_wmesh_int_p 	c2t_n_,
					   wmesh_int_p 		c2t_v_,
					   const_wmesh_int_p 	c2t_ld_,
					   
					   wmesh_int_t 		s_t2n_size_,
					   const_wmesh_int_p 	s_t2n_ptr_,
					   const_wmesh_int_p	s_t2n_m_,
					   const_wmesh_int_p	s_t2n_n_,
					   const_wmesh_int_p 	s_t2n_v_,
					   const_wmesh_int_p	s_t2n_ld_,
					   
					   wmesh_int_p		idx_,
					   wmesh_int_t		work_n_,
					   wmesh_int_p		work_)
  {
    
    bool match_mode = true;
    return bms_c2x<WMESH_ELEMENT_TRIANGLE>(c2n_size_,
					   
					   c2n_ptr_,
					   c2n_m_,
					   c2n_n_,
					   c2n_v_,
					   c2n_ld_,
					   
					   c2t_size_,
					   c2t_ptr_,
					   c2t_m_,
					   c2t_n_,
					   c2t_v_,
					   c2t_ld_,
					   
					   s_t2n_size_,
					   s_t2n_ptr_,
					   s_t2n_m_,
					   s_t2n_n_,
					   s_t2n_v_,
					   s_t2n_ld_,
					   
					   match_mode,
					   idx_,
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
struct wmesh_compare_triangle
{

public: 
  static inline bool are_same(const_wmesh_int_p this_t2n_,
			      const_wmesh_int_p that_t2n_) noexcept
  {
    //    ++count;
    if  ( ( (this_t2n_[0] == that_t2n_[2]) && (this_t2n_[1] == that_t2n_[1]) && (this_t2n_[2] == that_t2n_[0]) )
	  ||
	  ( (this_t2n_[0] == that_t2n_[1]) && (this_t2n_[1] == that_t2n_[0]) && (this_t2n_[2] == that_t2n_[2]) )
	  ||
	  ( (this_t2n_[0] == that_t2n_[0]) && (this_t2n_[1] == that_t2n_[2]) && (this_t2n_[2] == that_t2n_[1]) ) )
      {
	return true;
      }
#if 1
    if  ( ( (this_t2n_[0] == that_t2n_[0]) && (this_t2n_[1] == that_t2n_[1]) && (this_t2n_[2] == that_t2n_[2]) )
	  ||
	  ( (this_t2n_[0] == that_t2n_[1]) && (this_t2n_[1] == that_t2n_[2]) && (this_t2n_[2] == that_t2n_[0]) )
	  ||
	  ( (this_t2n_[0] == that_t2n_[2]) && (this_t2n_[1] == that_t2n_[0]) && (this_t2n_[2] == that_t2n_[1]) ) )
      {
	fprintf(stderr,"bad orientation\n");
	exit(1);
      }
#endif    
    return false;

  };
  
public:
  
  static inline wmesh_int_t hash(const_wmesh_int_p this_t2n_,
				 wmesh_int_t bound_) noexcept
  {
    wmesh_int_t a = (this_t2n_[0] < this_t2n_[1]) ? this_t2n_[0] : this_t2n_[1];
    wmesh_int_t b = (this_t2n_[0] < this_t2n_[1]) ? this_t2n_[1] : this_t2n_[0];    
    a = (a < this_t2n_[2]) ? a : this_t2n_[2];
    b = (b < this_t2n_[2]) ? this_t2n_[2] : b;
    wmesh_int_t sum = this_t2n_[0]+this_t2n_[1]+this_t2n_[2];    
    const wmesh_int_t h = (31 * a + 57*b + 79*sum) % bound_;
    return (h<0)
      ? -h
      : h;
  };  
  
};

static inline void half_triangle_decod(wmesh_int_t halfTriangleIndex_,
				       wmesh_int_t num_triangles_per_cell_,
				       wmesh_int_p cell_idx_,
				       wmesh_int_p ltriangle_idx_)
{
  *cell_idx_  = halfTriangleIndex_ / num_triangles_per_cell_;
  *ltriangle_idx_ = halfTriangleIndex_ % num_triangles_per_cell_;
};




wmesh_status_t wmesh_indexing_triangles_hash(wmesh_int_t 		cell_type_,
					     wmesh_int_t 		c2n_m_,
					     wmesh_int_t 		c2n_n_,
					     const_wmesh_int_p 		c2n_v_,
					     wmesh_int_t 		c2n_ld_,
					     wmesh_int_t 		c2e_m_,
					     wmesh_int_t 		c2e_n_,
					     wmesh_int_p 		c2e_v_,
					     wmesh_int_t 		c2e_ld_,				     
					     wmesh_int_t 		s_t2n_m_,
					     wmesh_int_t 		s_t2n_n_,
					     const_wmesh_int_p 		s_t2n_v_,
					     wmesh_int_t 		s_t2n_ld_,
					     wmesh_int_t		work_n_,
					     wmesh_int_p		work_)
{  
  const wmesh_int_t hash_size = work_[0];
  wmesh_int_t * hash_link = work_+1;

  wmesh_int_t t2n[3];	  
  wmesh_int_t c2n[8];  
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
      
      wmesh_int_t at = c2e_ld_ * cell_idx;
      for (wmesh_int_t t_lidx=c2e_m_-1;t_lidx>=0;--t_lidx)
	{ 
	  //
	  // Extract.
	  //
	  get_t2n(c2n,
		  t_lidx,
		  t2n,
		  s_t2n_m_,
		  s_t2n_n_,
		  s_t2n_v_,
		  s_t2n_ld_);
	  
	  //
	  // Compute hash value.
	  //
	  const wmesh_int_t hashValue = wmesh_compare_triangle::hash(t2n,
								     hash_size);
	  
	  //
	  // Assign the last half triangle index with the same hash value.
	  //
	  c2e_v_[at + t_lidx] = hash_link[hashValue];
	  
	  //
	  // Update the last half triangle index with respect to the hash value.
	  //
	  hash_link[hashValue] = - GenericEncoding<wmesh_int_t,2>::Encod(at + t_lidx + 1, cell_type_);
	}	
    }
  return WMESH_STATUS_SUCCESS;
};


wmesh_status_t wmesh_indexing_triangles_calculate(wmesh_int_t 				cell_type_,
						  
						  const_wmesh_int_p 			c2n_ptr_,
						  const_wmesh_int_p 			c2n_m_,
						  const_wmesh_int_p 			c2n_n_,
						  const_wmesh_int_p 			c2n_v_,
						  const_wmesh_int_p 			c2n_ld_,
						   
						  const_wmesh_int_p 			c2t_ptr_,
						  const_wmesh_int_p 			c2t_m_,
						  const_wmesh_int_p 			c2t_n_,
						  wmesh_int_p 				c2t_v_,
						  const_wmesh_int_p 			c2t_ld_,
						   
						  const_wmesh_int_p 			s_t2n_ptr_,
						  const_wmesh_int_p 			s_t2n_m_,
						  const_wmesh_int_p 			s_t2n_n_,
						  const_wmesh_int_p 			s_t2n_v_,
						  const_wmesh_int_p 			s_t2n_ld_,
						   
						  wmesh_int_p				t_idx_,
						  wmesh_int_t				work_n_,
						  wmesh_int_p				work_)
{

  WMESH_CHECK_POINTER(c2n_m_);
  WMESH_CHECK_POINTER(c2n_n_);
  WMESH_CHECK_POINTER(c2n_v_);
  WMESH_CHECK_POINTER(c2n_ld_);
  
  WMESH_CHECK_POINTER(c2t_m_);
  WMESH_CHECK_POINTER(c2t_n_);
  WMESH_CHECK_POINTER(c2t_v_);
  WMESH_CHECK_POINTER(c2t_ld_);
  
  WMESH_CHECK_POINTER(s_t2n_m_);
  WMESH_CHECK_POINTER(s_t2n_n_);
  WMESH_CHECK_POINTER(s_t2n_v_);
  WMESH_CHECK_POINTER(s_t2n_ld_);
  
  WMESH_CHECK_POINTER(t_idx_);
  WMESH_CHECK_POINTER(work_);
  
  wmesh_int_t
    t2n[3],
    c2n[8],
    t_idx = t_idx_[0];

  wmesh_int_t tested_cell_idx;
  wmesh_int_t tested_t_lidx;
  wmesh_int_t tested_t2n[3];
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
      for (int t_lidx = 0;t_lidx<c2t_m_[cell_type_];++t_lidx)
	{
	  //	  const wmesh_int_t start_half = -(c2t_ld_*cellIndex + 1*t_lidx + 0);
	  const wmesh_int_t start_half = - GenericEncoding<wmesh_int_t,2>::Encod(c2t_ld_[cell_type_] * cellIndex + t_lidx + 1, cell_type_);
	  wmesh_int_t next_half = c2t_v_[c2t_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-start_half)] + GenericEncoding<wmesh_int_t,2>::Up(-start_half)-1];
	  if (next_half < 0)
	    {
	      get_t2n(c2n,
		      t_lidx,
		      t2n,
		      s_t2n_m_[cell_type_],
		      s_t2n_n_[cell_type_],
		      s_t2n_v_ + s_t2n_ptr_[cell_type_],
		      s_t2n_ld_[cell_type_]);

	      //   std::cout << "starting with " << t2n[0] << " " << t2n[1] << std::endl;
	      ++t_idx;
	      
	      //
	      // Assign the interior dofs to the current new triangle
	      //
	      c2t_v_[c2t_ptr_[cell_type_] + c2t_ld_[cell_type_] * cellIndex  + t_lidx] = t_idx-1;
	      
	      wmesh_int_t last_half_non_matched = start_half;
	      while (next_half < 0 && next_half != default_hash)
		{		    
		  wmesh_int_t tested_cell_type = GenericEncoding<wmesh_int_t,2>::Low(-next_half);
		  wmesh_int_t tested_idx = GenericEncoding<wmesh_int_t,2>::Up(-next_half)-1;		  

		  half_triangle_decod(tested_idx,
				      c2t_ld_[tested_cell_type],
				      &tested_cell_idx,
				      &tested_t_lidx);
		  
		  //
		  // We need to extract the next half triangle before we overwrite it.
		  //
		  const wmesh_int_t next_halfBackup = next_half;
		  next_half = c2t_v_[c2t_ptr_[tested_cell_type] + tested_idx];

		  get_c2n(c2n_m_[tested_cell_type],
			  c2n_v_+c2n_ptr_[tested_cell_type],
			  c2n_ld_[tested_cell_type],
			  tested_cell_idx,
			  tested_c2n);
		  
		  get_t2n(tested_c2n,
			  tested_t_lidx,
			  tested_t2n,
			  s_t2n_m_[tested_cell_type],
			  s_t2n_n_[tested_cell_type],
			  s_t2n_v_ + s_t2n_ptr_[tested_cell_type],
			  s_t2n_ld_[tested_cell_type]);

		  if (wmesh_compare_triangle::are_same(t2n,
						       tested_t2n))
		    {
		      //
		      // We need to rebuild the linked list.
		      //
		      if (last_half_non_matched != start_half)
			{
			  const wmesh_int_t shift = c2t_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-last_half_non_matched)];
			  c2t_v_[ shift + ( GenericEncoding<wmesh_int_t,2>::Up(-last_half_non_matched) - 1) ] = next_half;
			}
		      
		      const wmesh_int_t at = c2t_ld_[tested_cell_type] * tested_cell_idx + tested_t_lidx;
		      c2t_v_[c2t_ptr_[tested_cell_type] + at] = t_idx-1;
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
	    }
	}
    }

  *t_idx_ = t_idx;  
  return WMESH_STATUS_SUCCESS;
};


extern "C"
{
  wmesh_status_t  wmesh_indexing_triangles_buffer_size(wmesh_int_t		c2n_size_,
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
  
  wmesh_status_t  wmesh_indexing_triangles(wmesh_int_t		c2n_size_,					   
					   const_wmesh_int_p 	c2n_ptr_,
					   const_wmesh_int_p 	c2n_m_,
					   const_wmesh_int_p 	c2n_n_,
					   const_wmesh_int_p 	c2n_v_,
					   const_wmesh_int_p 	c2n_ld_,
					   
					   wmesh_int_t 		c2t_size_,
					   const_wmesh_int_p 	c2t_ptr_,
					   const_wmesh_int_p 	c2t_m_,
					   const_wmesh_int_p 	c2t_n_,
					   wmesh_int_p 		c2t_v_,
					   const_wmesh_int_p 	c2t_ld_,
					   
					   wmesh_int_t 		s_t2n_size_,
					   const_wmesh_int_p 	s_t2n_ptr_,
					   const_wmesh_int_p	s_t2n_m_,
					   const_wmesh_int_p	s_t2n_n_,
					   const_wmesh_int_p 	s_t2n_v_,
					   const_wmesh_int_p	s_t2n_ld_,
					   
					   wmesh_int_p		triangle_idx_,
					   wmesh_int_t		work_n_,
					   wmesh_int_p		work_)
  {
    WMESH_CHECK_POINTER(c2n_m_);
    WMESH_CHECK_POINTER(c2n_n_);
    WMESH_CHECK_POINTER(c2n_v_);
    WMESH_CHECK_POINTER(c2n_ld_);
  
    WMESH_CHECK_POINTER(c2t_m_);
    WMESH_CHECK_POINTER(c2t_n_);
    WMESH_CHECK_POINTER(c2t_v_);
    WMESH_CHECK_POINTER(c2t_ld_);

    WMESH_CHECK_POINTER(s_t2n_m_);
    WMESH_CHECK_POINTER(s_t2n_n_);
    WMESH_CHECK_POINTER(s_t2n_v_);
    WMESH_CHECK_POINTER(s_t2n_ld_);
    WMESH_CHECK_POINTER(triangle_idx_);
    WMESH_CHECK_POINTER(work_);


    wmesh_status_t status;
    wmesh_int_t required_work_n;
    status = wmesh_indexing_triangles_buffer_size(c2n_size_,
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
    
    //
    //
    //
    //    std::cout << "triangle hashing ..." << std::endl;
    for (wmesh_int_t i=4-1;i>=0;--i)
      {
	if (c2n_n_[i]>0)
	  {
	    status = wmesh_indexing_triangles_hash(i,
								  c2n_m_[i],
								  c2n_n_[i],
								  c2n_v_ + c2n_ptr_[i],
								  c2n_ld_[i],
								  c2t_m_[i],
								  c2t_n_[i],
								  c2t_v_ + c2t_ptr_[i],
								  c2t_ld_[i],
								  s_t2n_m_[i],
								  s_t2n_n_[i],
								  s_t2n_v_ + s_t2n_ptr_[i],
								  s_t2n_ld_[i],
								  work_n_,
								  work_);	    
	    WMESH_STATUS_CHECK(status);
	  }
      }
    
#if 0
    for (int k=3;k>=0;--k)
      {
	if (c2t_n_[k]>0)
	  {
	    for (int j=c2t_n_[k]-1;j>=0;--j)
	      {
		for (int i=c2t_m_[k]-1;i>=0;--i)
		  {
		    std::cout << "triangle ("<< k << "," << j << "," << i << ")" <<  std::endl;
		    auto s = c2t_v_[c2t_ptr_[k]+j*c2t_ld_[k]+i];		    
		    static constexpr const wmesh_int_t d = - std::numeric_limits<wmesh_int_t>::max();
		    while (s < 0 && s != d)
		      {
			auto ltype = GenericEncoding<wmesh_int_t,2>::Low(-s);
			auto lidx = (GenericEncoding<wmesh_int_t,2>::Up(-s)-1)/c2t_ld_[ltype];
			auto llidx = (GenericEncoding<wmesh_int_t,2>::Up(-s)-1)%c2t_ld_[ltype];
			std::cout << "          linked triangle ("<< ltype << "," << lidx << "," << llidx << ")" <<  std::endl;
			s = c2t_v_[c2t_ptr_[ltype]+lidx*c2t_ld_[ltype]+llidx];		    
		      }		    
		  }
	      }
	  }
      }
#endif

    
#if 0
    //    for (int i=0;i<c2t_ptr_[4];++i) std::cout << " " << c2t_v_[i]<<std::endl;
    for (int i=c2t_ptr_[4]-1;i>=0;--i)
      {
	static constexpr const wmesh_int_t s = - std::numeric_limits<wmesh_int_t>::max();
	wmesh_int_t h = c2t_v_[i];
	if (h!=s)
	  {
	    wmesh_int_t cell_type = GenericEncoding<wmesh_int_t,2>::Low(-h-1);
	    wmesh_int_t idx  = GenericEncoding<wmesh_int_t,2>::Up(-h-1);
	    std::cout << " " << cell_type << " " << idx % c2t_ld_[cell_type] << std::endl;
	  }
	else
	  {
	    std::cout << " " << h << std::endl;
	    
	  }
      }
#endif
    
    //    std::cout << "triangle hashing done." << std::endl;
    //    exit(1);
#ifndef NDEBUG
    wmesh_int_t triangle_idx_bak = triangle_idx_[0];
#endif
    for (int i=0;i<4;++i)
      {
	wmesh_status_t status = wmesh_indexing_triangles_calculate(i,
								   c2n_ptr_,
								   c2n_m_,
								   c2n_n_,
								   c2n_v_,
								   c2n_ld_,
							   
								   c2t_ptr_,
								   c2t_m_,
								   c2t_n_,
								   c2t_v_,
								   c2t_ld_,
									  
								   s_t2n_ptr_,
								   s_t2n_m_,
								   s_t2n_n_,
								   s_t2n_v_,
								   s_t2n_ld_,
									  
								   triangle_idx_,
								   work_n_,
								   work_);
	WMESH_STATUS_CHECK(status);
      }

#ifndef NDEBUG
    std::cerr << "(FILE="
	      << __FILE__
	      << ",Line="
	      << __LINE__
	      << ") num triangles "
	      << triangle_idx_[0] - triangle_idx_bak
	      << std::endl;
#endif
    
    return WMESH_STATUS_SUCCESS;
  }
};
#endif
