
#include <limits>
#include <iostream>
#include <array>
#include "wmesh.hpp"

#include "GenericEncoding.hpp"

//!
//!
//!
static wmesh_int_t count=0;
struct wmesh_compare_triangle
{

public: 
  static inline bool are_same(const_wmesh_int_p this_t2n_,
			      const_wmesh_int_p that_t2n_) noexcept
  {
    ++count;
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


static inline void get_c2n(wmesh_int_t 		numCellNodes_,
			   const_wmesh_int_p 	cellsToNodes_,
			   wmesh_int_t 			cellsToNodesLd_,
			   wmesh_int_t 			cellIndex_,
			   wmesh_int_p 	cnc_)
{
  for (wmesh_int_t localNodeIndex=0;localNodeIndex<numCellNodes_;++localNodeIndex)
    {
      cnc_[localNodeIndex] = cellsToNodes_[cellIndex_*cellsToNodesLd_+localNodeIndex];      
    }
}


static  inline void get_t2n(const_wmesh_int_p	c2n_,
			    const wmesh_int_t 	t_lidx_,
			    wmesh_int_p		t2n_,
			    wmesh_int_t 		s_t2n_m_,
			    wmesh_int_t 		s_t2n_n_,
			    const_wmesh_int_p 	s_t2n_v_,
			    wmesh_int_t 		s_t2n_ld_)		    
{
  t2n_[0] = c2n_[s_t2n_v_[s_t2n_ld_ * t_lidx_ + 0]];
  t2n_[1] = c2n_[s_t2n_v_[s_t2n_ld_ * t_lidx_ + 1]];
  t2n_[2] = c2n_[s_t2n_v_[s_t2n_ld_ * t_lidx_ + 2]];
};

static wmesh_status_t wmesh_c2c_t_hash(wmesh_int_t 		cell_type_,
				       wmesh_int_t 		c2n_m_,
				       wmesh_int_t 		c2n_n_,
				       const_wmesh_int_p 	c2n_v_,
				       wmesh_int_t 		c2n_ld_,
				       wmesh_int_t 		c2e_m_,
				       wmesh_int_t 		c2e_n_,
				       wmesh_int_p 		c2e_v_,
				       wmesh_int_t 		c2e_ld_,				     
				       wmesh_int_t 		s_t2n_m_,
				       wmesh_int_t 		s_t2n_n_,
				       const_wmesh_int_p 	s_t2n_v_,
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


wmesh_status_t wmesh_c2c_t_search(wmesh_int_t 				cell_type_,
				     
				  const_wmesh_int_p 			c2n_ptr_,
				  const_wmesh_int_p 			c2n_m_,
				  const_wmesh_int_p 			c2n_n_,
				  const_wmesh_int_p 			c2n_v_,
				  const_wmesh_int_p 			c2n_ld_,
				     
				  const_wmesh_int_p 			c2c_t_ptr_,
				  const_wmesh_int_p 			c2c_t_m_,
				  const_wmesh_int_p 			c2c_t_n_,
				  wmesh_int_p 				c2c_t_v_,
				  const_wmesh_int_p 			c2c_t_ld_,
				     
				  const_wmesh_int_p 			s_t2n_ptr_,
				  const_wmesh_int_p 			s_t2n_m_,
				  const_wmesh_int_p 			s_t2n_n_,
				  const_wmesh_int_p 			s_t2n_v_,
				  const_wmesh_int_p 			s_t2n_ld_,
				     
				  wmesh_int_p				num_boundary_triangles_,
				  wmesh_int_t				work_n_,
				  wmesh_int_p				work_)
{

  WMESH_POINTER_CHECK(c2n_m_);
  WMESH_POINTER_CHECK(c2n_n_);
  WMESH_POINTER_CHECK(c2n_v_);
  WMESH_POINTER_CHECK(c2n_ld_);
  
  WMESH_POINTER_CHECK(c2c_t_m_);
  WMESH_POINTER_CHECK(c2c_t_n_);
  WMESH_POINTER_CHECK(c2c_t_v_);
  WMESH_POINTER_CHECK(c2c_t_ld_);
  
  WMESH_POINTER_CHECK(s_t2n_m_);
  WMESH_POINTER_CHECK(s_t2n_n_);
  WMESH_POINTER_CHECK(s_t2n_v_);
  WMESH_POINTER_CHECK(s_t2n_ld_);
  
  WMESH_POINTER_CHECK(num_boundary_triangles_);
  WMESH_POINTER_CHECK(work_);
  
  wmesh_int_t
    t2n[3],
    c2n[8],
    num_boundary_triangles = num_boundary_triangles_[0];

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
      for (int t_lidx = 0;t_lidx<c2c_t_m_[cell_type_];++t_lidx)
	{
	  //	  const wmesh_int_t start_half = -(c2c_t_ld_*cellIndex + 1*t_lidx + 0);
	  const wmesh_int_t start_half = - GenericEncoding<wmesh_int_t,2>::Encod(c2c_t_ld_[cell_type_] * cellIndex + t_lidx + 1, cell_type_);
	  wmesh_int_t next_half = c2c_t_v_[c2c_t_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-start_half)] + GenericEncoding<wmesh_int_t,2>::Up(-start_half)-1];
	  if (next_half < 0)
	    {
	      get_t2n(c2n,
		      t_lidx,
		      t2n,
		      s_t2n_m_[cell_type_],
		      s_t2n_n_[cell_type_],
		      s_t2n_v_ + s_t2n_ptr_[cell_type_],
		      s_t2n_ld_[cell_type_]);

	      bool matched = false;
	      //   std::cout << "starting with " << t2n[0] << " " << t2n[1] << std::endl;
	      //
	      // Assign the interior dofs to the current new triangle
	      //
	      c2c_t_v_[c2c_t_ptr_[cell_type_] + c2c_t_ld_[cell_type_] * cellIndex  + t_lidx] = 0;	      
	      wmesh_int_t last_half_non_matched = start_half;
	      while (next_half < 0 && next_half != default_hash)
		{		    
		  wmesh_int_t tested_cell_type = GenericEncoding<wmesh_int_t,2>::Low(-next_half);
		  wmesh_int_t tested_idx = GenericEncoding<wmesh_int_t,2>::Up(-next_half)-1;		  

		  half_triangle_decod(tested_idx,
				      c2c_t_ld_[tested_cell_type],
				      &tested_cell_idx,
				      &tested_t_lidx);
		  
		  //
		  // We need to extract the next half triangle before we overwrite it.
		  //
		  const wmesh_int_t next_halfBackup = next_half;
		  next_half = c2c_t_v_[c2c_t_ptr_[tested_cell_type] + tested_idx];

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
		      matched = true;
		      //
		      // We need to rebuild the linked list.
		      //
		      if (last_half_non_matched != start_half)
			{
			  const wmesh_int_t shift = c2c_t_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-last_half_non_matched)];
			  c2c_t_v_[ shift + ( GenericEncoding<wmesh_int_t,2>::Up(-last_half_non_matched) - 1) ] = next_half;
			}
		      
		      const wmesh_int_t at = c2c_t_ld_[tested_cell_type] * tested_cell_idx + tested_t_lidx;
		      
		      //		      c2c_t_v_[c2c_t_ptr_[tested_cell_type] + at] = t_idx-1;
		      //		      c2c_t_v_[c2c_t_ptr_[cell_type_]       + c2c_t_ld_[cell_type_] * cellIndex  + t_lidx] = t_idx-1;


		      c2c_t_v_[c2c_t_ptr_[tested_cell_type] + at] = GenericEncoding<wmesh_int_t,2>::Encod(cellIndex + 1,cell_type_);
		      c2c_t_v_[c2c_t_ptr_[cell_type_] + c2c_t_ld_[cell_type_] * cellIndex  + t_lidx] = GenericEncoding<wmesh_int_t,2>::Encod(tested_cell_idx + 1,tested_cell_type);
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
		  ++num_boundary_triangles;
		} 
	    }
	}
    }

  *num_boundary_triangles_ = num_boundary_triangles;  
  return WMESH_STATUS_SUCCESS;
};


extern "C"
{
  wmesh_status_t  wbms_c2c_t_calculate_buffer_size(wmesh_int_t		c2n_size_,
						   const_wmesh_int_p		c2n_n_,
						   wmesh_int_p 	work_n_)
  {
    WMESH_POINTER_CHECK(c2n_n_);
    WMESH_POINTER_CHECK(work_n_);
    wmesh_int_t mx_num_cells = 0;
    for (wmesh_int_t i=0;i<c2n_size_;++i)
      {
	//
	// higher it is, better it is.
	//
	mx_num_cells = std::max(mx_num_cells, c2n_n_[i] );
      }
    work_n_[0] = 1 + std::max(mx_num_cells,(wmesh_int_t)127);
    return WMESH_STATUS_SUCCESS;
  }
  
  wmesh_status_t  wbms_c2c_t_calculate(wmesh_int_t		c2n_size_,
							      
					const_wmesh_int_p 	c2n_ptr_,
					const_wmesh_int_p 	c2n_m_,
					const_wmesh_int_p 	c2n_n_,
					const_wmesh_int_p 	c2n_v_,
					const_wmesh_int_p 	c2n_ld_,
							      
					const_wmesh_int_p 	c2c_t_ptr_,
					const_wmesh_int_p 	c2c_t_m_,
					const_wmesh_int_p 	c2c_t_n_,
					wmesh_int_p 		c2c_t_v_,
					const_wmesh_int_p 	c2c_t_ld_,
							      
					const_wmesh_int_p 	s_t2n_ptr_,
					const_wmesh_int_p	s_t2n_m_,
					const_wmesh_int_p	s_t2n_n_,
					const_wmesh_int_p 	s_t2n_v_,
					const_wmesh_int_p	s_t2n_ld_,
							      
					wmesh_int_p		num_boundary_triangles_,
					wmesh_int_t		work_n_,
					wmesh_int_p		work_)
  {

    WMESH_POINTER_CHECK(c2n_n_);    
    WMESH_POINTER_CHECK(c2n_m_);
    WMESH_POINTER_CHECK(c2n_n_);
    WMESH_POINTER_CHECK(c2n_v_);
    WMESH_POINTER_CHECK(c2n_ld_);
  
    WMESH_POINTER_CHECK(c2c_t_m_);
    WMESH_POINTER_CHECK(c2c_t_n_);
    WMESH_POINTER_CHECK(c2c_t_v_);
    WMESH_POINTER_CHECK(c2c_t_ld_);

    WMESH_POINTER_CHECK(s_t2n_m_);
    WMESH_POINTER_CHECK(s_t2n_n_);
    WMESH_POINTER_CHECK(s_t2n_v_);
    WMESH_POINTER_CHECK(s_t2n_ld_);

    WMESH_POINTER_CHECK(num_boundary_triangles_);
    WMESH_POINTER_CHECK(work_);


    wmesh_status_t status;
    wmesh_int_t required_work_n;
    status = wbms_c2c_t_calculate_buffer_size(c2n_size_,
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
    for (wmesh_int_t i=c2n_size_-1;i>=0;--i)
      {
	if (c2n_n_[i]>0)
	  {
	    status = wmesh_c2c_t_hash(i,
				      c2n_m_[i],
				      c2n_n_[i],
				      c2n_v_ + c2n_ptr_[i],
				      c2n_ld_[i],
				      c2c_t_m_[i],
				      c2c_t_n_[i],
				      c2c_t_v_ + c2c_t_ptr_[i],
				      c2c_t_ld_[i],
				      s_t2n_m_[i],
				      s_t2n_n_[i],
				      s_t2n_v_ + s_t2n_ptr_[i],
				      s_t2n_ld_[i],
				      work_n_,
				      work_);	    
	    WMESH_STATUS_CHECK(status);
	  }
      }
    
    for (wmesh_int_t i=0;i<c2n_size_;++i)
      {
	status = wmesh_c2c_t_search(i,
				    c2n_ptr_,
				    c2n_m_,
				    c2n_n_,
				    c2n_v_,
				    c2n_ld_,
						      
				    c2c_t_ptr_,
				    c2c_t_m_,
				    c2c_t_n_,
				    c2c_t_v_,
				    c2c_t_ld_,
						      
				    s_t2n_ptr_,
				    s_t2n_m_,
				    s_t2n_n_,
				    s_t2n_v_,
				    s_t2n_ld_,
						      
				    num_boundary_triangles_,
				    work_n_,
				    work_);
	WMESH_STATUS_CHECK(status);
      }

    return WMESH_STATUS_SUCCESS;
  }
};
