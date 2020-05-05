
#include <limits>
#include <iostream>
#include <array>
#include "wmesh.hpp"

#include "GenericEncoding.hpp"

struct wmesh_compare_quadrilateral
{
public:  
  static inline void hash_dna(const_wmesh_int_p	this_q2n_,
			      wmesh_int_t 	dna[2]) noexcept
  {
    wmesh_int_t a = this_q2n_[0];
    int ia = 0;
    for (int i=1;i<4;++i)
      {
	if (a > this_q2n_[i])
	  {
	    a = this_q2n_[i];
	    ia = i;
	  }
      }
    dna[0] = a;
    dna[1] = this_q2n_[(ia+2)%4];
  };  
  
  static inline wmesh_int_t hash(const_wmesh_int_p this_q2n_,
				 wmesh_int_t bound_) noexcept
  {
    wmesh_int_t dna[2];
    
    hash_dna(this_q2n_,
	     dna);

    const wmesh_int_t h = (31 * dna[0] + 57*dna[1]) % bound_;
    return (h<0)
      ? -h
      : h;
  };  
  
};

static inline void half_quadrilateral_decod(wmesh_int_t halfQuadrilateralIndex_,
					    wmesh_int_t num_quadrilaterals_per_cell_,
					    wmesh_int_p cell_idx_,
					    wmesh_int_p lquadrilateral_idx_)
{
  *cell_idx_  = halfQuadrilateralIndex_ / num_quadrilaterals_per_cell_;
  *lquadrilateral_idx_ = halfQuadrilateralIndex_ % num_quadrilaterals_per_cell_;
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

static  inline void get_q2n(const_wmesh_int_p	c2n_,
			    const wmesh_int_t 	t_lidx_,
			    wmesh_int_p		q2n_,
			    wmesh_int_t 		s_q2n_m_,
			    wmesh_int_t 		s_q2n_n_,
			    const_wmesh_int_p 	s_q2n_v_,
			    wmesh_int_t 		s_q2n_ld_)		    
{
  q2n_[0] = c2n_[s_q2n_v_[s_q2n_ld_ * t_lidx_ + 0]];
  q2n_[1] = c2n_[s_q2n_v_[s_q2n_ld_ * t_lidx_ + 1]];
  q2n_[2] = c2n_[s_q2n_v_[s_q2n_ld_ * t_lidx_ + 2]];
  q2n_[3] = c2n_[s_q2n_v_[s_q2n_ld_ * t_lidx_ + 3]];
};

wmesh_status_t wmesh_c2c_q_hash(wmesh_int_t 		cell_type_,
				wmesh_int_t 		c2n_m_,
				wmesh_int_t 		c2n_n_,
				const_wmesh_int_p 	c2n_v_,
				wmesh_int_t 		c2n_ld_,
				wmesh_int_t 		c2q_m_,
				wmesh_int_t 		c2q_n_,
				wmesh_int_p 		c2q_v_,
				wmesh_int_t 		c2q_ld_,				     
				wmesh_int_t 		s_q2n_m_,
				wmesh_int_t 		s_q2n_n_,
				const_wmesh_int_p 	s_q2n_v_,
				wmesh_int_t 		s_q2n_ld_,
				wmesh_int_t		work_n_,
				wmesh_int_p		work_)
{  
  const wmesh_int_t hash_size = work_[0];
  wmesh_int_t * hash_link = work_+1;
  
  wmesh_int_t q2n[4];	  
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
      
      wmesh_int_t at = c2q_ld_ * cell_idx;
      for (wmesh_int_t t_lidx=c2q_m_-1;t_lidx>=0;--t_lidx)
	{ 
	  //
	  // Extract.
	  //
	  get_q2n(c2n,
		  t_lidx,
		  q2n,
		  s_q2n_m_,
		  s_q2n_n_,
		  s_q2n_v_,
		  s_q2n_ld_);
	  
	  //
	  // Compute hash value.
	  //
	  const wmesh_int_t hashValue = wmesh_compare_quadrilateral::hash(q2n,
									  hash_size);
	  
	  //
	  // Assign the last half quadrilateral index with the same hash value.
	  //
	  c2q_v_[at + t_lidx] = hash_link[hashValue];
	  
	  //
	  // Update the last half quadrilateral index with respect to the hash value.
	  //
	  hash_link[hashValue] = - GenericEncoding<wmesh_int_t,2>::Encod(at + t_lidx + 1, cell_type_);
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


wmesh_status_t wmesh_c2c_q_search(wmesh_int_t 				cell_type_,
				  
				  const_wmesh_int_p 			c2n_ptr_,
				  const_wmesh_int_p 			c2n_m_,
				  const_wmesh_int_p 			c2n_n_,
				  const_wmesh_int_p 			c2n_v_,
				  const_wmesh_int_p 			c2n_ld_,
				  
				  const_wmesh_int_p 			c2c_q_ptr_,
				  const_wmesh_int_p 			c2c_q_m_,
				  const_wmesh_int_p 			c2c_q_n_,
				  wmesh_int_p 				c2c_q_v_,
				  const_wmesh_int_p 			c2c_q_ld_,
						   
				  const_wmesh_int_p 			s_q2n_ptr_,
				  const_wmesh_int_p 			s_q2n_m_,
				  const_wmesh_int_p 			s_q2n_n_,
				  const_wmesh_int_p 			s_q2n_v_,
				  const_wmesh_int_p 			s_q2n_ld_,
						   
				  wmesh_int_p				t_idx_,
				  wmesh_int_t				work_n_,
				  wmesh_int_p				work_)
{

  WMESH_POINTER_CHECK(c2n_m_);
  WMESH_POINTER_CHECK(c2n_n_);
  WMESH_POINTER_CHECK(c2n_v_);
  WMESH_POINTER_CHECK(c2n_ld_);
  
  WMESH_POINTER_CHECK(c2c_q_m_);
  WMESH_POINTER_CHECK(c2c_q_n_);
  WMESH_POINTER_CHECK(c2c_q_v_);
  WMESH_POINTER_CHECK(c2c_q_ld_);
  
  WMESH_POINTER_CHECK(s_q2n_m_);
  WMESH_POINTER_CHECK(s_q2n_n_);
  WMESH_POINTER_CHECK(s_q2n_v_);
  WMESH_POINTER_CHECK(s_q2n_ld_);
  
  WMESH_POINTER_CHECK(t_idx_);
  WMESH_POINTER_CHECK(work_);
  
  wmesh_int_t
    q2n[4],
    c2n[8],
    t_idx = t_idx_[0];

  wmesh_int_t tested_cell_idx;
  wmesh_int_t tested_t_lidx;
  wmesh_int_t tested_q2n[4];
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
      for (int t_lidx = 0;t_lidx < c2c_q_m_[cell_type_];++t_lidx)
	{
	  //	  const wmesh_int_t start_half = -(c2c_q_ld_*cellIndex + 1*t_lidx + 0);
	  const wmesh_int_t start_half = - GenericEncoding<wmesh_int_t,2>::Encod(c2c_q_ld_[cell_type_] * cellIndex + t_lidx + 1, cell_type_);
	  wmesh_int_t next_half = c2c_q_v_[ c2c_q_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-start_half)] + GenericEncoding<wmesh_int_t,2>::Up(-start_half) - 1];
	  if (next_half < 0)
	    {
	      get_q2n(c2n,
		      t_lidx,
		      q2n,
		      s_q2n_m_[cell_type_],
		      s_q2n_n_[cell_type_],
		      s_q2n_v_ + s_q2n_ptr_[cell_type_],
		      s_q2n_ld_[cell_type_]);


	      
	      //
	      // Assign the interior dofs to the current new quadrilateral
	      //
	      c2c_q_v_[c2c_q_ptr_[cell_type_] + c2c_q_ld_[cell_type_] * cellIndex  + t_lidx] = 0;
	      
	      wmesh_int_t last_half_non_matched = start_half;

	      wmesh_int_t q2n_dna[2];
	      if (next_half < 0 && next_half != default_hash)
		{
		  wmesh_compare_quadrilateral::hash_dna(q2n,
							q2n_dna);
		}
	      bool matched = false;
	      while (next_half < 0 && next_half != default_hash)
		{		    
		  wmesh_int_t tested_cell_type = GenericEncoding<wmesh_int_t,2>::Low(-next_half);
		  wmesh_int_t tested_idx = GenericEncoding<wmesh_int_t,2>::Up(-next_half)-1;		  

		  half_quadrilateral_decod(tested_idx,
					   c2c_q_ld_[tested_cell_type],
					   &tested_cell_idx,
					   &tested_t_lidx);
		  
		  //
		  // We need to extract the next half quadrilateral before we overwrite it.
		  //
		  const wmesh_int_t next_halfBackup = next_half;
		  next_half = c2c_q_v_[c2c_q_ptr_[tested_cell_type] + tested_idx];

		  get_c2n(c2n_m_[tested_cell_type],
			  c2n_v_+c2n_ptr_[tested_cell_type],
			  c2n_ld_[tested_cell_type],
			  tested_cell_idx,
			  tested_c2n);
		  
		  get_q2n(tested_c2n,
			  tested_t_lidx,
			  tested_q2n,
			  s_q2n_m_[tested_cell_type],
			  s_q2n_n_[tested_cell_type],
			  s_q2n_v_ + s_q2n_ptr_[tested_cell_type],
			  s_q2n_ld_[tested_cell_type]);

		  
		  wmesh_int_t tested_q2n_dna[2];
		  wmesh_compare_quadrilateral::hash_dna(tested_q2n,
							tested_q2n_dna);

		  if (q2n_dna[0]==tested_q2n_dna[0] && q2n_dna[1] == tested_q2n_dna[1])
		    {
		      matched = true;
		      //
		      // We need to rebuild the linked list.
		      //
		      if (last_half_non_matched != start_half)
			{
			  const wmesh_int_t shift = c2c_q_ptr_[GenericEncoding<wmesh_int_t,2>::Low(-last_half_non_matched)];
			  c2c_q_v_[ shift + ( GenericEncoding<wmesh_int_t,2>::Up(-last_half_non_matched) - 1) ] = next_half;
			}
		      
		      const wmesh_int_t at = c2c_q_ld_[tested_cell_type] * tested_cell_idx + tested_t_lidx;


		      c2c_q_v_[c2c_q_ptr_[tested_cell_type] + at] = GenericEncoding<wmesh_int_t,2>::Encod(cellIndex + 1,cell_type_);
		      c2c_q_v_[c2c_q_ptr_[cell_type_] + c2c_q_ld_[cell_type_] * cellIndex  + t_lidx] = GenericEncoding<wmesh_int_t,2>::Encod(tested_cell_idx + 1,tested_cell_type);
		      break;
		    }
		  else
		    {
		      //
		      // this is NOT the same quadrilateral.
		      //
		      last_half_non_matched = next_halfBackup;
		    }
		}

	      if (false == matched)
		{
		  ++t_idx;
		}
	      
	    }
	}
    }

  *t_idx_ = t_idx;  
  return WMESH_STATUS_SUCCESS;
};


extern "C"
{

  wmesh_status_t  wbms_c2c_q_calculate_buffer_size(wmesh_int_t		c2n_size_,
						    const_wmesh_int_p	c2n_n_,
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
    work_n_[0] = 1 + std::max(mx_num_cells,(wmesh_int_t)128);
    return WMESH_STATUS_SUCCESS;
  }
  
  wmesh_status_t  wbms_c2c_q_calculate(wmesh_int_t		c2n_size_,
					
					const_wmesh_int_p 	c2n_ptr_,
					const_wmesh_int_p 	c2n_m_,
					const_wmesh_int_p 	c2n_n_,
					const_wmesh_int_p 	c2n_v_,
					const_wmesh_int_p 	c2n_ld_,
					  
					const_wmesh_int_p 	c2c_q_ptr_,
					const_wmesh_int_p 	c2c_q_m_,
					const_wmesh_int_p 	c2c_q_n_,
					wmesh_int_p 		c2c_q_v_,
					const_wmesh_int_p 	c2c_q_ld_,
					  
					const_wmesh_int_p 	s_q2n_ptr_,
					const_wmesh_int_p	s_q2n_m_,
					const_wmesh_int_p	s_q2n_n_,
					const_wmesh_int_p 	s_q2n_v_,
					const_wmesh_int_p	s_q2n_ld_,
				     
					wmesh_int_p		num_boundary_quadrilaterals_,
					wmesh_int_t		work_n_,
					wmesh_int_p		work_)
  {

    WMESH_POINTER_CHECK(c2n_m_);
    WMESH_POINTER_CHECK(c2n_n_);
    WMESH_POINTER_CHECK(c2n_v_);
    WMESH_POINTER_CHECK(c2n_ld_);
  
    WMESH_POINTER_CHECK(c2c_q_m_);
    WMESH_POINTER_CHECK(c2c_q_n_);
    WMESH_POINTER_CHECK(c2c_q_v_);
    WMESH_POINTER_CHECK(c2c_q_ld_);

    WMESH_POINTER_CHECK(s_q2n_m_);
    WMESH_POINTER_CHECK(s_q2n_n_);
    WMESH_POINTER_CHECK(s_q2n_v_);
    WMESH_POINTER_CHECK(s_q2n_ld_);

    WMESH_POINTER_CHECK(num_boundary_quadrilaterals_);
    WMESH_POINTER_CHECK(work_);

    wmesh_status_t status;
    wmesh_int_t required_work_n;
    status = wbms_c2c_q_calculate_buffer_size(c2n_size_,
					       c2n_n_,
					       &required_work_n);
    
    if (work_n_ < required_work_n)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
      }
    work_[0] = work_n_ - 1;
    
    //    std::cout << "quadrilateral index before " << num_boundary_quadrilaterals_[0] << std::endl;
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
    //    std::cout << "quadrilateral hashing ..." << std::endl;
    for (wmesh_int_t cell_type=4-1;cell_type>=0;--cell_type)
      {
	if (c2n_n_[cell_type]>0)
	  {
	     status = wmesh_c2c_q_hash(cell_type,
						     c2n_m_[cell_type],
						     c2n_n_[cell_type],
						     c2n_v_ + c2n_ptr_[cell_type],
						     c2n_ld_[cell_type],
						     c2c_q_m_[cell_type],
						     c2c_q_n_[cell_type],
						     c2c_q_v_ + c2c_q_ptr_[cell_type],
						     c2c_q_ld_[cell_type],
						     s_q2n_m_[cell_type],
						     s_q2n_n_[cell_type],
						     s_q2n_v_ + s_q2n_ptr_[cell_type],
						     s_q2n_ld_[cell_type],
						     work_n_,
						     work_);	    
	    WMESH_STATUS_CHECK(status);
	  }
      }
    
    
    //    std::cout << "quadrilateral hashing done." << std::endl;
    //    exit(1);
    for (int i=0;i<4;++i)
      {
	status = wmesh_c2c_q_search(i,
						   c2n_ptr_,
						   c2n_m_,
						   c2n_n_,
						   c2n_v_,
						   c2n_ld_,
							   
						   c2c_q_ptr_,
						   c2c_q_m_,
						   c2c_q_n_,
						   c2c_q_v_,
						   c2c_q_ld_,
									  
						   s_q2n_ptr_,
						   s_q2n_m_,
						   s_q2n_n_,
						   s_q2n_v_,
						   s_q2n_ld_,
									  
						   num_boundary_quadrilaterals_,
						   work_n_,
						   work_);
	WMESH_STATUS_CHECK(status);
      }
    //    std::cout << "count " << count << std::endl;
    return WMESH_STATUS_SUCCESS;
  }
};
