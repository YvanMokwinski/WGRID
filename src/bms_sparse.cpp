
#include "bms.h"
#include "wmesh-utils.hpp"
#include <iostream>
#include <algorithm>


template<typename T>
wmesh_status_t bms_sparse_add(wmesh_int_t 		idofs_n_,
			      const_wmesh_int_p 	idofs_,
			      wmesh_int_t 		idofs_inc_,

			      wmesh_int_t 		jdofs_n_,
			      const_wmesh_int_p 	jdofs_,
			      wmesh_int_t 		jdofs_inc_,

			      const T * __restrict__ 	lmat_,
			      wmesh_int_t 		lmat_ld_,
			      
			      wmesh_int_t 		csr_size_,
			      const_wmesh_int_p 	csr_ptr_,
			      const_wmesh_int_p 	csr_ind_,
			      T * __restrict__		csr_val_)
{
  WMESH_CHECK_POINTER(idofs_);
  WMESH_CHECK_POINTER(jdofs_);
  WMESH_CHECK_POINTER(lmat_);
  WMESH_CHECK(idofs_n_ <= lmat_ld_);
  WMESH_CHECK_POINTER(csr_ptr_);
  WMESH_CHECK_POINTER(csr_ind_);
  WMESH_CHECK_POINTER(csr_val_);
  //
  // Assembly.
  //
  for (wmesh_int_t i=0;i<idofs_n_;++i)
    {
      const wmesh_int_t idof = idofs_[i * idofs_inc_] - 1;
      for (wmesh_int_t j=0;j<jdofs_n_;++j)
	{
	  const wmesh_int_t jdof = jdofs_[j * jdofs_inc_] - 1;
	  bool found = false;
	  const wmesh_int_t bound = csr_ptr_[idof+1];
	  for (wmesh_int_t at = csr_ptr_[idof];at < bound;++at)
	    {
#ifndef NDEBUG
	      if (csr_ind_[at]==0)
		{
		  std::cerr << "pblm aaaaaaaassembly " << std::endl;
		  exit(1);
		  
		}
#endif
	      if (csr_ind_[at] - 1 == jdof)
		{
		  csr_val_[at] += lmat_[j*lmat_ld_ + i];
		  found = true;
		  break;
		}
	    }
	  if (!found)
	    {
	      std::cout << "Problem assembly (" << idof << "," << jdof << ")" << std::endl;
	      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
	    }
	}	  
    }
  return WMESH_STATUS_SUCCESS;
}

template
wmesh_status_t bms_sparse_add(wmesh_int_t 		idofs_n_,
			      const_wmesh_int_p 	idofs_,
			      wmesh_int_t 		idofs_inc_,

			      wmesh_int_t 		jdofs_n_,
			      const_wmesh_int_p 	jdofs_,
			      wmesh_int_t 		jdofs_inc_,

			      const float * __restrict__ 	lmat_,
			      wmesh_int_t 		lmat_ld_,
			      
			      wmesh_int_t 		csr_size_,
			      const_wmesh_int_p 	csr_ptr_,
			      const_wmesh_int_p 	csr_ind_,
			      float * __restrict__		csr_val_);

template
wmesh_status_t bms_sparse_add(wmesh_int_t 		idofs_n_,
			      const_wmesh_int_p 	idofs_,
			      wmesh_int_t 		idofs_inc_,

			      wmesh_int_t 		jdofs_n_,
			      const_wmesh_int_p 	jdofs_,
			      wmesh_int_t 		jdofs_inc_,

			      const double * __restrict__ 	lmat_,
			      wmesh_int_t 		lmat_ld_,
			      
			      wmesh_int_t 		csr_size_,
			      const_wmesh_int_p 	csr_ptr_,
			      const_wmesh_int_p 	csr_ind_,
			      double * __restrict__		csr_val_);



static inline void insert(wmesh_int_t 	jdof_,
			  wmesh_int_t& 	select_n_,
			  wmesh_int_p 	select_,
			  wmesh_int_p 	blank_)
{
  if (select_n_  > 0)
    {
      if (jdof_ > select_[select_n_-1])
	{
	  select_[select_n_] = jdof_;
	  blank_[jdof_] = ++select_n_;
	}
      else
	{
	  wmesh_int_t i;
	  for (i=0;i<select_n_;++i)
	    {
	      if (jdof_ < select_[i])
		{
		  for (wmesh_int_t j = select_n_;j>i;--j)
		    {
		      select_[j] = select_[j-1];
		      blank_[select_[j]] = j+1;
		    }
		  select_[i] = jdof_;
		  blank_[jdof_] = i+1;
		  break;
		}
	    }
	  ++select_n_;
	}
    }
  else
    {
      select_[select_n_] = jdof_;
      blank_[jdof_] = ++select_n_;
    }
}

extern "C"
{

  wmesh_status_t bms_sparse_buffer_size	(wmesh_int_t 		num_dofs_,
					 wmesh_int_t 		c2d_size_,
					 const_wmesh_int_p	c2d_m_,
					 const_wmesh_int_p	c2d_n_,
					 wmesh_int_p		iw_n_,
					 wmesh_int_p		num_table_coeffs_)
  {
    WMESH_CHECK_POINTER(c2d_m_);
    WMESH_CHECK_POINTER(c2d_n_);
    WMESH_CHECK_POINTER(iw_n_);
    WMESH_CHECK_POINTER(num_table_coeffs_);
    wmesh_int_t num_table_coeffs 	= 0;    
    for (wmesh_int_t i=0;i<c2d_size_;++i)
      {
	num_table_coeffs += c2d_m_[i] * c2d_n_[i];
      }
    iw_n_[0] = ( ((num_dofs_ + 1) + num_table_coeffs) + (2 * num_dofs_) );
    num_table_coeffs_[0] = num_table_coeffs;
    return WMESH_STATUS_SUCCESS;
  }
  

  wmesh_status_t bms_sparse_ptr	(wmesh_int_t 		num_dofs_,
				 
				 wmesh_int_t 		c2d_size_,
				 const_wmesh_int_p	c2d_ptr_,
				 const_wmesh_int_p	c2d_m_,
				 const_wmesh_int_p	c2d_n_,
				 const_wmesh_int_p	c2d_v_,
				 const_wmesh_int_p	c2d_ld_,
				 
				 wmesh_int_p 		csr_ptr_,				 
				 wmesh_int_t		iw_n_,
				 wmesh_int_p		iw_)
  {
    wmesh_status_t 	status;

    wmesh_int_t required_iw_n, num_table_coeffs;

    status = bms_sparse_buffer_size(num_dofs_,
				    c2d_size_,
				    c2d_m_,
				    c2d_n_,
				    &required_iw_n,
				    &num_table_coeffs);    
    WMESH_STATUS_CHECK(status);
    if (iw_n_ < required_iw_n)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
      }
    
    wmesh_int_t 	n2d_m 	= num_dofs_;
    wmesh_int_p 	n2d_ptr = iw_;
    wmesh_int_p 	n2d_v 	= n2d_ptr + (num_dofs_+1);
    wmesh_int_p 	blank 	= n2d_v + num_table_coeffs;
    wmesh_int_p 	select 	= blank + num_dofs_;
    
    for (wmesh_int_t i=0;i<num_dofs_;++i)
      {
	blank[i] = 0;
      }

    std::cout << "compute dofs-to-cells"  << std::endl;

    status = bms_n2c(c2d_size_,
		     c2d_ptr_,
		     c2d_m_,
		     c2d_n_,
		     c2d_v_,
		     c2d_ld_,		     
		     n2d_ptr,
		     n2d_m,
		     n2d_v);
    
    std::cout << "compute dofs-to-cells done."  << std::endl;
    
    WMESH_STATUS_CHECK(status);

    std::cout << "analyze dofs-to-cells"  << std::endl;
    
    csr_ptr_[0] = 0;    
    for (wmesh_int_t idof = 0;idof < n2d_m;++idof)
      {
	wmesh_int_t select_n = 0;
	//
	// Union 
	//
	for (wmesh_int_t s = n2d_ptr[idof];s<n2d_ptr[idof+1];++s)
	  {
	    wmesh_int_t c = n2d_v[s];
	    wmesh_int_t cindex, ctype;
	    
	    status = bms_n2c_cindex(c,
				    &cindex);
	    WMESH_STATUS_CHECK(status);
	    
	    status = bms_n2c_ctype	(c,
					 &ctype);
	    WMESH_STATUS_CHECK(status);
	    
	    for (wmesh_int_t i=0;i<c2d_m_[ctype];++i)
	      {
		wmesh_int_t jdof = c2d_v_[c2d_ptr_[ctype] + cindex * c2d_ld_[ctype] + i] - 1;
		if (0 == blank[jdof])
		  {
		    select[select_n] = jdof;
		    blank[jdof] = ++select_n;
		  }
	      }
	  }
	
	//
	// Reset.
	//
	for (wmesh_int_t s = 0;s<select_n;++s)
	  {
	    blank[select[s]] = 0;
	  }
	
	csr_ptr_[idof+1] = csr_ptr_[idof] + select_n;	
      }
    std::cout << "analyze done"  << std::endl;
    
    return WMESH_STATUS_SUCCESS;
  }


  
  
  wmesh_status_t bms_sparse	(wmesh_int_t 		num_dofs_,

				 wmesh_int_t 		c2d_size_,
				 const_wmesh_int_p	c2d_ptr_,
				 const_wmesh_int_p	c2d_m_,
				 const_wmesh_int_p	c2d_n_,
				 const_wmesh_int_p	c2d_v_,
				 const_wmesh_int_p	c2d_ld_,
				 
				 const_wmesh_int_p 	csr_ptr_,
				 wmesh_int_p 		csr_ind_,
				 wmesh_int_t		iw_n_,
				 wmesh_int_p		iw_)
  {
    wmesh_status_t 	status;
    wmesh_int_t
      required_iw_n,
      num_table_coeffs;
    status = bms_sparse_buffer_size(num_dofs_,
				    c2d_size_,
				    c2d_m_,
				    c2d_n_,
				    &required_iw_n,
				    &num_table_coeffs);    
    WMESH_STATUS_CHECK(status);
    if (iw_n_ < required_iw_n)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);
      }
    
    wmesh_int_t 	n2d_m 	= num_dofs_;
    wmesh_int_p 	n2d_ptr = iw_;
    wmesh_int_p 	n2d_v 	= n2d_ptr + (num_dofs_+1);
    wmesh_int_p 	blank 	= n2d_v + num_table_coeffs;
    wmesh_int_p 	select 	= blank + num_dofs_;
    
    for (wmesh_int_t idof=0;idof<n2d_m;++idof)
      {
	wmesh_int_t select_n = 0;
	//
	// Union 
	//
	for (wmesh_int_t s = n2d_ptr[idof];s<n2d_ptr[idof+1];++s)
	  {
	    wmesh_int_t c = n2d_v[s];
	    wmesh_int_t cindex, ctype;
	    
	    status = bms_n2c_cindex	(c,
					 &cindex);
	    WMESH_STATUS_CHECK(status);
	    
	    status = bms_n2c_ctype	(c,
					 &ctype);
	    WMESH_STATUS_CHECK(status);
	    for (wmesh_int_t i=0;i<c2d_m_[ctype];++i)
	      {
		wmesh_int_t jdof = c2d_v_[c2d_ptr_[ctype] + cindex * c2d_ld_[ctype] + i] - 1;	       
		if (0==blank[jdof])
		  {
#if 0
		    insert(jdof,
			   select_n,
			   select,
			   blank);
#endif
		    		    select[select_n] = jdof;
		    		    blank[jdof] = ++select_n;
		  }
	      }
	  }

	//
	// Reset.
	//
	for (wmesh_int_t s = 0;s<select_n;++s)
	  {
	    blank[select[s]] = 0;
	  }
	
	//
	// Sort.
	//

	std::sort(select,
		  select + select_n,
		  [](const wmesh_int_t& a,const wmesh_int_t& b){return a < b;});

	//
	// Copy back.
	//
	for (wmesh_int_t s = 0;s<select_n;++s)
	  {
	    csr_ind_[csr_ptr_[idof] + s] = select[s]+1;
	  }
      }

    return WMESH_STATUS_SUCCESS;
  }


#if 0
  struct wmesh_pairint_t {wmesh_int_t a;wmesh_int_t b;};

  wmesh_status_t bms_sparse_dg	(wmesh_int_t 		num_dofs_,
				 
				 wmesh_int_t 		c2c_size_,
				 const_wmesh_int_p	c2c_ptr_,
				 const_wmesh_int_p	c2c_m_,
				 const_wmesh_int_p	c2c_n_,
				 const_wmesh_int_p	c2c_v_,
				 const_wmesh_int_p	c2c_ld_,
				 
				 const_wmesh_int_p	row_size_blocks_,
				 const_wmesh_int_p	col_size_blocks_,
				 
				 wmesh_int_p 		csr_ptr__,
				 wmesh_int_p 		csr_ind__,
				 wmesh_int_t		iw_n_,
				 wmesh_int_p		iw_)
{
  
  WMESH_POINTER_CHECK(c2n_ptr_);
  WMESH_POINTER_CHECK(c2n_m_);
  WMESH_POINTER_CHECK(c2n_n_);
  WMESH_POINTER_CHECK(c2n_v_);
  WMESH_POINTER_CHECK(c2n_ld_);

  WMESH_POINTER_CHECK(c2c_ptr_);
  WMESH_POINTER_CHECK(c2c_m_);
  WMESH_POINTER_CHECK(c2c_n_);
  WMESH_POINTER_CHECK(c2c_v_);
  WMESH_POINTER_CHECK(c2c_ld_);

  WMESH_POINTER_CHECK(row_size_blocks_);
  WMESH_POINTER_CHECK(col_size_blocks_);

  WMESH_POINTER_CHECK(csr_m__);
  WMESH_POINTER_CHECK(csr_n__);
  WMESH_POINTER_CHECK(csr_nnz__);
  WMESH_POINTER_CHECK(csr_row_ptr__);
  WMESH_POINTER_CHECK(csr_row_ind__);


  //
  // Compute the number of rows and columns.
  //
  wmesh_int_t num_rows 		= 0;
  wmesh_int_t num_columns 	= 0;
  {
    for (wmesh_int_t cell_type = 0;cell_type < c2n_size_;++cell_type)
      {
	num_rows += c2c_n_[cell_type] * row_size_blocks_[cell_type];
	num_columns += c2c_n_[cell_type] * col_size_blocks_[cell_type];
      }
  }

  //
  // Allocate array csr_row_ptr.
  //
  csr_ptr__[0] = (wmesh_int_p)malloc(num_rows + 1,sizeof(wmesh_int_t) );
  if (!csr_ptr__[0])
    {
      WMESH_CHECK_STATUS(WMESH_CHECK_ERROR_MEMORY);
    }
  wmesh_int_p csr_ptr_ = csr_ptr__[0];
  
  //
  // Allocate array csr_ind.
  //
  wmesh_int_t global_shifts[5];
  global_shifts[0] = 0;
  for (wmesh_int_t cell_type = 0;cell_type < c2c_size_;++cell_type)
    {
      global_shifts[cell_type+1] = global_shifts[cell_type] + c2c_n_[cell_type];
    }

  wmesh_int_t global_dof_shifts[5];
  global_dof_shifts[0] = 0;
  for (wmesh_int_t cell_type = 0;cell_type < c2c_size_;++cell_type)
    {
      global_dof_shifts[cell_type+1] = global_dof_shifts[cell_type] + c2c_n_[cell_type] * row_size_blocks_[cell_type];
    }

  wmesh_int_t nnz = 0;  
  //
  // Get the number of interior faces.
  //
  csr_ind__[0] = (wmesh_int_p)malloc(nnz * sizeof(wmesh_int_t) );
  if (!csr_ind__[0])
    {
      WMESH_CHECK_STATUS(WMESH_CHECK_ERROR_MEMORY);
    }
  wmesh_int_p csr_ind_ = csr_ind__[0];
  wmesh_int_t at = 0;
  wmesh_int_t nei_ids[8*2];

  wmesh_int_t nblocks;
  wmesh_int_t blocks[8];
  wmesh_int_t row_idx = 0;
  
  wmesh_int_t acc=0;

  
  *csr_ptr_++=0;
  for (wmesh_int_t cell_type = 0;cell_type < c2c_size_;++cell_type)
    {
      wmesh_int_t row_idx_start = global_dof_shifts[cell_type];
      wmesh_int_t ib_size 	= row_size_blocks_[cell_type];      
      const wmesh_int_t c2c_n 	= c2c_n_[cell_type];
      //
      // Traverse the adjacency graph.
      //
      for (wmesh_int_t cell_idx=0;cell_idx < c2n_n;++cell_idx)
	{

	  //
	  // Diagonal
	  //
	  acc += col_size_blocks[cell_type];
	  for (wmesh_int_t lidx=0;lidx<c2c_m;++lidx)
	    {
	      wmesh_int_t info = c2c_v_[ c2c_ptr_[cell_type] + c2c_ld_[cell_type] * cell_idx + lidx];
	      if (info > 0)
		{
		  auto jtype = GenericEncoding<wmesh_int_t,2>::Low(info,tested_cell_type);
		  acc += col_size_blocks[jtype];
		}
	    }	  
	  *csr_ptr_++ = acc;
	}      
    } 
  csr_ptr_ = csr_ptr__[0];
  
  
  //
  // For each cell type.
  //
  for (wmesh_int_t cell_type = 0;cell_type < c2c_size_;++cell_type)
    {
      wmesh_int_t row_idx_start = global_dof_shifts[cell_type];
      wmesh_int_t ib_size 	= row_size_blocks_[cell_type];      
      const wmesh_int_t c2c_n 	= c2c_n_[cell_type];
      const wmesh_int_t c2c_m 	= c2c_m_[cell_type];
      
      //
      // Traverse the adjacency graph.
      //
      for (wmesh_int_t cell_idx=0;cell_idx < c2n_n;++cell_idx)
	{
	  wmesh_int_t row_idx = row_idx_start + cell_idx * ib_size;

	  //
	  // Get the global cell index.
	  //
	  wmesh_int_t global_cell_idx = global_shifts[cell_type] + cell_idx;

	  //
	  // Load nei_ids
	  //
	  nblocks = 0;
	  for (wmesh_int_t lidx=0;lidx<c2c_m;++lidx)
	    {
	      //
	      // Not exatly.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      //
	      wmesh_int_t info = c2c_v_[ c2c_ptr_[cell_type] + c2c_ld_[cell_type] * cell_idx + lidx];
	      if (info > 0)
		{
		  nei_ids[2*nblocks+0] = GenericEncoding<wmesh_int_t,2>::Up(info,tested_cell_type);
		  nei_ids[2*nblocks+1] = GenericEncoding<wmesh_int_t,2>::Low(info,tested_cell_type);
		  nblocks++;
		}
	    }
	  
	  std::sort((wmesh_pairint_t*)nei_ids,
		    ((wmesh_pairint_t*)nei_ids) + nblocks,
		    [](const wmesh_pairint_t& a,const wmesh_pairint_t& b){return a.a < b.a;});
	  
	  for (wmesh_int_t i=0;i<ib_size;++i)
	    {
	      for (wmesh_int_t li=0;li<nblocks;++li)
		{
		  if (nei_ids[2*li+0] > 0)
		    {
		      wmesh_int_t cell_j_idx  = nei_ids[2*li+0] - 1;
		      wmesh_int_t cell_j_type = nei_ids[2*li+1];
		      
		      wmesh_int_t global_cell_j_idx = global_shifts[cell_j_type] + cell_j_idx * col_size_blocks_[];
		      for (wmesh_int_t j=0;j<col_size_blocks_[cell_j_type];++j)
			{
			  *csr_ind_++ = global_cell_j_idx + j;
			}	      
		    }
		}
	    }
	}      
    }

  bsr_m__[0] = ncells_;
  bsr_n__[0] = ncells_;

}
#endif



};
