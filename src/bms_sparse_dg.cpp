#include "wmesh.h"
#include <iostream>
#include "wmesh-blas.hpp"
#include <algorithm>
#include "bms.h"

extern "C"
{
  wmesh_status_t bms_sparse_dg_nnz(const_wmesh_int_p 	size_blocks_,
				   wmesh_int_t 		c2c_size_,
				   const_wmesh_int_p 	c2c_ptr_,
				   const_wmesh_int_p 	c2c_m_,
				   const_wmesh_int_p 	c2c_n_,
				   const_wmesh_int_p 	c2c_v_,
				   const_wmesh_int_p 	c2c_ld_,
				   wmesh_int_p		num_dofs_,
				   wmesh_int_p		nnz_)
  {
    wmesh_int_t num_dofs = 0;
    wmesh_int_t nnz_diagonal = 0;
    wmesh_int_t nnz_extra_diagonal = 0;
    for (wmesh_int_t l=0;l<c2c_size_;++l)
      {
	const wmesh_int_t c2c_n  	= c2c_n_[l];
	if (c2c_n > 0)
	  {
	  
	    const wmesh_int_t isize_block	= size_blocks_[l];
	    num_dofs += isize_block * c2c_n;
	    nnz_diagonal += isize_block*isize_block * c2c_n;
	  
	    const wmesh_int_t c2c_m	= c2c_m_[l];
	    const wmesh_int_t c2c_ld	= c2c_ld_[l];
	    const_wmesh_int_p c2c_v 	= c2c_v_ + c2c_ptr_[l];  	  
	    for (wmesh_int_t j=0;j<c2c_n;++j)
	      {	      
		for (wmesh_int_t i=0;i<c2c_m;++i)
		  {
		    wmesh_int_t c2c_info = c2c_v[j*c2c_ld + i];
		    if (c2c_info != 0)
		      {
			wmesh_int_t jtype;
			// status = bms_c2c_cindex	(c2c_info,
			//				 &jelm);
			//		      WMESH_STATUS_CHECK(status);
			{
			  wmesh_status_t status;
			  status = bms_c2c_ctype	(c2c_info,
							 &jtype);
			  WMESH_STATUS_CHECK(status);
			}
		      
			const wmesh_int_t jsize_block = size_blocks_[jtype];
			nnz_extra_diagonal += jsize_block*isize_block;
		      }
		  }
	      }
	  }
      }
    num_dofs_[0] = num_dofs;
    nnz_[0] = nnz_extra_diagonal + nnz_diagonal;
    return WMESH_STATUS_SUCCESS;
  }

  struct wmesh_int_pair_t
  {
    wmesh_int_t first;
    wmesh_int_t second;
  };
  
  wmesh_status_t bms_sparse_dg(const_wmesh_int_p 	size_blocks_,
			       wmesh_int_t 		c2c_size_,
			       const_wmesh_int_p 	c2c_ptr_,
			       const_wmesh_int_p 	c2c_m_,
			       const_wmesh_int_p 	c2c_n_,
			       const_wmesh_int_p 	c2c_v_,
			       const_wmesh_int_p 	c2c_ld_,
			       wmesh_int_t		csr_n_,
			       wmesh_int_p		csr_ptr_,
			       wmesh_int_p		csr_ind_)
  {
    wmesh_int_t iw[(6+1)*2];
    wmesh_int_t idof = 0;
    wmesh_int_t shifts_dofs[4+1];

    shifts_dofs[0]=0;  
    for (wmesh_int_t l=0;l<c2c_size_;++l)
      {
	const wmesh_int_t c2c_n  	= c2c_n_[l];
	shifts_dofs[l+1] = shifts_dofs[l] + c2c_n * size_blocks_[l];
	std::cout << "SSSSSSSSSSSSSSSSs " << shifts_dofs[l] << std::endl;
      }
    
    wmesh_int_t shifts_cellidx[4+1];
    shifts_cellidx[0]=0;  
    for (wmesh_int_t l=0;l<c2c_size_;++l)
      {
	const wmesh_int_t c2c_n  	= c2c_n_[l];
	shifts_cellidx[l+1] = shifts_cellidx[l] + c2c_n;
      }
    
    csr_ptr_[0] = 0;
    for (wmesh_int_t l=0;l<c2c_size_;++l)
      {
	const wmesh_int_t c2c_n  	= c2c_n_[l];
	if (c2c_n > 0)
	  {
	  
	    const wmesh_int_t isize_block	= size_blocks_[l];	  
	    const wmesh_int_t c2c_m	= c2c_m_[l];
	    const wmesh_int_t c2c_ld	= c2c_ld_[l];
	    const_wmesh_int_p c2c_v 	= c2c_v_ + c2c_ptr_[l];
	    for (wmesh_int_t jidx=0;jidx<c2c_n;++jidx)
	      {
		//
		// Load
		//
		wmesh_int_t len = 0;
		for (wmesh_int_t i=0;i<c2c_m;++i)
		  {
		    wmesh_int_t c2c_info = c2c_v[jidx*c2c_ld + i];
		    if (c2c_info != 0)
		      {
			wmesh_int_t neielmidx, neitype;
			wmesh_status_t status;
			status = bms_c2c_cindex	(c2c_info,
		      				 &neielmidx);
			WMESH_STATUS_CHECK(status);
		      
			status = bms_c2c_ctype	(c2c_info,
						 &neitype);
			WMESH_STATUS_CHECK(status);
			iw[2*len+0] = shifts_cellidx[neitype] + neielmidx;
			iw[2*len+1] = neitype;
			++len;
		      }		      
		  }
	      
		iw[2*len+0] = shifts_cellidx[l] + jidx;
		iw[2*len+1] = l;
		++len;
	      
		//
		// sort.
		//
		std::sort((wmesh_int_pair_t*)iw,
			  (wmesh_int_pair_t*)(iw + len*2),
			  [](const wmesh_int_pair_t& a,const wmesh_int_pair_t& b){return a.first < b.first;});
	      
		for (wmesh_int_t i=0;i<isize_block;++i)
		  {
		    wmesh_int_t size_row=0;
		    for (wmesh_int_t k=0;k<len;++k)
		      {
			wmesh_int_t
			  jtype = iw[2*k + 1];
			wmesh_int_t  idx = iw[2*k + 0] - shifts_cellidx[jtype];
			
			const wmesh_int_t jsize_block = size_blocks_[jtype];		      
			for (wmesh_int_t j=0;j<jsize_block;++j)
			  {
			    const wmesh_int_t jdof = shifts_dofs[jtype] + (idx * jsize_block + j);
			    csr_ind_[csr_ptr_[idof] + size_row + j] = jdof + 1;
			  }
			size_row += jsize_block;		      
		      }
		    
		    csr_ptr_[idof+1] = csr_ptr_[idof] + size_row;
		    ++idof;
		  }	      
	      }
	  }

      }
    
    return WMESH_STATUS_SUCCESS;
  }

};

