#include <limits>
#include <iostream>
#include <array>
#include "wmesh_t.hpp"

#include "GenericEncoding.hpp"

extern "C"
{
  //
  //  x-----x-----x
  //  |     |     |
  //  x-----x-----x
  //  |     |     |
  //  x-----x-----x
  //
  //

  wmesh_status_t bms_n2c_cindex(wmesh_int_t 	c_,
				wmesh_int_p 	cindex_)
  {
    cindex_[0] = GenericEncoding<wmesh_int_t,2>::Up(c_) - 1;
    return WMESH_STATUS_SUCCESS;
  }

  wmesh_status_t bms_n2c_ctype(wmesh_int_t 	c_,
			       wmesh_int_p 	ctype_)
  {
    ctype_[0] = GenericEncoding<wmesh_int_t,2>::Low(c_);
    return WMESH_STATUS_SUCCESS;
  }
  
  wmesh_status_t bms_n2c(wmesh_int_t 		c2n_size_,
			 const_wmesh_int_p 	c2n_ptr_,
			 const_wmesh_int_p 	c2n_m_,
			 const_wmesh_int_p 	c2n_n_,
			 const_wmesh_int_p 	c2n_v_,
			 const_wmesh_int_p 	c2n_ld_,
			 
			 wmesh_int_p 		n2c_ptr_,
			 wmesh_int_t 		n2c_m_,
			 wmesh_int_p 		n2c_v_)
  {
    //
    // Blank n2c_ptr_.
    //
    for (wmesh_int_t i=0;i<=n2c_m_;++i)
      {
	n2c_ptr_[i] = 0;
      }

    //
    // Blank n2c_ptr_.
    //    
    for (wmesh_int_t cell_type=0;cell_type<c2n_size_;++cell_type)
      {
	for (wmesh_int_t j=0;j<c2n_n_[cell_type];++j)
	  {
	    for (wmesh_int_t i=0;i<c2n_m_[cell_type];++i)
	      {
		wmesh_int_t node_index = c2n_v_[c2n_ptr_[cell_type] + j * c2n_ld_[cell_type] + i] -1;
		n2c_ptr_[node_index + 1] += 1;
	      }	  
	  }
      }

    //
    // Set n2c_ptr_.
    //
    for (wmesh_int_t i=2;i<=n2c_m_;++i)
      {
	n2c_ptr_[i] += n2c_ptr_[i-1];
      }

    //
    // Blank n2c_ptr_.
    //    
    for (wmesh_int_t cell_type=0;cell_type<c2n_size_;++cell_type)
      {
	for (wmesh_int_t j=0;j<c2n_n_[cell_type];++j)
	  {
	    for (wmesh_int_t i=0;i<c2n_m_[cell_type];++i)
	      {
		wmesh_int_t node_index = c2n_v_[c2n_ptr_[cell_type] + j * c2n_ld_[cell_type] + i] -1;
		n2c_v_[n2c_ptr_[node_index]++] = GenericEncoding<wmesh_int_t,2>::Encod(j + 1, cell_type);
	      }	  
	  }
      }

    //
    // Set n2c_ptr_.
    //
    for (wmesh_int_t i = n2c_m_ - 1;i > 0;--i)
      {
	n2c_ptr_[i] = n2c_ptr_[i-1];
      }
    n2c_ptr_[0] = 0;
    
    return WMESH_STATUS_SUCCESS;
  };


  
#if 0
  
#if 0
  wmesh_status_t  wmesh_indexing_edges_buffer_size(wmesh_int_t		c2n_size_,
						   const_wmesh_int_p	c2n_n_,
						   wmesh_int_p 		work_n_)
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
  };
#endif  
  static wmesh_status_t wmesh_n2c_hash(wmesh_int_t 		c2n_size_,
				       const_wmesh_int_p 	c2n_ptr_,
				       const_wmesh_int_p 	c2n_m_,
				       const_wmesh_int_p 	c2n_n_,
				       const_wmesh_int_p 	c2n_v_,
				       const_wmesh_int_p 	c2n_ld_,
				       wmesh_int_p 		n2c_v_,
			  	       wmesh_int_t		work_n_,
				       wmesh_int_p		work_)
  {  
    
    const wmesh_int_t 	hash_size = work_[0];
    wmesh_int_t * 	hash_link = work_ + 1;

    //
    // Compute the 
    //
    
    {
      wmesh_int_t k = 0;
      for (wmesh_int_t cell_type=0;cell_type<c2n_size_;++cell_type)
	{
	  for (wmesh_int_t j=0;j<c2n_n_[cell_type];++j)
	    {
	      for (wmesh_int_t i=0;i<c2n_m_[cell_type];++i)
		{
		  auto v = c2n_v_[c2n_ptr_[cell_type] + j * c2n_ld_[cell_type] + i] - 1;
		  
		  auto h = ( (v < 0) ? -v : v) % hash_size; 
		  
		  n2c_v_[k] = hash_link[h];
		  
		  //
		  // Update the last half triangle index with respect to the hash value.
		  //
		  hash_link[h] = - GenericEncoding<wmesh_int_t,2>::Encod(j * c2n_ld_[cell_type] + i + 1, cell_type);
		  ++k;
		}	  
	    }
	}
    }
    
    return WMESH_STATUS_SUCCESS;
  };
  
  wmesh_status_t bms_n2c2(wmesh_int_t 		c2n_size_,
			  const_wmesh_int_p 	c2n_ptr_,
			  const_wmesh_int_p 	c2n_m_,
			  const_wmesh_int_p 	c2n_n_,
			  const_wmesh_int_p 	c2n_v_,
			  const_wmesh_int_p 	c2n_ld_,
			  
			  wmesh_int_p 		n2c_ptr_,
			  wmesh_int_t 		n2c_m_,
			  wmesh_int_p 		n2c_v_)
  {

    static constexpr const wmesh_int_t default_hash = - std::numeric_limits<wmesh_int_t>::max();    
    wmesh_int_t work_n 	= n2c_m_ / 2;
    wmesh_int_p work 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t) * work_n);
    for (wmesh_int_t i=0;i<n2c_m_;++i)
      {
	work[i] = default_hash;
      }

    wmesh_status_t status = wmesh_n2c_hash(c2n_size_,
					   c2n_ptr_,
					   c2n_m_,
					   c2n_n_,
					   c2n_v_,
					   c2n_ld_,
					   work_n,
					   work);

    free(work);
    WMESH_STATUS_CHECK(status);
    
    //
    // Blank n2c_ptr_.
    //
    for (wmesh_int_t i=0;i<=n2c_m_;++i)
      {
	n2c_ptr_[i] = 0;
      }

    //
    // Blank n2c_ptr_.
    //    
    for (wmesh_int_t cell_type=0;cell_type<c2n_size_;++cell_type)
      {
	for (wmesh_int_t j=0;j<c2n_n_[cell_type];++j)
	  {
	    for (wmesh_int_t i=0;i<c2n_m_[cell_type];++i)
	      {
		wmesh_int_t node_index = c2n_v_[c2n_ptr_[cell_type] + j * c2n_ld_[cell_type] + i] -1;
		n2c_ptr_[node_index + 1] += 1;
	      }	  
	  }
      }

    //
    // Set n2c_ptr_.
    //
    for (wmesh_int_t i=2;i<=n2c_m_;++i)
      {
	n2c_ptr_[i] += n2c_ptr_[i-1];
      }

    //
    // Blank n2c_ptr_.
    //    
    for (wmesh_int_t cell_type=0;cell_type<c2n_size_;++cell_type)
      {
	for (wmesh_int_t j=0;j<c2n_n_[cell_type];++j)
	  {
	    for (wmesh_int_t i=0;i<c2n_m_[cell_type];++i)
	      {
		wmesh_int_t node_index = c2n_v_[c2n_ptr_[cell_type] + j * c2n_ld_[cell_type] + i] -1;
		n2c_v_[n2c_ptr_[node_index]++] = GenericEncoding<wmesh_int_t,2>::Encod(j,cell_type);
	      }	  
	  }
      }

    //
    // Set n2c_ptr_.
    //
    for (wmesh_int_t i = n2c_m_ - 1;i > 0;--i)
      {
	n2c_ptr_[i] = n2c_ptr_[i-1];
      }
    n2c_ptr_[0] = 0;
    
    return WMESH_STATUS_SUCCESS;
  };
#endif
};
