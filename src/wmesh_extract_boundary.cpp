
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_t.hpp"
#include "bms.h"

extern "C"
{
  
    
  wmesh_status_t wmesh_extract_boundary(wmesh_t*    self_)
  {
    wmesh_status_t 	status;
    wmesh_int_t    	num_faces[2];
    wmesh_int_t    	num_bfaces[2];

    
    
    //
    // 2/ Create the cells to cells.
    //

    //
    // 2.a/ Use c2f since this is the same layout.
    //
    status = wmesh_init_c2f	(&self_->m_c2n,
				 &self_->m_c2c);
    WMESH_STATUS_CHECK(status);
      
    //
    // 2.a/ Define the part through triangles.
    //
    status = wmesh_init_c2f_t	(&self_->m_c2c,
				 &self_->m_c2c_t);
    WMESH_STATUS_CHECK(status);
      
    //
    // 2.a/ Define the part through quadrilaterals.
    //
    status = wmesh_init_c2f_q	(&self_->m_c2c,
				 &self_->m_c2c_q);
    WMESH_STATUS_CHECK(status);      
    
    //
    // 3/ Create the cells to edges.
    //    
    status = wmesh_init_c2e(&self_->m_c2n,
			    &self_->m_c2e,
			    self_->m_topology_dimension);
    WMESH_STATUS_CHECK(status);

    
    //
    // 4/ Create the cells to faces.
    //    
    status = wmesh_init_c2f(&self_->m_c2n,
			    &self_->m_c2f);
    WMESH_STATUS_CHECK(status);
      
    //
    // 4.a/ Define the part of triangles.
    //      
    status = wmesh_init_c2f_t	(&self_->m_c2f,
				 &self_->m_c2f_t);
    WMESH_STATUS_CHECK(status);
      
    //
    // 4.b/ Define the part of quadrilaterals.
    //      
    status = wmesh_init_c2f_q	(&self_->m_c2f,
				 &self_->m_c2f_q);
    WMESH_STATUS_CHECK(status);      
    
    //
    // Initialize the reference shape edges to nodes.
    //

    status = wmesh_create_c2c(self_,
			      num_faces,
			      num_bfaces);
    
    self_->m_num_triangles      = num_faces[0];
    self_->m_num_quadrilaterals = num_faces[1];
    self_->m_num_bfaces[0]      = num_bfaces[0];
    self_->m_num_bfaces[1]      = num_bfaces[1];


    WMESH_STATUS_CHECK(status);
#if 0
    std::cout << "boundary faces (triangles, quadrilaterals):  "
	      << num_bfaces[0]
	      << ", "
	      << num_bfaces[1]
	      << std::endl;
    
    std::cout << "total num faces (triangles, quadrilaterals):  "
	      << num_faces[0]
	      << ", "
	      << num_faces[1] << std::endl;
    
    std::cout << "total num edges:  " << self_->m_num_edges  << std::endl;
#endif

    wmesh_init_bf2n(num_bfaces[0],
		    num_bfaces[1],
		    &self_->m_bf2n);

    //    wmesh_int_p bf2n_m_ 	= self_->m_bf2n.m_m;
    wmesh_int_p bf2n_n_ 	= self_->m_bf2n.m_n;
    wmesh_int_p bf2n_ptr_ 	= self_->m_bf2n.m_ptr;
    wmesh_int_p bf2n_v_ 	= self_->m_bf2n.m_data;
    wmesh_int_p bf2n_ld_ 	= self_->m_bf2n.m_ld;


    wmesh_int_t c2n_size_ 	= self_->m_c2n.m_size;
    wmesh_int_p c2n_m_ 		= self_->m_c2n.m_m;
    wmesh_int_p c2n_n_ 		= self_->m_c2n.m_n;
    wmesh_int_p c2n_ptr_ 	= self_->m_c2n.m_ptr;
    wmesh_int_p c2n_v_ 		= self_->m_c2n.m_data;
    wmesh_int_p c2n_ld_ 	= self_->m_c2n.m_ld;
    
    // const_wmesh_int_p c2c_t_n_ 		= self_->m_c2c_t.m_n;
    const_wmesh_int_p c2c_t_m_ 		= self_->m_c2c_t.m_m;
    const_wmesh_int_p c2c_t_ptr_ 	= self_->m_c2c_t.m_ptr;
    const_wmesh_int_p c2c_t_v_ 		= self_->m_c2c_t.m_data;
    const_wmesh_int_p c2c_t_ld_ 	= self_->m_c2c_t.m_ld;

    // const_wmesh_int_p c2c_q_n_ 		= self_->m_c2c_q.m_n;
    const_wmesh_int_p c2c_q_m_ 		= self_->m_c2c_q.m_m;
    const_wmesh_int_p c2c_q_ptr_ 	= self_->m_c2c_q.m_ptr;
    const_wmesh_int_p c2c_q_v_ 		= self_->m_c2c_q.m_data;
    const_wmesh_int_p c2c_q_ld_ 	= self_->m_c2c_q.m_ld;

    const_wmesh_int_p s_t2n_n_ 		= self_->m_s_t2n.m_n;
    const_wmesh_int_p s_t2n_m_ 		= self_->m_s_t2n.m_m;
    const_wmesh_int_p s_t2n_ptr_ 	= self_->m_s_t2n.m_ptr;
    const_wmesh_int_p s_t2n_v_ 		= self_->m_s_t2n.m_data;
    const_wmesh_int_p s_t2n_ld_ 	= self_->m_s_t2n.m_ld;
    
    const_wmesh_int_p s_q2n_n_ 		= self_->m_s_q2n.m_n;
    const_wmesh_int_p s_q2n_m_ 		= self_->m_s_q2n.m_m;
    const_wmesh_int_p s_q2n_ptr_ 	= self_->m_s_q2n.m_ptr;
    const_wmesh_int_p s_q2n_v_ 		= self_->m_s_q2n.m_data;
    const_wmesh_int_p s_q2n_ld_ 	= self_->m_s_q2n.m_ld;

    wmesh_int_t c2n[8];
    wmesh_int_t t2n[3];
    wmesh_int_t q2n[4];
    wmesh_int_t itria=0;
    wmesh_int_t iquad=0;
    for (wmesh_int_t cell_type_=0;cell_type_<c2n_size_;++cell_type_)
      {
	for (wmesh_int_t cellIndex = 0;cellIndex < c2n_n_[cell_type_];++cellIndex)
	  {
	    bool extracted = false;
	    
	    //
	    // Reverse
	    //
	    //	    std::cout << " " << c2c_t_m_[cell_type_] << std::endl;
	    for (int t_lidx = 0;t_lidx<c2c_t_m_[cell_type_];++t_lidx)
	      {
		//std::cout << c2c_t_v_[c2c_t_ptr_[cell_type_] + c2c_t_ld_[cell_type_] * cellIndex + t_lidx] << std::endl;
		if (0 == c2c_t_v_[c2c_t_ptr_[cell_type_] + c2c_t_ld_[cell_type_] * cellIndex + t_lidx])
		  {

		    if (!extracted)
		      {

			extracted = true;
			//
			// Extract nodes.
			//
			get_c2n(c2n_m_[cell_type_],
				c2n_v_ + c2n_ptr_[cell_type_],
				c2n_ld_[cell_type_],
				cellIndex,
				c2n);

		      }
		    
		    get_t2n(c2n,
			    t_lidx,
			    t2n,
			    s_t2n_m_[cell_type_],
			    s_t2n_n_[cell_type_],
			    s_t2n_v_ + s_t2n_ptr_[cell_type_],
			    s_t2n_ld_[cell_type_]);

		    //
		    // Copy t2n.
		    //
		    bf2n_v_[bf2n_ptr_[0] + bf2n_ld_[0] * itria + 0] = t2n[0];
		    bf2n_v_[bf2n_ptr_[0] + bf2n_ld_[0] * itria + 1] = t2n[1];
		    bf2n_v_[bf2n_ptr_[0] + bf2n_ld_[0] * itria + 2] = t2n[2];
		    ++itria;
		  }
	      }
	  }
      }

    for (wmesh_int_t cell_type_=0;cell_type_<c2n_size_;++cell_type_)
      {
	for (wmesh_int_t cellIndex = 0;cellIndex < c2n_n_[cell_type_];++cellIndex)
	  {
	    bool extracted = false;
	    
	    //
	    // Reverse
	    //
	    for (int t_lidx = 0;t_lidx<c2c_q_m_[cell_type_];++t_lidx)
	      {
		if (0 == c2c_q_v_[c2c_q_ptr_[cell_type_] + c2c_q_ld_[cell_type_] * cellIndex + t_lidx ])
		  {		
		    if (!extracted)
		      {
			extracted = true;
			//
			// Extract nodes.
			//
			get_c2n(c2n_m_[cell_type_],
				c2n_v_ + c2n_ptr_[cell_type_],
				c2n_ld_[cell_type_],
				cellIndex,
				c2n);			
		      }
		    
		    get_q2n(c2n,
			    t_lidx,
			    q2n,
			    s_q2n_m_[cell_type_],
			    s_q2n_n_[cell_type_],
			    s_q2n_v_ + s_q2n_ptr_[cell_type_],
			    s_q2n_ld_[cell_type_]);

		    //
		    // Copy q2n.
		    //
		    bf2n_v_[bf2n_ptr_[1] + bf2n_ld_[1] * iquad + 0] = q2n[0];
		    bf2n_v_[bf2n_ptr_[1] + bf2n_ld_[1] * iquad + 1] = q2n[1];
		    bf2n_v_[bf2n_ptr_[1] + bf2n_ld_[1] * iquad + 2] = q2n[2];
		    bf2n_v_[bf2n_ptr_[1] + bf2n_ld_[1] * iquad + 3] = q2n[3];
		    ++iquad;

		  }
	      }
	  }

	
      }

    
    {
      wmesh_int_t bf_c_size 	= self_->m_bf2n.m_size;
      wmesh_int_t bf_c_ld[4] 	= {1,1,1,1};
      wmesh_int_t bf_c_ptr[5];
      wmesh_int_t bf_c_m[4] 	= {1,1,1,1};
      wmesh_int_t bf_c_n[4] 	= {0,0,0,0};
      for (wmesh_int_t i=0;i<self_->m_bf2n.m_size;++i)
	{
	  bf_c_n[i] = bf2n_n_[i];
	}
	
      bf_c_ptr[0] = 0;
      for (wmesh_int_t i=0;i<self_->m_bf2n.m_size;++i)
	{
	  bf_c_ptr[i+1] = bf_c_ptr[i] + bf_c_n[i] * bf_c_ld[i];
	}
      wmesh_int_p bf_c_v 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*bf_c_ptr[self_->m_bf2n.m_size]);

      status = wmesh_int_sparsemat_new(&self_->m_bf_c,
				       bf_c_size,
				       bf_c_ptr,
				       bf_c_m,
				       bf_c_n,
				       bf_c_v,
				       bf_c_ld);
      WMESH_STATUS_CHECK(status);
      wmesh_int_t N = bf_c_ptr[self_->m_bf2n.m_size];
      for (wmesh_int_t i=0;i<N;++i)
	{
	  bf_c_v[i] = 1000 + i + 1;
	}
    }
    
#if 0
    wmesh_int_sparsemat_fprintf(&bf2n,
				stdout);
#endif
    return WMESH_STATUS_SUCCESS;
  }
  


};
