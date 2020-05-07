#if 0
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "bms.h"

#include <chrono>
#include <iostream>

#define WINT_SPARSE_MAT_PARAMS(_f)				\
    _f.m_ptr,							\
      _f.m_m,							\
      _f.m_n,							\
      _f.m_data,						\
      _f.m_ld

using namespace std::chrono;
extern "C"
{
  
  static wmesh_status_t wmesh_init_c2d(const wmesh_int_sparsemat_t*	c2n_,
				       wmesh_int_sparsemat_t*		c2d_,
				       const_wmesh_int_p		c2d_m_,
				       const_wmesh_int_p		c2d_ld_)
  {   
    wmesh_int_t c2d_ptr[4+1];
    wmesh_int_t c2d_n[4];      
    for (wmesh_int_t i=0;i<4;++i) c2d_n[i] = c2n_->m_n[i];
    c2d_ptr[0] = 0;
    c2d_ptr[1] = c2d_ptr[0] + c2d_n[0] * c2d_ld_[0];
    c2d_ptr[2] = c2d_ptr[1] + c2d_n[1] * c2d_ld_[1];
    c2d_ptr[3] = c2d_ptr[2] + c2d_n[2] * c2d_ld_[2];
    c2d_ptr[4] = c2d_ptr[3] + c2d_n[3] * c2d_ld_[3];
    wmesh_int_t * __restrict__ c2d_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*c2d_ptr[4]);
    wmesh_int_sparsemat_def(c2d_,
			    4);
    
    wmesh_int_sparsemat_set(c2d_,
			    c2d_m_,
			    c2d_n,
			    c2d_ld_,
			    c2d_ptr,
			    c2d_v);
    return 0;
  };


  static wmesh_status_t wmesh_init_c2d_n(wmesh_int_sparsemat_t*	c2d_,
					 wmesh_int_sparsemat_t*	c2d_n_,
					 const_wmesh_int_p 	num_ldofs_,
					 wmesh_int_p 		shifts_)
  {
    wmesh_int_t c2d_n_m[4]  {4,5,6,8};
    wmesh_int_t c2d_n_ptr[4+1];
    
    c2d_n_m[0] *= num_ldofs_[0];
    c2d_n_m[1] *= num_ldofs_[1];
    c2d_n_m[2] *= num_ldofs_[2];
    c2d_n_m[3] *= num_ldofs_[3];

    c2d_n_ptr[0] = c2d_->m_ptr[0] + shifts_[0]; 
    c2d_n_ptr[1] = c2d_->m_ptr[1] + shifts_[1]; 
    c2d_n_ptr[2] = c2d_->m_ptr[2] + shifts_[2]; 
    c2d_n_ptr[3] = c2d_->m_ptr[3] + shifts_[3];
    c2d_n_ptr[4] = c2d_->m_ptr[4];
    
    shifts_[0] += c2d_n_m[0];
    shifts_[1] += c2d_n_m[1];
    shifts_[2] += c2d_n_m[2];
    shifts_[3] += c2d_n_m[3];
    
    wmesh_int_sparsemat_def(c2d_n_,
			    4);
    wmesh_int_sparsemat_set(c2d_n_,
			    c2d_n_m,
			    c2d_->m_n,
			    c2d_->m_ld,
			    c2d_n_ptr,
			    c2d_->m_data);
    return WMESH_STATUS_SUCCESS;
  }
  
  static wmesh_status_t wmesh_init_c2d_e(wmesh_int_sparsemat_t*	c2d_,
					 wmesh_int_sparsemat_t*	c2d_e_,
					 const_wmesh_int_p 	num_ldofs_,
					 wmesh_int_p 		shifts_)
  {
    wmesh_int_t c2d_e_m[4]  {6,8,9,12};
    wmesh_int_t c2d_e_ptr[4+1];
    
    c2d_e_m[0] *= num_ldofs_[0];
    c2d_e_m[1] *= num_ldofs_[1];
    c2d_e_m[2] *= num_ldofs_[2];
    c2d_e_m[3] *= num_ldofs_[3];

    c2d_e_ptr[0] = c2d_->m_ptr[0] + shifts_[0]; 
    c2d_e_ptr[1] = c2d_->m_ptr[1] + shifts_[1]; 
    c2d_e_ptr[2] = c2d_->m_ptr[2] + shifts_[2]; 
    c2d_e_ptr[3] = c2d_->m_ptr[3] + shifts_[3];
    c2d_e_ptr[4] = c2d_->m_ptr[4];
    
    shifts_[0] += c2d_e_m[0];
    shifts_[1] += c2d_e_m[1];
    shifts_[2] += c2d_e_m[2];
    shifts_[3] += c2d_e_m[3];

    wmesh_int_sparsemat_def(c2d_e_,
			    4);
    wmesh_int_sparsemat_set(c2d_e_,
			    c2d_e_m,
			    c2d_->m_n,
			    c2d_->m_ld,
			    c2d_e_ptr,
			    c2d_->m_data);
    return WMESH_STATUS_SUCCESS;
  }
  
  static wmesh_status_t wmesh_init_c2d_t(wmesh_int_sparsemat_t*	c2d_,
					 wmesh_int_sparsemat_t*	c2d_t_,
					 const_wmesh_int_p num_ldofs_,
					 wmesh_int_p 		shifts_)
  {
    wmesh_int_t c2d_t_m[4]  {4,4,2,0};
    wmesh_int_t c2d_t_ptr[4+1];
    
    c2d_t_m[0] *= num_ldofs_[0];
    c2d_t_m[1] *= num_ldofs_[1];
    c2d_t_m[2] *= num_ldofs_[2];
    c2d_t_m[3] *= num_ldofs_[3];

    c2d_t_ptr[0] = c2d_->m_ptr[0] + shifts_[0]; 
    c2d_t_ptr[1] = c2d_->m_ptr[1] + shifts_[1]; 
    c2d_t_ptr[2] = c2d_->m_ptr[2] + shifts_[2]; 
    c2d_t_ptr[3] = c2d_->m_ptr[3] + shifts_[3];
    c2d_t_ptr[4] = c2d_->m_ptr[4];
    
    shifts_[0] += c2d_t_m[0];
    shifts_[1] += c2d_t_m[1];
    shifts_[2] += c2d_t_m[2];
    shifts_[3] += c2d_t_m[3];

    wmesh_int_sparsemat_def(c2d_t_,
			    4);
    wmesh_int_sparsemat_set(c2d_t_,
			    c2d_t_m,
			    c2d_->m_n,
			    c2d_->m_ld,
			    c2d_t_ptr,
			    c2d_->m_data);
    return WMESH_STATUS_SUCCESS;
  }

  static wmesh_status_t wmesh_init_c2d_q(wmesh_int_sparsemat_t*	c2d_,
					 wmesh_int_sparsemat_t*		c2d_q_,
					 const_wmesh_int_p num_ldofs_,
					 wmesh_int_p shifts_)
  {
    wmesh_int_t c2d_q_m[4]  {0,1,3,6};
    wmesh_int_t c2d_q_ptr[4+1];
    
    c2d_q_m[0] *= num_ldofs_[0];
    c2d_q_m[1] *= num_ldofs_[1];
    c2d_q_m[2] *= num_ldofs_[2];
    c2d_q_m[3] *= num_ldofs_[3];

    c2d_q_ptr[0] = c2d_->m_ptr[0] + shifts_[0]; 
    c2d_q_ptr[1] = c2d_->m_ptr[1] + shifts_[1]; 
    c2d_q_ptr[2] = c2d_->m_ptr[2] + shifts_[2]; 
    c2d_q_ptr[3] = c2d_->m_ptr[3] + shifts_[3];
    c2d_q_ptr[4] = c2d_->m_ptr[4];
    
    shifts_[0] += c2d_q_m[0];
    shifts_[1] += c2d_q_m[1];
    shifts_[2] += c2d_q_m[2];
    shifts_[3] += c2d_q_m[3];

    wmesh_int_sparsemat_def(c2d_q_,
			    4);
    
    wmesh_int_sparsemat_set(c2d_q_,
			    c2d_q_m,
			    c2d_->m_n,
			    c2d_->m_ld,
			    c2d_q_ptr,
			    c2d_->m_data);
    return WMESH_STATUS_SUCCESS;
  }


  
  static wmesh_status_t wmesh_init_c2d_i(wmesh_int_sparsemat_t*		c2d_,
					 wmesh_int_sparsemat_t*		c2d_i_,
					 const_wmesh_int_p 		num_ldofs_,
					 wmesh_int_p 			shifts_)
  {
    wmesh_int_t c2d_i_m[4]  {1,1,1,1};
    wmesh_int_t c2d_i_ptr[4+1];
    
    c2d_i_m[0] *= num_ldofs_[0];
    c2d_i_m[1] *= num_ldofs_[1];
    c2d_i_m[2] *= num_ldofs_[2];
    c2d_i_m[3] *= num_ldofs_[3];

    c2d_i_ptr[0] = c2d_->m_ptr[0] + shifts_[0]; 
    c2d_i_ptr[1] = c2d_->m_ptr[1] + shifts_[1]; 
    c2d_i_ptr[2] = c2d_->m_ptr[2] + shifts_[2]; 
    c2d_i_ptr[3] = c2d_->m_ptr[3] + shifts_[3];
    c2d_i_ptr[4] = c2d_->m_ptr[4];
    
    shifts_[0] += c2d_i_m[0];
    shifts_[1] += c2d_i_m[1];
    shifts_[2] += c2d_i_m[2];
    shifts_[3] += c2d_i_m[3];

    wmesh_int_sparsemat_def(c2d_i_,
			    4);
    
    wmesh_int_sparsemat_set(c2d_i_,
			    c2d_i_m,
			    c2d_->m_n,
			    c2d_->m_ld,
			    c2d_i_ptr,
			    c2d_->m_data);
    return WMESH_STATUS_SUCCESS;
  }
  
  static wmesh_status_t wmesh_init_bf2n(wmesh_int_t 			num_triangles_,
					wmesh_int_t 			num_quadrilaterals_,
					wmesh_int_sparsemat_t*		bf2n_)
  {
    wmesh_int_t bf2n_m[2] {3,4};
    wmesh_int_t bf2n_ld[2] {3,4};
    wmesh_int_t bf2n_ptr[2+1];
    wmesh_int_t bf2n_n[2];

    bf2n_n[0] = num_triangles_;
    bf2n_n[1] = num_quadrilaterals_;
    
    bf2n_ptr[0] = 0;
    bf2n_ptr[1] = bf2n_ptr[0] + bf2n_n[0] * bf2n_ld[0];
    bf2n_ptr[2] = bf2n_ptr[1] + bf2n_n[1] * bf2n_ld[1];
    wmesh_int_t * __restrict__ bf2n_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*bf2n_ptr[2]);
    wmesh_int_sparsemat_def(bf2n_,
			    2);
    
    wmesh_int_sparsemat_set(bf2n_,
			    bf2n_m,
			    bf2n_n,
			    bf2n_ld,
			    bf2n_ptr,
			    bf2n_v);
    return 0;
  }

  
  static wmesh_status_t wmesh_init_c2f(const wmesh_int_sparsemat_t*	c2n_,
				       wmesh_int_sparsemat_t*		c2f_)
  {
    wmesh_int_t c2f_m[4] {4,5,5,6};
    wmesh_int_t c2f_ld[4] {4,5,5,6};
    wmesh_int_t c2f_ptr[4+1];
    wmesh_int_t c2f_n[4];      
    for (wmesh_int_t i=0;i<4;++i) c2f_n[i] = c2n_->m_n[i];
    c2f_ptr[0] = 0;
    c2f_ptr[1] = c2f_ptr[0] + c2f_n[0] * c2f_ld[0];
    c2f_ptr[2] = c2f_ptr[1] + c2f_n[1] * c2f_ld[1];
    c2f_ptr[3] = c2f_ptr[2] + c2f_n[2] * c2f_ld[2];
    c2f_ptr[4] = c2f_ptr[3] + c2f_n[3] * c2f_ld[3];
    wmesh_int_t * __restrict__ c2f_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*c2f_ptr[4]);
    
    wmesh_int_sparsemat_def(c2f_,
			    4);
    
    wmesh_int_sparsemat_set(c2f_,
			    c2f_m,
			    c2f_n,
			    c2f_ld,
			    c2f_ptr,
			    c2f_v);
    return 0;
  }

  static wmesh_status_t wmesh_init_c2f_q(wmesh_int_sparsemat_t*	c2f_,
					 wmesh_int_sparsemat_t*		c2f_q_)
  {
    wmesh_int_t c2f_q_m[4]  {0,1,3,6};
    wmesh_int_t c2f_q_ptr[4+1];        
    c2f_q_ptr[0] = c2f_->m_ptr[0] + 4; 
    c2f_q_ptr[1] = c2f_->m_ptr[1] + 4; 
    c2f_q_ptr[2] = c2f_->m_ptr[2] + 2; 
    c2f_q_ptr[3] = c2f_->m_ptr[3] + 0;
    c2f_q_ptr[4] = c2f_->m_ptr[4];
    
    wmesh_int_sparsemat_def(c2f_q_,
			    4);
    wmesh_int_sparsemat_set(c2f_q_,
			    c2f_q_m,
			    c2f_->m_n,
			    c2f_->m_ld,			    
			    c2f_q_ptr,
			    c2f_->m_data);
    return WMESH_STATUS_SUCCESS;
  }

  static wmesh_status_t wmesh_init_c2f_t(wmesh_int_sparsemat_t*	c2f_,
				       wmesh_int_sparsemat_t*	c2f_t_)
  {
    wmesh_int_t c2f_t_m[4]  {4,4,2,0};
    wmesh_int_sparsemat_def(c2f_t_,
			    4);
    wmesh_int_sparsemat_set(c2f_t_,
			    c2f_t_m,
			    c2f_->m_n,
			    c2f_->m_ld,
			    c2f_->m_ptr,
			    c2f_->m_data);
    return WMESH_STATUS_SUCCESS;
  }

  static wmesh_status_t wmesh_init_c2e(const wmesh_int_sparsemat_t*	c2n_,
				       wmesh_int_sparsemat_t*		c2e_)
  {
    wmesh_int_t c2e_m[4] {6,8,9,12};
    wmesh_int_t c2e_ptr[4+1];
    wmesh_int_t c2e_n[4];      
    for (wmesh_int_t i=0;i<4;++i) c2e_n[i] = c2n_->m_n[i];
    c2e_ptr[0] = 0;
    c2e_ptr[1] = c2e_ptr[0] + c2e_n[0] * c2e_m[0];
    c2e_ptr[2] = c2e_ptr[1] + c2e_n[1] * c2e_m[1];
    c2e_ptr[3] = c2e_ptr[2] + c2e_n[2] * c2e_m[2];
    c2e_ptr[4] = c2e_ptr[3] + c2e_n[3] * c2e_m[3];
    wmesh_int_t * __restrict__ c2e_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*c2e_ptr[4]);
    wmesh_int_sparsemat_def(c2e_,
		       4);
    wmesh_int_sparsemat_set(c2e_,
			    c2e_m,
			    c2e_n,
			    c2e_m,
			    c2e_ptr,
			    c2e_v);
    return 0;
  }

  wmesh_status_t wmesh_shape_eval_basis(wmesh_int_t     cell_type_,
					wmesh_int_t 	coo_m_,
					wmesh_int_t 	coo_n_,
					double * 	coo_x_,
					wmesh_int_t 	coo_ld_,
					double * 	x_,
					wmesh_int_t 	x_ld_)
  {
    
    static const double s_zero = ((double)0.0);
    static const double s_one = ((double)1.0);
    static constexpr double one = 1.0;
    if (cell_type_==0)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	    
	    x_[k*x_ld_+0] = s_one - (r+s+t);
	    x_[k*x_ld_+1] = r;
	    x_[k*x_ld_+2] = s;
	    x_[k*x_ld_+3] = t;
	  }		  
      }
    else if (cell_type_==1)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	    
	    const auto
	      tscale = (t < s_one) ? s_one / (s_one - t) : s_zero,
	      lts = s_one - t - s,
	      ltr = s_one - t - r;
	    
	    x_[k*x_ld_+0] = ltr  * lts * tscale;
	    x_[k*x_ld_+1] = r * lts * tscale;
	    x_[k*x_ld_+2] = r * s * tscale;
	    x_[k*x_ld_+3] = ltr * s  * tscale;
	    x_[k*x_ld_+4] = t;
	  }		  
      }
    else if (cell_type_==2)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	   
	    const auto
	      l = one - (r+s),
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = l * lt;
	    x_[k*x_ld_+1] = r * lt;
	    x_[k*x_ld_+2] = s * lt;
	    x_[k*x_ld_+3] = l * t;
	    x_[k*x_ld_+4] = r * t;
	    x_[k*x_ld_+5] = s * t;
	  }		  
      }
    else if (cell_type_==3)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];

	    const auto
	      lr = one - r,
	      ls = one - s,
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = lr* ls * lt;
	    x_[k*x_ld_+1] = r * ls * lt;
	    x_[k*x_ld_+2] = r *  s * lt;
	    x_[k*x_ld_+3] = lr*  s * lt;
	    x_[k*x_ld_+4] = lr* ls * t;
	    x_[k*x_ld_+5] =  r* ls * t;
	    x_[k*x_ld_+6] =  r*  s * t;
	    x_[k*x_ld_+7] = lr*  s * t;
	  }
      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
    return WMESH_STATUS_SUCCESS;
  }

  

  wmesh_status_t wmesh_refine(wmesh_t*		self_,
			      wmesh_int_t	degree_,
			      wmesh_t**		refined_mesh_)
  {
    wmesh_status_t status;
#if 0
    wmeshspace_t *meshspace;
    wmeshspace_def(&meshspace,
		   degree_,
		   self_);
#endif
    
#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space ..."
	<< std::endl;
#endif
      const wmesh_int_t num_dofs_per_node 		= (degree_ > 0) ? 1 : 0;      
      const wmesh_int_t num_dofs_per_edge 		= (degree_>0) ? degree_-1 : 0;
      const wmesh_int_t num_dofs_per_triangle		= (degree_>0) ? ((degree_-1)*(degree_-2))/2 : 0;
      const wmesh_int_t num_dofs_per_quadrilateral 	= (degree_>0) ? (degree_-1)*(degree_-1) : 0;
#if 0
      const wmesh_int_t ndofs_per_vol[4] = { (degree_ > 0) ? ((degree_-1)*(degree_-2)*(degree_-3)) / 6 : 1,
					     (degree_ > 0) ? ( (degree_-2)*(degree_-1)*(2*degree_-3) ) / 6 : 1,
					     (degree_ > 0) ? (degree_-2)*(degree_-1) : 1,
					     (degree_ > 0) ? (degree_-1)*(degree_-1)*(degree_-1) : 1};    
#endif      
      
      wmesh_int_t	dof_idx = 1;

#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space nodes [ndofs="
	<< num_dofs_per_node
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      if (num_dofs_per_node > 0)
	{
	  status = wmesh_space_indexing_nodes(4,
					      WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					      WINT_SPARSE_MAT_PARAMS(self_->m_c2d_n),
					      self_->m_num_nodes,
					      num_dofs_per_node,
					      dof_idx);
	  WMESH_STATUS_CHECK(status);
	  dof_idx += self_->m_num_nodes * num_dofs_per_node;	  
	}

#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space edges [ndofs="
	<< num_dofs_per_edge
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      if (num_dofs_per_edge > 0)
	{
	  status =  wmesh_space_indexing_edges(4,
					       WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					       WINT_SPARSE_MAT_PARAMS(self_->m_c2e),
					       WINT_SPARSE_MAT_PARAMS(self_->m_c2d_e),
					       WINT_SPARSE_MAT_PARAMS(self_->m_s_e2n),
					       self_->m_num_edges,
					       num_dofs_per_edge,
					       dof_idx);
	  WMESH_STATUS_CHECK(status);
	  dof_idx += self_->m_num_edges * num_dofs_per_edge;
	}
      //wmesh_int_sparsemat_fprintf(&self_->m_c2e,stdout);
//wmesh_int_sparsemat_fprintf(&self_->m_c2d_e,stdout);
//      exit(1);
#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space triangles [ndofs="
	<< num_dofs_per_triangle
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      if (num_dofs_per_triangle > 0)
	{
	  status =  wmesh_space_indexing_triangles(4,
						   WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
						   WINT_SPARSE_MAT_PARAMS(self_->m_c2f_t),
						   WINT_SPARSE_MAT_PARAMS(self_->m_c2d_t),
						   WINT_SPARSE_MAT_PARAMS(self_->m_s_t2n),
						   self_->m_num_triangles,
						   degree_,
						   num_dofs_per_triangle,
						   dof_idx);
	  WMESH_STATUS_CHECK(status);
	  dof_idx += self_->m_num_triangles * num_dofs_per_triangle;
	}

#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space quadrilateral [ndofs="
	<< num_dofs_per_quadrilateral
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif

      if (num_dofs_per_quadrilateral > 0)
	{
	  status = wmesh_space_indexing_quadrilaterals(4,
						       WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
						       WINT_SPARSE_MAT_PARAMS(self_->m_c2f_q),
						       WINT_SPARSE_MAT_PARAMS(self_->m_c2d_q),
						       WINT_SPARSE_MAT_PARAMS(self_->m_s_q2n),
						       self_->m_num_quadrilaterals,
						       degree_,
						       num_dofs_per_quadrilateral,
						       dof_idx);
	  WMESH_STATUS_CHECK(status);
	  dof_idx += self_->m_num_quadrilaterals * num_dofs_per_quadrilateral;
	}
      
#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space interior [ndofs="
	<< "x"
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      status =  wmesh_space_indexing_interior(4,
      					      WINT_SPARSE_MAT_PARAMS(self_->m_c2d_i),
      					      &dof_idx);
      WMESH_STATUS_CHECK(status);
      
#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space done [idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      
      double * refevals[4];
      wmesh_t * refmeshes[4];
      for (wmesh_int_t l=0;l<4;++l)
	{
	  refmeshes[l]	= nullptr;
	  refevals[l]	= nullptr;
	  if (self_->m_c2n.m_n[l]>0)
	    {
#ifndef NDEBUG
	      std::cerr << "// wmesh.ndebug.verbose: wmesh_analysis, refmesh celltype_=" << l << ", degree_ " << degree_ << std::endl;
#endif
	      printf("prepare treilli #########################\n.");      	      
	      wmesh_int_t work_n;
	      wmesh_int_t ref_num_entities[WMESH_ELEMENT_ALL];
	      wmesh_int_p work;	  
	      status = wmesh_treilli_buffer_size(l,
						 degree_,
						 &work_n,
						 ref_num_entities);
	      WMESH_STATUS_CHECK(status);
	      
	      work = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*work_n);
	      if (!work)
		{
		  WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
		}
	      
	      status = wmesh_treilli(&refmeshes[l],
				     l,
				     degree_,
				     work_n,
				     work);
	      free(work);
	      WMESH_STATUS_CHECK(status);


	      //wmesh_write_medit(refmeshes[l],"ahaha.mesh");
	      
	      if (l==0)
		{
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*4);
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      b[k*4+0] = ((double)1.0)-(r+s+t);
		      b[k*4+1] = r;
		      b[k*4+2] = s;
		      b[k*4+3] = t;
		    }		  
		  refevals[l] = b;		  
		}

	      if (l==1)
		{
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*5);
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      // std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;

#if 0
		      double one = 1.0;
		      b[k*5+0] = ( one - r )* ( one - s ) * ( one - t );
		      b[k*5+1] = ( ( r )* ( one - s ) * ( one - t ) );
		      b[k*5+2] = ( ( r )* ( s ) * ( one - t ) );
		      b[k*5+3] = ( one - r )* ( s ) * ( one - t );
		      b[k*5+4] = ( one - r )* ( one - s ) * ( t ) +
			( r )* ( one - s ) * ( t ) +
			( r )* ( s ) * ( t ) +
			( one - r )* ( s ) * ( t );
#endif

#if 0
		      r = 2.0*r-1.0;
		      s = 2.0*s-1.0;
		      t = 2.0*t-1.0;
		      double hr = (1.0+r)*( (1.0-t) / 2.0 ) -1.0;
		      double hs = (1.0+s)*( (1.0-t) / 2.0 ) -1.0;
		      double ht = t;
#endif
		      if (t<1.0)
			{
			  
			  double h = 1.0-t;
			  
			  double phir0 = (h-r)/h;
			  double phir1 = r/h;
		      
			  double phis0 = (h-s)/h;
			  double phis1 = s/h;
			  
			  b[k*5+0] = phir0 * phis0 * (1.0-t);
			  b[k*5+1] = phir1 * phis0 * (1.0-t);
			  b[k*5+2] = phir1 * phis1 * (1.0-t);
			  b[k*5+3] = phir0 * phis1 * (1.0-t);
			  b[k*5+4] = t;

			  b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			  b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			  b[k*5+2] = (r * s) / (1.0-t);
			  b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			  b[k*5+4] = t;
			  
#if 0
			  r = 2.0*r-1.0;
			  s = 2.0*s-1.0;
			  t = 2.0*t-1.0;
			  double c = (1.0-t)/2.0;
			  b[k*5+0] = ( (c-t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			  b[k*5+1] = ( (c+t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			  b[k*5+2] = ( (c+t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			  b[k*5+3] = ( (c-t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			  b[k*5+4] = (t+1.0)/2.0;
			  

#else
			  b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			  b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			  b[k*5+2] = (r * s) / (1.0-t);
			  b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			  b[k*5+4] = t;
#endif
			}
		      else
			{
			  b[k*5+0] = 0.0;
			  b[k*5+1] = 0.0;
			  b[k*5+2] = 0.0;
			  b[k*5+3] = 0.0;
			  b[k*5+4] = 1.0;			  
			}
		      
		      //
		      // Split into tetrahedra.
		      //
		      // 4
		      //
		      // 3 2
		      // 0 1
		      //
#if 0
		      if (r > s)
			{
			  b[k*5+0] = (1.0-r-s-t);
			  b[k*5+1] = r;
			  b[k*5+2] = s;
			  b[k*5+3] = 0.0;
			  b[k*5+4] = t;
			}
		      else
			{
			  b[k*5+0] = (1.0-r-s-t);
			  b[k*5+1] = 0.0;
			  b[k*5+2] = r;
			  b[k*5+3] = s;
			  b[k*5+4] = t;
			}
#endif
		      
		      
#if 0
		      b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		      b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		      b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		      b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		      b[k*5+4] = (1.0+t)/2.0;
		      r=hr;
		      s=hs;
		      t=ht;
		      b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		      b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		      b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		      b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		      b[k*5+4] = (1.0+t)/2.0;
		      b[k*5+0] = (1.0-r) * (1.0-s) * (1.0-t);
		      b[k*5+1] = (r) * (1.0-s) * (1.0-t);
		      b[k*5+2] = (r) * s * (1.0-t);
		      b[k*5+3] = (1.0-r) * s * (1.0-t);
		      b[k*5+4] = t;
#endif

#if 0
		      b[k*5+5] = ( one - r -t  )* ( one - s -t ) / (4.0 * (1-t));
		      b[k*5+1] = 
		      b[k*5+2] = ( ( r )* ( s ) * ( one - t ) );
		      b[k*5+3] = ( one - r )* ( s ) * ( one - t );
		      b[k*5+4] = ( one - r )* ( one - s ) * ( t );
#endif
		    }

		  refevals[l] = b;		  
#if 0
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*6);
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      double one=1.0;
		      b[k*8+0] = ( one - r - s) * ( one - t );
		      b[k*8+1] = r * ( one - t );
		      b[k*8+2] = s * ( one - t );
		      b[k*8+4] = ( one - r - s) * ( t );
		      b[k*8+5] = r * ( t );
		    }		  
		  refevals[l] = b;
#endif
		}
	      
	      if (l==2)
		{
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*6);
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      double one=1.0;
		      b[k*6+0] = ( one - r - s) * ( one - t );
		      b[k*6+1] = r * ( one - t );
		      b[k*6+2] = s * ( one - t );
		      b[k*6+3] = ( one - r - s) * ( t );
		      b[k*6+4] = r * ( t );
		      b[k*6+5] = s * ( t );
		    }		  
		  refevals[l] = b;		  
		}

	      if (l==3)
		{
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*8);
		  //		  std::cout << " " << std::endl;
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      //		      std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;
		      double one=1.0;
		      b[k*8+0] = ( one - r )* ( one - s ) * ( one - t );
		      b[k*8+1] = ( ( r )* ( one - s ) * ( one - t ) );
		      b[k*8+2] = ( ( r )* ( s ) * ( one - t ) );
		      b[k*8+3] = ( one - r )* ( s ) * ( one - t );
		      b[k*8+4] = ( one - r )* ( one - s ) * ( t );
		      b[k*8+5] = ( ( r )* ( one - s ) * ( t ) );
		      b[k*8+6] = ( ( r )* ( s ) * ( t ) );
		      b[k*8+7] = ( one - r )* ( s ) * ( t );
		    }

		  refevals[l] = b;		  
		}
	      
	    }

	}

      wmesh_int_t ndofs = dof_idx - 1;
      self_->m_ndofs 	= ndofs;
      printf("generate dofs " WMESH_INT_FORMAT "\n",self_->m_ndofs);      
      self_->m_coo_dofs = (double*)malloc(3*sizeof(double)*self_->m_ndofs);
      
      {
	double cell_xyz[32];
	wmesh_int_t cell_ld=3;
	// double * 	c_xyz 	= (double*)malloc(sizeof(double)*2048);
	// wmesh_int_p 	c_dofs 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*2048);
	for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
	  {
	    

	    double * 	refeval = refevals[l];

	    //
	    // Local c2d.
	    //	    
	    // wmesh_int_t c2d_n	= self_->m_c2d.m_n[l];
	    wmesh_int_t c2d_m 		= self_->m_c2d.m_m[l];
	    wmesh_int_t c2d_ld 		= self_->m_c2d.m_ld[l];
	    wmesh_int_p c2d_v 		= self_->m_c2d.m_data + self_->m_c2d.m_ptr[l];

	    //
	    // Local c2n.
	    //
	    wmesh_int_t c2n_n 		= self_->m_c2n.m_n[l];
	    wmesh_int_t c2n_m 		= self_->m_c2n.m_m[l];
	    wmesh_int_t c2n_ld		= self_->m_c2n.m_ld[l];
	    wmesh_int_p c2n_v 		= self_->m_c2n.m_data + self_->m_c2n.m_ptr[l];
	    
	    for (wmesh_int_t j=0;j<c2n_n;++j)
	      {
		
		//
		// Get the coordinates of the cell.
		//		
		for (wmesh_int_t i=0;i<c2n_m;++i)
		  {
		    wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
		    for (wmesh_int_t k=0;k<3;++k)
		      {
			cell_xyz[cell_ld * i + k] = self_->m_coo[3 * idx + k];
		      }
		  }
		
#if 0
		std::cout << "global coo " << std::endl;
		for (wmesh_int_t i=0;i<c2n_m;++i)
		  {
		    wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
		    for (wmesh_int_t k=0;k<3;++k)
		      {
			std::cout <<  " " << cell_xyz[cell_ld * i + k];
		      }
		    std::cout << std::endl;
		  }
#endif
		
		//
		// Get the physical coordinates of the dofs.
		//
		//		std::cout << "c2d_m " << c2d_m << std::endl;
		for (wmesh_int_t i=0;i<c2d_m;++i)
		  {
		    double x = 0.0, y =0.0, z = 0.0;
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      x += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 0];		    
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      y += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 1];		    
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      z += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 2];		    

		    wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		    // std::cout << "idx " << idx << std::endl;
		    self_->m_coo_dofs[3 * idx + 0] = x;
		    self_->m_coo_dofs[3 * idx + 1] = y;
		    self_->m_coo_dofs[3 * idx + 2] = z;
#if 0
		    std::cout << "x " << x << std::endl;
		    std::cout << "y " << y << std::endl;
		    std::cout << "z " << z << std::endl;
		    std::cout << " " << std::endl;
#endif
		  }
		//		std::cout  << std::endl;
	      }
	  }
      }


      //      std::cout << self_->m_coo_dofs[0] << " " << self_->m_coo_dofs[1] << std::endl;
      //      exit(1);

      printf("generate codes.\n");

      wmesh_int_p dofs_cod = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*self_->m_ndofs);      
      for (wmesh_int_t i=0;i<self_->m_num_nodes;++i)
	{
	  dofs_cod[i] = self_->m_n_c.v[i];
	}
      for (wmesh_int_t i=self_->m_num_nodes;i<self_->m_ndofs;++i)
	{
	  dofs_cod[i] = 1001;
	}

      //
      // Flag boundary nodes.
      //

      //
      // Flag boundary edges.
      //
      
      //
      // Flag boundary triangles.
      //
      
      //
      // Flag boundary quadrilaterals.
      //
      
      printf("generate connectivity.\n");
      {
	
	wmesh_int_t c2n_size = 4;
	wmesh_int_t c2n_m[4]{4,5,6,8};
	wmesh_int_t c2n_ld[4]{4,5,6,8};
	wmesh_int_t c2n_n[4]{0,0,0,0};	

	for (int i=0;i<4;++i)
	  {
	    if (refmeshes[i]!=nullptr)
	      {
		if (i==1)
		  {
		    c2n_n[0] += self_->m_c2n.m_n[i] * refmeshes[i]->m_c2n.m_n[0];
		    c2n_n[1] += self_->m_c2n.m_n[i] * refmeshes[i]->m_c2n.m_n[1];
		  }
		else
		  {
		    c2n_n[i] = self_->m_c2n.m_n[i] * refmeshes[i]->m_num_cells;
		  }
	      }
	  }
	
#if 0
	for (int i=0;i<4;++i)
	  {
	    std::cout << "dddddd " << c2n_n[i] << " " << refmeshes[i]->m_num_cells << std::endl;
	  }
#endif
	wmesh_int_t c2n_ptr[5];
	c2n_ptr[0]=0;
	for (int i=0;i<4;++i)
	  {
	    //	    std::cout << c2n_ptr[i] << std::endl;
	    c2n_ptr[i+1] = c2n_ptr[i] + c2n_ld[i] * c2n_n[i];
	  }
	//	    std::cout << c2n_ptr[4] << std::endl;
	
	wmesh_int_p c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[4]);
	
	wmesh_int_t mxdofs = 0;
	for (int i=0;i<4;++i)
	  //	  if (c2n_n[i] > 0) mxdofs = (refmeshes[i]->m_num_nodes > mxdofs) ? refmeshes[i]->m_num_nodes : mxdofs;
	  if (refmeshes[i])
	  mxdofs = (refmeshes[i]->m_num_nodes > mxdofs) ? refmeshes[i]->m_num_nodes : mxdofs;
	wmesh_int_p dofs = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);
	wmesh_int_p lidx = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);

	printf("generate connectivity.\n");
	wmesh_int_t idx[4]{0,0,0,0};
	for (wmesh_int_t l=0;l<4;++l)
	  {
	    auto ref_c2n = &refmeshes[l]->m_c2n;

	    wmesh_int_t ncells = self_->m_c2d.m_n[l];
	    for (wmesh_int_t j=0;j<ncells;++j)
	      {
		
		//
		// extract dofs.
		//
		for (wmesh_int_t i=0;i<self_->m_c2d.m_m[l];++i)
		  {
		    dofs[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[l] + j * self_->m_c2d.m_ld[l] + i];
		  }

		
		//
		// For each.
		//
		for (wmesh_int_t ref_cell_type=0;ref_cell_type<ref_c2n->m_size;++ref_cell_type)
		  {
		    wmesh_int_t ref_ncells = ref_c2n->m_n[ref_cell_type];
		    wmesh_int_t ref_nnodes_per_cell = ref_c2n->m_m[ref_cell_type];
		    for (wmesh_int_t sj=0;sj<ref_ncells;++sj)
		      {
			
			for (wmesh_int_t si=0;si<ref_nnodes_per_cell;++si)
			  {
			    lidx[si] = ref_c2n->m_data[ref_c2n->m_ptr[ref_cell_type] + sj * ref_c2n->m_ld[ref_cell_type]  + si ] - 1;
			  }
			
			for (wmesh_int_t si=0;si<ref_c2n->m_m[ref_cell_type];++si)
			  {
			    c2n_v[ c2n_ptr[ref_cell_type] + c2n_ld[ref_cell_type] * idx[ref_cell_type] + si] = dofs[lidx[si]];
			  }
			++idx[ref_cell_type];		       
		      }
		  }
	      }
	  }
      printf("generate connectivity done.\n");

#if 0	
	for (wmesh_int_t l=0;l<c2n_size;++l)
	  {
	    wmesh_int_t idx=0;
	    for (wmesh_int_t j=0;j<c2n_n[l];++j)
	      {
		for (wmesh_int_t i=0;i<c2n_m[l];++i)
		  {
		    std::cout << " " <<c2n_v[c2n_ptr[l] + j * c2n_ld[l] + i];
		  }
		std::cout << std::endl;
	      }
	  }
	exit(1);
#endif
	printf("define mesh.\n");
	//	wmesh_t * space_sublinear;
	//	c2n_n[3]=0;
	status =  wmesh_def	(refined_mesh_,
				 self_->m_topology_dimension,
				 self_->m_ndofs,
				 c2n_size,
				 c2n_ptr,
				 c2n_m,
				 c2n_n,
				 c2n_v,
				 c2n_ld,
				 self_->m_coo_dofs,
				 3);
	
      for (wmesh_int_t i=0;i<self_->m_ndofs;++i)
	{
	  refined_mesh_[0]->m_n_c.v[i] = dofs_cod[i];
	}
	
	WMESH_STATUS_CHECK(status);
	//	printf("write mesh.\n");
	
	//	wmesh_write(space_sublinear,"roger.mesh");
	//	printf("write mesh done.\n");
      }
      return WMESH_STATUS_SUCCESS;
  }


  wmesh_status_t wmeshspace_sublinearmesh(wmeshspace_t * 	self_,
					  wmesh_t ** 		mesh__)
  {
    std::cout << "gggggggggggggggggg   " << self_->m_ndofs << std::endl;
    wmesh_status_t status;
    double * refevals[4];
    for (wmesh_int_t l=0;l<4;++l)
      {

	refevals[l]	= nullptr;
	if (self_->m_mesh->m_c2n.m_n[l]>0)
	  {
#ifndef NDEBUG
	    //	    std::cerr << "// wmesh.ndebug.verbose: wmesh_analysis, refmesh celltype_=" << l << ", degree_ " << degree_ << std::endl;
#endif
	    //wmesh_write_medit(self_->m_patterns[l],"ahaha.mesh");
	    
	    if (l==0)
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*4);
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    b[k*4+0] = ((double)1.0)-(r+s+t);
		    b[k*4+1] = r;
		    b[k*4+2] = s;
		    b[k*4+3] = t;
		  }		  
		refevals[l] = b;		  
	      }

	    if (l==1)
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*5);
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    // std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;

#if 0
		    double one = 1.0;
		    b[k*5+0] = ( one - r )* ( one - s ) * ( one - t );
		    b[k*5+1] = ( ( r )* ( one - s ) * ( one - t ) );
		    b[k*5+2] = ( ( r )* ( s ) * ( one - t ) );
		    b[k*5+3] = ( one - r )* ( s ) * ( one - t );
		    b[k*5+4] = ( one - r )* ( one - s ) * ( t ) +
		      ( r )* ( one - s ) * ( t ) +
		      ( r )* ( s ) * ( t ) +
		      ( one - r )* ( s ) * ( t );
#endif

#if 0
		    r = 2.0*r-1.0;
		    s = 2.0*s-1.0;
		    t = 2.0*t-1.0;
		    double hr = (1.0+r)*( (1.0-t) / 2.0 ) -1.0;
		    double hs = (1.0+s)*( (1.0-t) / 2.0 ) -1.0;
		    double ht = t;
#endif
		    if (t<1.0)
		      {
			  
			double h = 1.0-t;
			  
			double phir0 = (h-r)/h;
			double phir1 = r/h;
		      
			double phis0 = (h-s)/h;
			double phis1 = s/h;
			  
			b[k*5+0] = phir0 * phis0 * (1.0-t);
			b[k*5+1] = phir1 * phis0 * (1.0-t);
			b[k*5+2] = phir1 * phis1 * (1.0-t);
			b[k*5+3] = phir0 * phis1 * (1.0-t);
			b[k*5+4] = t;

			b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			b[k*5+2] = (r * s) / (1.0-t);
			b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			b[k*5+4] = t;
			  
#if 0
			r = 2.0*r-1.0;
			s = 2.0*s-1.0;
			t = 2.0*t-1.0;
			double c = (1.0-t)/2.0;
			b[k*5+0] = ( (c-t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+1] = ( (c+t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+2] = ( (c+t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+3] = ( (c-t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			b[k*5+4] = (t+1.0)/2.0;
			  

#else
			b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			b[k*5+2] = (r * s) / (1.0-t);
			b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			b[k*5+4] = t;
#endif
		      }
		    else
		      {
			b[k*5+0] = 0.0;
			b[k*5+1] = 0.0;
			b[k*5+2] = 0.0;
			b[k*5+3] = 0.0;
			b[k*5+4] = 1.0;			  
		      }
		      
		    //
		    // Split into tetrahedra.
		    //
		    // 4
		    //
		    // 3 2
		    // 0 1
		    //
#if 0
		    if (r > s)
		      {
			b[k*5+0] = (1.0-r-s-t);
			b[k*5+1] = r;
			b[k*5+2] = s;
			b[k*5+3] = 0.0;
			b[k*5+4] = t;
		      }
		    else
		      {
			b[k*5+0] = (1.0-r-s-t);
			b[k*5+1] = 0.0;
			b[k*5+2] = r;
			b[k*5+3] = s;
			b[k*5+4] = t;
		      }
#endif
		      
		      
#if 0
		    b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+4] = (1.0+t)/2.0;
		    r=hr;
		    s=hs;
		    t=ht;
		    b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		    b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		    b[k*5+4] = (1.0+t)/2.0;
		    b[k*5+0] = (1.0-r) * (1.0-s) * (1.0-t);
		    b[k*5+1] = (r) * (1.0-s) * (1.0-t);
		    b[k*5+2] = (r) * s * (1.0-t);
		    b[k*5+3] = (1.0-r) * s * (1.0-t);
		    b[k*5+4] = t;
#endif

#if 0
		    b[k*5+5] = ( one - r -t  )* ( one - s -t ) / (4.0 * (1-t));
		    b[k*5+1] = 
		      b[k*5+2] = ( ( r )* ( s ) * ( one - t ) );
		    b[k*5+3] = ( one - r )* ( s ) * ( one - t );
		    b[k*5+4] = ( one - r )* ( one - s ) * ( t );
#endif
		  }

		refevals[l] = b;		  
#if 0
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*6);
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    double one=1.0;
		    b[k*8+0] = ( one - r - s) * ( one - t );
		    b[k*8+1] = r * ( one - t );
		    b[k*8+2] = s * ( one - t );
		    b[k*8+4] = ( one - r - s) * ( t );
		    b[k*8+5] = r * ( t );
		  }		  
		refevals[l] = b;
#endif
	      }
	      
	    if (l==2)
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*6);
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    double one=1.0;
		    b[k*6+0] = ( one - r - s) * ( one - t );
		    b[k*6+1] = r * ( one - t );
		    b[k*6+2] = s * ( one - t );
		    b[k*6+3] = ( one - r - s) * ( t );
		    b[k*6+4] = r * ( t );
		    b[k*6+5] = s * ( t );
		  }		  
		refevals[l] = b;		  
	      }

	    if (l==3)
	      {
		//
		// Build the shape functions.
		//
		double * b = (double*)malloc(sizeof(double)*self_->m_patterns[l]->m_num_nodes*8);
		//		  std::cout << " " << std::endl;
		for (wmesh_int_t k=0;k<self_->m_patterns[l]->m_num_nodes;++k)
		  {
		    double r = self_->m_patterns[l]->m_coo[3*k+0];
		    double s = self_->m_patterns[l]->m_coo[3*k+1];
		    double t = self_->m_patterns[l]->m_coo[3*k+2];
		    //		      std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;
		    double one=1.0;
		    b[k*8+0] = ( one - r )* ( one - s ) * ( one - t );
		    b[k*8+1] = ( ( r )* ( one - s ) * ( one - t ) );
		    b[k*8+2] = ( ( r )* ( s ) * ( one - t ) );
		    b[k*8+3] = ( one - r )* ( s ) * ( one - t );
		    b[k*8+4] = ( one - r )* ( one - s ) * ( t );
		    b[k*8+5] = ( ( r )* ( one - s ) * ( t ) );
		    b[k*8+6] = ( ( r )* ( s ) * ( t ) );
		    b[k*8+7] = ( one - r )* ( s ) * ( t );
		  }

		refevals[l] = b;		  
	      }
	      
	  }

      }

#if 1
    
    double * coo_dofs = (double*)malloc(3*sizeof(double)*self_->m_ndofs);
    
    {
      double cell_xyz[32];
      wmesh_int_t cell_ld=3;
      // double * 	c_xyz 	= (double*)malloc(sizeof(double)*2048);
      // wmesh_int_p 	c_dofs 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*2048);
      for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
	{
	  double * 	refeval = refevals[l];

	    //
	    // Local c2d.
	    //	    
	    // wmesh_int_t c2d_n	= self_->m_c2d.m_n[l];
	    wmesh_int_t c2d_m 		= self_->m_c2d.m_m[l];
	    wmesh_int_t c2d_ld 		= self_->m_c2d.m_ld[l];
	    wmesh_int_p c2d_v 		= self_->m_c2d.m_data + self_->m_c2d.m_ptr[l];

	    //
	    // Local c2n.
	    //
	    wmesh_int_t c2n_n 		= self_->m_mesh->m_c2n.m_n[l];
	    wmesh_int_t c2n_m 		= self_->m_mesh->m_c2n.m_m[l];
	    wmesh_int_t c2n_ld		= self_->m_mesh->m_c2n.m_ld[l];
	    wmesh_int_p c2n_v 		= self_->m_mesh->m_c2n.m_data + self_->m_mesh->m_c2n.m_ptr[l];
	    
	    for (wmesh_int_t j=0;j<c2n_n;++j)
	      {
		
		//
		// Get the coordinates of the cell.
		//		
		for (wmesh_int_t i=0;i<c2n_m;++i)
		  {
		    wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
		    for (wmesh_int_t k=0;k<3;++k)
		      {
			cell_xyz[cell_ld * i + k] = self_->m_mesh->m_coo[3 * idx + k];
		      }
		  }
		
		
		//
		// Get the physical coordinates of the dofs.
		//
		//		std::cout << "c2d_m " << c2d_m << std::endl;
		for (wmesh_int_t i=0;i<c2d_m;++i)
		  {
		    double x = 0.0, y =0.0, z = 0.0;
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      x += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 0];		    
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      y += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 1];		    
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      z += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 2];		    

		    wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		    // std::cout << "idx " << idx << std::endl;
		    coo_dofs[3 * idx + 0] = x;
		    coo_dofs[3 * idx + 1] = y;
		    coo_dofs[3 * idx + 2] = z;
		  }
	      }
	  }
      }

      printf("generate codes.\n");

      wmesh_int_p dofs_cod = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*self_->m_ndofs);      
      for (wmesh_int_t i=0;i<self_->m_mesh->m_num_nodes;++i)
	{
	  dofs_cod[i] = self_->m_mesh->m_n_c.v[i];
	}
      for (wmesh_int_t i=self_->m_mesh->m_num_nodes;i<self_->m_ndofs;++i)
	{
	  dofs_cod[i] = 1001;
	}



      printf("generate sublinear connectivity.\n");
      {
	
	wmesh_int_t c2n_size = 4;
	wmesh_int_t c2n_m[4]{4,5,6,8};
	wmesh_int_t c2n_ld[4]{4,5,6,8};
	wmesh_int_t c2n_n[4]{0,0,0,0};	

	for (int i=0;i<4;++i)
	  {
	    if (self_->m_patterns[i]!=nullptr)
	      {
		if (i==1)
		  {
		    c2n_n[0] += self_->m_mesh->m_c2n.m_n[i] * self_->m_patterns[i]->m_c2n.m_n[0];
		    c2n_n[1] += self_->m_mesh->m_c2n.m_n[i] * self_->m_patterns[i]->m_c2n.m_n[1];
		  }
		else
		  {
		    c2n_n[i] = self_->m_mesh->m_c2n.m_n[i] * self_->m_patterns[i]->m_num_cells;
		  }
	      }
	  }
	

	wmesh_int_t c2n_ptr[5];
	c2n_ptr[0]=0;
	for (int i=0;i<4;++i)
	  {
	    //	    std::cout << c2n_ptr[i] << std::endl;
	    c2n_ptr[i+1] = c2n_ptr[i] + c2n_ld[i] * c2n_n[i];
	  }
	//	    std::cout << c2n_ptr[4] << std::endl;
	
	wmesh_int_p c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[4]);
	
	wmesh_int_t mxdofs = 0;
	for (int i=0;i<4;++i)
	  //	  if (c2n_n[i] > 0) mxdofs = (self_->m_patterns[i]->m_num_nodes > mxdofs) ? self_->m_patterns[i]->m_num_nodes : mxdofs;
	  if (self_->m_patterns[i])
	  mxdofs = (self_->m_patterns[i]->m_num_nodes > mxdofs) ? self_->m_patterns[i]->m_num_nodes : mxdofs;
	wmesh_int_p dofs = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);
	wmesh_int_p lidx = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);

	printf("generate connectivity.\n");
	wmesh_int_t idx[4]{0,0,0,0};
	for (wmesh_int_t l=0;l<4;++l)
	  {
	    auto ref_c2n = &self_->m_patterns[l]->m_c2n;

	    wmesh_int_t ncells = self_->m_c2d.m_n[l];
	    for (wmesh_int_t j=0;j<ncells;++j)
	      {
		
		//
		// extract dofs.
		//
		for (wmesh_int_t i=0;i<self_->m_c2d.m_m[l];++i)
		  {
		    dofs[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[l] + j * self_->m_c2d.m_ld[l] + i];
		  }

		
		//
		// For each.
		//
		for (wmesh_int_t ref_cell_type=0;ref_cell_type<ref_c2n->m_size;++ref_cell_type)
		  {
		    wmesh_int_t ref_ncells = ref_c2n->m_n[ref_cell_type];
		    wmesh_int_t ref_nnodes_per_cell = ref_c2n->m_m[ref_cell_type];
		    for (wmesh_int_t sj=0;sj<ref_ncells;++sj)
		      {
			
			for (wmesh_int_t si=0;si<ref_nnodes_per_cell;++si)
			  {
			    lidx[si] = ref_c2n->m_data[ref_c2n->m_ptr[ref_cell_type] + sj * ref_c2n->m_ld[ref_cell_type]  + si ] - 1;
			  }
			
			for (wmesh_int_t si=0;si<ref_c2n->m_m[ref_cell_type];++si)
			  {
			    c2n_v[ c2n_ptr[ref_cell_type] + c2n_ld[ref_cell_type] * idx[ref_cell_type] + si] = dofs[lidx[si]];
			  }
			++idx[ref_cell_type];		       
		      }
		  }
	      }
	  }
      printf("generate connectivity done.\n");

      status =  wmesh_def(mesh__,
			  self_->m_mesh->m_topology_dimension,				 
			  self_->m_ndofs,
			  c2n_size,
			  c2n_ptr,
			  c2n_m,
			  c2n_n,
			  c2n_v,
			  c2n_ld,
			  coo_dofs,
			  3);
      
      for (wmesh_int_t i=0;i<self_->m_ndofs;++i)
	{
	  mesh__[0]->m_n_c.v[i] = dofs_cod[i];
	}
      
      WMESH_STATUS_CHECK(status);
      //	printf("write mesh.\n");
      
      //	wmesh_write(space_sublinear,"roger.mesh");
      //	printf("write mesh done.\n");
      }
#endif
      return WMESH_STATUS_SUCCESS;
  }

  
  
  
  static wmesh_status_t wmeshspace_compute(wmeshspace_t * 	space_)
  {
    wmesh_t*		self_ 	= space_->m_mesh;    
    wmesh_int_t		degree_ = space_->m_degree;    
    wmesh_status_t 	status;
    const wmesh_int_t num_dofs_per_node 		= (degree_ > 0) ? 1 : 0;      
    const wmesh_int_t num_dofs_per_edge 		= (degree_>0) ? degree_-1 : 0;
    const wmesh_int_t num_dofs_per_triangle		= (degree_>0) ? ((degree_-1)*(degree_-2))/2 : 0;
    const wmesh_int_t num_dofs_per_quadrilateral 	= (degree_>0) ? (degree_-1)*(degree_-1) : 0;    
    wmesh_int_t	dof_idx = 1;
    if (num_dofs_per_node > 0)
      {
	status = wmesh_space_indexing_nodes(self_->m_c2n.m_size,
					    WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					    WINT_SPARSE_MAT_PARAMS(space_->m_c2d_n),
					    self_->m_num_nodes,
					    num_dofs_per_node,
					    dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += self_->m_num_nodes * num_dofs_per_node;	  
      }

    if (num_dofs_per_edge > 0)
      {
	status =  wmesh_space_indexing_edges(self_->m_c2n.m_size,
					     WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					     WINT_SPARSE_MAT_PARAMS(self_->m_c2e),
					     WINT_SPARSE_MAT_PARAMS(space_->m_c2d_e),
					     WINT_SPARSE_MAT_PARAMS(self_->m_s_e2n),
					     self_->m_num_edges,
					     num_dofs_per_edge,
					     dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += self_->m_num_edges * num_dofs_per_edge;
      }

    if (num_dofs_per_triangle > 0)
      {
	status =  wmesh_space_indexing_triangles(self_->m_c2n.m_size,
						 WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
						 WINT_SPARSE_MAT_PARAMS(self_->m_c2f_t),
						 WINT_SPARSE_MAT_PARAMS(space_->m_c2d_t),
						 WINT_SPARSE_MAT_PARAMS(self_->m_s_t2n),
						 self_->m_num_triangles,
						 degree_,
						 num_dofs_per_triangle,
						 dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += self_->m_num_triangles * num_dofs_per_triangle;
      }

    if (num_dofs_per_quadrilateral > 0)
      {
	status = wmesh_space_indexing_quadrilaterals(self_->m_c2n.m_size,
						     WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
						     WINT_SPARSE_MAT_PARAMS(self_->m_c2f_q),
						     WINT_SPARSE_MAT_PARAMS(space_->m_c2d_q),
						     WINT_SPARSE_MAT_PARAMS(self_->m_s_q2n),
						     self_->m_num_quadrilaterals,
						     degree_,
						     num_dofs_per_quadrilateral,
						     dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += self_->m_num_quadrilaterals * num_dofs_per_quadrilateral;
      }

    status =  wmesh_space_indexing_interior(self_->m_c2n.m_size,
					    WINT_SPARSE_MAT_PARAMS(space_->m_c2d_i),
					    &dof_idx);
    WMESH_STATUS_CHECK(status);

    space_->m_ndofs 	= dof_idx - 1;
    printf("generate dofs " WMESH_INT_FORMAT "\n",
	   space_->m_ndofs);      
    
    return WMESH_STATUS_SUCCESS;
  }
  
  wmesh_status_t wmeshspace_def(wmeshspace_t ** self__,
				wmesh_int_t 	degree_,
				wmesh_t * 	mesh_)
  {
    self__[0] = (wmeshspace_t*)calloc(1,sizeof(wmeshspace_t));
    if (!self__[0])
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    
    wmeshspace_t * self_ 	= self__[0];
    self_->m_mesh 		= mesh_;
    self_->m_degree 		= degree_;
    
    wmesh_int_t work_n;
    wmesh_int_t ref_num_entities[32];
    wmesh_int_p work;
    
    wmesh_status_t status;
    if (mesh_->m_c2n.m_size == 4)
      {
	for (wmesh_int_t l=0;l<mesh_->m_c2n.m_size;++l)
	  {
	    status = wmesh_treilli_buffer_size(l,
					       degree_,
					       &work_n,
					       ref_num_entities);
	    WMESH_STATUS_CHECK(status);
	    
	    work = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*work_n);
	    if (!work)
	      {
		WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	      }	    
	    status = wmesh_treilli(&self_->m_patterns[l],
				   l,
				   degree_,
				   work_n,
				   work);
	    free(work);	    
	    WMESH_STATUS_CHECK(status);
	  }
      }
    else if (mesh_->m_c2n.m_size == 2)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
      }
    WMESH_STATUS_CHECK(status);


    wmesh_int_t degree = degree_;
    wmesh_int_t ndofs[4]= { ((degree+1)*(degree+2)*(degree+3)) / 6,
			    ((degree+2)*(degree+1)*(2*degree+3)) / 6,
			    ((degree+1)*(degree+1)*(degree+2)) / 2,
			    (degree+1)*(degree+1)*(degree+1)};
      
      
    wmesh_int_t ndofs_n[4] = {1,1,1,1};
    wmesh_int_t ndofs_e[4] = {degree-1,degree-1,degree-1,degree-1};

    wmesh_int_t ndofs_t[4] = {(degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0};

    wmesh_int_t ndofs_q[4] = {(degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0};
            
    wmesh_int_t ndofs_i[4] = { (degree>0) ? ((degree-1)*(degree-2)*(degree-3)) / 6 : 1,
			       (degree>0) ? ( (degree-2)*(degree-1)*(2*degree-3) ) / 6 : 1,
			       (degree>0) ? ((degree-1)*(degree-1)*(degree-2)) / 2 : 1, //				 (degree>2) ? (degree-2)*(degree-1) : 0,
			       (degree>0) ? (degree-1)*(degree-1)*(degree-1) : 1 };

    //
    // cells to dofs.
    //
    {
      wmesh_status_t status;
      
      //
      // Use c2f since this is the same layout.
      //
      wmesh_int_t shifts[4] = {0,0,0,0};
      
      status = wmesh_init_c2d	(&self_->m_mesh->m_c2n,
				 &self_->m_c2d,
				 ndofs,
				 ndofs);
      WMESH_STATUS_CHECK(status);

      //
      // Node-based dofs.
      //
      status = wmesh_init_c2d_n	(&self_->m_c2d,
				 &self_->m_c2d_n,
				 ndofs_n,
				 shifts);
      WMESH_STATUS_CHECK(status);

      //
      // edge-based dofs.
      //
      status = wmesh_init_c2d_e	(&self_->m_c2d,
				 &self_->m_c2d_e,
				 ndofs_e,
				 shifts);
      WMESH_STATUS_CHECK(status);

      //
      // triangle-based dofs.
      //
      status = wmesh_init_c2d_t	(&self_->m_c2d,
				 &self_->m_c2d_t,
				 ndofs_t,
				 shifts);
      WMESH_STATUS_CHECK(status);
      
      //
      // quadrilateral-based dofs.
      //
      status = wmesh_init_c2d_q	(&self_->m_c2d,
				 &self_->m_c2d_q,
				 ndofs_q,
				 shifts);
      WMESH_STATUS_CHECK(status);      
      
      //
      // cell-based dofs.
      //
      status = wmesh_init_c2d_i(&self_->m_c2d,
				&self_->m_c2d_i,
				ndofs_i,
				shifts);
      WMESH_STATUS_CHECK(status);
    }


    std::cout << "gggggggggggggggggg   " << self_->m_ndofs << std::endl;
    status = wmeshspace_compute(self_);
    std::cout << "gggggggggggggggggg after   " << self_->m_ndofs << std::endl;
    WMESH_STATUS_CHECK(status);
    return WMESH_STATUS_SUCCESS;
  };

  wmesh_status_t wmesh_refine_calculate(wmesh_t*		self_,
					wmesh_int_t 		degree_,
					wmesh_t**		refined_mesh_)
  {
    wmesh_status_t status;

#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space ..."
	<< std::endl;
#endif
      const wmesh_int_t num_dofs_per_node 		= (degree_ > 0) ? 1 : 0;      
      const wmesh_int_t num_dofs_per_edge 		= (degree_>0) ? degree_-1 : 0;
      const wmesh_int_t num_dofs_per_triangle		= (degree_>0) ? ((degree_-1)*(degree_-2))/2 : 0;
      const wmesh_int_t num_dofs_per_quadrilateral 	= (degree_>0) ? (degree_-1)*(degree_-1) : 0;
#if 0
      const wmesh_int_t ndofs_per_vol[4] = { (degree_ > 0) ? ((degree_-1)*(degree_-2)*(degree_-3)) / 6 : 1,
					     (degree_ > 0) ? ( (degree_-2)*(degree_-1)*(2*degree_-3) ) / 6 : 1,
					     (degree_ > 0) ? (degree_-2)*(degree_-1) : 1,
					     (degree_ > 0) ? (degree_-1)*(degree_-1)*(degree_-1) : 1};    
#endif      
      
      wmesh_int_t	dof_idx = 1;

#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space nodes [ndofs="
	<< num_dofs_per_node
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      if (num_dofs_per_node > 0)
	{
	  status = wmesh_space_indexing_nodes(4,
					      WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					      WINT_SPARSE_MAT_PARAMS(self_->m_c2d_n),
					      self_->m_num_nodes,
					      num_dofs_per_node,
					      dof_idx);
	  WMESH_STATUS_CHECK(status);
	  dof_idx += self_->m_num_nodes * num_dofs_per_node;	  
	}

#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space edges [ndofs="
	<< num_dofs_per_edge
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      if (num_dofs_per_edge > 0)
	{
	  status =  wmesh_space_indexing_edges(4,
					       WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					       WINT_SPARSE_MAT_PARAMS(self_->m_c2e),
					       WINT_SPARSE_MAT_PARAMS(self_->m_c2d_e),
					       WINT_SPARSE_MAT_PARAMS(self_->m_s_e2n),
					       self_->m_num_edges,
					       num_dofs_per_edge,
					       dof_idx);
	  WMESH_STATUS_CHECK(status);
	  dof_idx += self_->m_num_edges * num_dofs_per_edge;
	}
      //wmesh_int_sparsemat_fprintf(&self_->m_c2e,stdout);
//wmesh_int_sparsemat_fprintf(&self_->m_c2d_e,stdout);
//      exit(1);
#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space triangles [ndofs="
	<< num_dofs_per_triangle
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      if (num_dofs_per_triangle > 0)
	{
	  status =  wmesh_space_indexing_triangles(4,
						   WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
						   WINT_SPARSE_MAT_PARAMS(self_->m_c2f_t),
						   WINT_SPARSE_MAT_PARAMS(self_->m_c2d_t),
						   WINT_SPARSE_MAT_PARAMS(self_->m_s_t2n),
						   self_->m_num_triangles,
						   degree_,
						   num_dofs_per_triangle,
						   dof_idx);
	  WMESH_STATUS_CHECK(status);
	  dof_idx += self_->m_num_triangles * num_dofs_per_triangle;
	}

#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space quadrilateral [ndofs="
	<< num_dofs_per_quadrilateral
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif

      if (num_dofs_per_quadrilateral > 0)
	{
	  status = wmesh_space_indexing_quadrilaterals(4,
						       WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
						       WINT_SPARSE_MAT_PARAMS(self_->m_c2f_q),
						       WINT_SPARSE_MAT_PARAMS(self_->m_c2d_q),
						       WINT_SPARSE_MAT_PARAMS(self_->m_s_q2n),
						       self_->m_num_quadrilaterals,
						       degree_,
						       num_dofs_per_quadrilateral,
						       dof_idx);
	  WMESH_STATUS_CHECK(status);
	  dof_idx += self_->m_num_quadrilaterals * num_dofs_per_quadrilateral;
	}
      
#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space interior [ndofs="
	<< "x"
	<< ", idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      status =  wmesh_space_indexing_interior(4,
      					      WINT_SPARSE_MAT_PARAMS(self_->m_c2d_i),
      					      &dof_idx);
      WMESH_STATUS_CHECK(status);
      
#ifndef NDEBUG
      std::cerr
	<< "//wmesh.verbose: wmesh_analysis_space done [idx="
	<< dof_idx
	<< "]"
	<< std::endl;
#endif
      
      double * refevals[4];
      wmesh_t * refmeshes[4];
      for (wmesh_int_t l=0;l<4;++l)
	{
	  refmeshes[l]	= nullptr;
	  refevals[l]	= nullptr;
	  if (self_->m_c2n.m_n[l]>0)
	    {
#ifndef NDEBUG
	      std::cerr << "// wmesh.ndebug.verbose: wmesh_analysis, refmesh celltype_=" << l << ", degree_ " << degree_ << std::endl;
#endif
	      printf("prepare treilli #########################\n.");      	      
	      wmesh_int_t work_n;
	      wmesh_int_t ref_num_entities[WMESH_ELEMENT_ALL];
	      wmesh_int_p work;	  
	      status = wmesh_treilli_buffer_size(l,
						 degree_,
						 &work_n,
						 ref_num_entities);
	      WMESH_STATUS_CHECK(status);
	      
	      work = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*work_n);
	      if (!work)
		{
		  WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
		}
	      
	      status = wmesh_treilli(&refmeshes[l],
				     l,
				     degree_,
				     work_n,
				     work);
	      free(work);
	      WMESH_STATUS_CHECK(status);


	      wmesh_write_medit(refmeshes[l],"ahaha.mesh");
	      
	      if (l==0)
		{
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*4);
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      b[k*4+0] = ((double)1.0)-(r+s+t);
		      b[k*4+1] = r;
		      b[k*4+2] = s;
		      b[k*4+3] = t;
		    }		  
		  refevals[l] = b;		  
		}

	      if (l==1)
		{
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*5);
		  //		  std::cout << " " << std::endl;
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      // std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;

#if 0
		      double one = 1.0;
		      b[k*5+0] = ( one - r )* ( one - s ) * ( one - t );
		      b[k*5+1] = ( ( r )* ( one - s ) * ( one - t ) );
		      b[k*5+2] = ( ( r )* ( s ) * ( one - t ) );
		      b[k*5+3] = ( one - r )* ( s ) * ( one - t );
		      b[k*5+4] = ( one - r )* ( one - s ) * ( t ) +
			( r )* ( one - s ) * ( t ) +
			( r )* ( s ) * ( t ) +
			( one - r )* ( s ) * ( t );
#endif

#if 0
		      r = 2.0*r-1.0;
		      s = 2.0*s-1.0;
		      t = 2.0*t-1.0;
		      double hr = (1.0+r)*( (1.0-t) / 2.0 ) -1.0;
		      double hs = (1.0+s)*( (1.0-t) / 2.0 ) -1.0;
		      double ht = t;
#endif
		      if (t<1.0)
			{
			  
			  double h = 1.0-t;
			  
			  double phir0 = (h-r)/h;
			  double phir1 = r/h;
		      
			  double phis0 = (h-s)/h;
			  double phis1 = s/h;
			  
			  b[k*5+0] = phir0 * phis0 * (1.0-t);
			  b[k*5+1] = phir1 * phis0 * (1.0-t);
			  b[k*5+2] = phir1 * phis1 * (1.0-t);
			  b[k*5+3] = phir0 * phis1 * (1.0-t);
			  b[k*5+4] = t;

			  b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			  b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			  b[k*5+2] = (r * s) / (1.0-t);
			  b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			  b[k*5+4] = t;
			  
#if 0
			  r = 2.0*r-1.0;
			  s = 2.0*s-1.0;
			  t = 2.0*t-1.0;
			  double c = (1.0-t)/2.0;
			  b[k*5+0] = ( (c-t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			  b[k*5+1] = ( (c+t) * (c-s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			  b[k*5+2] = ( (c+t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			  b[k*5+3] = ( (c-t) * (c+s) ) / 4.0 / c/c * (1.0-t) / (2.0);
			  b[k*5+4] = (t+1.0)/2.0;
			  

#else
			  b[k*5+0] = (1.0 - t  - r) * (1.0-t-s) / (1.0-t);
			  b[k*5+1] = (r * (1.0-t-s)) / (1.0-t);
			  b[k*5+2] = (r * s) / (1.0-t);
			  b[k*5+3] = ((1.0 - t - r) * s) / (1.0-t);
			  b[k*5+4] = t;
#endif
			}
		      else
			{
			  b[k*5+0] = 0.0;
			  b[k*5+1] = 0.0;
			  b[k*5+2] = 0.0;
			  b[k*5+3] = 0.0;
			  b[k*5+4] = 1.0;			  
			}
		      
		      //
		      // Split into tetrahedra.
		      //
		      // 4
		      //
		      // 3 2
		      // 0 1
		      //
#if 0
		      if (r > s)
			{
			  b[k*5+0] = (1.0-r-s-t);
			  b[k*5+1] = r;
			  b[k*5+2] = s;
			  b[k*5+3] = 0.0;
			  b[k*5+4] = t;
			}
		      else
			{
			  b[k*5+0] = (1.0-r-s-t);
			  b[k*5+1] = 0.0;
			  b[k*5+2] = r;
			  b[k*5+3] = s;
			  b[k*5+4] = t;
			}
#endif
		      
		      
#if 0
		      b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		      b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		      b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		      b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		      b[k*5+4] = (1.0+t)/2.0;
		      r=hr;
		      s=hs;
		      t=ht;
		      b[k*5+0] = ((r+t)*(s+t))/(2.0*(1.0-t));
		      b[k*5+1] = -((r+t)*(s+1.0))/(2.0*(1.0-t));
		      b[k*5+2] = -((r+t)*(s+t))/(2.0*(1.0-t));
		      b[k*5+3] = ((r+1.0)*(s+1.0))/(2.0*(1.0-t));
		      b[k*5+4] = (1.0+t)/2.0;
		      b[k*5+0] = (1.0-r) * (1.0-s) * (1.0-t);
		      b[k*5+1] = (r) * (1.0-s) * (1.0-t);
		      b[k*5+2] = (r) * s * (1.0-t);
		      b[k*5+3] = (1.0-r) * s * (1.0-t);
		      b[k*5+4] = t;
#endif

#if 0
		      b[k*5+5] = ( one - r -t  )* ( one - s -t ) / (4.0 * (1-t));
		      b[k*5+1] = 
		      b[k*5+2] = ( ( r )* ( s ) * ( one - t ) );
		      b[k*5+3] = ( one - r )* ( s ) * ( one - t );
		      b[k*5+4] = ( one - r )* ( one - s ) * ( t );
#endif
		    }

		  refevals[l] = b;		  
#if 0
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*6);
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      double one=1.0;
		      b[k*8+0] = ( one - r - s) * ( one - t );
		      b[k*8+1] = r * ( one - t );
		      b[k*8+2] = s * ( one - t );
		      b[k*8+4] = ( one - r - s) * ( t );
		      b[k*8+5] = r * ( t );
		    }		  
		  refevals[l] = b;
#endif
		}
	      
	      if (l==2)
		{
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*6);
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      double one=1.0;
		      b[k*6+0] = ( one - r - s) * ( one - t );
		      b[k*6+1] = r * ( one - t );
		      b[k*6+2] = s * ( one - t );
		      b[k*6+3] = ( one - r - s) * ( t );
		      b[k*6+4] = r * ( t );
		      b[k*6+5] = s * ( t );
		    }		  
		  refevals[l] = b;		  
		}

	      if (l==3)
		{
		  //
		  // Build the shape functions.
		  //
		  double * b = (double*)malloc(sizeof(double)*refmeshes[l]->m_num_nodes*8);
		  //		  std::cout << " " << std::endl;
		  for (wmesh_int_t k=0;k<refmeshes[l]->m_num_nodes;++k)
		    {
		      double r = refmeshes[l]->m_coo[3*k+0];
		      double s = refmeshes[l]->m_coo[3*k+1];
		      double t = refmeshes[l]->m_coo[3*k+2];
		      //		      std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;
		      double one=1.0;
		      b[k*8+0] = ( one - r )* ( one - s ) * ( one - t );
		      b[k*8+1] = ( ( r )* ( one - s ) * ( one - t ) );
		      b[k*8+2] = ( ( r )* ( s ) * ( one - t ) );
		      b[k*8+3] = ( one - r )* ( s ) * ( one - t );
		      b[k*8+4] = ( one - r )* ( one - s ) * ( t );
		      b[k*8+5] = ( ( r )* ( one - s ) * ( t ) );
		      b[k*8+6] = ( ( r )* ( s ) * ( t ) );
		      b[k*8+7] = ( one - r )* ( s ) * ( t );
		    }

		  refevals[l] = b;		  
		}
	      
	    }

	}

      wmesh_int_t ndofs = dof_idx - 1;
      self_->m_ndofs 	= ndofs;
      printf("generate dofs " WMESH_INT_FORMAT "\n",self_->m_ndofs);      
      self_->m_coo_dofs = (double*)malloc(3*sizeof(double)*self_->m_ndofs);
      
      {
	double cell_xyz[32];
	wmesh_int_t cell_ld=3;
	// double * 	c_xyz 	= (double*)malloc(sizeof(double)*2048);
	// wmesh_int_p 	c_dofs 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*2048);
	for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
	  {
	    

	    double * 	refeval = refevals[l];

	    //
	    // Local c2d.
	    //	    
	    // wmesh_int_t c2d_n	= self_->m_c2d.m_n[l];
	    wmesh_int_t c2d_m 		= self_->m_c2d.m_m[l];
	    wmesh_int_t c2d_ld 		= self_->m_c2d.m_ld[l];
	    wmesh_int_p c2d_v 		= self_->m_c2d.m_data + self_->m_c2d.m_ptr[l];

	    //
	    // Local c2n.
	    //
	    wmesh_int_t c2n_n 		= self_->m_c2n.m_n[l];
	    wmesh_int_t c2n_m 		= self_->m_c2n.m_m[l];
	    wmesh_int_t c2n_ld		= self_->m_c2n.m_ld[l];
	    wmesh_int_p c2n_v 		= self_->m_c2n.m_data + self_->m_c2n.m_ptr[l];
	    
	    for (wmesh_int_t j=0;j<c2n_n;++j)
	      {
		
		//
		// Get the coordinates of the cell.
		//		
		for (wmesh_int_t i=0;i<c2n_m;++i)
		  {
		    wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
		    for (wmesh_int_t k=0;k<3;++k)
		      {
			cell_xyz[cell_ld * i + k] = self_->m_coo[3 * idx + k];
		      }
		  }
		
#if 0
		std::cout << "global coo " << std::endl;
		for (wmesh_int_t i=0;i<c2n_m;++i)
		  {
		    wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
		    for (wmesh_int_t k=0;k<3;++k)
		      {
			std::cout <<  " " << cell_xyz[cell_ld * i + k];
		      }
		    std::cout << std::endl;
		  }
#endif
		
		//
		// Get the physical coordinates of the dofs.
		//
		//		std::cout << "c2d_m " << c2d_m << std::endl;
		for (wmesh_int_t i=0;i<c2d_m;++i)
		  {
		    double x = 0.0, y =0.0, z = 0.0;
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      x += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 0];		    
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      y += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 1];		    
		    for (wmesh_int_t k=0;k<c2n_m;++k)
		      z += refeval[c2n_m * i + k] * cell_xyz[cell_ld * k + 2];		    

		    wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		    // std::cout << "idx " << idx << std::endl;
		    self_->m_coo_dofs[3 * idx + 0] = x;
		    self_->m_coo_dofs[3 * idx + 1] = y;
		    self_->m_coo_dofs[3 * idx + 2] = z;
#if 0
		    std::cout << "x " << x << std::endl;
		    std::cout << "y " << y << std::endl;
		    std::cout << "z " << z << std::endl;
		    std::cout << " " << std::endl;
#endif
		  }
		//		std::cout  << std::endl;
	      }
	  }
      }


      //      std::cout << self_->m_coo_dofs[0] << " " << self_->m_coo_dofs[1] << std::endl;
      //      exit(1);

      printf("generate codes.\n");

      wmesh_int_p dofs_cod = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*self_->m_ndofs);      
      for (wmesh_int_t i=0;i<self_->m_num_nodes;++i)
	{
	  dofs_cod[i] = self_->m_n_c.v[i];
	}
      for (wmesh_int_t i=self_->m_num_nodes;i<self_->m_ndofs;++i)
	{
	  dofs_cod[i] = 1001;
	}

      //
      // Flag boundary nodes.
      //

      //
      // Flag boundary edges.
      //
      
      //
      // Flag boundary triangles.
      //
      
      //
      // Flag boundary quadrilaterals.
      //
      
      printf("generate connectivity.\n");
      {
	

	wmesh_int_t c2n_ld[4]{4,5,6,8};
	wmesh_int_t c2n_n[4]{0,0,0,0};	

	for (int i=0;i<4;++i)
	  {
	    if (refmeshes[i]!=nullptr)
	      {
		if (i==1)
		  {
		    c2n_n[0] += self_->m_c2n.m_n[i] * refmeshes[i]->m_c2n.m_n[0];
		    c2n_n[1] += self_->m_c2n.m_n[i] * refmeshes[i]->m_c2n.m_n[1];
		  }
		else
		  {
		    c2n_n[i] = self_->m_c2n.m_n[i] * refmeshes[i]->m_num_cells;
		  }
	      }
	  }
	
#if 0
	for (int i=0;i<4;++i)
	  {
	    std::cout << "dddddd " << c2n_n[i] << " " << refmeshes[i]->m_num_cells << std::endl;
	  }
#endif
	wmesh_int_t c2n_ptr[5];
	c2n_ptr[0]=0;
	for (int i=0;i<4;++i)
	  {
	    //	    std::cout << c2n_ptr[i] << std::endl;
	    c2n_ptr[i+1] = c2n_ptr[i] + c2n_ld[i] * c2n_n[i];
	  }
	//	    std::cout << c2n_ptr[4] << std::endl;
	
	wmesh_int_p c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[4]);
	
	wmesh_int_t mxdofs = 0;
	for (int i=0;i<4;++i)
	  //	  if (c2n_n[i] > 0) mxdofs = (refmeshes[i]->m_num_nodes > mxdofs) ? refmeshes[i]->m_num_nodes : mxdofs;
	  if (refmeshes[i])
	  mxdofs = (refmeshes[i]->m_num_nodes > mxdofs) ? refmeshes[i]->m_num_nodes : mxdofs;
	wmesh_int_p dofs = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);
	wmesh_int_p lidx = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);

	printf("generate connectivity.\n");
	wmesh_int_t idx[4]{0,0,0,0};
	for (wmesh_int_t l=0;l<4;++l)
	  {
	    auto ref_c2n = &refmeshes[l]->m_c2n;

	    wmesh_int_t ncells = self_->m_c2d.m_n[l];
	    for (wmesh_int_t j=0;j<ncells;++j)
	      {
		
		//
		// extract dofs.
		//
		for (wmesh_int_t i=0;i<self_->m_c2d.m_m[l];++i)
		  {
		    dofs[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[l] + j * self_->m_c2d.m_ld[l] + i];
		  }

		
		//
		// For each.
		//
		for (wmesh_int_t ref_cell_type=0;ref_cell_type<ref_c2n->m_size;++ref_cell_type)
		  {
		    wmesh_int_t ref_ncells = ref_c2n->m_n[ref_cell_type];
		    wmesh_int_t ref_nnodes_per_cell = ref_c2n->m_m[ref_cell_type];
		    for (wmesh_int_t sj=0;sj<ref_ncells;++sj)
		      {
			
			for (wmesh_int_t si=0;si<ref_nnodes_per_cell;++si)
			  {
			    lidx[si] = ref_c2n->m_data[ref_c2n->m_ptr[ref_cell_type] + sj * ref_c2n->m_ld[ref_cell_type]  + si ] - 1;
			  }
			
			for (wmesh_int_t si=0;si<ref_c2n->m_m[ref_cell_type];++si)
			  {
			    c2n_v[ c2n_ptr[ref_cell_type] + c2n_ld[ref_cell_type] * idx[ref_cell_type] + si] = dofs[lidx[si]];
			  }
			++idx[ref_cell_type];		       
		      }
		  }
	      }
	  }
      printf("generate connectivity done.\n");

#if 0	
	for (wmesh_int_t l=0;l<c2n_size;++l)
	  {
	    wmesh_int_t idx=0;
	    for (wmesh_int_t j=0;j<c2n_n[l];++j)
	      {
		for (wmesh_int_t i=0;i<c2n_m[l];++i)
		  {
		    std::cout << " " <<c2n_v[c2n_ptr[l] + j * c2n_ld[l] + i];
		  }
		std::cout << std::endl;
	      }
	  }
	exit(1);
#endif
	printf("define mesh.\n");
	//	wmesh_t * space_sublinear;
	//	c2n_n[3]=0;
	status =  wmesh_def	(refined_mesh_,
				 self_->m_topology_dimension,
				 self_->m_ndofs,
				 self_->m_c2d.m_size,
				 WINT_SPARSE_MAT_PARAMS(self_->m_c2d),
				 self_->m_coo_dofs,
				 3);
	
      for (wmesh_int_t i=0;i<self_->m_ndofs;++i)
	{
	  refined_mesh_[0]->m_n_c.v[i] = dofs_cod[i];
	}
	
	WMESH_STATUS_CHECK(status);
	//	printf("write mesh.\n");
	
	//	wmesh_write(space_sublinear,"roger.mesh");
	//	printf("write mesh done.\n");
      }
      return WMESH_STATUS_SUCCESS;
  }

  
  wmesh_status_t wmesh_analysis_edges(wmesh_t* 		self_,
				      wmesh_int_p 	num_edges_)
  {
    wmesh_int_t edge_idx = 0;
    
    wmesh_int_t work_n;
    wmesh_status_t status =  wmesh_indexing_edges_buffer_size(4,
							      self_->m_c2n.m_n,
							      &work_n);      
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_t * work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }      
    status =  wmesh_indexing_edges(4,
				   WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
				   WINT_SPARSE_MAT_PARAMS(self_->m_c2e),
				   WINT_SPARSE_MAT_PARAMS(self_->m_s_e2n),				     
				   &edge_idx,
				   work_n,
				   work);
    free(work);
    WMESH_STATUS_CHECK(status);
    num_edges_[0] = edge_idx;
    return WMESH_STATUS_SUCCESS;
  }
  

  wmesh_status_t wmesh_analysis_faces(wmesh_t* 		self_)
  {
    wmesh_status_t status;

#define WINT_SPARSE_MAT_PARAMS(_f)				\
    _f.m_ptr,							\
      _f.m_m,							\
      _f.m_n,							\
      _f.m_data,						\
      _f.m_ld
    
    wmesh_int_t work_n;
    status =  wmesh_indexing_triangles_buffer_size(self_->m_c2n.m_size,
						   self_->m_c2n.m_n,
						   &work_n);
    WMESH_STATUS_CHECK(status);
    wmesh_int_p work = (wmesh_int_p)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }      
    
    //
    // Enumerates triangles faces.
    //
    wmesh_int_t face_idx = 0;
    status =  wmesh_indexing_triangles(self_->m_c2n.m_size,
				       WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
				       WINT_SPARSE_MAT_PARAMS(self_->m_c2f_t),
				       WINT_SPARSE_MAT_PARAMS(self_->m_s_t2n),				     
				       &face_idx,
				       work_n,
				       work);      
    
    WMESH_STATUS_CHECK(status);      
    face_idx = 0;
    status = wmesh_indexing_quadrilaterals(4,
					   WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					   WINT_SPARSE_MAT_PARAMS(self_->m_c2f_q),
					   WINT_SPARSE_MAT_PARAMS(self_->m_s_q2n),				     
					   &face_idx,
					   work_n,
					   work);      
    WMESH_STATUS_CHECK(status);
    return WMESH_STATUS_SUCCESS;
  }

  
  wmesh_status_t wmesh_analysis_neighbors(wmesh_t* 		self_,
					  wmesh_int_p 		num_faces_,
					  wmesh_int_p 		num_boundary_faces_)
  {    
    
#define WINT_SPARSE_MAT_PARAMS2(_f)				\
    _f.m_size,							\
      _f.m_ptr,							\
      _f.m_m,							\
      _f.m_n,							\
      _f.m_data,						\
      _f.m_ld
    
    wmesh_int_t work_n;
    wmesh_int_t * work  = nullptr;
    wmesh_status_t status;
    
    status =  wbms_c2c_calculate_buffer_size(self_->m_c2n.m_size,
					     self_->m_c2n.m_n,
					     &work_n);    
    WMESH_STATUS_CHECK(status);    
    work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }      

    num_boundary_faces_[0] 	= 0;
    num_boundary_faces_[1] 	= 0;
    num_faces_[0] 		= 0;
    num_faces_[1] 		= 0;
    
    status = wbms_c2c_calculate(WINT_SPARSE_MAT_PARAMS2(self_->m_c2n),
				WINT_SPARSE_MAT_PARAMS2(self_->m_c2c_t),
				WINT_SPARSE_MAT_PARAMS2(self_->m_c2c_q),
				WINT_SPARSE_MAT_PARAMS2(self_->m_s_t2n),
				WINT_SPARSE_MAT_PARAMS2(self_->m_s_q2n),
				work_n,
				work,
				num_faces_,
				num_boundary_faces_);

    free(work);      
    WMESH_STATUS_CHECK(status);    
    return WMESH_STATUS_SUCCESS;
    
#undef WINT_SPARSE_MAT_PARAMS2    
  }



  
  using timing_t = high_resolution_clock::time_point;
  inline timing_t timing_stop(){ return high_resolution_clock::now(); }
  inline timing_t timing_start(){ return high_resolution_clock::now(); }
  inline double timing_seconds(timing_t&t,timing_t&t1){ return duration_cast<duration<double>>(t1-t).count(); }




  
  wmesh_status_t wmesh_analysis(wmesh_t*    self_,
				wmesh_int_t degree)
  {
    const wmesh_int_t topology_dimension = self_->m_topology_dimension;
    
    wmesh_status_t 	status;
    wmesh_int_t 	num_faces[2];
    wmesh_int_t 	num_bfaces[2];
    
#ifndef NDEBUG
    std::cerr << "wmesh_analysis: topology_dimension " << topology_dimension << std::endl;
#endif

    
#if 0
    status    = wmesh_analysis_neighbors(self_,
					 num_faces,
					 num_bfaces);
    
    self_->m_num_triangles      = num_faces[0];
    self_->m_num_quadrilaterals = num_faces[1];
    self_->m_num_bfaces[0]      = num_bfaces[0];
    self_->m_num_bfaces[1]      = num_bfaces[1];
    
    //
    //
    //
    status = wmesh_adjacencies(self_);

#endif

    //
    // 1/ Create the cells to cells.
    //
    //
    // 1.a/ Use c2f since this is the same layout.
    //
    status = wmesh_init_c2f	(&self_->m_c2n,
				 &self_->m_c2c);
    WMESH_STATUS_CHECK(status);
	
    if (topology_dimension == 3)
      {
	status = wmesh_init_c2f_t(&self_->m_c2c,
				  &self_->m_c2c_t);
	WMESH_STATUS_CHECK(status);
	
	//
	// 1.a/ Define the part through quadrilaterals.
	//
	status = wmesh_init_c2f_q	(&self_->m_c2c,
					 &self_->m_c2c_q);
	WMESH_STATUS_CHECK(status);      
      }
  


    //
    // 3/ Create the cells to edges.
    //    
    {
      wmesh_status_t status;
      status = wmesh_init_c2e(&self_->m_c2n,
			      &self_->m_c2e);
      WMESH_STATUS_CHECK(status);
    }
    
    //
    // 4/ Create the cells to faces.
    //    
    {
      wmesh_status_t status;
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
    }
    
    //
    // Initialize the reference shape edges to nodes.
    //
    {
      wmesh_status_t status;
#if 1
      {
	auto start = timing_start();      
	status    = wmesh_analysis_neighbors(self_,
					     num_faces,
					     num_bfaces);
	
	self_->m_num_triangles      = num_faces[0];
	self_->m_num_quadrilaterals = num_faces[1];
	self_->m_num_bfaces[0]      = num_bfaces[0];
	self_->m_num_bfaces[1]      = num_bfaces[1];
	auto stop = timing_stop();
	std::cout << "neighbors detection, elapsed " << timing_seconds(start,stop) << std::endl;
      }
#endif      
      //      WMESH_STATUS_CHECK(status);

      {
	auto start = timing_start();      
	status = wmesh_analysis_edges(self_, &self_->m_num_edges);
	auto stop = timing_stop();
	std::cout << "edge analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	WMESH_STATUS_CHECK(status);
      }
      
      {
	std::cout << "// wmesh_analysis ... faces indexing ..."  << std::endl;
	auto start = timing_start();      
	status    = wmesh_analysis_faces(self_);
	auto stop = timing_stop();
	std::cout << "faces analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	WMESH_STATUS_CHECK(status);
      }

      WMESH_STATUS_CHECK(status);
    }

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




#if 0
    
    wmesh_int_t ndofs[4]= { ((degree+1)*(degree+2)*(degree+3)) / 6,
			    ((degree+2)*(degree+1)*(2*degree+3)) / 6,
			    ((degree+1)*(degree+1)*(degree+2)) / 2,
			    (degree+1)*(degree+1)*(degree+1)};
      
      
    wmesh_int_t ndofs_n[4] = {1,1,1,1};
    wmesh_int_t ndofs_e[4] = {degree-1,degree-1,degree-1,degree-1};

    wmesh_int_t ndofs_t[4] = {(degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0};

    wmesh_int_t ndofs_q[4] = {(degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0};
            
    wmesh_int_t ndofs_i[4] = { (degree>0) ? ((degree-1)*(degree-2)*(degree-3)) / 6 : 1,
			       (degree>0) ? ( (degree-2)*(degree-1)*(2*degree-3) ) / 6 : 1,
			       (degree>0) ? ((degree-1)*(degree-1)*(degree-2)) / 2 : 1, //				 (degree>2) ? (degree-2)*(degree-1) : 0,
			       (degree>0) ? (degree-1)*(degree-1)*(degree-1) : 1 };
#if 1
    //
    // cells to dofs.
    //
    {
      wmesh_status_t status;
      
      //
      // Use c2f since this is the same layout.
      //
      wmesh_int_t shifts[4] = {0,0,0,0};
      
      status = wmesh_init_c2d	(&self_->m_c2n,
				 &self_->m_c2d,
				 ndofs,
				 ndofs);
      WMESH_STATUS_CHECK(status);

      //
      // Node-based dofs.
      //
      status = wmesh_init_c2d_n	(&self_->m_c2d,
				 &self_->m_c2d_n,
				 ndofs_n,
				 shifts);
      WMESH_STATUS_CHECK(status);

      //
      // edge-based dofs.
      //
      status = wmesh_init_c2d_e	(&self_->m_c2d,
				 &self_->m_c2d_e,
				 ndofs_e,
				 shifts);
      WMESH_STATUS_CHECK(status);

      //
      // triangle-based dofs.
      //
      status = wmesh_init_c2d_t	(&self_->m_c2d,
				 &self_->m_c2d_t,
				 ndofs_t,
				 shifts);
      WMESH_STATUS_CHECK(status);
      
      //
      // quadrilateral-based dofs.
      //
      status = wmesh_init_c2d_q	(&self_->m_c2d,
				 &self_->m_c2d_q,
				 ndofs_q,
				 shifts);
      WMESH_STATUS_CHECK(status);      
      
      //
      // cell-based dofs.
      //
      status = wmesh_init_c2d_i(&self_->m_c2d,
				&self_->m_c2d_i,
				ndofs_i,
				shifts);
      WMESH_STATUS_CHECK(status);
    }
#endif
#if 0
    {
      auto start = timing_start();      
      wmesh_status_t status = wmeshspace_analysis(self_,
						  degree);
      auto stop = timing_stop();
      std::cout << "analysis space, elapsed " << timing_seconds(start,stop) << std::endl;
      WMESH_STATUS_CHECK(status);
    }
#endif
    {
#if 0
      
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_n,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_e,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_t,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_q,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_i,stdout);
#endif
    }


#if 0
    {
      for (wmesh_int_t degree = 2;degree<=11;++degree)
	{
	  
	  wmesh_int_t ndofs[4]= { ((degree+1)*(degree+2)*(degree+3)) / 6,
				  ((degree+2)*(degree+1)*(2*degree+3)) / 6,
				  ((degree+1)*(degree+1)*(degree+2)) / 2,
				  (degree+1)*(degree+1)*(degree+1)};
	  wmesh_int_t ndofs_per_node = 1;
	  wmesh_int_t ndofs_per_edge = degree-1;
	  wmesh_int_t ndofs_per_faces[2] = {(degree>2) ? ((degree-1)*(degree-2))/2 : 0,
					    (degree>1) ? (degree-1)*(degree-1) : 0};
	  
	  wmesh_int_t ndofs_per_vol[4] = { (degree>2) ? ((degree+1)*(degree+2)*(degree+3)) / 6 - ndofs_per_faces[0]*4-4-6*ndofs_per_edge : 0,
					   (degree>0) ? ( (degree-2)*(degree-1)*(2*degree-3) ) / 6 : 0,
					   (degree>2) ? (degree-2)*(degree-1) : 0,
					   (degree>1) ? (degree-1)*(degree-1)*(degree-1) : 0 };
#if 0
	  std::cout << "edge part "<<ndofs_per_edge * num_edges[0] << std::endl;
	  std::cout << "face part "<<num_faces[0] * ndofs_per_faces[0] +
	    num_faces[1] * ndofs_per_faces[1]  << std::endl;
	  std::cout << "volume part "<<self_->m_c2n.m_n[0] * ndofs_per_vol[0] +
	    self_->m_c2n.m_n[1] * ndofs_per_vol[1] +
	    self_->m_c2n.m_n[2] * ndofs_per_vol[2] +
	    self_->m_c2n.m_n[3] * ndofs_per_vol[3] << std::endl;
#endif
	  std::cout << self_->m_num_nodes << std::endl;
	  std::cout << "P1 " <<  ndofs_per_node * self_->m_num_nodes << std::endl;
	  std::cout << "P2 " <<  ndofs_per_node * self_->m_num_nodes + ndofs_per_edge * self_->m_num_edges << std::endl;
	  std::cout << "P3 " <<
	    num_faces[0] * ndofs_per_faces[0] +
	    num_faces[1] * ndofs_per_faces[1] +
	    ndofs_per_node * self_->m_num_nodes +
	    ndofs_per_edge * self_->m_num_edges << std::endl;

	  for (int i=0;i<4;++i) std::cout << "$#$# " << ndofs_per_vol[i]<< std::endl;
	  unsigned long long int total_ndofs =
	    ndofs_per_node * self_->m_num_nodes +
	    ndofs_per_edge * self_->m_num_edges +
	    num_faces[0] * ndofs_per_faces[0] +
	    num_faces[1] * ndofs_per_faces[1] +
	    self_->m_c2n.m_n[0] * ndofs_per_vol[0] +
	    self_->m_c2n.m_n[1] * ndofs_per_vol[1] +
	    self_->m_c2n.m_n[2] * ndofs_per_vol[2] +
	    self_->m_c2n.m_n[3] * ndofs_per_vol[3];

	  unsigned long long int total_dndofs =
	    self_->m_c2n.m_n[0] * ndofs[0] +
	    self_->m_c2n.m_n[1] * ndofs[1] +
	    self_->m_c2n.m_n[2] * ndofs[2] +
	    self_->m_c2n.m_n[3] * ndofs[3];
	  
	  std::cout << "total ndofs[" << degree << "]:  " << total_ndofs  << " " << total_dndofs  << std::endl;
	}
      
    }
#endif
#endif
    return WMESH_STATUS_SUCCESS;
  }




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



  static wmesh_status_t wmesh_refine_c2c(wmesh_t* 		self_,
					 wmesh_int_p 		num_faces_,
					 wmesh_int_p 		num_boundary_faces_)
  {    
    
#define WINT_SPARSE_MAT_PARAMS2(_f)				\
    _f.m_size,							\
      _f.m_ptr,							\
      _f.m_m,							\
      _f.m_n,							\
      _f.m_data,						\
      _f.m_ld
    
    wmesh_int_t work_n;
    wmesh_int_t * work  = nullptr;
    wmesh_status_t status;
    
    status =  wbms_c2c_calculate_buffer_size(self_->m_c2n.m_size,
					     self_->m_c2n.m_n,
					     &work_n);    
    WMESH_STATUS_CHECK(status);    
    work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }      
    num_boundary_faces_[0] = 0;
    num_boundary_faces_[1] = 0;
    num_faces_[0] = 0;
    num_faces_[1] = 0;
    status = wbms_c2c_calculate(WINT_SPARSE_MAT_PARAMS2(self_->m_c2n),
				WINT_SPARSE_MAT_PARAMS2(self_->m_c2c_t),
				WINT_SPARSE_MAT_PARAMS2(self_->m_c2c_q),
				WINT_SPARSE_MAT_PARAMS2(self_->m_s_t2n),
				WINT_SPARSE_MAT_PARAMS2(self_->m_s_q2n),
				work_n,
				work,
				num_faces_,
				num_boundary_faces_);
    free(work);      
    WMESH_STATUS_CHECK(status);    
    return WMESH_STATUS_SUCCESS;



    
#undef WINT_SPARSE_MAT_PARAMS2
    
  }
    
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
			    &self_->m_c2e);
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
    auto start = timing_start();      
    status = wmesh_refine_c2c(self_,
			      num_faces,
			      num_bfaces);
    
    self_->m_num_triangles      = num_faces[0];
    self_->m_num_quadrilaterals = num_faces[1];
    self_->m_num_bfaces[0]      = num_bfaces[0];
    self_->m_num_bfaces[1]      = num_bfaces[1];
    auto stop = timing_stop();
    std::cout << "neighbors detection, elapsed " << timing_seconds(start,stop) << std::endl;
    WMESH_STATUS_CHECK(status);
    
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

      wmesh_int_sparsemat_def(&self_->m_bf_c,
			      bf_c_size);
    
      wmesh_int_sparsemat_set(&self_->m_bf_c,
			      bf_c_m,
			      bf_c_n,
			      bf_c_ld,
			      bf_c_ptr,
			      bf_c_v);
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
  
  wmesh_status_t wmesh_refine_dev(wmesh_t*    self_,			      
			      wmesh_int_t degree,
			      wmesh_t**   refined_mesh_)
  {
    //    wmesh_extract_boundary(self_);
    // wmesh_int_t degree = 4;
    wmesh_int_t num_faces[2];
    wmesh_int_t num_bfaces[2];
    
    //
    // 2/ Create the cells to cells.
    //
    {
      wmesh_status_t status;

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
    }


    //
    // 3/ Create the cells to edges.
    //    
    {
      wmesh_status_t status;
      status = wmesh_init_c2e(&self_->m_c2n,
			      &self_->m_c2e);
      WMESH_STATUS_CHECK(status);
    }
    
    //
    // 4/ Create the cells to faces.
    //    
    {
      wmesh_status_t status;
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
    }
    
    //
    // Initialize the reference shape edges to nodes.
    //
    {
      wmesh_status_t status;
#if 1
      {
	auto start = timing_start();      
	status = wmesh_analysis_neighbors(self_,
					  num_faces,
					  num_bfaces);
	
	self_->m_num_triangles      = num_faces[0];
	self_->m_num_quadrilaterals = num_faces[1];
	self_->m_num_bfaces[0]      = num_bfaces[0];
	self_->m_num_bfaces[1]      = num_bfaces[1];
	auto stop = timing_stop();
	std::cout << "neighbors detection, elapsed " << timing_seconds(start,stop) << std::endl;
      }
#endif      
      //      WMESH_STATUS_CHECK(status);

      {
	auto start = timing_start();      
	status = wmesh_analysis_edges(self_, &self_->m_num_edges);
	auto stop = timing_stop();
	std::cout << "edge analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	WMESH_STATUS_CHECK(status);
      }
      
      {
	std::cout << "// wmesh_analysis ... faces indexing ..."  << std::endl;
	auto start = timing_start();      
	status    = wmesh_analysis_faces(self_);
	auto stop = timing_stop();
	std::cout << "faces analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	WMESH_STATUS_CHECK(status);
      }

      WMESH_STATUS_CHECK(status);
    }

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

    wmesh_int_t ndofs[4]= { ((degree+1)*(degree+2)*(degree+3)) / 6,
			    ((degree+2)*(degree+1)*(2*degree+3)) / 6,
			    ((degree+1)*(degree+1)*(degree+2)) / 2,
			    (degree+1)*(degree+1)*(degree+1)};
      
      
    wmesh_int_t ndofs_n[4] = {1,1,1,1};
    wmesh_int_t ndofs_e[4] = {degree-1,degree-1,degree-1,degree-1};

    wmesh_int_t ndofs_t[4] = {(degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0,
			      (degree>2) ? ((degree-1)*(degree-2))/2 : 0};

    wmesh_int_t ndofs_q[4] = {(degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0,
			      (degree>1) ? (degree-1)*(degree-1) : 0};
            
    wmesh_int_t ndofs_i[4] = { (degree>0) ? ((degree-1)*(degree-2)*(degree-3)) / 6 : 1,
			       (degree>0) ? ( (degree-2)*(degree-1)*(2*degree-3) ) / 6 : 1,
			       (degree>0) ? ((degree-1)*(degree-1)*(degree-2)) / 2 : 1, //				 (degree>2) ? (degree-2)*(degree-1) : 0,
			       (degree>0) ? (degree-1)*(degree-1)*(degree-1) : 1 };

#if 1
    //
    // cells to dofs.
    //
    {
      wmesh_status_t status;
      
      //
      // Use c2f since this is the same layout.
      //
      wmesh_int_t shifts[4] = {0,0,0,0};
      
      status = wmesh_init_c2d	(&self_->m_c2n,
				 &self_->m_c2d,
				 ndofs,
				 ndofs);
      WMESH_STATUS_CHECK(status);

      //
      // Node-based dofs.
      //
      status = wmesh_init_c2d_n	(&self_->m_c2d,
				 &self_->m_c2d_n,
				 ndofs_n,
				 shifts);
      WMESH_STATUS_CHECK(status);

      //
      // edge-based dofs.
      //
      status = wmesh_init_c2d_e	(&self_->m_c2d,
				 &self_->m_c2d_e,
				 ndofs_e,
				 shifts);
      WMESH_STATUS_CHECK(status);

      //
      // triangle-based dofs.
      //
      status = wmesh_init_c2d_t	(&self_->m_c2d,
				 &self_->m_c2d_t,
				 ndofs_t,
				 shifts);
      WMESH_STATUS_CHECK(status);
      
      //
      // quadrilateral-based dofs.
      //
      status = wmesh_init_c2d_q	(&self_->m_c2d,
				 &self_->m_c2d_q,
				 ndofs_q,
				 shifts);
      WMESH_STATUS_CHECK(status);      
      
      //
      // cell-based dofs.
      //
      status = wmesh_init_c2d_i(&self_->m_c2d,
				&self_->m_c2d_i,
				ndofs_i,
				shifts);
      WMESH_STATUS_CHECK(status);
    }
#endif
#if 1
    {
      auto start = timing_start();      
      wmesh_status_t status = wmesh_refine_calculate(self_,
						     degree,
						     refined_mesh_);
      auto stop = timing_stop();
      std::cout << "analysis space, elapsed " << timing_seconds(start,stop) << std::endl;
      WMESH_STATUS_CHECK(status);
    }
#endif
    {
#if 0
      
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_n,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_e,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_t,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_q,stdout);
      fprintf(stdout,"------\n");
      wmesh_int_sparsemat_fprintf(&self_->m_c2d_i,stdout);
#endif
    }



    {
      for (wmesh_int_t degree = 2;degree<=11;++degree)
	{
	  
	  wmesh_int_t ndofs[4]= { ((degree+1)*(degree+2)*(degree+3)) / 6,
				  ((degree+2)*(degree+1)*(2*degree+3)) / 6,
				  ((degree+1)*(degree+1)*(degree+2)) / 2,
				  (degree+1)*(degree+1)*(degree+1)};
	  wmesh_int_t ndofs_per_node = 1;
	  wmesh_int_t ndofs_per_edge = degree-1;
	  wmesh_int_t ndofs_per_faces[2] = {(degree>2) ? ((degree-1)*(degree-2))/2 : 0,
					    (degree>1) ? (degree-1)*(degree-1) : 0};
	  
	  wmesh_int_t ndofs_per_vol[4] = { (degree>2) ? ((degree+1)*(degree+2)*(degree+3)) / 6 - ndofs_per_faces[0]*4-4-6*ndofs_per_edge : 0,
					   (degree>0) ? ( (degree-2)*(degree-1)*(2*degree-3) ) / 6 : 0,
					   (degree>2) ? (degree-2)*(degree-1) : 0,
					   (degree>1) ? (degree-1)*(degree-1)*(degree-1) : 0 };
#if 0
	  std::cout << "edge part "<<ndofs_per_edge * num_edges[0] << std::endl;
	  std::cout << "face part "<<num_faces[0] * ndofs_per_faces[0] +
	    num_faces[1] * ndofs_per_faces[1]  << std::endl;
	  std::cout << "volume part "<<self_->m_c2n.m_n[0] * ndofs_per_vol[0] +
	    self_->m_c2n.m_n[1] * ndofs_per_vol[1] +
	    self_->m_c2n.m_n[2] * ndofs_per_vol[2] +
	    self_->m_c2n.m_n[3] * ndofs_per_vol[3] << std::endl;
#endif
	  std::cout << self_->m_num_nodes << std::endl;
	  std::cout << "P1 " <<  ndofs_per_node * self_->m_num_nodes << std::endl;
	  std::cout << "P2 " <<  ndofs_per_node * self_->m_num_nodes + ndofs_per_edge * self_->m_num_edges << std::endl;
	  std::cout << "P3 " <<
	    num_faces[0] * ndofs_per_faces[0] +
	    num_faces[1] * ndofs_per_faces[1] +
	    ndofs_per_node * self_->m_num_nodes +
	    ndofs_per_edge * self_->m_num_edges << std::endl;

	  for (int i=0;i<4;++i) std::cout << "$#$# " << ndofs_per_vol[i]<< std::endl;
	  unsigned long long int total_ndofs =
	    ndofs_per_node * self_->m_num_nodes +
	    ndofs_per_edge * self_->m_num_edges +
	    num_faces[0] * ndofs_per_faces[0] +
	    num_faces[1] * ndofs_per_faces[1] +
	    self_->m_c2n.m_n[0] * ndofs_per_vol[0] +
	    self_->m_c2n.m_n[1] * ndofs_per_vol[1] +
	    self_->m_c2n.m_n[2] * ndofs_per_vol[2] +
	    self_->m_c2n.m_n[3] * ndofs_per_vol[3];

	  unsigned long long int total_dndofs =
	    self_->m_c2n.m_n[0] * ndofs[0] +
	    self_->m_c2n.m_n[1] * ndofs[1] +
	    self_->m_c2n.m_n[2] * ndofs[2] +
	    self_->m_c2n.m_n[3] * ndofs[3];
	  
	  std::cout << "total ndofs[" << degree << "]:  " << total_ndofs  << " " << total_dndofs  << std::endl;
	}
      
    }
    return WMESH_STATUS_SUCCESS;
  }


};
#endif
