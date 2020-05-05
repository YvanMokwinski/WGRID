#if 0
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "wfe_element.hpp"
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
    wmesh_int_p 	c2d_v    = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*c2d_ptr[4]);
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



  
  static wmesh_status_t wmesh_init_c2f(const wmesh_int_sparsemat_t*	c2n_,
				       wmesh_int_sparsemat_t*		c2f_,
				       wmesh_int_t 			ix_ptr_,
				       wmesh_int_t * 			ix_v_)
  {
    wmesh_int_t c2f_m[4] {4,5,5,6};
    wmesh_int_t c2f_ld[4] {4,5,5,6};
    wmesh_int_t c2f_ptr[4+1];
    wmesh_int_t c2f_n[4];      
    for (wmesh_int_t i=0;i<4;++i) c2f_n[i] = c2n_->m_n[i];

    c2f_ptr[0] = ix_ptr_;
    c2f_ptr[1] = c2f_ptr[0] + c2f_n[0] * c2f_ld[0];
    c2f_ptr[2] = c2f_ptr[1] + c2f_n[1] * c2f_ld[1];
    c2f_ptr[3] = c2f_ptr[2] + c2f_n[2] * c2f_ld[2];
    c2f_ptr[4] = c2f_ptr[3] + c2f_n[3] * c2f_ld[3];
    wmesh_int_sparsemat_def(c2f_,
			    4);
    
    wmesh_int_sparsemat_set(c2f_,
			    c2f_m,
			    c2f_n,
			    c2f_ld,
			    c2f_ptr,
			    ix_v_);
    return WMESH_STATUS_SUCCESS;
  }

  static wmesh_status_t wmesh_init_c2f_q(wmesh_int_sparsemat_t*	c2f_,
					 wmesh_int_sparsemat_t*	c2f_q_)
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
				       wmesh_int_sparsemat_t*		c2e_,
				       wmesh_int_t *			ix_ptr_,
				       wmesh_int_t * 			ix_v_)
  {
    wmesh_int_t c2e_m[4] {6,8,9,12};
    wmesh_int_t c2e_ptr[4+1];
    wmesh_int_t c2e_n[4];      
    for (wmesh_int_t i=0;i<4;++i) c2e_n[i] = c2n_->m_n[i];
    c2e_ptr[0] = ix_ptr_[0];
    c2e_ptr[1] = c2e_ptr[0] + c2e_n[0] * c2e_m[0];
    c2e_ptr[2] = c2e_ptr[1] + c2e_n[1] * c2e_m[1];
    c2e_ptr[3] = c2e_ptr[2] + c2e_n[2] * c2e_m[2];
    c2e_ptr[4] = c2e_ptr[3] + c2e_n[3] * c2e_m[3];
    ix_ptr_[0] = c2e_ptr[4];
    wmesh_int_sparsemat_def(c2e_,
			    4);
    wmesh_int_sparsemat_set(c2e_,
			    c2e_m,
			    c2e_n,
			    c2e_m,
			    c2e_ptr,
			    ix_v_);
    return 0;
  }


  static wmesh_status_t wmesh_space_num_interior_dofs(wmesh_int_p 	num_dofs_,
						      wmesh_int_t	degree_)
  {    
    WMESH_POINTER_CHECK(num_dofs_);
    if (degree_ < 0)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
    for (wmesh_int_t i=0;i<WFE_ELEMENT_ALL;++i) num_dofs_[i] = 0;
    num_dofs_[WFE_ELEMENT_POINT] 		= (degree_>0) ? 1 : 0; 
    num_dofs_[WFE_ELEMENT_EDGE]  		= (degree_>0) ? degree_-1 : 0;      
    num_dofs_[WFE_ELEMENT_TRIANGLE]  		= (degree_>0) ? ((degree_-1)*(degree_-2))/2 : 0; 
    num_dofs_[WFE_ELEMENT_QUADRILATERAL]  	= (degree_>0) ? (degree_-1)*(degree_-1) : 0;      
    num_dofs_[WFE_ELEMENT_TETRAHEDRON]  	= (degree_>0) ? ((degree_-1)*(degree_-2)*(degree_-3)) / 6 : 1;
    num_dofs_[WFE_ELEMENT_PYRAMID]  		= (degree_>0) ? ( (degree_-1)*(degree_-2)*(2*degree_-3) ) / 6 : 1;
    num_dofs_[WFE_ELEMENT_WEDGE] 		= (degree_>0) ? (( (degree_-1)*(degree_-2) ) / 2) * (degree_-1) : 1;
    num_dofs_[WFE_ELEMENT_HEXAHEDRON] 		= (degree_>0) ? (degree_-1)*(degree_-1)*(degree_-1) : 1;
    return WMESH_STATUS_SUCCESS;
  };

  static wmesh_status_t wmesh_space_num_dofs(wmesh_int_p 	num_dofs_,
						 wmesh_int_t	degree_)
  {    
    WMESH_POINTER_CHECK(num_dofs_);
    if (degree_ < 0)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
    for (wmesh_int_t i=0;i<WFE_ELEMENT_ALL;++i) num_dofs_[i] = 0;
    num_dofs_[WFE_ELEMENT_POINT] 		= (degree_>0) ? 1 : 0; 
    num_dofs_[WFE_ELEMENT_EDGE]  		= (degree_>0) ? degree_+1 : 0;      
    num_dofs_[WFE_ELEMENT_TRIANGLE]  		= (degree_>0) ? ((degree_+1)*(degree_+2))/2 : 0; 
    num_dofs_[WFE_ELEMENT_QUADRILATERAL]  	= (degree_>0) ? (degree_+1)*(degree_+1) : 0;      
    num_dofs_[WFE_ELEMENT_TETRAHEDRON]  	= (degree_>0) ? ((degree_+1)*(degree_+2)*(degree_+3)) / 6 : 1;
    num_dofs_[WFE_ELEMENT_PYRAMID]  		= (degree_>0) ? ( (degree_+1)*(degree_+2)*(2*degree_+3) ) / 6 : 1;
    num_dofs_[WFE_ELEMENT_WEDGE] 		= (degree_>0) ? ( ( (degree_+1)*(degree_+2) ) /2 ) * (degree_+1) : 1;
    num_dofs_[WFE_ELEMENT_HEXAHEDRON] 		= (degree_>0) ? (degree_+1)*(degree_+1)*(degree_+1) : 1;
    return WMESH_STATUS_SUCCESS;
  };
  
  wmesh_status_t wmesh_space_def(wmesh_t**		self_,
				 const wmesh_t*		mesh_linear_,
				 wmesh_int_t 		degree_)
  {
    wmesh_status_t status;

#ifndef NDEBUG
    std::cerr
      << "//wmesh.verbose: wmesh_space_def ..."
      << std::endl;
#endif
    
    wmesh_int_sparsemat_t 	c2d;
    wmesh_int_sparsemat_t 	c2d_n;
    wmesh_int_sparsemat_t 	c2d_e;
    wmesh_int_sparsemat_t 	c2d_t;
    wmesh_int_sparsemat_t 	c2d_q;
    wmesh_int_sparsemat_t 	c2d_i;    
    double * 			coo_dofs;    


    wmesh_int_t num_interior_dofs[WFE_ELEMENT_ALL];
    wmesh_int_t num_dofs[WFE_ELEMENT_ALL];
    
    //
    // Get the number of dofs per entities.
    //
    status = wmesh_space_num_interior_dofs(num_interior_dofs,
					   degree_);
    WMESH_STATUS_CHECK(status);
    
    const wmesh_int_t num_dofs_per_node 		= (degree_ > 0) ? 1 : 0;      
    const wmesh_int_t num_dofs_per_edge 		= (degree_>0) ? degree_-1 : 0;
    const wmesh_int_t num_dofs_per_triangle		= (degree_>0) ? ((degree_-1)*(degree_-2))/2 : 0;
    const wmesh_int_t num_dofs_per_quadrilateral 	= (degree_>0) ? (degree_-1)*(degree_-1) : 0;      
    const wmesh_int_t ndofs_per_vol[4] 			= { (degree_ > 0) ? ((degree_-1)*(degree_-2)*(degree_-3)) / 6 : 1,
							    (degree_ > 0) ? ( (degree_-2)*(degree_-1)*(2*degree_-3) ) / 6 : 1,
							    (degree_ > 0) ? (degree_-2)*(degree_-1) : 1,
							    (degree_ > 0) ? (degree_-1)*(degree_-1)*(degree_-1) : 1};    

    //
    // Start with the base.
    //
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
					    WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
					    WINT_SPARSE_MAT_PARAMS(c2d_n),
					    mesh_linear_->m_num_nodes,
					    num_dofs_per_node,
					    dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += mesh_linear_->m_num_nodes * num_dofs_per_node;	  
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
					     WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
					     WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2e),
					     WINT_SPARSE_MAT_PARAMS(c2d_e),
					     WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_e2n),
					     mesh_linear_->m_num_edges,
					     num_dofs_per_edge,
					     dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += mesh_linear_->m_num_edges * num_dofs_per_edge;
      }
    //wmesh_int_sparsemat_fprintf(&mesh_linear_->m_c2e,stdout);
    //wmesh_int_sparsemat_fprintf(&&c2d_e,stdout);
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
						 WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
						 WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_t),
						 WINT_SPARSE_MAT_PARAMS(c2d_t),
						 WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_t2n),
						 mesh_linear_->m_num_triangles,
						 degree_,
						 num_dofs_per_triangle,
						 dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += mesh_linear_->m_num_triangles * num_dofs_per_triangle;
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
						     WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
						     WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_q),
						     WINT_SPARSE_MAT_PARAMS(c2d_q),
						     WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_q2n),
						     mesh_linear_->m_num_quadrilaterals,
						     degree_,
						     num_dofs_per_quadrilateral,
						     dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += mesh_linear_->m_num_quadrilaterals * num_dofs_per_quadrilateral;
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
					    WINT_SPARSE_MAT_PARAMS(c2d_i),
					    &dof_idx);
    WMESH_STATUS_CHECK(status);
#if 0
    printf("----\n");
    wmesh_int_sparsemat_fprintf(&c2d,stdout);
    printf("----\n");
    exit(1);
#endif
      
#if 0      
    wmesh_int_t ndofs_i[4] = { (degree_>2) ? ((degree_-1)*(degree_-2)*(degree_-3)) / 6 : 0,
			       (degree_>0) ? ( (degree_-2)*(degree_-1)*(2*degree_-3) ) / 6 : 0,
			       (degree_>2) ? (degree_-2)*(degree_-1) : 0,
			       (degree_>1) ? (degree_-1)*(degree_-1)*(degree_-1) : 0 };
    for (wmesh_int_t i=0;i<4;++i)
      {
	return wmesh_space_indexing_cell_interior(0,
						  WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
						  WINT_SPARSE_MAT_PARAMS(c2d_i),
						  mesh_linear_->m_num_tetrahedra,
						  ndofs_i[0],
						  mesh_linear_->m_num_nodes * num_dofs_per_node + mesh_linear_->m_num_edges * num_dofs_per_edge + mesh_linear_->m_num_triangles * num_dofs_per_triangle );
      }
#endif
      
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
	refmeshes[l]=nullptr;
	refevals[l]=nullptr;
	if (mesh_linear_->m_c2n.m_n[l]>0)
	  {
#ifndef NDEBUG
	    std::cerr << "// wmesh.ndebug.verbose: wmesh_analysis, refmesh celltype_=" << l << ", degree_ " << degree_ << std::endl;
#endif
	    printf("prepare treilli\n.");      	      
	    wmesh_int_t work_n;
	    wmesh_int_t ref_num_entities[WFE_ELEMENT_ALL];
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
		      std::cout << "rst = " << r << " " << s << " " << t << " " << std::endl;
		      
#if 0
		      double one=1.0;
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

      printf("generate dofs.\n");      
      wmesh_int_t ndofs = dof_idx - 1;
      self_[0]->m_num_nodes 	= ndofs;
      self_[0]->m_coo = (double*)malloc(3*sizeof(double)*self_[0]->m_num_nodes);

      {
	double cell_xyz[32];
	wmesh_int_t cell_ld=3;
	double * 	c_xyz 	= (double*)malloc(sizeof(double)*2048);
	wmesh_int_p 	c_dofs 	= (wmesh_int_p)malloc(sizeof(wmesh_int_t)*2048);
	for (wmesh_int_t l=0;l<c2d.m_size;++l)
	  {
	    
	    wmesh_t * 	refmesh = refmeshes[l];
	    double * 	refeval = refevals[l];

	    //
	    // Local c2d.
	    //	    
	    wmesh_int_t c2d_n 		= c2d.m_n[l];
	    wmesh_int_t c2d_m 		= c2d.m_m[l];
	    wmesh_int_t c2d_ld 		= c2d.m_ld[l];
	    wmesh_int_p c2d_v 		= c2d.m_data + c2d.m_ptr[l];

	    //
	    // Local c2n.
	    //
	    wmesh_int_t c2n_n 		= mesh_linear_->m_c2n.m_n[l];
	    wmesh_int_t c2n_m 		= mesh_linear_->m_c2n.m_m[l];
	    wmesh_int_t c2n_ld		= mesh_linear_->m_c2n.m_ld[l];
	    wmesh_int_p c2n_v 		= mesh_linear_->m_c2n.m_data + mesh_linear_->m_c2n.m_ptr[l];
	    
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
			cell_xyz[cell_ld * i + k] = mesh_linear_->m_coo[3 * idx + k];
		      }
		  }
		
		//
		// Get the physical coordinates of the dofs.
		//
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
		    coo_dofs[3 * idx + 0] = x;
		    coo_dofs[3 * idx + 1] = y;
		    coo_dofs[3 * idx + 2] = z;
		  }
	      }
	  }
      }

      {
	wmesh_int_t c2n_size = 4;
	wmesh_int_t c2n_m[4]{4,5,6,8};
	wmesh_int_t c2n_ld[4]{4,5,6,8};
	wmesh_int_t c2n_n[4]{0,0,0,0};
	for (int i=0;i<4;++i)
	  {
	    if (refmeshes[i]!=nullptr)
	      {
		for (int j=0;j<4;++j)
		  {
		    c2n_n[j] += mesh_linear_->m_c2n.m_n[i] * refmeshes[i]->m_c2n.m_n[j];
		  }
#if 0		
		if (i==1)
		  {
		    c2n_n[0] += mesh_linear_->m_c2n.m_n[i] * refmeshes[i]->m_c2n.m_n[0];
		    c2n_n[1] += mesh_linear_->m_c2n.m_n[i] * refmeshes[i]->m_c2n.m_n[1];
		  }
		else
		  {
		    c2n_n[i] = mesh_linear_->m_c2n.m_n[i] * refmeshes[i]->m_num_cells;
		  }
#endif
	      }
	  }

	
	wmesh_int_t c2n_ptr[5];
	c2n_ptr[0]=0;
	for (int i=0;i<4;++i)
	  {
	    c2n_ptr[i+1] = c2n_ptr[i] + c2n_ld[i] * c2n_n[i];
	  }
	
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

	    wmesh_int_t ncells = c2d.m_n[l];
	    for (wmesh_int_t j=0;j<ncells;++j)
	      {
		//
		// extract dofs.
		//
		for (wmesh_int_t i=0;i<c2d.m_m[l];++i)
		  {
		    dofs[i] = c2d.m_data[c2d.m_ptr[l] + j * c2d.m_ld[l] + i];
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
      printf("define mesh.\n");
      
      wmesh_t * space_sublinear;
      //	c2n_n[3]=0;
      status =  wmesh_def	(&space_sublinear,
				 mesh_linear_->m_ndofs,
				 c2n_size,
				 c2n_ptr,
				 c2n_m,
				 c2n_n,
				 c2n_v,
				 c2n_ld,
				 coo_dofs,
				 3);
      
      WMESH_STATUS_CHECK(status);
      }
      return WMESH_STATUS_SUCCESS;
  };

  
  wmesh_status_t wmesh_analysis_edges(wmesh_t* 		mesh_linear_,
				      wmesh_int_p 	num_edges_)
  {
    //    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    wmesh_int_t edge_idx = 0;
    
    {
      wmesh_int_t work_n = 0;
      wmesh_int_t * work  = nullptr;
      wmesh_status_t status =  wmesh_indexing_edges_hybrid(4,
							   WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
							   WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2e),
							   WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_e2n),				     
							   &edge_idx,
							   &work_n,
							   work);      
      WMESH_STATUS_CHECK(status);
      
      work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
      if (!work)
	{
	  return WMESH_STATUS_ERROR_MEMORY;
	}      

      status =  wmesh_indexing_edges_hybrid(4,
					    WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
					    WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2e),
					    WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_e2n),				     
					    &edge_idx,
					    &work_n,
					    work);
      free(work);      
      WMESH_STATUS_CHECK(status);

    }
    
    num_edges_[0] = edge_idx;
  //    wmesh_int_sparsemat_fprintf(&mesh_linear_->m_c2e,stdout);
    return 0;
  }


  wmesh_status_t wmesh_analysis_faces(wmesh_t* 		mesh_linear_)
  {    
    wmesh_int_t face_idx = 0;

#define WINT_SPARSE_MAT_PARAMS(_f)				\
    _f.m_ptr,							\
      _f.m_m,							\
      _f.m_n,							\
      _f.m_data,						\
      _f.m_ld
    
    //
    // Enumerates triangles faces.
    //
    {
      wmesh_int_t work_n=0;
      wmesh_int_t * work  = nullptr;
      wmesh_status_t status =  wmesh_indexing_triangles_hybrid(4,
							       WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
							       WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_t),
							       WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_t2n),				     
							       &face_idx,
							       &work_n,
							       work);      
      WMESH_STATUS_CHECK(status);
      
      work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
      if (!work)
	{
	  return WMESH_STATUS_ERROR_MEMORY;
	}      

      status =  wmesh_indexing_triangles_hybrid(4,
						WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
						WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_t),
						WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_t2n),				     
						&face_idx,
						&work_n,
						work);
      free(work);      
      WMESH_STATUS_CHECK(status);

    }

    //
    // Enumerates quadrilateral faces.
    //
    {
      face_idx = 0;
      wmesh_int_t work_n;
      wmesh_int_t * work  = nullptr;
      wmesh_status_t status =  wmesh_indexing_quadrilaterals_hybrid(4,
								    WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
								    WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_q),
								    WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_q2n),				     
								    &face_idx,
								    &work_n,
								    work);      
      WMESH_STATUS_CHECK(status);
      
      work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
      if (!work)
	{
	  return WMESH_STATUS_ERROR_MEMORY;
	}      

      status =  wmesh_indexing_quadrilaterals_hybrid(4,
						     WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
						     WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_q),
						     WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_q2n),				     
						     &face_idx,
						     &work_n,
						     work);
      free(work);      
      WMESH_STATUS_CHECK(status);
    }
    

    return 0;
  }


  

  using timing_t = high_resolution_clock::time_point;
  inline timing_t timing_stop(){ return high_resolution_clock::now(); }
  inline timing_t timing_start(){ return high_resolution_clock::now(); }
  inline double timing_seconds(timing_t&t,timing_t&t1){ return duration_cast<duration<double>>(t1-t).count(); }
  
  


  
  wmesh_status_t wmesh_space_prepare(wmesh_t* 		self_,
				     wmesh_int_t 	degree_,
				     wmesh_int_t *      c2d_ptr_,
				     wmesh_int_t * 	c2d_v_)
  {
    //
    // Get the number of dofs per entities.
    //
    static constexpr wmesh_int_t dimension = 3;
    
    //    wmesh_int_t c2d_ptr = c2d_ptr_[0];    
    wmesh_int_t 	num_dofs[WFE_ELEMENT_ALL];
    wmesh_int_t 	num_faces[2];
    wmesh_int_t 	num_bfaces[2];
    wmesh_status_t 	status;

    //
    // 1.a/ Local edges to nodes.
    //
    status  = wmesh_build_s_e2n(&self_->m_s_e2n, dimension);
    WMESH_STATUS_CHECK(status);
    
    //
    // 1.b/ Local triangles to nodes.
    //
    status  = wmesh_build_s_t2n(&self_->m_s_t2n, dimension);
    WMESH_STATUS_CHECK(status);
    
    //
    // 1.c/ Local quadrilaterals to nodes.
    //
    status = wmesh_build_s_q2n(&self_->m_s_q2n, dimension);
    WMESH_STATUS_CHECK(status);
    
    {      
      wmesh_int_t c2n_ptr[4+1];
      wmesh_int_t c2e_ptr[4+1];
      wmesh_int_t c2f_ptr[4+1];
      
      wmesh_int_t c2n_m[4]  { 4, 5, 6, 8};
      wmesh_int_t c2e_m[4]  { 6, 8, 9, 12};
      wmesh_int_t c2f_m[4]  { 4, 5, 5, 6};

      wmesh_int_t c2n_n[4];
      wmesh_int_t c2e_n[4];      
      wmesh_int_t c2f_n[4];      

      wmesh_int_t c2n_ld[4];
      wmesh_int_t c2e_ld[4];      
      wmesh_int_t c2f_ld[4];      
      
      {
	for (int i=0;i<4;++i)
	  c2n_ld[i] = c2n_m[i] + c2e_m[i] + c2f_m[i];
	
	for (int i=0;i<4;++i)
	  c2e_ld[i] = c2n_ld[i];
	for (int i=0;i<4;++i)
	  c2f_ld[i] = c2n_ld[i];
      }
      
      {
	for (int i=0;i<4;++i)
	  c2n_n[i] = c2n_->m_n[i];
	for (int i=0;i<4;++i)
	  c2e_n[i] = c2n_n[i];
	for (int i=0;i<4;++i)
	  c2f_n[i] = c2n_n[i];
      }
      
      c2n_ptr[0] = 0;
      for (int i=0;i<4;++i)
	c2n_ptr[i+1] = c2n_ptr[i] + c2n_ld[i] * c2n_n[i];
      
      //
      //
      //
      for (int i=0;i<4;++i)
	c2e_ptr[i] = c2n_ptr[i] + c2n_m[i];
      c2e_ptr[4] = c2n_ptr[4];
      
      for (int i=0;i<4;++i)
	c2f_ptr[i] = c2e_ptr[i] + c2e_m[i];
      c2f_ptr[4] = c2f_ptr[4];

      wmesh_int_t work_n;
      wmesh_int_t * work  = nullptr;
      
      //
      // Initialize the reference shape edges to nodes.
      //
      {
	auto start = timing_start();      

#define WINT_SPARSE_MAT_PARAMS(_f)		\
	  _f.m_ptr,				\
	  _f.m_m,				\
	  _f.m_n,				\
	  _f.m_data,				\
	  _f.m_ld
	
	//
	// Enumerates triangles faces.
	//
	{
	  wmesh_status_t status;
	  wmesh_int_t num_boundary_triangles = 0;

	  work_n = 0;
	  status = wmesh_c2c_t_calculate(4,
					 WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					 WINT_SPARSE_MAT_PARAMS(self_->m_c2c_t),
					 WINT_SPARSE_MAT_PARAMS(self_->m_s_t2n),
					 &num_boundary_triangles,
					 &work_n,
					 work);    
	  WMESH_STATUS_CHECK(status);

	  work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
	  if (!work)
	    {
	      return WMESH_STATUS_ERROR_MEMORY;
	    }
	  
	  status =  wmesh_c2c_t_calculate(4,
					  WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					  WINT_SPARSE_MAT_PARAMS(self_->m_c2c_t),
					  WINT_SPARSE_MAT_PARAMS(self_->m_s_t2n),
					  &num_boundary_triangles,
					  &work_n,
					  work);    
	  free(work);      
	  WMESH_STATUS_CHECK(status);
	  num_boundary_faces_[0] = num_boundary_triangles;
	  num_faces_[0] = ((self_->m_c2n.m_n[0]*4+self_->m_c2n.m_n[1]*4+self_->m_c2n.m_n[2]*2) + num_boundary_triangles)/2 ;
	}
    
	//
	// Enumerates quadrilaterals faces.
	//
	{
	  wmesh_int_t work_n;
	  wmesh_int_t * work  = nullptr;
	  wmesh_status_t status;
	  wmesh_int_t num_boundary_quadrilaterals = 0;
	  status =  wmesh_c2c_q_calculate(4,
					  WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					  WINT_SPARSE_MAT_PARAMS(self_->m_c2c_q),
					  WINT_SPARSE_MAT_PARAMS(self_->m_s_q2n),
					  &num_boundary_quadrilaterals,
					  &work_n,
					  work);    
	  WMESH_STATUS_CHECK(status);
      
	  work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
	  if (!work)
	    {
	      return WMESH_STATUS_ERROR_MEMORY;
	    }      
      
	  status =  wmesh_c2c_q_calculate(4,
					  WINT_SPARSE_MAT_PARAMS(self_->m_c2n),
					  WINT_SPARSE_MAT_PARAMS(self_->m_c2c_q),
					  WINT_SPARSE_MAT_PARAMS(self_->m_s_q2n),
					  &num_boundary_quadrilaterals,
					  &work_n,
					  work);    
	  free(work);      
	  WMESH_STATUS_CHECK(status);
	  num_boundary_faces_[1] = num_boundary_quadrilaterals;
	  num_faces_[1] = ((self_->m_c2n.m_n[1]*1+self_->m_c2n.m_n[2]*3+self_->m_c2n.m_n[3]*6) + num_boundary_quadrilaterals)/2;
	}
	//
	//	status = wmesh_analysis_neighbors(self_,
	//					  num_faces,
	//					  num_bfaces);
//
    
	self_->m_num_triangles      = num_faces[0];
	self_->m_num_quadrilaterals = num_faces[1];
	self_->m_num_bfaces[0]      = num_bfaces[0];
	self_->m_num_bfaces[1]      = num_bfaces[1];
	auto stop = timing_stop();
	std::cout << "neighbors detection, elapsed " << timing_seconds(start,stop) << std::endl;
      }
      
      WMESH_STATUS_CHECK(status);

      wmesh_int_t num_edges = 0;
      {
	//    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	wmesh_int_t edge_idx = 0;
	
	{
	  wmesh_int_t work_n = 0;
	  wmesh_int_t * work  = nullptr;
	  status =  wmesh_indexing_edges_hybrid(4,
						WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
						c2e_ptr,
						c2e_m,
						c2e_n,
						c2e_v_,
						c2e_ld,							       
						WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_e2n),				     
						&edge_idx,
						&work_n,
						work);      
	  WMESH_STATUS_CHECK(status);	  

	  work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
	  if (!work)
	    {
	      return WMESH_STATUS_ERROR_MEMORY;
	    }
	  
	  status =  wmesh_indexing_edges_hybrid(4,
						WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
						c2e_ptr,
						c2e_m,
						c2e_n,
						c2e_v_,
						c2e_ld,							       
						WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_e2n),				     
						&edge_idx,
						&work_n,
						work);
	  free(work);
	  WMESH_STATUS_CHECK(status);	  
	}
	
	num_edges = edge_idx;
      }
      
      //
      // Enumerates triangles faces.
      //
      {
	triangle_face_idx 	= 0;
	work_n 			= 0;
	work  			= nullptr;
	wmesh_int_t work_n	= 0;
	wmesh_int_t * work  	= nullptr;
	wmesh_status_t status 	= wmesh_indexing_triangles_hybrid(4,
								  WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
								  WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_t),
								  WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_t2n),				     
								  &triangle_face_idx,
								  &work_n,
								  work);      
	WMESH_STATUS_CHECK(status);
	
	work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
	if (!work)
	  {
	    return WMESH_STATUS_ERROR_MEMORY;
	  }      
	
	status =  wmesh_indexing_triangles_hybrid(4,
						  WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
						  WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_t),
						  WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_t2n),				     
						  &triangle_face_idx,
						  &work_n,
						  work);
	free(work);      
	WMESH_STATUS_CHECK(status);	
      }
      
      //
      // Enumerates quadrilateral faces.
      //
      {

	face_idx 	= 0;
	work_n 		= 0;
	work  		= nullptr;
	
	wmesh_status_t status =  wmesh_indexing_quadrilaterals_hybrid(4,
								      WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
								      WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_q),
								      WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_q2n),				     
								      &face_idx,
								      &work_n,
								      work);      
	WMESH_STATUS_CHECK(status);
	
	work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
	if (!work)
	  {
	    return WMESH_STATUS_ERROR_MEMORY;
	  }      
	
	status =  wmesh_indexing_quadrilaterals_hybrid(4,
						       WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
						       WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_q),
						       WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_q2n),				     
						       &face_idx,
						       &work_n,
						       work);
	free(work);      
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
    return WMESH_STATUS_SUCCESS;
  }

  
  wmesh_status_t wmesh_indexing_entities(wmesh_int_t 		c2n_size_,
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

					 wmesh_int_t 		c2t_size_,
					 const_wmesh_int_p 	c2t_ptr_,
					 const_wmesh_int_p 	c2t_m_,
					 const_wmesh_int_p 	c2t_n_,
					 wmesh_int_p 		c2t_v_,
					 const_wmesh_int_p 	c2t_ld_,

					 wmesh_int_t 		c2q_size_,
					 const_wmesh_int_p 	c2q_ptr_,
					 const_wmesh_int_p 	c2q_m_,
					 const_wmesh_int_p 	c2q_n_,
					 wmesh_int_p 		c2q_v_,
					 const_wmesh_int_p 	c2q_ld_,

					 wmesh_int_t 		s_e2n_size_,
					 const_wmesh_int_p 	s_e2n_ptr_,
					 const_wmesh_int_p 	s_e2n_m_,
					 const_wmesh_int_p 	s_e2n_n_,
					 const_wmesh_int_p 	s_e2n_v_,
					 const_wmesh_int_p 	s_e2n_ld_,
					 
					 wmesh_int_t 		s_t2n_size_,
					 const_wmesh_int_p 	s_t2n_ptr_,
					 const_wmesh_int_p 	s_t2n_m_,
					 const_wmesh_int_p 	s_t2n_n_,
					 const_wmesh_int_p 	s_t2n_v_,
					 const_wmesh_int_p 	s_t2n_ld_,

					 wmesh_int_t 		s_q2n_size_,
					 const_wmesh_int_p 	s_q2n_ptr_,
					 const_wmesh_int_p 	s_q2n_m_,
					 const_wmesh_int_p 	s_q2n_n_,
					 const_wmesh_int_p 	s_q2n_v_,
					 const_wmesh_int_p 	s_q2n_ld_)
  {
    //
    // Get the number of dofs per entities.
    //
    static constexpr wmesh_int_t dimension = 3;
    
    wmesh_status_t 	status;

#if 0
    wmesh_int_t 	num_dofs[WFE_ELEMENT_ALL];
    wmesh_int_t 	num_faces[2];
    wmesh_int_t 	num_bfaces[2];
    
    wmesh_int_t 	c2n_ptr[4+1];
    wmesh_int_t 	c2e_ptr[4+1];
    wmesh_int_t 	c2f_ptr[4+1];
    
    wmesh_int_t 	c2n_m[4]  { 4, 5, 6, 8};
    wmesh_int_t 	c2e_m[4]  { 6, 8, 9, 12};
    wmesh_int_t 	c2f_m[4]  { 4, 5, 5, 6};
    
    wmesh_int_t 	c2n_n[4];
    wmesh_int_t 	c2e_n[4];      
    wmesh_int_t 	c2f_n[4];      
    
    wmesh_int_t 	c2n_ld[4];
    wmesh_int_t 	c2e_ld[4];      
    wmesh_int_t 	c2f_ld[4];      
      
    {
      for (int i=0;i<4;++i)
	c2n_ld[i] = c2n_m[i] + c2e_m[i] + c2f_m[i];
      
      for (int i=0;i<4;++i)
	c2e_ld[i] = c2n_ld[i];
      for (int i=0;i<4;++i)
	c2f_ld[i] = c2n_ld[i];

      for (int i=0;i<4;++i)
	c2n_n[i] = c2n_->m_n[i];
      for (int i=0;i<4;++i)
	c2e_n[i] = c2n_n[i];
      for (int i=0;i<4;++i)
	c2f_n[i] = c2n_n[i];

      c2n_ptr[0] = 0;
      for (int i=0;i<4;++i)
	c2n_ptr[i+1] = c2n_ptr[i] + c2n_ld[i] * c2n_n[i];
    }
    
    for (int i=0;i<4;++i)
      c2e_ptr[i] = c2n_ptr[i] + c2n_m[i];
    c2e_ptr[4] = c2n_ptr[4];
      
    for (int i=0;i<4;++i)
      c2f_ptr[i] = c2e_ptr[i] + c2e_m[i];
    c2f_ptr[4] = c2f_ptr[4];
    
    wmesh_int_t work_n;
    wmesh_int_t * work  = nullptr;
#endif
    
    //
    // Initialize the reference shape edges to nodes.
    //
    auto start = timing_start();
    
#define WINT_SPARSE_MAT_PARAMS(_f)		\
      _f.m_ptr,					\
	_f.m_m,					\
	_f.m_n,					\
	_f.m_data,				\
	_f.m_ld
    
    auto stop = timing_stop();
    std::cout << "neighbors detection, elapsed " << timing_seconds(start,stop) << std::endl;
    
    WMESH_STATUS_CHECK(status);
    
    wmesh_int_t num_edges = 0;
    //    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    wmesh_int_t edge_idx = 0;
    triangle_face_idx 		= 0;    
    quadrilateral_face_idx 	= 0;    
    
    status =  wmesh_indexing_edges(c2n_size_,
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
				   
				   &edge_idx,
				   work_n_,
				   work);
    
    WMESH_STATUS_CHECK(status);	  
    
    num_edges = edge_idx;


    status =  wmesh_indexing_triangles(c2n_size_,
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

				       &triangle_face_idx,
				       work_n_,
				       work_);
    WMESH_STATUS_CHECK(status);	  
    
    status = wmesh_indexing_quadrilaterals(4,
					   WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2n),
					   WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_c2f_q),
					   WINT_SPARSE_MAT_PARAMS(mesh_linear_->m_s_q2n),				     
					   &quadrilateral_face_idx,
					   work_n_,
					   work_);    
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
    return WMESH_STATUS_SUCCESS;
  };    

};

#endif
