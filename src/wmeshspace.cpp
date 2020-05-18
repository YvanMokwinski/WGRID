
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "wmesh_utils.hpp"
#include "bms.h"

#include <chrono>
#include <iostream>

extern "C"
{
  static wmesh_status_t wmesh_init_c2d(const wmesh_int_sparsemat_t*	c2n_,
				       wmesh_int_sparsemat_t*		c2d_,
				       const_wmesh_int_p		c2d_m_,
				       const_wmesh_int_p		c2d_ld_)
  {

    wmesh_int_t c2d_ptr[4+1];
    wmesh_int_t c2d_n[4];      
    wmesh_int_t c2d_size = c2n_->m_size;
    for (wmesh_int_t i=0;i<c2d_size;++i) c2d_n[i] = c2n_->m_n[i];
    c2d_ptr[0]=0;
    for (wmesh_int_t i=0;i<c2d_size;++i) c2d_ptr[i+1] = c2d_ptr[i] + c2d_n[i] * c2d_ld_[i];
    
    wmesh_int_t * __restrict__ c2d_v = (wmesh_int_t * __restrict__)malloc(sizeof(wmesh_int_t)*c2d_ptr[c2d_size]);
    if (!c2d_v)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    return wmesh_int_sparsemat_new(c2d_,
				   c2d_size,
				   c2d_ptr,
				   c2d_m_,
				   c2d_n,
				   c2d_v,
				   c2d_ld_);
  };


  static wmesh_status_t wmesh_init_c2d_n(wmesh_int_t 		topodim_,
					 wmesh_int_sparsemat_t*	c2d_,
					 wmesh_int_sparsemat_t*	c2d_n_,
					 const_wmesh_int_p 	num_ldofs_,
					 wmesh_int_p 		shifts_)
  {
    wmesh_status_t status;
    wmesh_int_t c2d_size;
    wmesh_int_t elements[4];
    wmesh_int_t c2d_n_m[4];//  {4,5,6,8};
    wmesh_int_t c2d_n_ptr[4+1];
    status = wmesh_topodim2elements(topodim_,
				    &c2d_size,
				    elements);
    WMESH_STATUS_CHECK(status);    
    
    status = wmesh_elements_num_nodes(c2d_size,
				      elements,
				      c2d_n_m);
    WMESH_STATUS_CHECK(status);
    
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
    
    return wmesh_int_sparsemat_new(c2d_n_,
				   c2d_size,
				   c2d_n_ptr,
				   c2d_n_m,
				   c2d_->m_n,
				   c2d_->m_data,
				   c2d_->m_ld);
  }
  
  static wmesh_status_t wmesh_init_c2d_e(wmesh_int_t 		topodim_,
					 wmesh_int_sparsemat_t*	c2d_,
					 wmesh_int_sparsemat_t*	c2d_e_,
					 const_wmesh_int_p 	num_ldofs_,
					 wmesh_int_p 		shifts_)
  {
    wmesh_status_t status;
    wmesh_int_t c2d_e_size;
    wmesh_int_t elements[4];
    wmesh_int_t c2d_e_m[4];
    wmesh_int_t c2d_e_ptr[4+1];
    status = wmesh_topodim2elements(topodim_,
				    &c2d_e_size,
				    elements);
    WMESH_STATUS_CHECK(status);    
    
    status = wmesh_elements_num_edges(c2d_e_size,
				      elements,
				      c2d_e_m);
    WMESH_STATUS_CHECK(status);

    for (wmesh_int_t i=0;i<c2d_e_size;++i)
      {	
	c2d_e_m[i] *= num_ldofs_[i];
	std::cout << c2d_e_m[i] << std::endl;
      }
    
    for (wmesh_int_t i=0;i<c2d_e_size;++i)
      c2d_e_ptr[i] = c2d_->m_ptr[i] + shifts_[i]; 

    c2d_e_ptr[c2d_e_size] = c2d_->m_ptr[c2d_e_size];

    for (wmesh_int_t i=0;i<c2d_e_size;++i)
      shifts_[i] += c2d_e_m[i];

    return wmesh_int_sparsemat_new(c2d_e_,
				   c2d_e_size,
				   c2d_e_ptr,
				   c2d_e_m,
				   c2d_->m_n,
				   c2d_->m_data,
				   c2d_->m_ld);
  }
  
  static wmesh_status_t wmesh_init_c2d_t(wmesh_int_sparsemat_t*	c2d_,
					 wmesh_int_sparsemat_t*	c2d_t_,
					 const_wmesh_int_p 	num_ldofs_,
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

    return wmesh_int_sparsemat_new(c2d_t_,
				   4,
				   c2d_t_ptr,
				   c2d_t_m,
				   c2d_->m_n,
				   c2d_->m_data,
				   c2d_->m_ld);
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

    return wmesh_int_sparsemat_new(c2d_q_,
				   4,
				   c2d_q_ptr,
				   c2d_q_m,
				   c2d_->m_n,
				   c2d_->m_data,
				   c2d_->m_ld);
    return WMESH_STATUS_SUCCESS;
  }

  static wmesh_status_t wmesh_init_c2d_i(wmesh_int_t 			topodim_,
					 wmesh_int_sparsemat_t*		c2d_,
					 wmesh_int_sparsemat_t*		c2d_i_,
					 const_wmesh_int_p 		num_ldofs_,
					 wmesh_int_p 			shifts_)
  {
    wmesh_status_t status;
    wmesh_int_t c2d_i_m[4]  {1,1,1,1};
    wmesh_int_t c2d_i_ptr[4+1];
    wmesh_int_t c2d_i_size;
    
    status = wmesh_topodim2numtypes(topodim_,
				    &c2d_i_size);
    WMESH_STATUS_CHECK(status);
    
    for (wmesh_int_t i=0;i<c2d_i_size;++i)
      {
	c2d_i_m[i] *= num_ldofs_[i];
      }	

    for (wmesh_int_t i=0;i<c2d_i_size;++i)
      {
	c2d_i_ptr[i] = c2d_->m_ptr[i] + shifts_[i];
      }
    c2d_i_ptr[c2d_i_size] = c2d_->m_ptr[c2d_i_size];
    
    for (wmesh_int_t i=0;i<c2d_i_size;++i)
      {
	shifts_[i] += c2d_i_m[i];
      }
    
    return wmesh_int_sparsemat_new(c2d_i_,
				   c2d_i_size,
				   c2d_i_ptr,
				   c2d_i_m,
				   c2d_->m_n,
				   c2d_->m_data,
				   c2d_->m_ld);
  }

  static wmesh_status_t wmeshspace_compute(wmeshspace_t * 	space_)
  {
    wmesh_t*		self_ 	= space_->m_mesh;
    const wmesh_int_t topodim = self_->m_topology_dimension;
    wmesh_int_t		degree_ = space_->m_degree;    
    wmesh_status_t 	status;
    const wmesh_int_t num_dofs_per_node 		= (degree_ > 0) ? 1 : 0;      
    const wmesh_int_t num_dofs_per_edge 		= (degree_>0) ? degree_-1 : 0;
    const wmesh_int_t num_dofs_per_triangle		= (degree_>0) ? ((degree_-1)*(degree_-2))/2 : 0;
    const wmesh_int_t num_dofs_per_quadrilateral 	= (degree_>0) ? (degree_-1)*(degree_-1) : 0;    

    wmesh_int_t	dof_idx = 1;
    
    if (num_dofs_per_node > 0)
      {
	status = wmesh_space_indexing_nodes(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
					    WMESH_INT_SPARSEMAT_FORWARD(space_->m_c2d_n),
					    self_->m_num_nodes,
					    num_dofs_per_node,
					    dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += self_->m_num_nodes * num_dofs_per_node;	  
      }

    if (num_dofs_per_edge > 0)
      {
	status =  wmesh_space_indexing_edges(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
					     WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2e),
					     WMESH_INT_SPARSEMAT_FORWARD(space_->m_c2d_e),
					     WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_e2n),
					     self_->m_num_edges,
					     num_dofs_per_edge,
					     dof_idx);
	WMESH_STATUS_CHECK(status);
	dof_idx += self_->m_num_edges * num_dofs_per_edge;
      }

    if (topodim > 2)
      {
	if (num_dofs_per_triangle > 0)
	  {
	    status =  wmesh_space_indexing_triangles(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
						     WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2f_t),
						     WMESH_INT_SPARSEMAT_FORWARD(space_->m_c2d_t),
						     WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_t2n),
						     self_->m_num_triangles,
						     degree_,
						     num_dofs_per_triangle,
						     dof_idx);
	    WMESH_STATUS_CHECK(status);
	    dof_idx += self_->m_num_triangles * num_dofs_per_triangle;
	  }
	
	if (num_dofs_per_quadrilateral > 0)
	  {
	    status = wmesh_space_indexing_quadrilaterals(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
							 WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2f_q),
							 WMESH_INT_SPARSEMAT_FORWARD(space_->m_c2d_q),
							 WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_q2n),
							 self_->m_num_quadrilaterals,
							 degree_,
							 num_dofs_per_quadrilateral,
							 dof_idx);
	    WMESH_STATUS_CHECK(status);
	    dof_idx += self_->m_num_quadrilaterals * num_dofs_per_quadrilateral;
	  }
      }

    status =  wmesh_space_indexing_interior(WMESH_INT_SPARSEMAT_FORWARD(space_->m_c2d_i),
					    &dof_idx);
    WMESH_STATUS_CHECK(status);
#if 0
    printf("c2d_i \n");
    wmesh_int_sparsemat_fprintf(&space_->m_c2d_i,
				stdout);
    printf("c2d_f_t \n");
    wmesh_int_sparsemat_fprintf(&space_->m_c2d_t,
				stdout);
   printf("c2d_e \n");
    wmesh_int_sparsemat_fprintf(&space_->m_c2d_e,
				stdout);
   printf("c2d_n \n");
    wmesh_int_sparsemat_fprintf(&space_->m_c2d_n,
				stdout);
    wmesh_int_sparsemat_fprintf(&space_->m_c2d,
				stdout);
#endif
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
    
    wmesh_status_t status;
    if (mesh_->m_c2n.m_size == 4)
      {
	for (wmesh_int_t l=0;l<mesh_->m_c2n.m_size;++l)
	  {
	    status = wmesh_rmacro_def(&self_->m_patterns[l],
				      4+l,
				      degree_);
	    WMESH_STATUS_CHECK(status);
	  }
      }
    else if (mesh_->m_c2n.m_size == 2)
      {
	for (wmesh_int_t l=0;l<mesh_->m_c2n.m_size;++l)
	  {
	    status = wmesh_rmacro_def(&self_->m_patterns[l],
				      2+l,
				      degree_);
	    WMESH_STATUS_CHECK(status);
	  }
      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
      }
    WMESH_STATUS_CHECK(status);


    
    wmesh_int_t topodim = mesh_->m_topology_dimension;


    wmesh_int_t degree = degree_;

    
    if (topodim==3)
      {
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
	status = wmesh_init_c2d_n	(topodim,
					 &self_->m_c2d,
					 &self_->m_c2d_n,
					 ndofs_n,
					 shifts);
	WMESH_STATUS_CHECK(status);

	//
	// edge-based dofs.
	//
	status = wmesh_init_c2d_e	(topodim,
					 &self_->m_c2d,
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
	status = wmesh_init_c2d_i(topodim,
				  &self_->m_c2d,
				  &self_->m_c2d_i,
				  ndofs_i,
				  shifts);
	WMESH_STATUS_CHECK(status);
      }
    else if (topodim==2)
      {
	wmesh_int_t ndofs[2]= { ((degree+1)*(degree+2)) / 2,
				((degree+1)*(degree+1))};

	wmesh_int_t ndofs_n[2] = {1,1};
	wmesh_int_t ndofs_e[2] = {degree-1,degree-1};
	wmesh_int_t ndofs_i[2] = { (degree>0) ? ((degree-1)*(degree-2)) / 2 : 1,
				   (degree>0) ? (degree-1)*(degree-1) : 1 };	
	wmesh_int_t shifts[2] = {0,0};
	
	status = wmesh_init_c2d	(&self_->m_mesh->m_c2n,
				 &self_->m_c2d,
				 ndofs,
				 ndofs);
	WMESH_STATUS_CHECK(status);

	
	//
	// Node-based dofs.
	//
	status = wmesh_init_c2d_n(topodim,
				  &self_->m_c2d,
				  &self_->m_c2d_n,
				  ndofs_n,
				  shifts);
	WMESH_STATUS_CHECK(status);

	//
	// edge-based dofs.
	//
	status = wmesh_init_c2d_e	(topodim,
					 &self_->m_c2d,
					 &self_->m_c2d_e,
					 ndofs_e,
					 shifts);
	WMESH_STATUS_CHECK(status);
	
	//
	// cell-based dofs.
	//
	status = wmesh_init_c2d_i(topodim,
				  &self_->m_c2d,
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
  
};
