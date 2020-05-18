#ifndef WMESH_FUNCTIONS_H
#define WMESH_FUNCTIONS_H

#include "wmesh-status.h"

#ifdef __cplusplus
extern "C"
{
#endif
  
  
  wmesh_status_t  wmesh_space_indexing_nodes(wmesh_int_t 	c2n_size_,
					     const_wmesh_int_p 	c2n_ptr_,
					     const_wmesh_int_p 	c2n_m_,
					     const_wmesh_int_p 	c2n_n_,
					     const_wmesh_int_p 	c2n_v_,
					     const_wmesh_int_p 	c2n_ld_,
					     
					     wmesh_int_t 	c2d_n_size_,
					     const_wmesh_int_p 	c2d_n_ptr_,
					     const_wmesh_int_p 	c2d_n_m_,
					     const_wmesh_int_p 	c2d_n_n_,
					     wmesh_int_p 	c2d_n_v_,
					     const_wmesh_int_p 	c2d_n_ld_,
					     
					     wmesh_int_t 	num_nodes_,
					     wmesh_int_t 	num_dofs_per_node_,
					     wmesh_int_t 	dof_idx_origin_);

  wmesh_status_t  wmesh_space_indexing_edges(wmesh_int_t 	c2n_size_,
					     const_wmesh_int_p 	c2n_ptr_,
					     const_wmesh_int_p 	c2n_m_,
					     const_wmesh_int_p 	c2n_n_,
					     const_wmesh_int_p 	c2n_v_,
					     const_wmesh_int_p 	c2n_ld_,
					     
					     wmesh_int_t 	c2e_size_,
					     const_wmesh_int_p 	c2e_ptr_,
					     const_wmesh_int_p 	c2e_m_,
					     const_wmesh_int_p 	c2e_n_,
					     const_wmesh_int_p 	c2e_v_,
					     const_wmesh_int_p 	c2e_ld_,
					     
					     wmesh_int_t 	c2d_e_size_,
					     const_wmesh_int_p 	c2d_e_ptr_,
					     const_wmesh_int_p 	c2d_e_m_,
					     const_wmesh_int_p 	c2d_e_n_,
					     wmesh_int_p 	c2d_e_v_,
					     const_wmesh_int_p 	c2d_e_ld_,

					     wmesh_int_t 	s_e2n_size_,
					     const_wmesh_int_p 	s_e2n_ptr_,
					     const_wmesh_int_p	s_e2n_m_,
					     const_wmesh_int_p	s_e2n_n_,
					     const_wmesh_int_p 	s_e2n_v_,
					     const_wmesh_int_p	s_e2n_ld_,
					     
					     wmesh_int_t 	num_edges_,
					     wmesh_int_t 	num_dofs_per_edge_,
					       wmesh_int_t 	dof_idx_origin_);
    wmesh_status_t  wmesh_space_indexing_triangles(wmesh_int_t 		c2n_size_,
						 const_wmesh_int_p 	c2n_ptr_,
						 const_wmesh_int_p 	c2n_m_,
						 const_wmesh_int_p 	c2n_n_,
						 const_wmesh_int_p 	c2n_v_,
						 const_wmesh_int_p 	c2n_ld_,
						 
						 wmesh_int_t 		c2f_t_size_,
						 const_wmesh_int_p 	c2f_t_ptr_,
						 const_wmesh_int_p 	c2f_t_m_,
						 const_wmesh_int_p 	c2f_t_n_,
						 const_wmesh_int_p 	c2f_t_v_,
						 const_wmesh_int_p 	c2f_t_ld_,
						 
						 wmesh_int_t 		c2d_t_size_,
						 const_wmesh_int_p 	c2d_t_ptr_,
						 const_wmesh_int_p 	c2d_t_m_,
						 const_wmesh_int_p 	c2d_t_n_,
						 wmesh_int_p 		c2d_t_v_,
						 const_wmesh_int_p 	c2d_t_ld_,

						 wmesh_int_t 		s_t2n_size_,
						 const_wmesh_int_p 	s_t2n_ptr_,
						 const_wmesh_int_p	s_t2n_m_,
						 const_wmesh_int_p	s_t2n_n_,
						 const_wmesh_int_p 	s_t2n_v_,
						 const_wmesh_int_p	s_t2n_ld_,
						      
						 wmesh_int_t 		num_triangles_,
						 wmesh_int_t 		degree_,
						 wmesh_int_t 		num_dofs_per_triangle_,
						   wmesh_int_t 		dof_idx_origin_);

  wmesh_status_t  wmesh_space_indexing_quadrilaterals(wmesh_int_t 		c2n_size_,					     
						      const_wmesh_int_p 	c2n_ptr_,
						      const_wmesh_int_p 	c2n_m_,
						      const_wmesh_int_p 	c2n_n_,
						      const_wmesh_int_p 	c2n_v_,
						      const_wmesh_int_p 	c2n_ld_,
						 
						      wmesh_int_t 		c2f_q_size_,					     
						      const_wmesh_int_p 	c2f_q_ptr_,
						      const_wmesh_int_p 	c2f_q_m_,
						      const_wmesh_int_p 	c2f_q_n_,
						      const_wmesh_int_p 	c2f_q_v_,
						      const_wmesh_int_p 	c2f_q_ld_,
						 
						      wmesh_int_t 		c2d_q_size_,					     
						      const_wmesh_int_p 	c2d_q_ptr_,
						      const_wmesh_int_p 	c2d_q_m_,
						      const_wmesh_int_p 	c2d_q_n_,
						      wmesh_int_p 		c2d_q_v_,
						      const_wmesh_int_p 	c2d_q_ld_,
						 
						      wmesh_int_t 		s_q2n_size_,					     
						      const_wmesh_int_p 	s_q2n_ptr_,
						      const_wmesh_int_p		s_q2n_m_,
						      const_wmesh_int_p		s_q2n_n_,
						      const_wmesh_int_p 	s_q2n_v_,
						      const_wmesh_int_p		s_q2n_ld_,
						      
						      wmesh_int_t 		num_quadrilaterals_,
						      wmesh_int_t 		degree_,
						      wmesh_int_t 		num_dofs_per_quadrilateral_,
						      wmesh_int_t 		dof_idx_origin_);  
  wmesh_status_t
  wmesh_space_indexing_interior
  (wmesh_int_t 		c2d_i_size_,						
   const_wmesh_int_p 	c2d_i_ptr_,
   const_wmesh_int_p 	c2d_i_m_,
   const_wmesh_int_p 	c2d_i_n_,
   wmesh_int_p 		c2d_i_v_,
   const_wmesh_int_p 	c2d_i_ld_,						
   wmesh_int_p 		dof_idx_origin_);
  
#ifdef __cplusplus
}
#endif

#endif
