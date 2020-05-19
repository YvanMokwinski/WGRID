#ifndef BMS_H
#define BMS_H
#include <stdlib.h>
#include "wmesh-status.h"

#ifdef __cplusplus
extern "C"
{
#endif

  //!
  //! @brief Partitioning with space filling curve.
  //!
  wmesh_status_t bms_partitioning_sfc	(wmesh_int_t 				nparts_,
					 wmesh_int_t 				num_cells_,
					 wmesh_int_p 				p_v_,				      
					 wmesh_int_t  				p_ld_,				      
					 
					 wmesh_int_t 				c2n_size_,
					 const_wmesh_int_p 			c2n_ptr_,
					 const_wmesh_int_p 			c2n_m_,
					 const_wmesh_int_p 			c2n_n_,
					 const_wmesh_int_p 			c2n_v_,
					 const_wmesh_int_p 			c2n_ld_,
					 
					 wmesh_int_t 				coo_m_,
					 wmesh_int_t 				coo_n_,
					 const double *__restrict__ 		coo_v_,
					 wmesh_int_t 				coo_ld_,
					 
					 size_t*__restrict__                    work_size_,
					 void *__restrict__                     work_);
  
  
  //!
  //! @brief Compute the nodes-to-cells topological graph.
  //!
  wmesh_status_t bms_n2c		(wmesh_int_t 		c2n_size_,
					 const_wmesh_int_p 	c2n_ptr_,
					 const_wmesh_int_p 	c2n_m_,
					 const_wmesh_int_p 	c2n_n_,
					 const_wmesh_int_p 	c2n_v_,
					 const_wmesh_int_p 	c2n_ld_,				 
					 wmesh_int_p 		n2c_ptr_,
					 wmesh_int_t 		n2c_m_,
					 wmesh_int_p 		n2c_v_);
  
  wmesh_status_t bms_n2c_cindex		(wmesh_int_t 		c_,
					 wmesh_int_p 		cindex_);
  
  wmesh_status_t bms_n2c_ctype		(wmesh_int_t 		c_,
					 wmesh_int_p 		ctype_);
  
  //!
  //! @brief Get the size of the needed buffer.
  //!
  wmesh_status_t
  bms_c2c_e_buffer_size
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_n_,
   wmesh_int_p 		work_n_);
    
  //!
  //! @brief Build the cells-to-cells topological graph through edges.
  //! @remark The topology must be manifold.
  //!
  wmesh_status_t
  bms_c2c_e
  (wmesh_int_t		c2n_size_,   
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
		       
   wmesh_int_t 		c2c_size_,
   const_wmesh_int_p 	c2c_e_ptr_,
   const_wmesh_int_p 	c2c_e_m_,
   const_wmesh_int_p 	c2c_e_n_,
   wmesh_int_p 		c2c_e_v_,
   const_wmesh_int_p 	c2c_e_ld_,
		       
   wmesh_int_t 		s_e2n_size_,
   const_wmesh_int_p 	s_e2n_ptr_,
   const_wmesh_int_p	s_e2n_m_,
   const_wmesh_int_p	s_e2n_n_,
   const_wmesh_int_p 	s_e2n_v_,
   const_wmesh_int_p	s_e2n_ld_,
		       
   wmesh_int_p		num_boundary_edges_,
   wmesh_int_t		work_n_,
   wmesh_int_p		work_);

  //!
  //! @brief Get the size of the needed buffer.
  //!
  wmesh_status_t
  bms_c2c_t_buffer_size
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_n_,
   wmesh_int_p 		work_n_);
    
  //!
  //! @brief Build the cells-to-cells topological graph through triangle faces.
  //!
  wmesh_status_t
  bms_c2c_t
  (wmesh_int_t		c2n_size_,   
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
		       
   wmesh_int_t 		c2c_size_,
   const_wmesh_int_p 	c2c_t_ptr_,
   const_wmesh_int_p 	c2c_t_m_,
   const_wmesh_int_p 	c2c_t_n_,
   wmesh_int_p 		c2c_t_v_,
   const_wmesh_int_p 	c2c_t_ld_,
		       
   wmesh_int_t 		s_t2n_size_,
   const_wmesh_int_p 	s_t2n_ptr_,
   const_wmesh_int_p	s_t2n_m_,
   const_wmesh_int_p	s_t2n_n_,
   const_wmesh_int_p 	s_t2n_v_,
   const_wmesh_int_p	s_t2n_ld_,
		       
   wmesh_int_p		num_boundary_triangles_,
   wmesh_int_t		work_n_,
   wmesh_int_p		work_);
  
  //!
  //! @brief Get the size of the needed buffer.
  //!
  wmesh_status_t
  bms_c2c_q_buffer_size
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_n_,
   wmesh_int_p 		work_n_);
  
  //!
  //! @brief Build the cells-to-cells topological graph through triangle faces.
  //!
  wmesh_status_t  bms_c2c_q
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
   
   wmesh_int_t 		c2c_q_size_,
   const_wmesh_int_p 	c2c_q_ptr_,
   const_wmesh_int_p 	c2c_q_m_,
   const_wmesh_int_p 	c2c_q_n_,
   wmesh_int_p 		c2c_q_v_,
   const_wmesh_int_p 	c2c_q_ld_,
   
   wmesh_int_t 		s_q2n_size_,
   const_wmesh_int_p 	s_q2n_ptr_,
   const_wmesh_int_p	s_q2n_m_,
   const_wmesh_int_p	s_q2n_n_,
   const_wmesh_int_p 	s_q2n_v_,
   const_wmesh_int_p	s_q2n_ld_,
   
   wmesh_int_p		num_boundary_quadrilaterals_,
   wmesh_int_t		work_n_,
   wmesh_int_p		work_);
  
  //!
  //! @brief Build the cell-to-cell graph.
  //!
  wmesh_status_t
  bms_c2c_buffer_size
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_n_,
   wmesh_int_p 		work_n_);
  
  //!
  //! @brief Build the cells-to-cells topological graph.
  //!
  wmesh_status_t
  bms_c2c
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
   
   wmesh_int_t		c2c_t_size_,
   const_wmesh_int_p 	c2c_t_ptr_,
   const_wmesh_int_p 	c2c_t_m_,
   const_wmesh_int_p 	c2c_t_n_,
   wmesh_int_p		c2c_t_v_,
   const_wmesh_int_p 	c2c_t_ld_,
						 
   wmesh_int_t		c2c_q_size_,
   const_wmesh_int_p 	c2c_q_ptr_,
   const_wmesh_int_p 	c2c_q_m_,
   const_wmesh_int_p 	c2c_q_n_,
   wmesh_int_p		c2c_q_v_,
   const_wmesh_int_p 	c2c_q_ld_,
						 
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
   const_wmesh_int_p 	s_q2n_ld_,
						 
   wmesh_int_t 		work_n_,
   wmesh_int_p 		work_,
   wmesh_int_p 		num_faces_,
   wmesh_int_p 		num_bfaces_);




  
  
  //!
  //! @brief Get the size of the needed buffer.
  //!
  wmesh_status_t
  bms_c2e_buffer_size
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p	c2n_n_,
   wmesh_int_p 		work_n_);
  
  //!
  //! @brief Indexing of the edges.
  //!
  wmesh_status_t
  bms_c2e
  (wmesh_int_t		c2n_size_,   
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
				       
   wmesh_int_t		c2e_size_,   
   const_wmesh_int_p 	c2e_ptr_,
   const_wmesh_int_p 	c2e_m_,
   const_wmesh_int_p 	c2e_n_,
   wmesh_int_p 		c2e_v_,
   const_wmesh_int_p 	c2e_ld_,
				       
   wmesh_int_t		s_e2n_size_,   
   const_wmesh_int_p	s_e2n_ptr_,
   const_wmesh_int_p	s_e2n_m_,
   const_wmesh_int_p	s_e2n_n_,
   const_wmesh_int_p	s_e2n_v_,
   const_wmesh_int_p	s_e2n_ld_,

   wmesh_int_p		edge_idx_,
   wmesh_int_t		work_n_,
   wmesh_int_p		work_);


  //!
  //! @brief Get the size of the needed buffer.
  //!
  wmesh_status_t
  bms_c2f_t_buffer_size
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p	c2n_n_,
   wmesh_int_p 		work_n_);
  
  //!
  //! @brief Indexing of the triangle faces.
  //!
  wmesh_status_t
  bms_c2f_t
  (wmesh_int_t		c2n_size_,   
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
						 
   wmesh_int_t		c2t_size_,   
   const_wmesh_int_p 	c2t_ptr_,
   const_wmesh_int_p 	c2t_m_,
   const_wmesh_int_p 	c2t_n_,
   wmesh_int_p 		c2t_v_,
   const_wmesh_int_p 	c2t_ld_,
						 
   wmesh_int_t		s_t2n_size_,   
   const_wmesh_int_p	s_t2n_ptr_,
   const_wmesh_int_p	s_t2n_m_,
   const_wmesh_int_p	s_t2n_n_,
   const_wmesh_int_p	s_t2n_v_,
   const_wmesh_int_p	s_t2n_ld_,
						 
   wmesh_int_p		t_idx_,
   wmesh_int_t		work_n_,
   wmesh_int_p		work_);

  
  //!
  //! @brief Get the size of the needed buffer.
  //!
  wmesh_status_t
  bms_c2f_q_buffer_size
  (wmesh_int_t			c2n_size_,
   const_wmesh_int_p 		c2n_n_,						
   wmesh_int_t*__restrict__ 	size_);
  
  //!
  //! @brief Indexing of the quadrilateral faces.
  //!
  wmesh_status_t
  bms_c2f_q
  (wmesh_int_t		c2n_size_,
   const_wmesh_int_p 	c2n_ptr_,
   const_wmesh_int_p 	c2n_m_,
   const_wmesh_int_p 	c2n_n_,
   const_wmesh_int_p 	c2n_v_,
   const_wmesh_int_p 	c2n_ld_,
   
   wmesh_int_t		c2q_size_,   
   const_wmesh_int_p 	c2q_ptr_,
   const_wmesh_int_p 	c2q_m_,
   const_wmesh_int_p 	c2q_n_,
   wmesh_int_p 		c2q_v_,
   const_wmesh_int_p 	c2q_ld_,
   
   wmesh_int_t		s_q2n_size_,   
   const_wmesh_int_p 	s_q2n_ptr_,
   const_wmesh_int_p	s_q2n_m_,
   const_wmesh_int_p	s_q2n_n_,
   const_wmesh_int_p 	s_q2n_v_,
   const_wmesh_int_p	s_q2n_ld_,
						  
   wmesh_int_p		quadrilateral_idx_,
   wmesh_int_t		work_n_,
   wmesh_int_p		work_);


  
  wmesh_status_t
  bms_nodes_buffer_size(wmesh_int_t 	element_,
			wmesh_int_t 	family_,
			wmesh_int_t	degree_,
			wmesh_int_p	work_n_);
  
  wmesh_status_t
  bms_snodes(wmesh_int_t 	element_,
	     wmesh_int_t 	family_,
	     wmesh_int_t	degree_,
	     wmesh_int_t	b_storage_,
	     wmesh_int_t	b_m_,
	     wmesh_int_t	b_n_,
	     const_wmesh_int_p 	b_v_,
	     wmesh_int_t	b_ld_,
	     wmesh_int_t	c_storage_,
	     wmesh_int_t	c_m_,
	     wmesh_int_t	c_n_,
	     float* 		c_v_,
	     wmesh_int_t	c_ld_,
	     wmesh_int_t	work_n_,
	     float* 		work_);

  wmesh_status_t
  bms_dnodes(wmesh_int_t 	element_,
	     wmesh_int_t 	family_,
	     wmesh_int_t	degree_,
	     wmesh_int_t	b_storage_,
	     wmesh_int_t	b_m_,
	     wmesh_int_t	b_n_,
	     const_wmesh_int_p 	b_v_,
	     wmesh_int_t	b_ld_,
	     wmesh_int_t	c_storage_,
	     wmesh_int_t	c_m_,
	     wmesh_int_t	c_n_,
	     double* 		c_v_,
	     wmesh_int_t	c_ld_,
	     wmesh_int_t	work_n_,
	     double* 		work_);
  

  //
  // Polynomials.
  //
  
  wmesh_status_t bms_djacobip(wmesh_int_t 	alpha_,
			      wmesh_int_t 	beta_,
			      wmesh_int_t 	N_,
			      wmesh_int_t 	x_n_,
			      const double * 	x_,
			      wmesh_int_t  	x_ld_,
			      double * 		y_,
			      wmesh_int_t  	y_ld_,
			      wmesh_int_t 	work_n_,			   
			      double * 		work_);
  
  
  //!
  //! @brief I/O
  //!

  wmesh_status_t
  bms_read_medit_open
  (wmesh_str_t 		filename_,
   int64_t*		inm_,
   int32_t * 		version_,
   int32_t * 		dim_);
  
  wmesh_status_t
  bms_read_medit_stat
  (int64_t 		inm_,
   wmesh_int_p 		num_nodes_,
   wmesh_int_p 		num_edges_,
   wmesh_int_p 		num_triangles_,
   wmesh_int_p 		num_quadrilaterals_,
   wmesh_int_p 		num_tetrahedra_,
   wmesh_int_p 		num_pyramids_,
   wmesh_int_p 		num_wedges_,
   wmesh_int_p 		num_hexahedra_);

  wmesh_status_t
  bms_read_medit_topology
  (int64_t		inm_,
   
   wmesh_int_t 		c2n_size,
   const_wmesh_int_p 	c2n_ptr,
   const_wmesh_int_p 	c2n_m,
   const_wmesh_int_p 	c2n_n,
   wmesh_int_p 		c2n_v,
   const_wmesh_int_p 	c2n_ld,
			  
   wmesh_int_t 		c_c_size,
   const_wmesh_int_p 	c_c_ptr,
   const_wmesh_int_p 	c_c_m,
   const_wmesh_int_p 	c_c_n,
   wmesh_int_p 		c_c_v,
   const_wmesh_int_p 	c_c_ld);

  wmesh_status_t
  bms_read_medit_geometry
  (int64_t 		inm_,
   wmesh_int_t		coo_m_,
   wmesh_int_t		coo_n_,
   double *__restrict__ coo_,
   wmesh_int_t 		coo_ld_,
   wmesh_int_p  	nflags_,
   wmesh_int_t 		nflags_ld_);

  wmesh_status_t
  bms_write_medit_open
  (int64_t*		inm_,
   wmesh_str_t 		filename_,
   wmesh_int_t          precision_,
   wmesh_int_t          dimension_);
  
  wmesh_status_t
  bms_write_medit_topology
  (int64_t		inm_,
   
   wmesh_int_t 		c2n_size,
   const_wmesh_int_p 	c2n_ptr,
   const_wmesh_int_p 	c2n_m,
   const_wmesh_int_p 	c2n_n,
   const_wmesh_int_p 	c2n_v,
   const_wmesh_int_p 	c2n_ld,
   
   wmesh_int_t 		c_c_size,
   const_wmesh_int_p 	c_c_ptr,
   const_wmesh_int_p 	c_c_m,
   const_wmesh_int_p 	c_c_n,
   const_wmesh_int_p 	c_c_v,
   const_wmesh_int_p 	c_c_ld);
  
  
  wmesh_status_t
  bms_write_medit_geometry
  (int64_t 			inm_,
   wmesh_int_t			coo_m_,
   wmesh_int_t			coo_n_,
   const double *__restrict__ 	coo_,
   wmesh_int_t 			coo_ld_,
   const_wmesh_int_p  		nflags_,
   wmesh_int_t 			nflags_ld_);

  wmesh_status_t
  bms_medit_close
  (int64_t		inm_);
  

  //!
  //! @brief Get the required size for the needed buffer.
  //!
  wmesh_status_t
  bms_rmacro_buffer_size
  (wmesh_int_t 	element_,
   wmesh_int_t 	degree_,
   wmesh_int_p	work_n_,
   wmesh_int_p	num_entities_);
  
  wmesh_status_t
  bms_rmacro
  (wmesh_int_t 		element_,
   wmesh_int_t 		degree_,
   wmesh_int_t 		num_nodes_,
       
   wmesh_int_t 		c2n_size_,
   const_wmesh_int_p	c2n_ptr_,
   const_wmesh_int_p	c2n_m_,
   const_wmesh_int_p	c2n_n_,
   wmesh_int_p 		c2n_v_,
   const_wmesh_int_p	c2n_ld_,
				
   wmesh_int_t 		icoo_m_,
   wmesh_int_t 		icoo_n_,
   wmesh_int_p 		icoo_v_,
   wmesh_int_t 		icoo_ld_,
				
   wmesh_int_t		work_n_,
   wmesh_int_p		work_);


  wmesh_status_t bms_ordering_edge		(wmesh_int_t	degree_,
						 wmesh_int_t	c_storage_,
						 wmesh_int_t	c_m_,
						 wmesh_int_t	c_n_,
						 wmesh_int_p	c_v_,
						 wmesh_int_t	c_ld_);
  
  wmesh_status_t bms_ordering_triangle		(wmesh_int_t 	degree_,
						 wmesh_int_t 	c_storage_,
						 wmesh_int_t	c_m_,
						 wmesh_int_t	c_n_,
						 wmesh_int_p 	c_v_,
						 wmesh_int_t 	c_ld_,
						 wmesh_int_t 	shift_);
  
  wmesh_status_t bms_ordering_quadrilateral	(wmesh_int_t 	degree_,
						 wmesh_int_t 	c_storage_,
						 wmesh_int_t	c_m_,
						 wmesh_int_t	c_n_,
						 wmesh_int_p 	c_v_,
						 wmesh_int_t 	c_ld_,
						 wmesh_int_t 	shift_);
  
  wmesh_status_t bms_ordering_face	(wmesh_int_t		element_,
					 wmesh_int_t		degree_,
					 wmesh_int_t		c_storage_,
					 wmesh_int_t		c_m_,
					 wmesh_int_t		c_n_,
					 wmesh_int_p		c_v_,
					 wmesh_int_t		c_ld_);
  
  wmesh_status_t bms_ordering_volume	(wmesh_int_t		element_,
					 wmesh_int_t		degree_,
					 
					 wmesh_int_t		c_storage_,
					 wmesh_int_t		c_m_,
					 wmesh_int_t		c_n_,
					 wmesh_int_p		c_v_,
					 wmesh_int_t		c_ld_,
					 
					 wmesh_int_t 		s_e2n_m_,
					 wmesh_int_t 		s_e2n_n_,
					 const_wmesh_int_p 	s_e2n_v_,
					 wmesh_int_t 		s_e2n_ld_,
				     
					 wmesh_int_t 		s_t2n_m_,
					 wmesh_int_t 		s_t2n_n_,
					 const_wmesh_int_p 	s_t2n_v_,
					 wmesh_int_t 		s_t2n_ld_,
				     
					 wmesh_int_t 		s_q2n_m_,
					 wmesh_int_t 		s_q2n_n_,
					 const_wmesh_int_p 	s_q2n_v_,
					 wmesh_int_t 		s_q2n_ld_);

  wmesh_status_t bms_ordering_topoid(wmesh_int_t		element_,
				     wmesh_int_t		degree_,
				     wmesh_int_t		topoid_n_,
				     wmesh_int_p		topoid_v_,
				     wmesh_int_t		topoid_inc_);

  wmesh_status_t bms_ordering_flat(wmesh_int_t		degree_,
				   wmesh_int_t		b_storage_,
				   wmesh_int_t		b_m_,
				   wmesh_int_t		b_n_,
				   const_wmesh_int_p	b_v_,
				   wmesh_int_t		b_ld_,
				   wmesh_int_p		flat_);

  
  wmesh_status_t bms_ordering(wmesh_int_t		element_,
			      wmesh_int_t		degree_,
			      wmesh_int_t		c_storage_,
			      wmesh_int_t		c_m_,
			      wmesh_int_t		c_n_,
			      wmesh_int_p		c_v_,
			      wmesh_int_t		c_ld_);



  wmesh_status_t bms_s_e2n(wmesh_int_t 	dim_,
			   wmesh_int_p	s_size_,
			   wmesh_int_p	s_ptr_,
			   wmesh_int_p	s_m_,
			   wmesh_int_p	s_n_,
			   wmesh_int_p	s_v_,
			   wmesh_int_p	s_ld_);

  wmesh_status_t bms_s_q2n(wmesh_int_p	s_size_,
			   wmesh_int_p	s_ptr_,
			   wmesh_int_p	s_m_,
			   wmesh_int_p	s_n_,
			   wmesh_int_p	s_v_,
			   wmesh_int_p	s_ld_);
  
  wmesh_status_t bms_s_t2n(wmesh_int_p	s_size_,
			   wmesh_int_p	s_ptr_,
			   wmesh_int_p	s_m_,
			   wmesh_int_p	s_n_,
			   wmesh_int_p	s_v_,
			   wmesh_int_p	s_ld_);
    
#ifdef __cplusplus
}
#endif


#endif
