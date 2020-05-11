#ifndef WMESH_H
#define WMESH_H

#include "wmesh-functions.h"

#define WMESH_ELEMENT_NODE 		0
#define WMESH_ELEMENT_EDGE 		1
#define WMESH_ELEMENT_TRIANGLE 		2
#define WMESH_ELEMENT_QUADRILATERAL 	3
#define WMESH_ELEMENT_TETRAHEDRON 	4
#define WMESH_ELEMENT_PYRAMID 		5
#define WMESH_ELEMENT_WEDGE 		6
#define WMESH_ELEMENT_HEXAHEDRON	7
#define WMESH_ELEMENT_ALL		8


#ifdef __cplusplus
extern "C"
{
#endif
  
  struct wmesh_t;
  struct wmeshspace_t;
  struct wmesh_bspline_t;


  
  wmesh_status_t wmeshspace_def			(wmeshspace_t ** 	self__,
						 wmesh_int_t 		degree_,
						 wmesh_t * 		mesh_);
  
  wmesh_status_t wmeshspace_sublinearmesh	(wmeshspace_t * 	self_,
						 wmesh_t ** 		mesh__);
  

  wmesh_status_t wmesh_extract_boundary		(wmesh_t* self_);
  wmesh_status_t wmesh_kill			(wmesh_t* self_);
  wmesh_status_t wmesh_def_polar_extrusion	(wmesh_t ** self_,
						 wmesh_int_t 		dim_,
						 const_wmesh_int_p 	n_,
						 const double * 	x_,
						 const_wmesh_int_p  	nbRotations_);
  
  wmesh_status_t wmesh_extrusion		(wmesh_t ** 			self_,
						 const wmesh_t * 		surface_,
						 wmesh_int_t 			nz_,
						 wmesh_int_t 			ndz_,
						 const double*__restrict__	dz_,
						 const_wmesh_int_p		bfaces_ids);
  
  wmesh_status_t  wmesh_spline_extrusion(wmesh_t *x,
					 const wmesh_t *	surface_,
					 const wmesh_bspline_t * spline);
  
  wmesh_status_t  wmesh_curve_extrusion		(wmesh_t *		self_,
						 const wmesh_t *	surface_,
						 const wmesh_t * 	curve_);
  
  
  wmesh_status_t wmesh_treilli_buffer_size	(wmesh_int_t 	cell_type_,
						 wmesh_int_t 	degree_,
						 wmesh_int_p	work_n_,
						 wmesh_int_p	num_entities_);
  
  
  wmesh_status_t wmesh_treilli_calculate	(wmesh_int_t 		cell_type_,
						 wmesh_int_t 		degree_,
						 wmesh_int_t 		num_nodes_,
						     
						 const_wmesh_int_p 	c2n_ptr_,
						 const_wmesh_int_p 	c2n_m_,
						 const_wmesh_int_p 	c2n_n_,
						 wmesh_int_p 		c2n_v_,
						 const_wmesh_int_p 	c2n_ld_,
						 
						 wmesh_int_t 		icoo_m_,
						 wmesh_int_t 		icoo_n_,
						 wmesh_int_p 		icoo_v_,
						 wmesh_int_t 		icoo_ld_,
						 
						 wmesh_int_t		work_n_,
						 wmesh_int_p		work_);


    
  //!
  //! @brief Kill the mesh.
  //!
  wmesh_status_t wmesh_rmacro_def(wmesh_t ** 	mesh__,
				  wmesh_int_t 	element_,
				  wmesh_int_t 	degree_);

  wmesh_status_t wmesh_refine			(wmesh_t*    self_,			      
						 wmesh_int_t degree,
						 wmesh_t**   refined_mesh_);
  
  wmesh_status_t wmesh_fespace_endomorphism	(const wmesh_t*	mesh_,
						 wmesh_int_t 	degree_,
						 wmesh_int_p 	csr_size_,
						 wmesh_int_p* 	csr_ptr_,
						 wmesh_int_p* 	csr_ind_);
  
  wmesh_status_t wmesh_fespace			(wmesh_t**	self_,
						 const wmesh_t*	mesh_,
						 wmesh_int_t 	degree_);

  wmesh_status_t wmesh_partitioning		(wmesh_t*	self_,
						 wmesh_int_t 	nparts_);


  wmesh_status_t wmesh_treilli(wmesh_t ** 	mesh__,
			       wmesh_int_t 	cell_type_,
			       wmesh_int_t 	degree_,
			       wmesh_int_t	work_n_,
			       wmesh_int_p	work_);

  //!
  //! @brief Kill the mesh.
  //!
  wmesh_status_t wmesh_reorder(wmesh_t*self_);
  
  //!
  //! @brief Get the number of entities with a specific dimension.
  //!
  wmesh_status_t wmesh_analysis(wmesh_t*		self_);
  
  //!
  //! @brief Get the number of entities with a specific dimension.
  //!
  wmesh_status_t wmesh_info	(const wmesh_t*		self_,
				 FILE * out_);
  
  wmesh_status_t wmesh_factory	(wmesh_t** 		self__,
				 wmesh_int_t 		topology_dimension_,
				 
				 wmesh_int_t 		c2n_size_,
				 const_wmesh_int_p 	c2n_ptr_,
				 const_wmesh_int_p 	c2n_m_,
				 const_wmesh_int_p 	c2n_n_,
				 wmesh_int_p 		c2n_v_,
				 const_wmesh_int_p	c2n_ld_,

				 wmesh_int_t 		c_c_size_,
				 const_wmesh_int_p 	c_c_ptr_,
				 const_wmesh_int_p 	c_c_m_,
				 const_wmesh_int_p 	c_c_n_,
				 wmesh_int_p 		c_c_v_,
				 const_wmesh_int_p	c_c_ld_,
				 
				 wmesh_int_t 		bf2n_size_,
				 const_wmesh_int_p 	bf2n_ptr_,
				 const_wmesh_int_p 	bf2n_m_,
				 const_wmesh_int_p 	bf2n_n_,
				 wmesh_int_p 		bf2n_v_,
				 const_wmesh_int_p	bf2n_ld_,
				 
				 wmesh_int_t 		bf_c_size_,
				 const_wmesh_int_p 	bf_c_ptr_,
				 const_wmesh_int_p 	bf_c_m_,
				 const_wmesh_int_p 	bf_c_n_,
				 wmesh_int_p 		bf_c_v_,
				 const_wmesh_int_p	bf_c_ld_,
				 
				 wmesh_int_t		coo_m_,
				 wmesh_int_t		coo_n_,
				 double * 		coo_v_,
				 wmesh_int_t 		coo_ld_,
				 
				 wmesh_int_p 		n_c_v_,
				 wmesh_int_t		n_c_ld_);

  //!
  //! @brief Get the number of entities with a specific dimension.
  //!
  wmesh_status_t wmesh_def	(wmesh_t** 		self_,
				 wmesh_int_t 		topology_dimension_,
				 wmesh_int_t 		c2n_size_,
				 const_wmesh_int_p 	c2n_ptr_,
				 const_wmesh_int_p 	c2n_m_,
				 const_wmesh_int_p 	c2n_n_,
				 wmesh_int_p 		c2n_v_,
				 const_wmesh_int_p	c2n_ld_,

				 wmesh_int_t		coo_m_,
				 wmesh_int_t		coo_n_,
				 double * 		coo_,
				 wmesh_int_t 		coo_ld_);

  //!
  //! @brief Get the number of entities with a specific dimension.
  //!
  wmesh_status_t wmesh_read	(wmesh_t ** 		self_,
				 const char * 		filename_);
  
  //!
  //! @brief Get the number of entities with a specific dimension.
  //!
  wmesh_status_t wmesh_write	(const wmesh_t * 	self_,
				 const char * 	filename_);
  wmesh_status_t wmesh_write_vtk(const wmesh_t* 		self_,
				 const char * 		filename_,
				 ...);
  
  
  //!
  //! @brief Get the number of entities with a specific dimension.
  //!
  wmesh_status_t wmesh_num_entities		(const wmesh_t* 	self_,
						 wmesh_int_t*	out_num_) ;

  wmesh_status_t wmesh_entity_idx		(const wmesh_t* 	self_,
						 wmesh_int_t	entity_id_,
						 wmesh_int_t*	out_entity_idx_); 
  
  wmesh_status_t wmesh_entity_id		(const wmesh_t* 	self_,
						 wmesh_int_t	entity_idx_,
						 wmesh_int_t*	out_entity_id_); 
  
  wmesh_status_t wmesh_entity2nodes		(const wmesh_t* 	self_,
						 wmesh_int_t	entity_id_,
						 wmesh_int_t*	ent2n_n_,
						 wmesh_int_t*	ent2n_);
  
  wmesh_status_t wmesh_entity2edges		(const wmesh_t* 	self_,
						 wmesh_int_t	entity_id_,
						 wmesh_int_t*	ent2e_n_,
						 wmesh_int_t*	ent2e_);
  
  wmesh_status_t wmesh_entity2triangles	(const wmesh_t* 	self_,
						 wmesh_int_t	entity_id_,
						 wmesh_int_t*	ent2tr_n_,
						 wmesh_int_t*	ent2tr);
  
  wmesh_status_t wmesh_entity2quadrilaterals	(const wmesh_t* 	self_,
						 wmesh_int_t	entity_id_,
						 wmesh_int_t*	ent2tr_n_,
						 wmesh_int_t*	ent2tr);



#ifdef __cplusplus
}
#endif


#endif