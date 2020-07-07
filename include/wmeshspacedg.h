#ifndef WMESHSPACEDG_H
#define WMESHSPACEDG_H

#include "wmeshspace.h"

#ifdef __cplusplus
extern "C"
{
#endif

  struct wmeshspacedg_t;
  
  
  wmesh_status_t 	wmeshspacedg_sublinearmesh	(wmeshspacedg_t * 			self_,
							 wmesh_t ** 				mesh__);
  wmesh_status_t wmeshspacedg_def		(wmeshspacedg_t** 		self_,
						 wmesh_int_t 			nodes_family_,
						 wmesh_int_t 			degree_,
						 wmesh_t*			mesh_);
  
  wmesh_status_t wmeshspacedg_get_dofs_ids	(const wmeshspacedg_t * 	self_,
						 wmesh_int_t 			element_type_,
						 wmesh_int_t 			element_idx_,
						 wmesh_int_p 			dofs_,
						 wmesh_int_t 			dofs_inc_);
  
  wmesh_status_t wmeshspacedg_sparse	(const wmesh_t*__restrict__ 	self_,
					 wmesh_int_t 			degree_,
					 wmesh_int_p 			csr_size_,
					 wmesh_int_p*__restrict__ 	csr_ptr_,
					 wmesh_int_p*__restrict__ 	csr_ind_);

  wmesh_status_t wmeshspacedg_advection(const wmeshspacedg_t*__restrict__ 	self_,
					const wmesh_cubature_info_t* 		cubature_info_, 		
					const wmesh_shape_info_t* 		shape_info_element_, 		
					const wmesh_shape_info_t* 		shape_info_trial_,
					const wmesh_shape_info_t* 		shape_info_test_,
					const wmesh_shape_info_t* 		shape_info_velocity_,
					const wmeshspace_t *			velocity_space_,
					wmesh_int_t 				velocity_storage_,
					wmesh_int_t 				velocity_m_,
					wmesh_int_t 				velocity_n_,
					double *__restrict__ 				velocity_,
					wmesh_int_t 				velocity_ld_,
					wmesh_int_t				csr_size_,
					const_wmesh_int_p			csr_ptr_,
					const_wmesh_int_p			csr_ind_,
					double * 	__restrict__			csr_val_,
					double * 		__restrict__		rhs_);

#ifdef __cplusplus
}
#endif


#endif
