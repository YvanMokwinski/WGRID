#ifndef WMESHSPACE_H
#define WMESHSPACE_H
#include "wmesh.h"
#ifdef __cplusplus
extern "C"
{
#endif
  
  struct wmeshspace_t;
  wmesh_status_t wmeshspace_get_dofs_ids	(const wmeshspace_t * 	self_,
						 wmesh_int_t 		element_type_,
						 wmesh_int_t 		element_idx_,
						 wmesh_int_p 		dofs_,
						 wmesh_int_t 		dofs_inc_);

  wmesh_status_t wmeshspace_generate_dcoodofs	(const wmeshspace_t * 	self_,
						 wmesh_int_t 		coo_storage_,
						 wmesh_int_t 		coo_m_,
						 wmesh_int_t 		coo_n_,
						 double * 		coo_,
						 wmesh_int_t 		coo_ld_);

  wmesh_status_t wmeshspace_generate_scoodofs	(const wmeshspace_t * 	self_,
						 wmesh_int_t 		coo_storage_,
						 wmesh_int_t 		coo_m_,
						 wmesh_int_t 		coo_n_,
						 float * 		coo_,
						 wmesh_int_t 		coo_ld_);
  
  wmesh_status_t 	wmeshspace_get_ndofs(const wmeshspace_t*__restrict__ 	self_,
					     wmesh_int_p 				ndofs_);
  

    
  wmesh_status_t 	wmeshspace_def			(wmeshspace_t ** 			self__,
							 wmesh_int_t 				family_,
							 wmesh_int_t 				degree_,
							 wmesh_t * 				mesh_);

  
  wmesh_status_t 	wmeshspace_sublinearmesh	(wmeshspace_t * 			self_,
							 wmesh_t ** 				mesh__);

  wmesh_status_t 	wmeshspace_sparse		(const wmeshspace_t*__restrict__ 	self_,
							 wmesh_int_p 				csr_size_,
							 wmesh_int_p*__restrict__ 		csr_ptr_,
							 wmesh_int_p*__restrict__ 		csr_ind_);

  
  wmesh_status_t wmeshspace_laplace_old(const wmeshspace_t*__restrict__ 	self_,
				    wmesh_int_t				csr_size_,
				    const_wmesh_int_p			csr_ptr_,
				    const_wmesh_int_p			csr_ind_,
				    double * 				csr_val_,
					double * 				rhs_);
  
    wmesh_status_t wmeshspace_laplace(const wmeshspace_t*__restrict__ 	self_,
				      const wmesh_cubature_info_t* 	cubature_info_, 		
				      const wmesh_shape_info_t* 		shape_info_element_, 		
				      const wmesh_shape_info_t* 		shape_info_trial_,
				    const wmesh_shape_info_t* 		shape_info_test_,
				    const wmesh_shape_info_t* 		shape_info_a_,
				    wmesh_int_t				csr_size_,
				    const_wmesh_int_p			csr_ptr_,
				    const_wmesh_int_p			csr_ind_,
				    double * 				csr_val_,
				      double * 				rhs_);

#ifdef __cplusplus
}
#endif


#endif
