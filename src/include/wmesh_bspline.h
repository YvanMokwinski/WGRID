#ifndef WMESH_BSPLINE_H
#define WMESH_BSPLINE_H

#include "wmesh-functions.h"
#include "wmesh-status.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  struct wmesh_bspline_t;
  

  wmesh_status_t  wmesh_bspline_transform_coordinates(const wmesh_bspline_t*	self_,
						      wmesh_int_t		nv_,
						      double * 			coo_,
						      wmesh_int_t 		coo_ld_,
						      wmesh_int_p		ids_,
						      const_wmesh_int_p	 	offsetNode_);
  wmesh_status_t 	wmesh_bspline_seval			(const wmesh_bspline_t*		self_,
								 const float *__restrict__ 	parametric_position_,
								 wmesh_int_t			edge_idx_,
								 float*__restrict__ 		spatial_position_);

  wmesh_status_t 	wmesh_bspline_deval			(const wmesh_bspline_t*		self_,
								 const double *__restrict__ 	parametric_position_,
								 wmesh_int_t			edge_idx_,
								 double*__restrict__ 		spatial_position_);
  
  wmesh_status_t	wmesh_bspline_num_nodes			(const wmesh_bspline_t* self_,
								 wmesh_int_p 		out_num_nodes_);
  wmesh_status_t	wmesh_bspline_num_edges			(const wmesh_bspline_t* self_,
								 wmesh_int_p 		out_num_edges_);  
  wmesh_status_t	wmesh_bspline_dimension			(const wmesh_bspline_t* self_,
								 wmesh_int_p 		out_dimension_);  

  wmesh_status_t	wmesh_bspline_slength			(const wmesh_bspline_t* self_,
								 float*  		out_length_);

  wmesh_status_t	wmesh_bspline_dlength			(const wmesh_bspline_t* self_,
								 double* 		out_length_);

  wmesh_status_t 			wmesh_bspline_scoo			(const wmesh_bspline_t* 	self_,
								 wmesh_int_t		idx_,
								 float*__restrict__ 		coo_);
  wmesh_status_t 			wmesh_bspline_dcoo			(const wmesh_bspline_t* 	self_,
								 wmesh_int_t		idx_,
								 double*__restrict__ 		coo_);

  
  wmesh_status_t 			wmesh_bspline_stan			(const wmesh_bspline_t* 	self_,
										 wmesh_int_t			idx_,
										 float*__restrict__ 		tan_);
  wmesh_status_t 			wmesh_bspline_dtan			(const wmesh_bspline_t* 	self_,
										 wmesh_int_t			idx_,
										 double*__restrict__ 		tan_);
  
  wmesh_status_t 			wmesh_bspline_free			(wmesh_bspline_t*self_);
  wmesh_status_t  wmesh_bspline_write_medit(const wmesh_bspline_t*	bspline_,
					    const char * 			filename_);

  wmesh_status_t 			wmesh_bspline_ddef			(wmesh_bspline_t** 		self_,
										 wmesh_int_t  			dim_,
										 wmesh_int_t 			numPoints_,
										 const double *  		__restrict__ coo_,
										 wmesh_int_t 			coo_ld_);
  
  void 					wmesh_bspline_sdef			(wmesh_bspline_t* 		self_,
										 wmesh_int_t  			dim_,
										 wmesh_int_t 			numPoints_,
										 const double *  		coo_,
										 wmesh_int_t 			coo_ld_);
  
  wmesh_status_t 			wmesh_bspline_dnew			(wmesh_int_t  			dim_,
										 wmesh_int_t 			n_,
										 const double *__restrict__     coo_,
										 wmesh_int_t 			coo_ld_);

  wmesh_status_t 			wmesh_bspline_snew			(wmesh_int_t  			dim_,
										 wmesh_int_t 			n_,
										 const float *__restrict__     	coo_,
										 wmesh_int_t 			coo_ld_);

#if 0
  void 			wmesh_bspline_get_dofGeometry		(const wmesh_bspline_t* 	const 	self_,
								 cst_pI 		const	edgeIndex_,
								 pR 			const	dofGeometry_,
								 cst_pI			const 	dofGeometryOff_);

  void 			wmesh_bspline_transform_coordinates	(const wmesh_bspline_t* 	const	self_,
								 pMeshNode 		const 	meshNode_,
								 cst_pI			const 	offsetNode_);
#endif
  
#ifdef __cplusplus
}
#endif

#endif
