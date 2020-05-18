#ifndef WMESH_NODES_H
#define WMESH_NODES_H

#include "wmesh_nodes_family.h"

#ifdef __cplusplus
extern "C"
{
#endif
  typedef wmesh_int_t wmesh_nodes_t;
  
  wmesh_int_t wmesh_nodes_dim		(const wmesh_nodes_t * 	self_);
  wmesh_int_t wmesh_nodes_ndofs		(const wmesh_nodes_t *  self_);
  wmesh_int_t wmesh_nodes_degree	(const wmesh_nodes_t *  self_);
  wmesh_int_t wmesh_nodes_family	(const wmesh_nodes_t *  self_);
  wmesh_int_t wmesh_nodes_element	(const wmesh_nodes_t *  self_);
  wmesh_int_t wmesh_nodes_def		(wmesh_nodes_t *  	self_,
					 wmesh_int_t 		element_,
					 wmesh_int_t 		family_,
					 wmesh_int_t 		degree_);


#ifdef WMESH_NODES_INTERFACE_PARAMS
#error WMESH_NODES_INTERFACE_PARAMS already defined
#endif
  
#define WMESH_NODES_INTERFACE_PARAMS(_T)				\
  const wmesh_nodes_t*__restrict__ self_,					\
    wmesh_int_t 		storagec_,				\
    wmesh_int_t 		mc_,					\
    wmesh_int_t 		nc_,					\
    _T*__restrict__ 		c_,					\
    wmesh_int_t 		ldc_,					\
    wmesh_int_t 		work_n_,				\
    _T*__restrict__ 		work_
  
  
#ifdef WMESH_NODES_INTERFACE_FORWARD
#error WMESH_NODES_INTERFACE_FORWARD already defined
#endif
  
#define WMESH_NODES_INTERFACE_FORWARD self_,	\
    storagec_,					\
    mc_,					\
    nc_,					\
    c_,						\
    ldc_,					\
    work_n_,					\
    work_
  
  wmesh_status_t wmesh_snodes			(WMESH_NODES_INTERFACE_PARAMS(float));
  wmesh_status_t wmesh_dnodes			(WMESH_NODES_INTERFACE_PARAMS(double));
  
#ifdef __cplusplus
};
#endif

#endif

