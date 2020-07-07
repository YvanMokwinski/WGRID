#include "wmesh.hpp"
#include "bms.h"
#include <iostream>
#include <string.h>



wmesh_status_t wmeshspacedg_get_dofs_ids(const wmeshspacedg_t * 	self_,
					 wmesh_int_t 		element_type_,
					 wmesh_int_t 		element_idx_,
					 wmesh_int_p 		dofs_,
					 wmesh_int_t 		dofs_inc_)
{
  const wmesh_int_t num_dofs = self_->m_dofs_m[element_type_];
  const wmesh_int_t shift = self_->m_dofs_ptr[element_type_];
  for (wmesh_int_t i=0;i<num_dofs;++i)
    {
      dofs_[dofs_inc_*i] = 1 + shift  + element_idx_ * num_dofs + i;
    }
  return WMESH_STATUS_SUCCESS;
}


wmesh_status_t wmeshspacedg_def(wmeshspacedg_t ** 		self__,
				wmesh_int_t 			nodes_family_,
				wmesh_int_t 			degree_,
				wmesh_t*			mesh_)
{
  WMESH_CHECK_POINTER(self__);
  WMESH_CHECK_POINTER(mesh_);
  wmeshspacedg_t *self_ = (wmeshspacedg_t *)malloc(sizeof(wmeshspacedg_t));
  self__[0] = self_;
  memset(self_,0,sizeof(wmeshspacedg_t));
  self_->m_mesh 		= mesh_;
  const wmesh_int_t topodim 	= mesh_->m_topology_dimension;
  self_->m_nodes_family = nodes_family_;
  self_->m_degree = degree_;
  self_->m_dofs_ptr[0] 	= 0;
  for (wmesh_int_t l=0;l<mesh_->m_c2n.m_size;++l)
    {
      const wmesh_int_t element = l + ( (topodim==3) ? 4 : ( (topodim==2) ? 2 : 1 ) );
      wmesh_int_t num_dofs_element;
      wmesh_status_t status = bms_ndofs(element,degree_,&num_dofs_element);
      WMESH_STATUS_CHECK(status);
      self_->m_dofs_m[l] = num_dofs_element;
      self_->m_dofs_ptr[l+1] = self_->m_dofs_ptr[l] + num_dofs_element * mesh_->m_c2n.m_n[l];
    }

  
  wmesh_status_t status;
  for (wmesh_int_t l=0;l<mesh_->m_c2n.m_size;++l)
    {
      const wmesh_int_t element = l + ( (topodim==3) ? 4 : ( (topodim==2) ? 2 : 1 ) );
      status = wmesh_def_rmacro(&self_->m_patterns[l],
				element,
				nodes_family_,
				degree_);
      WMESH_STATUS_CHECK(status);
    }
  WMESH_STATUS_CHECK(status);
  self_->m_ndofs = self_->m_dofs_ptr[mesh_->m_c2n.m_size];
  
  return WMESH_STATUS_SUCCESS;
}
