#include "wmesh_t.hpp"
#include "bms.h"

extern "C"
{

  wmesh_status_t wmesh_partitioning(wmesh_t*	self_,
				    wmesh_int_t nparts_)
  {
    WMESH_CHECK_POINTER(self_);
    WMESH_CHECK( (nparts_ > 0) );
    
    size_t  work_size 	= 0;
    void *  work 	= nullptr;
    wmesh_status_t status;
    status = bms_partitioning_sfc(nparts_,
				  self_->m_num_cells,
				  self_->m_c_c.m_data,
				  1,
				  self_->m_c2n.m_size,
				  self_->m_c2n.m_ptr,
				  self_->m_c2n.m_m,
				  self_->m_c2n.m_n,
				  self_->m_c2n.m_data,
				  self_->m_c2n.m_ld,
				  
				  3,
				  self_->m_num_nodes,
				  self_->m_coo,
				  3,
				  
				  &work_size,
				  work);
    WMESH_STATUS_CHECK(status);

    work = (void*)malloc(work_size);
    status = bms_partitioning_sfc(nparts_,
				  self_->m_num_cells,
				  self_->m_c_c.m_data,
				  1,
				  self_->m_c2n.m_size,
				  self_->m_c2n.m_ptr,
				  self_->m_c2n.m_m,
				  self_->m_c2n.m_n,
				  self_->m_c2n.m_data,
				  self_->m_c2n.m_ld,
				  
				  3,
				  self_->m_num_nodes,
				  self_->m_coo,
				  3,
				  
				  &work_size,
				  work);
    free(work);
    work = nullptr;
    if (WMESH_STATUS_SUCCESS != status)
      {
	WMESH_STATUS_CHECK(status);
      }
    return WMESH_STATUS_SUCCESS;
  }
  
};
