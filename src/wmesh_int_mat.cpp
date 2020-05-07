#include "wmesh-types.hpp"

extern "C"
{

  wmesh_status_t wmesh_int_mat_def(wmesh_int_mat_t * self_,
			 wmesh_int_t m,
			 wmesh_int_t n,
			 wmesh_int_t*v,
			 wmesh_int_t ld)
  {
    self_->m = m;
    self_->n = n;
    self_->v = v;
    self_->ld = ld;
    return WMESH_STATUS_SUCCESS;
  }

}
