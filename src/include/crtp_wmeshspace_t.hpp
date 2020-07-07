#pragma once

template<typename impl_t>
struct crtp_wmeshspace_t
{
  
  inline const wmesh_t* get_mesh() const
  {
    return static_cast<impl_t&>(*this).get_mesh();
  }

  inline const wmesh_t* get_refinement_pattern(wmesh_int_t element_type_) const
  {
    return static_cast<impl_t&>(*this).get_pattern(element_type_);
  }
  
  inline wmesh_int_t get_ndofs() const
  {
    return static_cast<impl_t&>(*this).get_ndofs();
  };

  inline wmesh_int_t get_degree() const
  {
    return static_cast<impl_t&>(*this).get_degree();
  };
  
  inline wmesh_status_t get_dofs_ids(wmesh_int_t 		element_type_,
				     wmesh_int_t 		element_idx_,
				     wmesh_int_p 		dofs_,
				     wmesh_int_t 		dofs_inc_) const
  {
    return static_cast<impl_t&>(*this).get_dofs_ids(element_type_,
						    element_idx_,
						    dofs_,
						    dofs_inc_);    
  };
  
};
