#pragma once

#include "wmeshspace.h"
#include "wmesh_t.hpp"
#include "crtp_wmeshspace_t.hpp"

struct wmeshspace_t : public crtp_wmeshspace_t<wmeshspace_t>
{
private:

public:
  wmesh_t * 			m_mesh;
  wmesh_int_t 			m_degree;
  wmesh_t * 			m_patterns[4];
  wmesh_int_t 			m_ndofs;
  wmesh_int_sparsemat_t 	m_c2d;
  wmesh_int_sparsemat_t 	m_c2d_n;
  wmesh_int_sparsemat_t 	m_c2d_e;
  wmesh_int_sparsemat_t 	m_c2d_t;
  wmesh_int_sparsemat_t 	m_c2d_q;
  wmesh_int_sparsemat_t 	m_c2d_i;
  wmesh_int_p			m_dof_codes;
    
  inline const wmesh_t* 	get_mesh() const
  {
    return this->m_mesh;
  };
    
  inline const wmesh_t* 	get_refinement_pattern(wmesh_int_t element_type_) const
  {
    return this->m_patterns[element_type_];
  };
  inline wmesh_t* 	get_refinement_pattern(wmesh_int_t element_type_) 
  {
    return this->m_patterns[element_type_];
  };
    
  inline wmesh_int_t 		get_ndofs() const
  {
    return this->m_ndofs;
  };

  inline wmesh_int_t 		get_degree() const
  {
    return this->m_degree;
  };
    
  inline wmesh_status_t 	get_dofs_ids(wmesh_int_t 		element_type_,
					     wmesh_int_t 		element_idx_,
					     wmesh_int_p 		dofs_,
					     wmesh_int_t 		dofs_inc_) const
  {
    const wmesh_int_t num_dofs = this->m_c2d.m_m[element_type_];
    const wmesh_int_t shift = this->m_c2d.m_ptr[element_type_];
    const wmesh_int_t ld = this->m_c2d.m_ld[element_type_];
    for (wmesh_int_t i=0;i<num_dofs;++i)
      {
	dofs_[dofs_inc_*i] = this->m_c2d.m_data[ shift  + element_idx_ * ld + i ];
      }
    return WMESH_STATUS_SUCCESS;
  };
    
    
  //    double* 		m_coo_dofs;
  //    wmesh_int_t		m_coo_dofs_ld;
};


template<typename T>
wmesh_status_t wmeshspace_generate_coodofs(const wmeshspace_t * __restrict__ 	self_,
					   wmesh_int_t 				coo_storage_,
					   wmesh_int_t 				coo_m_,
					   wmesh_int_t 				coo_n_,
					   T *  __restrict__			coo_,
					   wmesh_int_t 				coo_ld_);




template<typename T>
wmesh_status_t wmeshspace_get_dof_values(const wmeshspace_t& 	self_,
					 wmesh_int_t 		itype_,
					 wmesh_int_t 		idx_elm_,
					 wmesh_int_t 		velocity_storage_,
					 const wmesh_mat_t<T>& 	velocity_,
					 wmesh_int_t 		velocity_dofs_storage_,
					 wmesh_mat_t<T>& 	velocity_dofs_);
