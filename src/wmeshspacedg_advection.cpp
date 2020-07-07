#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh.hpp"
#include <chrono>
#include <iostream>
#include "bms.h"
#include "wmesh-utils.hpp"
#include "wmesh-blas.h"
#include "bms_templates.hpp"

#include "wmesh_cubature_factory_t.hpp"
#include "wmesh_cubature_info_t.hpp"
#include "wmesh_shape_eval_factory_t.hpp"
#include "wmesh_integral_convection_t.hpp"
#include "wmesh_integral_flux_t.hpp"
#include "bms.hpp"
#include "wmesh_nodes_boundary_factory_t.hpp"
template<typename T>
static  std::ostream& operator<<(std::ostream&out_,
				 const wmesh_mat_t<T>&that_)
{
  for (wmesh_int_t i=0;i<that_.m;++i)
    {
      for (wmesh_int_t j=0;j<that_.n;++j)
	{
	  out_ << " " << that_.v[that_.ld * j + i];
	}
      out_ << std::endl;
    }
  return out_;
};



#if 0
template <typename T>
wmesh_status_t wmeshspace_get_dofelm(const wmeshspace_t * __restrict__	self_,
				     wmesh_int_t			itype_,
				     wmesh_int_t			ielm_,
				     wmesh_int_t			dofs_storage_,
				     wmesh_int_t			dofs_m_,
				     wmesh_int_t			dofs_n_,
				     const T * 				dofs_,		
				     wmesh_int_t			dofs_ld_,
				     wmesh_int_t			elmsdofs_storage_,
				     wmesh_int_t			elmsdofs_m_,
				     wmesh_int_t			elmsdofs_n_,
				     T * __restrict__			elmsdofs_,
				     wmesh_int_t			elmsdofs_ld_)
{
  WMESH_CHECK_POINTER(self_);
  WMESH_CHECK_POINTER(elmsdofs_);
  switch(elmsdofs_storage_)
    {
    case WMESH_STORAGE_INTERLEAVE:
      {
	switch(dofs_storage_)
	  {
	  case WMESH_STORAGE_INTERLEAVE:
	    {
	      for (wmesh_int_t j=0;j<elmsdofs_n_;++j)
		{
		  for (wmesh_int_t i=0;i<elmsdofs_m_;++i)
		    {
		      elmsdofs_[elmsdofs_ld_*j+i]
			= dofs_[dofs_ld_ * ( self_->m_c2d.m_data[self_->m_c2d.m_ptr[itype_] + self_->m_c2d.m_ld[itype_] * ielm_ + j] - 1) + i];
		    }
		}
	      return WMESH_STATUS_SUCCESS;
	    }
	  case WMESH_STORAGE_BLOCK:
	    {
	      for (wmesh_int_t j=0;j<elmsdofs_n_;++j)
		{
		  for (wmesh_int_t i=0;i<elmsdofs_m_;++i)
		    {
		      elmsdofs_[elmsdofs_ld_*j+i]
			= dofs_[dofs_ld_ * i + ( self_->m_c2d.m_data[self_->m_c2d.m_ptr[itype_] + self_->m_c2d.m_ld[itype_] * ielm_ + j] - 1)];
		    }
		}
	      return WMESH_STATUS_SUCCESS;
	    }
	  }
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
      }
      
    case WMESH_STORAGE_BLOCK:
      {
	switch(dofs_storage_)
	  {
	  case WMESH_STORAGE_INTERLEAVE:
	    {
	      for (wmesh_int_t j=0;j<elmsdofs_n_;++j)
		{
		  for (wmesh_int_t i=0;i<elmsdofs_m_;++i)
		    {
		      elmsdofs_[elmsdofs_ld_*i+j]
			= dofs_[dofs_ld_ * ( self_->m_c2d.m_data[self_->m_c2d.m_ptr[itype_] + self_->m_c2d.m_ld[itype_] * ielm_ + j] - 1) + i];
		    }
		}
	      return WMESH_STATUS_SUCCESS;  
	    }
	  case WMESH_STORAGE_BLOCK:
	    {
	      for (wmesh_int_t j=0;j<elmsdofs_n_;++j)
		{
		  for (wmesh_int_t i=0;i<elmsdofs_m_;++i)
		    {
		      elmsdofs_[elmsdofs_ld_*i+j]
			= dofs_[dofs_ld_ * i + ( self_->m_c2d.m_data[self_->m_c2d.m_ptr[itype_] + self_->m_c2d.m_ld[itype_] * ielm_ + j] - 1)];
		    }
		}
	      return WMESH_STATUS_SUCCESS;  
	    }
	  }
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ENUM);
      }

    }
  return WMESH_STATUS_INVALID_ENUM;  
}
#endif

//template<typename T>
//static  std::ostream& operator<<(std::ostream&out_,
//				 const wmesh_mat_t<T>&that_)
//{
//  for (wmesh_int_t i=0;i<that_.m;++i)
//    {
//      for (wmesh_int_t j=0;j<that_.n;++j)
//	{
//	  out_ << " " << that_.v[that_.ld * j + i];
//	}
//      out_ << std::endl;
//    }
//  return out_;
//};


template<typename IMPL, typename T>
struct wmesh_pde_t
{


  wmesh_status_t jacobian(wmesh_int_t			csr_size_,
			  const_wmesh_int_p		csr_ptr_,
			  const_wmesh_int_p		csr_ind_,
			  T * __restrict__		csr_val_)
  {
    return static_cast<IMPL&>(*this).jacobian(csr_size_,
					      csr_ptr_,
					      csr_ind_,
					      csr_val_);
  };
  
  wmesh_status_t residual(T * __restrict__    rhs_,
			  wmesh_int_t 	      rhs_inc_)
  {
    return static_cast<IMPL&>(*this).residual(rhs_,
					      rhs_inc_);
  };
  
};


template<typename T>
struct wmesh_pde_advection_t : public wmesh_pde_t< wmesh_pde_advection_t<T>, T >
{
  const wmesh_t * 		m_mesh;
  const wmeshspacedg_t * 	m_trial_space;
  const wmeshspacedg_t * 	m_test_space;
  const wmeshspace_t *   	m_velocity_space;
  wmesh_int_t 			m_velocity_storage;
  const wmesh_mat_t<T> *   	m_velocity;
  
  wmesh_cubature_info_t 	m_cubature_info;
  wmesh_shape_info_t 		m_shape_info_element;
  wmesh_shape_info_t 		m_shape_info_trial;
  wmesh_shape_info_t 		m_shape_info_test;
  wmesh_shape_info_t 		m_shape_info_velocity;
  
  wmesh_template_integral_convection_t<T> * m_integral_convection[WMESH_ELEMENT_ALL]{};
  wmesh_template_integral_flux_t<T> * m_integral_flux[WMESH_ELEMENT_ALL]{};


  const wmesh_nodes_boundary_t<T>* nodes_boundary_trial[4]{};
  const wmesh_nodes_boundary_t<T>* nodes_boundary_test[4]{};
  const wmesh_nodes_boundary_t<T>* nodes_boundary_velocity[4]{};
  const wmesh_nodes_boundary_t<T>* nodes_boundary_element[4]{};
    
  const wmesh_shape_eval_t<T> * shape_eval_boundary_element[6][8]{};
  const wmesh_shape_eval_t<T> * shape_eval_boundary_trial[6][8]{};
  const wmesh_shape_eval_t<T> * shape_eval_boundary_test[6][8]{};
  const wmesh_shape_eval_t<T> * shape_eval_boundary_velocity[6][8]{};

  
  wmesh_pde_advection_t(const wmeshspacedg_t * 		trial_space_,
			const wmeshspacedg_t * 		test_space_,
			const wmeshspace_t * 		velocity_space_,
			wmesh_int_t 			velocity_storage_,
			const wmesh_mat_t<T> *   	velocity_,
			const wmesh_cubature_info_t& 	cubature_info_, 		
			const wmesh_shape_info_t& 	shape_info_element_, 		
			const wmesh_shape_info_t& 	shape_info_trial_,
			const wmesh_shape_info_t& 	shape_info_test_,
			const wmesh_shape_info_t& 	shape_info_velocity_)
    : m_mesh(trial_space_->m_mesh),
      m_trial_space(trial_space_),
      m_test_space(test_space_),
      m_velocity_space(velocity_space_),
      m_velocity_storage(velocity_storage_),
      m_velocity(velocity_),
      m_cubature_info(cubature_info_),
      m_shape_info_element(shape_info_element_),
      m_shape_info_trial(shape_info_trial_),
      m_shape_info_test(shape_info_test_),
      m_shape_info_velocity(shape_info_velocity_)
  {
    const wmesh_int_t 	topodim = trial_space_->m_mesh->m_topology_dimension;      
    const wmesh_int_t 	ntypes 	= trial_space_->m_mesh->m_c2n.m_size;      
    wmesh_int_t facets_flags[WMESH_ELEMENT_ALL]{};
    for (wmesh_int_t itype=0;itype<ntypes;++itype)
      {
	wmesh_int_t facets[6];
	wmesh_int_t element = (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );
	const wmesh_int_t num_elements = trial_space_->m_mesh->m_c2n.m_n[itype];
	if (num_elements > 0)
	  {
	    this->m_integral_convection[element] = new wmesh_template_integral_convection_t<T>(element, 
											       cubature_info_, 
											       this->m_shape_info_element, 
											       this->m_shape_info_velocity, 
											       this->m_shape_info_test, 
											       this->m_shape_info_trial);
	    wmesh_int_t num_facets;
	    wmesh_status_t status = bms_element_facets(element,
						       &num_facets,
						       facets);
	    if (status!=0)
	      {
		
	      }

	    for (wmesh_int_t i=0;i<num_facets;++i)
	      {
		if (!facets_flags[i])
		  {
		    facets_flags[i] = true;
		    wmesh_cubature_info_t 	cubature_info_facet;
		    status = wmesh_cubature_info_def(&cubature_info_facet,
						     this->m_cubature_info.m_family,
						     this->m_cubature_info.m_degree);
		    if (!status)
		      {
			
		      }
		    
		    wmesh_shape_info_t 	shape_info_facet;
		    status = wmesh_shape_info_def(&shape_info_facet,
						  this->m_shape_info_element.m_family,
						  this->m_shape_info_element.m_degree);
		    if (!status)
		      {
			
		      }
		    
		    this->m_integral_flux[ facets[i] ] = new wmesh_template_integral_flux_t<T>(facets[i], 
											       cubature_info_facet, 
											       shape_info_facet, 
											       this->m_shape_info_velocity, 
											       this->m_shape_info_test, 
											       this->m_shape_info_trial);
		  }
	      }
	  }
      }


    const wmesh_int_sparsemat_t&trial_c2n = this->m_mesh->m_c2n;
    for (wmesh_int_t itype=0;itype<trial_c2n.m_size;++itype)
      {
	if (trial_c2n.m_n[itype] > 0)
	  {	    
	    const wmesh_int_t element = (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );	   
	    nodes_boundary_trial[itype] = wmesh_nodes_boundary_factory_t<T>::nodes_boundary_instance	(element,
													 this->m_integral_convection[element]->m_shape_trial.m_nodes_family,
													 this->m_integral_convection[element]->m_shape_trial.m_degree);
	    
	    nodes_boundary_test[itype] = wmesh_nodes_boundary_factory_t<T>::nodes_boundary_instance	(element,
													 this->m_integral_convection[element]->m_shape_test.m_nodes_family,
													 this->m_integral_convection[element]->m_shape_test.m_degree);
	    
	    nodes_boundary_velocity[itype] = wmesh_nodes_boundary_factory_t<T>::nodes_boundary_instance	(element,
													 this->m_integral_convection[element]->m_shape_velocity.m_nodes_family,
													 this->m_integral_convection[element]->m_shape_velocity.m_degree);
	    
	    nodes_boundary_element[itype] = wmesh_nodes_boundary_factory_t<T>::nodes_boundary_instance	(element,
													 this->m_integral_convection[element]->m_shape_element.m_nodes_family,
													 this->m_integral_convection[element]->m_shape_element.m_degree);
	    
	  }
      }

    for (wmesh_int_t itype=0;itype<trial_c2n.m_size;++itype)
      {
	if (nodes_boundary_trial[itype])
	  {
	const wmesh_int_t element = (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );	   
	    wmesh_int_t num_facets;
	    wmesh_int_t facets[6];
	    wmesh_status_t status = bms_element_facets(element,
						       &num_facets,
						       facets);
	    if (status!=0)
	      {
		
	      }

	    for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
	      {
		const wmesh_int_t facet		= facets[ifacet];
		const wmesh_int_t facet_type 	= (topodim==3) ? (facet - 2) : ((topodim==2) ? (facet-1) : ( (topodim==1) ? (facet-1) : -1));
		const wmesh_int_t n = (topodim==3) ? (facet_type + 3) : ( (topodim==2) ? 2 : 1);
		
		for (wmesh_int_t irot=0;irot < 2*n;++irot)
		  {
		    shape_eval_boundary_element[ifacet][irot] = wmesh_shape_eval_factory_t<T>::shape_eval_instance(element,
														   this->m_integral_convection[element]->m_shape_element.m_family,
														   this->m_integral_convection[element]->m_shape_element.m_degree,
														   nodes_boundary_element[itype]->m_facets_nodes_storage,
														   &nodes_boundary_element[itype]->m_facets_nodes[ifacet][irot]);
		  }

		for (wmesh_int_t irot=0;irot < 2*n;++irot)
		  {
		    shape_eval_boundary_trial[ifacet][irot] = wmesh_shape_eval_factory_t<T>::shape_eval_instance(element,
														 this->m_integral_convection[element]->m_shape_trial.m_family,
														 this->m_integral_convection[element]->m_shape_trial.m_degree,
														 nodes_boundary_trial[itype]->m_facets_nodes_storage,
														 &nodes_boundary_trial[itype]->m_facets_nodes[ifacet][irot]);
		  }
		
		for (wmesh_int_t irot=0;irot < 2*n;++irot)
		  {
		    shape_eval_boundary_velocity[ifacet][irot] = wmesh_shape_eval_factory_t<T>::shape_eval_instance(element,
														    this->m_integral_convection[element]->m_shape_velocity.m_family,
														    this->m_integral_convection[element]->m_shape_velocity.m_degree,
														    nodes_boundary_velocity[itype]->m_facets_nodes_storage,
														    &nodes_boundary_velocity[itype]->m_facets_nodes[ifacet][irot]);
		  }

		for (wmesh_int_t irot=0;irot < 2*n;++irot)
		  {
		    shape_eval_boundary_test[ifacet][irot] = wmesh_shape_eval_factory_t<T>::shape_eval_instance(element,
														 this->m_integral_convection[element]->m_shape_test.m_family,
														 this->m_integral_convection[element]->m_shape_test.m_degree,
														 nodes_boundary_test[itype]->m_facets_nodes_storage,
														 &nodes_boundary_test[itype]->m_facets_nodes[ifacet][irot]);
		  }
	      }
	  }
      }
    


    
  };
  
  wmesh_status_t jacobian(wmesh_int_t			csr_size_,
			  const_wmesh_int_p		csr_ptr_,
			  const_wmesh_int_p		csr_ind_,
			  T * __restrict__		csr_val_)
  {

  for (wmesh_int_t i=0;i<this->m_trial_space->m_ndofs;++i)
    {     
      for (wmesh_int_t k = csr_ptr_[i];k<csr_ptr_[i+1];++k)
	{
csr_val_[k] = static_cast<T>(0);
	}
    }

    const wmesh_int_t 	topodim = this->m_trial_space->m_mesh->m_topology_dimension;      
    //
    // Loop over the trial space.
    //
    using convection_data_t = typename wmesh_template_integral_convection_t<T>::data_t;
    convection_data_t * integral_convection_data[WMESH_ELEMENT_ALL]{};
    
    using flux_data_t = typename wmesh_template_integral_flux_t<T>::data_t;
    flux_data_t * integral_flux_data[WMESH_ELEMENT_ALL]{};
    
    const wmesh_int_t 	ntypes 	= this->m_trial_space->m_mesh->m_c2n.m_size;      
    for (wmesh_int_t itype=0;itype<ntypes;++itype)
      {
	wmesh_int_t element = (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );
	const wmesh_int_t num_elements = this->m_trial_space->m_mesh->m_c2n.m_n[itype];
	if (num_elements > 0)
	  {
	    integral_convection_data[element] = new convection_data_t(*this->m_integral_convection[element]);
	  }
      } 






    
    const wmesh_int_sparsemat_t& trial_c2n = this->m_trial_space->m_mesh->m_c2n;
    wmesh_int_t s_m[2];
    wmesh_int_t s_n[2];
    wmesh_int_t s_ld[2];
    wmesh_int_t s_data[4*3 + 6*4];
    wmesh_int_p s_v[2] = {&s_data[0], &s_data[4*3]};
    wmesh_status_t status;
    for (wmesh_int_t itype=0;itype<trial_c2n.m_size;++itype)
      {
	if (trial_c2n.m_n[itype] > 0)
	  {	    

	    if (topodim==2)
	      {
		integral_flux_data[WMESH_ELEMENT_EDGE] = new flux_data_t(*this->m_integral_flux[WMESH_ELEMENT_EDGE]);	
		status = bms_s_e2n_type(itype,
					topodim,
					&s_m[0],
					&s_n[0],
					s_v[0],
					&s_ld[0]);
	      }
	    else if (topodim==3)
	      {
		if ((this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_TETRAHEDRON] > 0)||(this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_PYRAMID] > 0)||(this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_WEDGE] > 0))
		  {
		    integral_flux_data[WMESH_ELEMENT_TRIANGLE] = new flux_data_t(*this->m_integral_flux[WMESH_ELEMENT_TRIANGLE]);
		    status = bms_s_t2n_type(itype,
					    &s_m[0],
					    &s_n[0],
					    s_v[0],
					    &s_ld[0]);
		    
		  }
		if ((this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_HEXAHEDRON] > 0)||(this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_PYRAMID] > 0)||(this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_WEDGE] > 0))
		  {
		    integral_flux_data[WMESH_ELEMENT_QUADRILATERAL] = new flux_data_t(*this->m_integral_flux[WMESH_ELEMENT_QUADRILATERAL]);
		    status = bms_s_q2n_type(itype,
					    &s_m[1],
					    &s_n[1],
					    s_v[1],
					    &s_ld[1]);
		  }
		
	      }
	  }
      }


    wmesh_int_t rw_n = 0;
    wmesh_mat_t<T> local_matrix_tmp;
    for (wmesh_int_t itype=0;itype<trial_c2n.m_size;++itype)
      {
	if (trial_c2n.m_n[itype] > 0)
	  {
	    const wmesh_int_t element 	= (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );
	    wmesh_int_t facets[6];
	    wmesh_int_t num_facets;
	     status = bms_element_facets(element,
					 &num_facets,
					 facets);
	     WMESH_STATUS_CHECK(status);
	     for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
	      {
		const wmesh_mat_t<T>& eval_btest= shape_eval_boundary_test[ifacet][0]->m_f;
		rw_n = (rw_n < eval_btest.m* eval_btest.n) ? (eval_btest.m* eval_btest.n) : rw_n;
	      }
	  }
      }
     T * __restrict__ rw = (T*__restrict__)malloc(sizeof(T)*(rw_n));
    
    wmesh_int_t rw2_n = 2048;
    wmesh_mat_t<T> local_matrix_tmp2;
     T * __restrict__ rw2 = (T*__restrict__)malloc(sizeof(T)*(rw2_n));
    const wmesh_int_t 	idofs_inc = 1;
    const wmesh_int_t 	jdofs_inc = 1;

    wmesh_int_t idofs_ids[1024];
    wmesh_int_t jdofs_ids[1024];
    wmesh_int_t nei_jdofs_ids[1024];

    wmesh_int_t udofs_ids[1024];


//    const wmesh_int_sparsemat_t& test_c2n = this->m_trial_space->m_mesh->m_c2n;
    for (wmesh_int_t itype=0;itype<trial_c2n.m_size;++itype)
      {
	if (trial_c2n.m_n[itype] > 0)
	  {
	    const wmesh_int_t element 	= (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );
#if 0
	    const wmesh_int_t trial_c2n_m  	= trial_c2n.m_m[itype];
	    const wmesh_int_t trial_c2n_n  	= trial_c2n.m_n[itype];
	    const wmesh_int_t trial_c2n_ld 	= trial_c2n.m_ld[itype];
	    const_wmesh_int_p trial_c2n_v 	= trial_c2n.m_data + trial_c2n.m_ptr[itype];
	    
	    const wmesh_int_t test_c2n_m  	= test_c2n.m_m[itype];
	    const wmesh_int_t test_c2n_n  	= test_c2n.m_n[itype];
	    const wmesh_int_t test_c2n_ld 	= test_c2n.m_ld[itype];
	    const_wmesh_int_p test_c2n_v 	= test_c2n.m_data + test_c2n.m_ptr[itype];
#ifndef NDEBUG
	    WMESH_CHECK(test_c2n_n == trial_c2n_n);
#endif
#endif

	    wmesh_int_t facets[6];
	    wmesh_int_t num_facets;
	    wmesh_status_t status = bms_element_facets(element,
						       &num_facets,
						       facets);
	    WMESH_STATUS_CHECK(status);
	    
	    auto * data = integral_convection_data[element];

	    const wmesh_int_t 	idofs_n 	= this->m_test_space->m_dofs_m[itype];
	    const wmesh_int_t 	jdofs_n 	= this->m_trial_space->m_dofs_m[itype];



	    
	    for (wmesh_int_t idx_elm=0;idx_elm<trial_c2n.m_n[itype];++idx_elm)
	      {
		status = wmeshspacedg_get_dofs_ids	(this->m_test_space,
							 itype,
							 idx_elm,
							 idofs_ids,
							 1);
		WMESH_STATUS_CHECK(status);
		
		status = wmeshspacedg_get_dofs_ids	(this->m_trial_space,
							 itype,
							 idx_elm,
							 jdofs_ids,
							 1);
		WMESH_STATUS_CHECK(status);
		
		status = wmeshspace_get_dofs_ids	(this->m_velocity_space,
							 itype,
							 idx_elm,
							 udofs_ids,
							 1);
		WMESH_STATUS_CHECK(status);

		
		//
		// Get idofs and jdofs from test and trial spaces.
		//
		// const_wmesh_int_p 	idofs 	= test_c2d_v + test_c2d_ld * idx_elm;
		// const_wmesh_int_p 	jdofs 	= trial_c2d_v + trial_c2d_ld * idx_elm;
		
		//
		// Extract coordinates.
		//
		status = wmesh_get_cooelm(this->m_mesh,
					  itype,
					  idx_elm,
					  data->m_dofs_element_storage,
					  WMESH_MAT_FORWARD(data->m_dofs_element));
		WMESH_STATUS_CHECK(status);
		

		//
		// Extract velocity dofs.
		//
		status = wmeshspace_get_dof_values(*this->m_velocity_space,
						   itype,
						   idx_elm,
						   this->m_velocity_storage,
						   *this->m_velocity,
						   data->m_dofs_velocity_storage,
						   data->m_dofs_velocity);
		WMESH_STATUS_CHECK(status);

		//		std::cout << data->m_dofs_velocity << std::endl;
		
		//
		// Compute local matrix.
		//
#if 1
		status = this->m_integral_convection[element]->eval(*data,
								    static_cast<T>(0),
								    data->m_local_matrix);
		WMESH_STATUS_CHECK(status);
#else
		wmesh_mat_t<T>::zero(data->m_local_matrix);
#endif
		//
		// 
		//
#if 1

		for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
		  {
		    const wmesh_int_t facet 	= facets[ifacet];
		    // const wmesh_int_t facet_type= (topodim==3) ? (facet - 2) : ((topodim==2) ? (facet-1) : ( (topodim==1) ? (facet-1) : -1));
		    auto * flux_data = integral_flux_data[facet];
		    //
		    // Now restrict the coordinates of the element.
		    //
		    //		    std::cout << "-------------------" << std::endl;
		    //		    std::cout << data->m_dofs_element << std::endl;
		    //		    std::cout << "-------------------" << std::endl;
		    //		    std::cout << shape_eval_boundary_element[ifacet][0]->m_f << std::endl;
		    //		    std::cout << "-------------------" << std::endl;
		    wmesh_mat_gemm(static_cast<T>(1),
				   data->m_dofs_element,
				   shape_eval_boundary_element[ifacet][0]->m_f,
				   static_cast<T>(0),
				   flux_data->m_dofs_element);
		    //		    std::cout << "####### " << flux_data->m_dofs_element.m << std::endl;
		    //		    std::cout << "####### " << flux_data->m_dofs_element.n << std::endl;
		    
		    //		    std::cout << "####### FACE" << std::endl;
		    //		    std::cout << flux_data->m_dofs_element << std::endl;

		    
		    //		    std::cout << "???" << std::endl;
		    
		    //
		    // Now restrict the velocity.
		    //
		    wmesh_mat_gemm(static_cast<T>(1),
				   data->m_dofs_velocity,
				   shape_eval_boundary_velocity[ifacet][0]->m_f,
				   static_cast<T>(0),				   
				   flux_data->m_dofs_velocity);
		    //		    std::cout << "???" << std::endl;
		    //		    std::cout << "####### VELOCITY" << std::endl;
		    //		    std::cout << flux_data->m_dofs_velocity << std::endl;
		    
		    
		    const wmesh_mat_t<T>& eval_btrial 		= shape_eval_boundary_trial[ifacet][0]->m_f;		      
		    const wmesh_mat_t<T>& eval_btest 		= shape_eval_boundary_test[ifacet][0]->m_f;		      
		    
		    
		    //
		    // Compute the flux on the facet.
		    //
		    status = this->m_integral_flux[facet]->eval(*flux_data,
								static_cast<T>(0),
								flux_data->m_local_matrix);
		    //		    std::cout << "flux_data->m_local_matrix " << std::endl;
		    //		    std::cout << flux_data->m_local_matrix << std::endl;

		    //
		    // Up the dimensions. As it is (INTERLEAVE MODE FOR SHAPE FUNCTIONS), it gives eval_btest * flux_data->m_local_matrix * transpose(eval_btrial)
		    //
		    //		    std::cout << "???1" << std::endl;
		    wmesh_mat_t<T>::define(&local_matrix_tmp,
					   eval_btest.m,
					   eval_btest.n,
					   rw,
					   eval_btest.m);
		    
		    wmesh_mat_gemm(static_cast<T>(1),
				   eval_btest,
				   flux_data->m_local_matrix,
				   static_cast<T>(0),				   
				   local_matrix_tmp);
		    
		    //		    std::cout << "???2" << std::endl;
		    ///		    std::cout << 				   flux_data->m_local_matrix_tmp << std::endl;
		    //  exit(1);
		    
		    //		    std::cout << "???2bbb" << std::endl;
		    //		    std::cout << eval_btrial << std::endl;

		    wmesh_mat_gemm("N",
				   "T",
				   static_cast<T>(1),
				   local_matrix_tmp,
				   eval_btrial,
				   static_cast<T>(1),
				   data->m_local_matrix);
		    //		    std::cout << "local flux matrix " << std::endl;
		    //		    std::cout << data->m_local_matrix << std::endl;
		    //    std::cout << "???3" << std::endl;

#if 0
		    std::cout << "-------------- ifacet " << ifacet << std::endl;
		    std::cout << "################################################## " << std::endl;
		    std::cout << btrial << std::endl;
		    std::cout << "################################################## " << std::endl;
		    std::cout << btest << std::endl;
		    std::cout << "################################################## " << std::endl;
		    std::cout << bvelocity << std::endl;
#endif
		    
		  }

#endif

		//
		// Add the local matrix.
		//
		status =  bms_sparse_add(idofs_n,
					 idofs_ids,
					 idofs_inc,
					 
					 jdofs_n,
					 jdofs_ids,
					 jdofs_inc,

					 data->m_local_matrix.v,
					 data->m_local_matrix.ld,
					 
					 csr_size_,
					 csr_ptr_,
					 csr_ind_,
					 csr_val_);		
		WMESH_STATUS_CHECK(status);
#if 0
		wmesh_int_t bnumtypes=topodim-1;
		for (wmesh_int_t ibtype=0;ibtype<bnumtypes;++ibtype)
		  {
		    wmesh_int_t nfacet_of_type              = s_n[ibtype];
		    wmesh_int_t num_nodes_in_nfacet_of_type = s_m[ibtype];
		  }
#endif

		//
		// Now extra-diagonal block.
		//
		wmesh_int_t cncfacet_n;
		wmesh_int_t cncneifacet_n;
		wmesh_int_t cncfacet[4]{},cncneifacet[4]{};
		wmesh_int_t sss[4] = {4,4,2,0};
		for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
		  {
		    wmesh_int_t nei_info = this->m_mesh->m_c2c.m_data[this->m_mesh->m_c2c.m_ptr[itype] + this->m_mesh->m_c2c.m_ld[itype]*idx_elm + ifacet];
		    if (!nei_info)  continue;

		    
		    const wmesh_int_t facet 	= facets[ifacet];
		    const wmesh_int_t facet_type= (topodim==3) ? (facet - 2) : ((topodim==2) ? (facet-1) : ( (topodim==1) ? (facet-1) : -1));
		    
		    status =  wmesh_get_facet_dofs_ids(this->m_mesh,
						       itype,
						       idx_elm,
						       facet_type,
						       (facet == WMESH_ELEMENT_QUADRILATERAL) ? (ifacet - sss[itype])  : ifacet,
						       &cncfacet_n,
						       cncfacet,
						       1);
		    
		    //
		    //
		    //
		    wmesh_int_t nei_facet_lidx = 0;
		    
		    wmesh_int_t nei_idx;
		    wmesh_int_t nei_type;
		    status = bms_c2c_cindex(nei_info,
					    &nei_idx);
		    WMESH_STATUS_CHECK(status);		    
		    status = bms_c2c_ctype(nei_info,
					   &nei_type);
		    WMESH_STATUS_CHECK(status);
		    for (;nei_facet_lidx< this->m_mesh->m_c2c.m_m[nei_type];++nei_facet_lidx)
		      {
			wmesh_int_t nei2_info = this->m_mesh->m_c2c.m_data[this->m_mesh->m_c2c.m_ptr[nei_type] + this->m_mesh->m_c2c.m_ld[nei_type]*nei_idx + nei_facet_lidx];			
			wmesh_int_t nei2_idx;
			wmesh_int_t nei2_type;
			status = bms_c2c_cindex(nei2_info,
						&nei2_idx);
			WMESH_STATUS_CHECK(status);		    
			status = bms_c2c_ctype(nei2_info,
					       &nei2_type);
			WMESH_STATUS_CHECK(status);
			if ((itype == nei2_type) && (idx_elm == nei2_idx))
			  {
			    break;
			  }			
		      }
		    WMESH_CHECK(nei_facet_lidx< this->m_mesh->m_c2c.m_m[nei_type]);

		    
		    status =  wmesh_get_facet_dofs_ids(this->m_mesh,
						       nei_type,
						       nei_idx,
						       facet_type,
						       (facet == WMESH_ELEMENT_QUADRILATERAL) ? (nei_facet_lidx - sss[nei_type])  : nei_facet_lidx,
						       &cncneifacet_n,
						       cncneifacet,
						       1);
		    

		    //
		    // Now calculate the rotation.
		    //
#if 0
		    std::cout << "bertrand ################### " << std::endl;
		    std::cout << " " << cncfacet[0] << std::endl;
		    std::cout << " " << cncfacet[1] << std::endl;

		    std::cout << "delanoe ################### " << std::endl;
		    std::cout << " " << cncneifacet[0] << std::endl;
		    std::cout << " " << cncneifacet[1] << std::endl;
#endif
		    wmesh_int_t signed_rotation = 0;
		    {
		      for (;signed_rotation<cncneifacet_n;++signed_rotation)
			{
			  if (cncneifacet[signed_rotation] == cncfacet[0])
			    {
			      break;
			    }
			}
		    }
		    //		    std::cout << signed_rotation << " " << s_m[facet_type] << std::endl;
		    WMESH_CHECK(signed_rotation<cncneifacet_n);
		    //		    signed_rotation += cncneifacet_n;
#if 0
		    std::cout << "nei idx  " << nei_idx  << std::endl;
		    std::cout << "nei type " << nei_type << std::endl;
		    //		    --signed_rotation;
		    std::cout << "yo signed rotation " << signed_rotation << std::endl;
#endif

		    wmesh_int_t nei_element 	= (topodim==3) ? (4 + nei_type) : ( (topodim==2) ? (2+nei_type) : (1+nei_type) );
		    auto * a = wmesh_shape_eval_factory_t<T>::shape_eval_instance(nei_element,
										  this->m_integral_convection[nei_element]->m_shape_trial.m_family,
										  this->m_integral_convection[nei_element]->m_shape_trial.m_degree,
										  nodes_boundary_trial[nei_type]->m_facets_nodes_storage,
										  &nodes_boundary_trial[nei_type]->m_facets_nodes[nei_facet_lidx][signed_rotation]);
		    
		    const wmesh_mat_t<T>& eval_nei_btrial = a->m_f;		    
		    auto * flux_data = integral_flux_data[facet];
		    
		    wmesh_mat_gemm(static_cast<T>(1),
				   data->m_dofs_element,
				   shape_eval_boundary_element[ifacet][0]->m_f,
				   static_cast<T>(0),
				   flux_data->m_dofs_element);
		    
		    wmesh_mat_gemm(static_cast<T>(1),
				   data->m_dofs_velocity,
				   shape_eval_boundary_velocity[ifacet][0]->m_f,
				   static_cast<T>(0),				   
				   flux_data->m_dofs_velocity);
		    

		    const wmesh_mat_t<T>& eval_btest 		= shape_eval_boundary_test[ifacet][0]->m_f;
				    
		    //
		    // Compute the flux on the facet.
		    //
		    status = this->m_integral_flux[facet]->eval(*flux_data,
								static_cast<T>(0),
								flux_data->m_local_matrix);

		    

		    //
		    // Up the dimensions. As it is (INTERLEAVE MODE FOR SHAPE FUNCTIONS), it gives eval_btest * flux_data->m_local_matrix * transpose(eval_btrial)
		    //
		    //		    std::cout << "???1" << std::endl;
		    wmesh_mat_t<T>::define(&local_matrix_tmp,
					   eval_btest.m,
					   eval_btest.n,
					   rw,
					   eval_btest.m);
		    
		    wmesh_mat_gemm(static_cast<T>(1),
				   eval_btest,
				   flux_data->m_local_matrix,
				   static_cast<T>(0),				   
				   local_matrix_tmp);
		    
		    //		    std::cout << "???2" << std::endl;
		    ///		    std::cout << 				   flux_data->m_local_matrix_tmp << std::endl;
		    //  exit(1);
		    
		    //		    std::cout << "???2bbb" << std::endl;
		    //		    std::cout << eval_btrial << std::endl;

		    wmesh_mat_t<T>::define(&local_matrix_tmp2,
					   this->m_test_space->m_dofs_m[itype],
					   this->m_trial_space->m_dofs_m[nei_type],
					   rw2,
					   this->m_test_space->m_dofs_m[itype]);

		    wmesh_mat_gemm("N",
				   "T",
				   static_cast<T>(-1),
				   local_matrix_tmp,
				   eval_nei_btrial,
				   static_cast<T>(0),
				   local_matrix_tmp2);
#if 0
		    std::cout << "eval_btest " << std::endl;
		    std::cout << eval_btest << std::endl;
		    std::cout << "GGGGGGGGGGGGGGGGGg " << std::endl;
		    std::cout << eval_nei_btrial << std::endl;
		    std::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT " << std::endl;
		    std::cout << local_matrix_tmp2 << std::endl;
#endif

		    {
		      const wmesh_int_t 	nei_jdofs_n 	= this->m_trial_space->m_dofs_m[nei_type];
		      status = wmeshspacedg_get_dofs_ids	(this->m_trial_space,
								 nei_type,
								 nei_idx,
								 nei_jdofs_ids,
								 1);
		      WMESH_STATUS_CHECK(status);
		      
		      status =  bms_sparse_add(idofs_n,
					       idofs_ids,
					       idofs_inc,
					       
					       nei_jdofs_n,
					       nei_jdofs_ids,
					       1,
					       
					       local_matrix_tmp2.v,
					       local_matrix_tmp2.ld,
					       
					       csr_size_,
					       csr_ptr_,
					       csr_ind_,
					       csr_val_);		
		      WMESH_STATUS_CHECK(status);
		    }

		    //		    std::cout << "aaa " << std::endl;
		  }

	      }

	  }
	
      }
    
    return WMESH_STATUS_SUCCESS;
  };


  
  wmesh_status_t residual(T * __restrict__		rhs_,
			  wmesh_int_t 			rhs_inc_)
  {


    
    const wmesh_int_t topodim = this->m_mesh->m_topology_dimension;
    using flux_data_t = typename wmesh_template_integral_flux_t<T>::data_t;
    flux_data_t * integral_flux_data[WMESH_ELEMENT_ALL]{};

    const wmesh_int_sparsemat_t& trial_c2n = this->m_trial_space->m_mesh->m_c2n;
    wmesh_int_t s_m[2];
    wmesh_int_t s_n[2];
    wmesh_int_t s_ld[2];
    wmesh_int_t s_data[4*3 + 6*4];
    wmesh_int_p s_v[2] = {&s_data[0], &s_data[4*3]};
    wmesh_status_t status;
    for (wmesh_int_t itype=0;itype<trial_c2n.m_size;++itype)
      {
	if (trial_c2n.m_n[itype] > 0)
	  {	    

	    if (topodim==2)
	      {
		integral_flux_data[WMESH_ELEMENT_EDGE] = new flux_data_t(*this->m_integral_flux[WMESH_ELEMENT_EDGE]);	
		status = bms_s_e2n_type(itype,
					topodim,
					&s_m[0],
					&s_n[0],
					s_v[0],
					&s_ld[0]);
		WMESH_STATUS_CHECK(status);
	      }
	    else if (topodim==3)
	      {
		if ((this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_TETRAHEDRON] > 0)||(this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_PYRAMID] > 0)||(this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_WEDGE] > 0))
		  {
		    integral_flux_data[WMESH_ELEMENT_TRIANGLE] = new flux_data_t(*this->m_integral_flux[WMESH_ELEMENT_TRIANGLE]);
		    status = bms_s_t2n_type(itype,
					    &s_m[0],
					    &s_n[0],
					    s_v[0],
					    &s_ld[0]);
		    WMESH_STATUS_CHECK(status);	    
		  }
		if ((this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_HEXAHEDRON] > 0)||(this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_PYRAMID] > 0)||(this->m_mesh->m_c2n.m_n[WMESH_ELEMENT_WEDGE] > 0))
		  {
		    integral_flux_data[WMESH_ELEMENT_QUADRILATERAL] = new flux_data_t(*this->m_integral_flux[WMESH_ELEMENT_QUADRILATERAL]);
		    status = bms_s_q2n_type(itype,
					    &s_m[1],
					    &s_n[1],
					    s_v[1],
					    &s_ld[1]);
		    WMESH_STATUS_CHECK(status);
		  }
		
	      }
	  }
      }

    wmesh_int_t idofs_ids[1024];

    wmesh_int_t udofs_ids[1024];

    wmesh_int_t ntypes = this->m_mesh->m_c2n.m_size;
    wmesh_int_t mx_ndofs_velocity = 0;    
    for (wmesh_int_t itype=0;itype<ntypes;++itype)
      {
	if (this->m_mesh->m_c2n.m_n[itype] > 0)
	  {
	    mx_ndofs_velocity = (mx_ndofs_velocity < this->m_velocity_space->m_c2d.m_m[itype]) ? this->m_velocity_space->m_c2d.m_m[itype] : mx_ndofs_velocity;
	  }
      }
    const wmesh_int_t dofs_velocity_storage  = WMESH_STORAGE_INTERLEAVE;
    const wmesh_int_t dofs_velocity_m = (this->m_velocity_storage == WMESH_STORAGE_INTERLEAVE) ? this->m_velocity->m : this->m_velocity->n;
    const wmesh_int_t dofs_velocity_ld = dofs_velocity_m;
    T*__restrict__ dofs_velocity_data = (T*__restrict__)malloc(sizeof(T)*dofs_velocity_m*mx_ndofs_velocity);
    wmesh_mat_t<T> dofs_velocity;


    wmesh_int_t mx_ndofs_test = 0;    
    for (wmesh_int_t itype=0;itype<ntypes;++itype)
      {
	if (this->m_mesh->m_c2n.m_n[itype] > 0)
	  {
	    mx_ndofs_test = (mx_ndofs_test < this->m_test_space->m_dofs_m[itype]) ? this->m_test_space->m_dofs_m[itype] : mx_ndofs_test;
	  }
      }

    T*__restrict__ local_rhs_data = (T*__restrict__)malloc(sizeof(T)*mx_ndofs_test);
    wmesh_mat_t<T> local_rhs;
    

    wmesh_int_t mx_ndofs_element = 0;    
    for (wmesh_int_t itype=0;itype<this->m_mesh->m_c2n.m_size;++itype)
      {
	if (this->m_mesh->m_c2n.m_n[itype] > 0)
	  {
	    mx_ndofs_element = (mx_ndofs_element < this->m_mesh->m_c2n.m_m[itype]) ? this->m_mesh->m_c2n.m_m[itype] : mx_ndofs_element;
	  }
      }
    const wmesh_int_t dofs_element_storage  = WMESH_STORAGE_INTERLEAVE;
    //    const wmesh_int_t dofs_element_m = (this->m_mesh->m_coo_storage == WMESH_STORAGE_INTERLEAVE) ? this->m_mesh->m_coo_m : this->m_mesh->m_coo_n;
    const wmesh_int_t dofs_element_m = this->m_mesh->m_coo_m;
    const wmesh_int_t dofs_element_ld = dofs_element_m;
    T*__restrict__ dofs_element_data = (T*__restrict__)malloc(sizeof(T)*dofs_element_m*mx_ndofs_element);
    wmesh_mat_t<T> dofs_element;

    for (wmesh_int_t i=0;i<this->m_test_space->m_ndofs;++i)
      {     
	rhs_[i*rhs_inc_] = 0.0;
      }
    
    for (wmesh_int_t itype=0;itype<this->m_mesh->m_c2n.m_size;++itype)
      {
	if (this->m_mesh->m_c2n.m_n[itype] > 0)
	  {
	    const wmesh_int_t element 	= (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );

	    wmesh_mat_t<T>::define(&local_rhs,
				   this->m_test_space->m_dofs_m[itype],
				   1,
				   local_rhs_data,
				   this->m_test_space->m_dofs_m[itype]);
	    
	    wmesh_mat_t<T>::define(&dofs_element,dofs_element_m,this->m_mesh->m_c2n.m_m[itype],dofs_element_data,dofs_element_ld);
	    wmesh_mat_t<T>::define(&dofs_velocity,dofs_velocity_m,this->m_velocity_space->m_c2d.m_m[itype],dofs_velocity_data,dofs_velocity_ld);

	    wmesh_int_t facets[6];
	    wmesh_int_t num_facets;
	    status = bms_element_facets(element,
						       &num_facets,
						       facets);
	    WMESH_STATUS_CHECK(status);

	    const wmesh_int_t 	idofs_n 	= this->m_test_space->m_dofs_m[itype];
	    for (wmesh_int_t idx_elm=0;idx_elm<this->m_mesh->m_c2n.m_n[itype];++idx_elm)
	      {
		status = wmeshspacedg_get_dofs_ids	(this->m_test_space,
							 itype,
							 idx_elm,
							 idofs_ids,
							 1);
		WMESH_STATUS_CHECK(status);
				
		status = wmeshspace_get_dofs_ids	(this->m_velocity_space,
							 itype,
							 idx_elm,
							 udofs_ids,
							 1);
		WMESH_STATUS_CHECK(status);
		
		//
		// Extract coordinates.
		//
		status = wmesh_get_cooelm(this->m_mesh,
					  itype,
					  idx_elm,
					  dofs_element_storage,
					  WMESH_MAT_FORWARD(dofs_element));
		WMESH_STATUS_CHECK(status);
		

		//
		// Extract velocity dofs.
		//
		status = wmeshspace_get_dof_values(*this->m_velocity_space,
						   itype,
						   idx_elm,
						   this->m_velocity_storage,
						   *this->m_velocity,
						   dofs_velocity_storage,
						   dofs_velocity);
		WMESH_STATUS_CHECK(status);
		//		std::cout << dofs_velocity << std::endl;

		for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
		  {
		    wmesh_int_t nei_info = this->m_mesh->m_c2c.m_data[this->m_mesh->m_c2c.m_ptr[itype] + this->m_mesh->m_c2c.m_ld[itype]*idx_elm + ifacet];
		    if (nei_info)  continue;
		    
		    const wmesh_int_t facet 	= facets[ifacet];
		    auto * flux_data = integral_flux_data[facet];
		    wmesh_mat_gemm(static_cast<T>(1),
				   dofs_element,
				   shape_eval_boundary_element[ifacet][0]->m_f,
				   static_cast<T>(0),
				   flux_data->m_dofs_element);

		    //
		    // Nasty skip 
		    //
		    bool v = true; 
		    for (wmesh_int_t i=0;i<flux_data->m_dofs_element.n;++i)
		      {
			if (flux_data->m_dofs_element.v[flux_data->m_dofs_element.ld*i+0] != 0)
			  {
			    v = false;
			    break;
			  }
		      }
		    if (!v)
		      {
			continue;
		      }
		    
		    wmesh_mat_gemm(static_cast<T>(1),
				   dofs_velocity,
				   shape_eval_boundary_velocity[ifacet][0]->m_f,
				   static_cast<T>(0),				   
				   flux_data->m_dofs_velocity);



		    //
		    // Get the coordinates of the quadrature.
		    // 
		    wmesh_mat_gemm(static_cast<T>(1),
				   flux_data->m_dofs_element,
				   this->m_integral_flux[facet]->m_shape_eval_element->m_f,
				   static_cast<T>(0),
				   flux_data->m_q_coo);
		    
		    for (wmesh_int_t i=0;i<flux_data->m_q_a.m;++i)
		      {
			T y = flux_data->m_q_coo.v[flux_data->m_q_coo.ld * i + 1];
			flux_data->m_q_a.v[i] = y;
		      }
		    
		    const wmesh_mat_t<T>& eval_btest 		= shape_eval_boundary_test[ifacet][0]->m_f;		      
		    status = this->m_integral_flux[facet]->eval_residual(*flux_data,
									 static_cast<T>(0),
									 flux_data->m_local_rhs,
									 flux_data->m_q_a.v,
									 1);
		    //		    std::cout << flux_data->m_local_rhs << std::endl;
		    const T r1 = static_cast<T>(1);
		    const T r0 = static_cast<T>(0);
		    const wmesh_int_t n1 = static_cast<wmesh_int_t>(1);
		    xgemv("N",
			  &eval_btest.m,
			  &eval_btest.n,

			  &r1,
			  eval_btest.v,
			  &eval_btest.ld,
			  flux_data->m_local_rhs.v,
			  &n1,

			  &r0,				   
			  local_rhs.v,
			  &n1);		    
		    //		    std::cout <<  "yo" <<  std::endl;
		    //		    std::cout << local_rhs << std::endl;
		    for (wmesh_int_t i=0;i<idofs_n;++i)
		      {
			rhs_[rhs_inc_*(idofs_ids[i]-1)] += local_rhs.v[i];
		      }
		  }
	      }
	  }
      }

    return WMESH_STATUS_SUCCESS;
  };

};


extern "C"
{

  wmesh_status_t wmeshspacedg_advection(const wmeshspacedg_t*__restrict__ 	self_,
					const wmesh_cubature_info_t* 		cubature_info_, 		
					const wmesh_shape_info_t* 		shape_info_element_, 		
					const wmesh_shape_info_t* 		shape_info_trial_,
					const wmesh_shape_info_t* 		shape_info_test_,
					const wmesh_shape_info_t* 		shape_info_velocity_,
					const wmeshspace_t *			velocity_space_,
					wmesh_int_t 				velocity_storage_,
					const wmesh_mat_t<double>* 		velocity_,
					wmesh_int_t				csr_size_,
					const_wmesh_int_p			csr_ptr_,
					const_wmesh_int_p			csr_ind_,
					double * 				csr_val_,
					double * 				rhs_)
  {
    std::cout << "define pde data ..." << std::endl;
    auto start_time0 = std::chrono::high_resolution_clock::now();
    wmesh_pde_advection_t<double> pde_advection(self_,
						self_,
						velocity_space_,
						velocity_storage_,
						velocity_,
						*cubature_info_,
						*shape_info_element_,
						*shape_info_trial_,
						*shape_info_test_,
						*shape_info_velocity_);
    
    std::cout << "define pde data done." << std::endl;
    wmesh_status_t status;

    auto end_time0 = std::chrono::high_resolution_clock::now();
    auto time0 = end_time0 - start_time0;
    std::cout << "elapsed  " << time0/std::chrono::milliseconds(1) << "ms.\n";      

    {
      auto start_time = std::chrono::high_resolution_clock::now();

      std::cout << "residual ..." << std::endl;
      status = pde_advection.residual(rhs_,
				    1);
      auto end_time = std::chrono::high_resolution_clock::now();
      auto time = end_time - start_time;
      std::cout << "elapsed  " << time/std::chrono::milliseconds(1) << "ms.\n";      
    }

    
    {
      auto start_time = std::chrono::high_resolution_clock::now();
      
      std::cout << "jacobian ..." << std::endl;
      status = pde_advection.jacobian(csr_size_,
				  csr_ptr_,
				  csr_ind_,
				  csr_val_);

    WMESH_STATUS_CHECK(status);
      auto end_time = std::chrono::high_resolution_clock::now();
      auto time = end_time - start_time;
      std::cout << "elapsed  " << time/std::chrono::milliseconds(1) << "ms.\n";      
    }


    return WMESH_STATUS_SUCCESS;
  };

};
