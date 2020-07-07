#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_t.hpp"
#include <chrono>
#include <iostream>
#include "bms.h"
#include "wmesh-utils.hpp"
#include "wmesh-blas.h"
#include "bms_templates.hpp"

#include "wmesh_cubature_factory_t.hpp"
#include "wmesh_cubature_info_t.hpp"
#include "wmesh_shape_eval_factory_t.hpp"
#include "wmesh_integral_diffusion_t.hpp"
#include "bms.hpp"
#include "wmeshspace_t.hpp"

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
struct wmesh_pde_laplace_t : public wmesh_pde_t< wmesh_pde_laplace_t<T>, T >
{
  const wmeshspace_t * 	m_trial_space;
  const wmeshspace_t * 	m_test_space;
  wmesh_cubature_info_t m_cubature_info;
  wmesh_shape_info_t 	m_shape_info_element;
  wmesh_shape_info_t 	m_shape_info_trial;
  wmesh_shape_info_t 	m_shape_info_test;
  wmesh_shape_info_t 	m_shape_info_a;
  
  wmesh_template_integral_diffusion_t<T> * m_integral_diffusion[WMESH_ELEMENT_ALL]{};
  
  wmesh_pde_laplace_t(const wmeshspace_t * 		trial_space_,
		      const wmeshspace_t * 		test_space_,
		      const wmesh_cubature_info_t& 	cubature_info_, 		
		      const wmesh_shape_info_t& 	shape_info_element_, 		
		      const wmesh_shape_info_t& 	shape_info_trial_,
		      const wmesh_shape_info_t& 	shape_info_test_,
		      const wmesh_shape_info_t& 	shape_info_a_)
    : m_trial_space(trial_space_),
      m_test_space(test_space_),
      m_cubature_info(cubature_info_),
      m_shape_info_element(shape_info_element_),
      m_shape_info_trial(shape_info_trial_),
      m_shape_info_test(shape_info_test_),
      m_shape_info_a(shape_info_a_)
  {
    const wmesh_int_t 	topodim = trial_space_->m_mesh->m_topology_dimension;      
    const wmesh_int_t 	ntypes 	= trial_space_->m_mesh->m_c2n.m_size;      
    for (wmesh_int_t itype=0;itype<ntypes;++itype)
      {
	wmesh_int_t element = (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );
	const wmesh_int_t num_elements = trial_space_->m_mesh->m_c2n.m_n[itype];
	if (num_elements > 0)
	  {
	    this->m_integral_diffusion[element] = new wmesh_template_integral_diffusion_t<T>(element, 
											     cubature_info_, 
											     this->m_shape_info_element, 
											     this->m_shape_info_a, 
											     this->m_shape_info_test, 
											     this->m_shape_info_trial);
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
    using data_t = typename wmesh_template_integral_diffusion_t<T>::data_t;
    data_t * integral_diffusion_data[WMESH_ELEMENT_ALL]{};
    
    const wmesh_int_t 	ntypes 	= this->m_trial_space->m_mesh->m_c2n.m_size;      
    for (wmesh_int_t itype=0;itype<ntypes;++itype)
      {
	wmesh_int_t element = (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );
	const wmesh_int_t num_elements = this->m_trial_space->m_mesh->m_c2n.m_n[itype];
	if (num_elements > 0)
	  {
	    integral_diffusion_data[element] = new data_t(*this->m_integral_diffusion[element]);
	  }
      } 

    const wmesh_int_t 	idofs_inc = 1;
    const wmesh_int_t 	jdofs_inc = 1;    
    const wmesh_int_sparsemat_t& trial_c2d = this->m_trial_space->m_c2d;
    const wmesh_int_sparsemat_t& test_c2d = this->m_trial_space->m_c2d;
    for (wmesh_int_t itype=0;itype<trial_c2d.m_size;++itype)
      {
	if (trial_c2d.m_n[itype] > 0)
	  {
	    const wmesh_int_t element 	= (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );

	    const wmesh_int_t trial_c2d_m  	= trial_c2d.m_m[itype];
	    const wmesh_int_t trial_c2d_n  	= trial_c2d.m_n[itype];
	    const wmesh_int_t trial_c2d_ld 	= trial_c2d.m_ld[itype];
	    const_wmesh_int_p trial_c2d_v 	= trial_c2d.m_data + trial_c2d.m_ptr[itype];
	    
	    const wmesh_int_t test_c2d_m  	= test_c2d.m_m[itype];
	    const wmesh_int_t test_c2d_ld 	= test_c2d.m_ld[itype];
	    const_wmesh_int_p test_c2d_v 	= test_c2d.m_data + test_c2d.m_ptr[itype];

#ifndef NDEBUG
	    const wmesh_int_t test_c2d_n  	= test_c2d.m_n[itype];
	    WMESH_CHECK(test_c2d_n == trial_c2d_n);
#endif

	    auto * data = integral_diffusion_data[element];
	    const wmesh_int_t 	idofs_n 	= test_c2d_m;
	    const wmesh_int_t 	jdofs_n 	= trial_c2d_m;		
	    
	    for (wmesh_int_t idx_elm=0;idx_elm<trial_c2d_n;++idx_elm)
	      {
		//
		// Get idofs and jdofs from test and trial spaces.
		//
		const_wmesh_int_p 	idofs 	= test_c2d_v + test_c2d_ld * idx_elm;
		const_wmesh_int_p 	jdofs 	= trial_c2d_v + trial_c2d_ld * idx_elm;
		
		//
		// Extract coordinates.
		//
		wmesh_status_t status = wmesh_get_cooelm(this->m_trial_space->m_mesh,
							 itype,
							 idx_elm,
							 data->m_dofs_element_storage,
							 WMESH_MAT_FORWARD(data->m_dofs_element));
		WMESH_STATUS_CHECK(status);
		
#if 0
		//
		// Extract coeffs.
		//
		status = wmeshspace_get_dofselm(this->m_a_space,
						itype,
						idx_ielm,
						coeffs_storage,
						WMESH_MAT_FORWARD(coeffs));
		WMESH_STATUS_CHECK(status);
#endif
		
		//
		// Compute local matrix.
		//
		status = this->m_integral_diffusion[element]->eval(*data,
								   static_cast<T>(0),
								   data->m_local_matrix);
#if 0
		std::cout << "local matrix" <<std::endl;
		std::cout << data->m_local_matrix << std::endl;
		exit(1);
#endif
		//
		// Add the local matrix.
		//

		status =  bms_sparse_add(idofs_n,
					 idofs,
					 idofs_inc,
					 
					 jdofs_n,
					 jdofs,
					 jdofs_inc,

					 data->m_local_matrix.v,
					 data->m_local_matrix.ld,
					 
					 csr_size_,
					 csr_ptr_,
					 csr_ind_,
					 csr_val_);

	      }

	  }
	
      }

    {
      auto start_time = std::chrono::high_resolution_clock::now();
      
      std::cout << "apply bc ..." << std::endl;
      
      //
      // Now apply bc.
      //
      for (wmesh_int_t i=0;i<this->m_trial_space->m_ndofs;++i)
	{
	  if (this->m_trial_space->m_dof_codes[i] == 1)
	    {
	      for (wmesh_int_t k = csr_ptr_[i];k<csr_ptr_[i+1];++k)
		{
		  wmesh_int_t j = csr_ind_[k]-1;
		  if (j==i)
		    {
		      csr_val_[k] = 1.0;
		    }
		  else
		    {
		      csr_val_[k] = 0.0;
		    }
		}
	      
	      //	  rhs_[i] = 0.0;
	      
	    }
	}

      auto end_time = std::chrono::high_resolution_clock::now();
      auto time = end_time - start_time;
      std::cout << "elapsed  " << time/std::chrono::milliseconds(1) << "ms.\n";      
    }

    
    return WMESH_STATUS_SUCCESS;
  };

  wmesh_status_t residual(T * __restrict__		rhs_,
			  wmesh_int_t 			rhs_inc_)
  {

  for (wmesh_int_t i=0;i<this->m_trial_space->m_ndofs;++i)
    {     
      rhs_[i*rhs_inc_] = 0.0;
    }

    const wmesh_int_t 	topodim = this->m_test_space->m_mesh->m_topology_dimension;      


    //
    // Loop over the trial space.
    //
    using data_t = typename wmesh_template_integral_diffusion_t<T>::data_t;
    data_t * integral_diffusion_data[WMESH_ELEMENT_ALL]{};
    
    const wmesh_int_t 	ntypes 	= this->m_test_space->m_mesh->m_c2n.m_size;      
    for (wmesh_int_t itype=0;itype<ntypes;++itype)
      {
	wmesh_int_t element = (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );
	const wmesh_int_t num_elements = this->m_test_space->m_mesh->m_c2n.m_n[itype];
	if (num_elements > 0)
	  {
	    integral_diffusion_data[element] = new data_t(*this->m_integral_diffusion[element]);
	  }
      } 

    const wmesh_int_t 	idofs_inc = 1;


    const wmesh_int_sparsemat_t& test_c2d = this->m_test_space->m_c2d;
    for (wmesh_int_t itype=0;itype<test_c2d.m_size;++itype)
      {
	if (test_c2d.m_n[itype] > 0)
	  {
	    const wmesh_int_t element 	= (topodim==3) ? (4+itype) : ( (topodim==2) ? 2+itype : 1+itype );
	    
	    const wmesh_int_t test_c2d_m  	= test_c2d.m_m[itype];
	    const wmesh_int_t test_c2d_n  	= test_c2d.m_n[itype];
	    const wmesh_int_t test_c2d_ld 	= test_c2d.m_ld[itype];
	    const_wmesh_int_p test_c2d_v 	= test_c2d.m_data + test_c2d.m_ptr[itype];

	    const wmesh_int_t test_ndofs  	= test_c2d_m;

	    auto * data = integral_diffusion_data[element];


	    const wmesh_mat_t<T>& eval_test = this->m_integral_diffusion[element]->m_shape_eval_test->m_f;

	    const wmesh_int_t 		q_n 	= this->m_integral_diffusion[element]->m_cubature->m_w.n;
	    const T * __restrict__ 	q_w 	= this->m_integral_diffusion[element]->m_cubature->m_w.v;
	    const wmesh_int_t q_w_inc 		= this->m_integral_diffusion[element]->m_cubature->m_w.ld;

	    
	    for (wmesh_int_t idx_elm=0;idx_elm<test_c2d_n;++idx_elm)
	      {
		//
		// Get idofs and jdofs from test and trial spaces.
		//
		const_wmesh_int_p 	idofs 	= test_c2d_v + test_c2d_ld * idx_elm;
		
		//
		// Extract coordinates.
		//
		wmesh_status_t status = wmesh_get_cooelm(this->m_test_space->m_mesh,
							 itype,
							 idx_elm,
							 data->m_dofs_element_storage,
							 WMESH_MAT_FORWARD(data->m_dofs_element));
		WMESH_STATUS_CHECK(status);
		
#if 0
		//
		// Extract coeffs.
		//
		status = wmeshspace_get_dofselm(this->m_a_space,
						itype,
						idx_ielm,
						coeffs_storage,
						WMESH_MAT_FORWARD(coeffs));
		WMESH_STATUS_CHECK(status);
		
#endif

#if 0		
		//
		// Compute local matrix.
		//
		status = this->m_integral_diffusion[element]->eval(*data,
								   static_cast<T>(0),
								   data->m_local_matrix);
#endif
		

		status =  bms_element_jacobians(this->m_integral_diffusion[element]->m_element,
						data->m_dofs_element_storage,
						data->m_dofs_element,
						this->m_integral_diffusion[element]->m_shape_eval_element->m_diff_storage,
						&this->m_integral_diffusion[element]->m_shape_eval_element->m_diff[0],
						data->m_q_jacobians,
						data->m_q_jacobians_det);
		WMESH_STATUS_CHECK(status);

		//
		// data->m_q_jacobians_det.v[data->m_q_jacobians_det.ld * k + 0]
		//
		
		//
		// Reset the local rhs.
		//
		for (wmesh_int_t i=0;i<test_ndofs;++i)
		  {
		    data->m_local_rhs[i] = static_cast<T>(0);	      
		  }
		
		for (wmesh_int_t k=0;k<q_n;++k)
		  {
		    for (wmesh_int_t i=0;i<test_ndofs;++i)
		      {
			data->m_local_rhs[i] += eval_test.v[k*eval_test.ld + i] * data->m_q_jacobians_det.v[data->m_q_jacobians_det.ld * k + 0] * q_w[q_w_inc * k];
		      }      
		  }

		for (wmesh_int_t i=0;i<test_ndofs;++i)
		  {
		    rhs_[ ( idofs[idofs_inc * i] - 1 ) ] += data->m_local_rhs[i];
		  }
		
	      }
	  }
      }

    //
    // apply boundary conditions.
    //
    for (wmesh_int_t i=0;i<this->m_trial_space->m_ndofs;++i)
      {
	if (this->m_trial_space->m_dof_codes[i] == 1)
	  {	    
	    rhs_[i] = 0.0;	    
	  }
      }
    
    return WMESH_STATUS_SUCCESS;
  };

};


extern "C"
{
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
				    double * 				rhs_)
  {
    std::cout << "define pde data ..." << std::endl;
      auto start_time0 = std::chrono::high_resolution_clock::now();
    wmesh_pde_laplace_t<double> pde_laplace(self_,
					    self_,
					    *cubature_info_,
					    *shape_info_element_,
					    *shape_info_trial_,
					    *shape_info_test_,
					    *shape_info_a_);
    
    std::cout << "define pde data done." << std::endl;
    wmesh_status_t status;

    auto end_time0 = std::chrono::high_resolution_clock::now();
    auto time0 = end_time0 - start_time0;
    std::cout << "elapsed  " << time0/std::chrono::milliseconds(1) << "ms.\n";      

    {
      auto start_time = std::chrono::high_resolution_clock::now();

      std::cout << "residual ..." << std::endl;
      status = pde_laplace.residual(rhs_,
				    1);
      auto end_time = std::chrono::high_resolution_clock::now();
      auto time = end_time - start_time;
      std::cout << "elapsed  " << time/std::chrono::milliseconds(1) << "ms.\n";      
    }

    
    {
      auto start_time = std::chrono::high_resolution_clock::now();
      
      std::cout << "jacobian ..." << std::endl;
      status = pde_laplace.jacobian(csr_size_,
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
