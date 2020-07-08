#include "bms.hpp"
#include "wmesh_shape_t.hpp"
#include "wmesh_cubature_t.hpp"
#include "wmesh_shape_eval_t.hpp"
#include "wmesh_cubature_info_t.hpp"
#include "wmesh_shape_info_t.hpp"
#include "wmesh_integral_convection_t.hpp"
#include "wmesh_shape_eval_factory_t.hpp"
#include "wmesh_cubature_factory_t.hpp"
#include "wmesh-blas.hpp"
#include "bms_templates.hpp"

template<typename T>
wmesh_template_integral_convection_t<T>::wmesh_template_integral_convection_t(wmesh_int_t 			element_,
									      const wmesh_cubature_info_t& 	cubature_info_,
									      const wmesh_shape_info_t& 		shape_info_element_,
									      const wmesh_shape_info_t& 		shape_info_velocity_,
									      const wmesh_shape_info_t& 		shape_info_trial_,
									      const wmesh_shape_info_t& 		shape_info_test_)
: m_shape_element(element_,shape_info_element_),
  m_shape_velocity(element_,shape_info_velocity_),
  m_shape_trial(element_,shape_info_trial_),
  m_shape_test(element_,shape_info_test_)
  
{
  wmesh_status_t status;

  wmesh_int_t topodim;
  status = bms_element2topodim(element_,
			       &topodim);  
  WMESH_STATUS_CHECK_EXIT(status);
  this->m_topodim = topodim;

  
  this->m_cubature			= wmesh_cubature_factory_t<T>::cubature_instance	(element_,
												 cubature_info_.m_family,
												 cubature_info_.m_degree);
  
  this->m_shape_eval_velocity  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(this->m_shape_velocity,
												 this->m_cubature);
  
  this->m_shape_eval_element  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(this->m_shape_element,
												 this->m_cubature);

  this->m_shape_eval_trial  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(this->m_shape_trial,
												 this->m_cubature);
  
  this->m_shape_eval_test  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(this->m_shape_test,
												 this->m_cubature);
  
  this->m_build_storage 		= WMESH_STORAGE_INTERLEAVE;
  const wmesh_mat_t<T>& q_weights 	= this->m_cubature->get_weights();	    
  const wmesh_int_t q_n 	= q_weights.n;
  const wmesh_int_t test_ndofs 	= this->m_shape_test.get_ndofs();
  const wmesh_int_t trial_ndofs = this->m_shape_trial.get_ndofs();

  wmesh_mat_t<T>::alloc(&this->m_build,
			test_ndofs * trial_ndofs,
			q_n * topodim);

  wmesh_mat_t<T>::zero(this->m_build);


  const wmesh_int_t inc_i = (this->m_shape_eval_test->m_diff_storage == WMESH_STORAGE_INTERLEAVE) ? 1 : this->m_shape_eval_test->m_diff[0].ld;
  const wmesh_int_t inc_j = (this->m_shape_eval_trial->m_diff_storage == WMESH_STORAGE_INTERLEAVE) ? 1 : this->m_shape_eval_trial->m_diff[0].ld;
  
  const wmesh_int_t ni = test_ndofs;
  const wmesh_int_t nj = trial_ndofs;

  for (wmesh_int_t k=0;k<q_n;++k)
    {
      const T wk = q_weights.v[q_weights.ld*k];
      const T * __restrict__  test_i = this->m_shape_eval_test->m_f.v + this->m_shape_eval_test->m_f.ld * k;
      for (wmesh_int_t idim=0;idim<topodim;++idim)
	{
	  xger(&ni,
	       &nj,
	       &wk,
	       test_i,
	       &inc_i,
	       this->m_shape_eval_trial->m_diff[idim].v + this->m_shape_eval_trial->m_diff[idim].ld * k,
	       &inc_j,
	       this->m_build.v + this->m_build.ld * ( k * topodim + idim),
	       &ni);
	}
    }      

};



template<typename T>
wmesh_template_integral_convection_t<T>::data_t::data_t(const wmesh_template_integral_convection_t<T>&  parent_)
{

  
  const wmesh_int_t q_n 		= parent_.m_cubature->get_num_points();
  const wmesh_int_t ndofs_velocity 	= parent_.m_shape_velocity.get_ndofs();
  const wmesh_int_t ndofs_element 	= parent_.m_shape_element.get_ndofs();
  const wmesh_int_t ndofs_test 		= parent_.m_shape_test.get_ndofs();
  const wmesh_int_t ndofs_trial 	= parent_.m_shape_trial.get_ndofs();
  const wmesh_int_t topodim 		= parent_.m_topodim;
  
  this->m_dofs_element_storage  = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&this->m_dofs_element, topodim, ndofs_element);

  this->m_dofs_velocity_storage  = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&this->m_dofs_velocity, topodim, ndofs_velocity);
  
  this->m_q_velocity_storage 		= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&this->m_q_velocity, topodim, q_n);

  this->m_q_jacobians_storage 	= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&this->m_q_jacobians,topodim*topodim, q_n);
  
  wmesh_mat_t<T>::alloc(&this->m_q_jacobians_det,1,q_n);
  wmesh_mat_t<T>::alloc(&this->m_build_rhs,1, q_n * (topodim*topodim));
  wmesh_mat_t<T>::alloc(&this->m_local_matrix,ndofs_test,ndofs_trial);

  this->m_local_rhs = (T*__restrict__)malloc(sizeof(T)*ndofs_test);
}

template<typename T>
wmesh_status_t wmesh_template_integral_convection_t<T>::eval(wmesh_template_integral_convection_t<T>::data_t&  		data_,
							    const T 							alpha_,
							    wmesh_mat_t<T>& 						local_matrix_) const
{
  
  //
  //  Evaluate the velocity.
  //
  wmesh_mat_gemm(static_cast<T>(1),
		 data_.m_dofs_velocity,
		 this->m_shape_eval_velocity->m_f,
		 static_cast<T>(0),
		 data_.m_q_velocity);
#if 0
  wmesh_mat_gemm(this->m_shape_eval_velocity->m_f_storage,
		 data_.m_dofs_velocity_storage,
		 data_.m_q_velocity_storage,
		 
		 static_cast<T>(1),
		 this->m_shape_eval_velocity->m_f,
		 data_.m_dofs_velocity,
		 
		 static_cast<T>(0),
		 data_.m_q_velocity);
#endif  
  //
  // Compute the jacobian.
  //
  wmesh_status_t status;
  status =  bms_element_jacobians(this->m_element,
				  data_.m_dofs_element_storage,
				  data_.m_dofs_element,
				  this->m_shape_eval_element->m_diff_storage,
				  &this->m_shape_eval_element->m_diff[0],
				  data_.m_q_jacobians,
				  data_.m_q_jacobians_det);
  WMESH_STATUS_CHECK(status);


  
  //
  // Transform the velocity.
  //

  const wmesh_int_t n1 		= static_cast<wmesh_int_t>(1);
  const T r1 			= static_cast<T>(1);
  const T r0 			= static_cast<T>(0);


  const wmesh_int_t q_n 	= this->m_cubature->get_num_points();
  const wmesh_int_t topodim 	= this->m_topodim;

  //
  // This should be 'batched' with other technologies.
  // 
  T tmpu[3];
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      const T * __restrict__ q_J    = data_.m_q_jacobians.v + data_.m_q_jacobians.ld * k;
      const T * __restrict__ det_q_J = data_.m_q_jacobians_det.v + data_.m_q_jacobians_det.ld * k;

      if (data_.m_q_velocity_storage==WMESH_STORAGE_BLOCK)
	{
	  for (wmesh_int_t i=0;i<topodim;++i)
	    {
	      tmpu[i] = *(data_.m_q_velocity.v + data_.m_q_velocity.ld*i+k);
	    }
	}
      else
	{
	  for (wmesh_int_t i=0;i<topodim;++i)
	    {
	      tmpu[i] = *(data_.m_q_velocity.v + data_.m_q_velocity.ld*k+i);
	    }
	}
      
      xgemv("N",
	    &topodim,
	    &topodim,
	    det_q_J,
	    q_J,
	    &topodim,
	    tmpu,
	    &n1,
	    &r0,
	    data_.m_build_rhs.v + k * topodim,
	    &n1);
      
    }

  //
  // Build the matrix.
  //
  if (local_matrix_.ld == local_matrix_.m)
    {
      xgemv("N",
	    &this->m_build.m,
	    &this->m_build.n,
	    &r1,
	    this->m_build.v,
	    &this->m_build.ld,
	    data_.m_build_rhs.v,
	    &n1,
	    &alpha_,
	    local_matrix_.v,
	    &n1);
    }
  else
    {
      std::cerr << "need a correction." << std::endl;
      exit(1);
    }

  return WMESH_STATUS_SUCCESS;
};

template
struct wmesh_template_integral_convection_t<float>;
template
struct wmesh_template_integral_convection_t<double>;
