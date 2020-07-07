#include "bms.hpp"
#include "wmesh_shape_t.hpp"
#include "wmesh_cubature_t.hpp"
#include "wmesh_shape_eval_t.hpp"
#include "wmesh_cubature_info_t.hpp"
#include "wmesh_shape_info_t.hpp"
#include "wmesh_integral_diffusion_t.hpp"
#include "wmesh_shape_eval_factory_t.hpp"
#include "wmesh_cubature_factory_t.hpp"
#include "wmesh-blas.hpp"
#include "bms_templates.hpp"
#if 0
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
#endif

template<typename T>
wmesh_template_integral_diffusion_t<T>::wmesh_template_integral_diffusion_t(wmesh_int_t 			element_,
									    const wmesh_cubature_info_t& 	cubature_info_,
									    const wmesh_shape_info_t& 		shape_info_element_,
									    const wmesh_shape_info_t& 		shape_info_a_,
									    const wmesh_shape_info_t& 		shape_info_trial_,
									    const wmesh_shape_info_t& 		shape_info_test_)
{
  wmesh_status_t status;
  
  status = wmesh_shape_def(&this->m_shape_element,
			   element_,
			   shape_info_element_.m_family,
			   shape_info_element_.m_degree);
  WMESH_STATUS_CHECK_FAIL(status);
  
  status = wmesh_shape_def(&this->m_shape_a,
			   element_,
			   shape_info_a_.m_family,
			   shape_info_a_.m_degree);
  WMESH_STATUS_CHECK_FAIL(status);
  
  status = wmesh_shape_def(&this->m_shape_trial,
			   element_,
			   shape_info_trial_.m_family,
			   shape_info_trial_.m_degree);
  WMESH_STATUS_CHECK_FAIL(status);

  status = wmesh_shape_def(&this->m_shape_test,
			   element_,
			   shape_info_test_.m_family,
			   shape_info_test_.m_degree);
  WMESH_STATUS_CHECK_FAIL(status);

  wmesh_int_t topodim;
  status = bms_element2topodim(element_,
			       &topodim);  
  WMESH_STATUS_CHECK_FAIL(status);
  this->m_topodim = topodim;

  
  this->m_cubature		= wmesh_cubature_factory_t<T>::cubature_instance	(element_,
											 cubature_info_.m_family,
											 cubature_info_.m_degree);
#if 0
  this->m_shape_eval_a  	= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(this->m_shape_a,
											 this->m_cubature);
#endif
  this->m_shape_eval_element  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(this->m_shape_element,
												 this->m_cubature);
  

  this->m_shape_eval_trial  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(this->m_shape_trial,
												 this->m_cubature);
  
  this->m_shape_eval_test  		= wmesh_shape_eval_factory_t<T>::shape_eval_instance	(this->m_shape_test,
												 this->m_cubature);
  
  this->m_build_storage 		= WMESH_STORAGE_INTERLEAVE;
  
  const wmesh_int_t q_n 	= this->m_cubature->m_w.n;
  const wmesh_int_t test_ndofs 	= this->m_shape_test.m_ndofs;
  const wmesh_int_t trial_ndofs = this->m_shape_trial.m_ndofs;

  wmesh_mat_t<T>::alloc(&this->m_build,
			test_ndofs * trial_ndofs,
			q_n * topodim * topodim);
  wmesh_mat_t<T>::zero(this->m_build);

  const wmesh_int_t inc_i = (this->m_shape_eval_test->m_diff_storage == WMESH_STORAGE_INTERLEAVE) ? 1 : this->m_shape_eval_test->m_diff[0].ld;
  const wmesh_int_t inc_j = (this->m_shape_eval_trial->m_diff_storage == WMESH_STORAGE_INTERLEAVE) ? 1 : this->m_shape_eval_trial->m_diff[0].ld;
  
  const wmesh_int_t ni = test_ndofs;
  const wmesh_int_t nj = trial_ndofs;
  
#if 1
  const T * __restrict__ trial_nabla_j[3]{};
  const T * __restrict__ test_nabla_i[3]{};
  const wmesh_int_t trial_nabla_ld = this->m_shape_eval_trial->m_diff[0].ld;
  const wmesh_int_t test_nabla_ld = this->m_shape_eval_test->m_diff[0].ld;

  for (wmesh_int_t idim=0;idim<topodim;++idim)
    {
      trial_nabla_j[idim] = this->m_shape_eval_trial->m_diff[idim].v;
    }
  for (wmesh_int_t idim=0;idim<topodim;++idim)
    {
      test_nabla_i[idim] = this->m_shape_eval_test->m_diff[idim].v;
    }
  
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      const T wk = this->m_cubature->m_w.v[this->m_cubature->m_w.ld*k]; // static_cast<T>(1);
      for (wmesh_int_t jdim=0;jdim<topodim;++jdim)
	{
	  const T * __restrict__  dj = trial_nabla_j[jdim] + k * trial_nabla_ld;
	  for (wmesh_int_t idim=0;idim<topodim;++idim)
	    {
	      const T * __restrict__  di = test_nabla_i[idim]  + k * test_nabla_ld;	      
	      xger(&ni, &nj, &wk, di, &inc_i, dj, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + (jdim*topodim + idim)) , &ni);
	    }
	}
    }
#endif
  
#if 0 
  if (topodim==2)
    {
      for (wmesh_int_t k=0;k<q_n;++k)
	{
      const T wk = this->m_cubature->m_w.v[this->m_cubature->m_w.ld*k]; // static_cast<T>(1);
      //	  const T wk = static_cast<T>(1);//this->m_cubature->m_w.v[this->m_cubature->m_w.ld*k];
	  const T * __restrict__  dr_j = this->m_shape_eval_trial->m_diff[0].v + this->m_shape_eval_trial->m_diff[0].ld * k;
	  const T * __restrict__  ds_j = this->m_shape_eval_trial->m_diff[1].v + this->m_shape_eval_trial->m_diff[1].ld * k;
	  const T * __restrict__  dr_i = this->m_shape_eval_test->m_diff[0].v + this->m_shape_eval_test->m_diff[0].ld * k;
	  const T * __restrict__  ds_i = this->m_shape_eval_test->m_diff[1].v + this->m_shape_eval_test->m_diff[1].ld * k;
	  
	  xger(&ni, &nj, &wk, dr_i, &inc_i, dr_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 0) , &ni);
	  xger(&ni, &nj, &wk, ds_i, &inc_i, dr_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 1) , &ni);

	  xger(&ni, &nj, &wk, dr_i, &inc_i, ds_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 2) , &ni);
	  xger(&ni, &nj, &wk, ds_i, &inc_i, ds_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 3) , &ni);
	}
    }
  else if (topodim==3)
    {
      for (wmesh_int_t k=0;k<q_n;++k)
	{
	  const T wk = static_cast<T>(1);//this->m_cubature->m_w.v[this->m_cubature->m_w.ld*k];

	  const T * __restrict__  dr_j = this->m_shape_eval_trial->m_diff[0].v + this->m_shape_eval_trial->m_diff[0].ld * k;
	  const T * __restrict__  ds_j = this->m_shape_eval_trial->m_diff[1].v + this->m_shape_eval_trial->m_diff[1].ld * k;
	  const T * __restrict__  dt_j = this->m_shape_eval_trial->m_diff[2].v + this->m_shape_eval_trial->m_diff[2].ld * k;
	  const T * __restrict__  dr_i = this->m_shape_eval_test->m_diff[0].v + this->m_shape_eval_test->m_diff[0].ld * k;
	  const T * __restrict__  ds_i = this->m_shape_eval_test->m_diff[1].v + this->m_shape_eval_test->m_diff[1].ld * k;
	  const T * __restrict__  dt_i = this->m_shape_eval_test->m_diff[2].v + this->m_shape_eval_test->m_diff[2].ld * k;
	  
	  xger(&ni, &nj, &wk, dr_i, &inc_i, dr_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 0) , &ni);
	  xger(&ni, &nj, &wk, ds_i, &inc_i, dr_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 1) , &ni);
	  xger(&ni, &nj, &wk, dt_i, &inc_i, dr_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 2) , &ni);
	  
	  xger(&ni, &nj, &wk, dr_i, &inc_i, ds_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 3) , &ni);
	  xger(&ni, &nj, &wk, ds_i, &inc_i, ds_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 4) , &ni);
	  xger(&ni, &nj, &wk, dt_i, &inc_i, ds_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 5) , &ni);
	  
	  xger(&ni, &nj, &wk, dr_i, &inc_i, dt_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 6) , &ni);
	  xger(&ni, &nj, &wk, ds_i, &inc_i, dt_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 7) , &ni);
	  xger(&ni, &nj, &wk, dt_i, &inc_i, dt_j, &inc_j, this->m_build.v + this->m_build.ld * ( k * (topodim*topodim) + 8) , &ni);	  
	}
    }
#endif
  
};



template<typename T>
wmesh_template_integral_diffusion_t<T>::data_t::data_t(const wmesh_template_integral_diffusion_t<T>&  parent_)
{

  
  const wmesh_int_t q_n 		= parent_.m_cubature->m_w.n;
  //  const wmesh_int_t ndofs_a 		= parent_.m_shape_a.m_ndofs;
  const wmesh_int_t ndofs_element 	= parent_.m_shape_element.m_ndofs;
  const wmesh_int_t ndofs_test 		= parent_.m_shape_test.m_ndofs;
  const wmesh_int_t ndofs_trial 	= parent_.m_shape_trial.m_ndofs;
  const wmesh_int_t topodim 		= parent_.m_topodim;

  //  this->m_dofs_a_storage    = WMESH_STORAGE_INTERLEAVE;
  //  wmesh_mat_t<T>::alloc(&this->m_dofs_a, 1, ndofs_a);


  
  this->m_dofs_element_storage  = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&this->m_dofs_element, topodim, ndofs_element);

  //  this->m_q_a_storage 		= WMESH_STORAGE_INTERLEAVE;
  //  wmesh_mat_t<T>::alloc(&this->m_q_a, 1, q_n);

  this->m_q_jacobians_storage 	= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&this->m_q_jacobians,topodim*topodim, q_n);
  
  wmesh_mat_t<T>::alloc(&this->m_q_jacobians_det,1,q_n);
  wmesh_mat_t<T>::alloc(&this->m_build_rhs,1, q_n * (topodim*topodim));
  wmesh_mat_t<T>::alloc(&this->m_local_matrix,ndofs_test,ndofs_trial);
  this->m_local_rhs = (T*__restrict__)malloc(sizeof(T)*ndofs_test);
}


template<typename T>
wmesh_status_t wmesh_template_integral_diffusion_t<T>::eval(wmesh_template_integral_diffusion_t<T>::data_t&  		data_,
							    const T 							alpha_,
							    wmesh_mat_t<T>& 						local_matrix_) const
{
  
  //
  //  Evaluate the coefficients.
  //
#if 0
  wmesh_mat_gemm(this->m_shape_eval_a->m_f_storage,
		 data_.m_dofs_a_storage,
		 data_.m_q_a_storage,
		 
		 static_cast<T>(1),
		 this->m_shape_eval_a->m_f,
		 data_.m_dofs_a,
		 
		 static_cast<T>(0),
		 data_.m_q_a);
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
  const wmesh_int_t n1 		= static_cast<wmesh_int_t>(1);
  const T r1 			= static_cast<T>(1);
  const T r0 			= static_cast<T>(0);
  const wmesh_int_t q_n 	= this->m_cubature->m_w.n;
  const wmesh_int_t topodim 	= this->m_topodim;
  T a[9], b[9];

  for (wmesh_int_t k=0;k<q_n;++k)
    {
      for (wmesh_int_t i=0;i<topodim*topodim;++i)
	{
	  a[i]=data_.m_q_jacobians.v[data_.m_q_jacobians.ld*k+i];	  
	}

      for (wmesh_int_t i=0;i<topodim*topodim;++i)
	{
	  b[i]=data_.m_q_jacobians.v[data_.m_q_jacobians.ld*k+i];	  
	}
      
      T beta = data_.m_q_jacobians_det.v[data_.m_q_jacobians_det.ld*k];
      xgemm("N",
	    "T",
	    &topodim,
	    &topodim,
	    &topodim,
	    &beta,
	    a,
	    &topodim,
	    b,
	    &topodim,
	    &r0,
	    data_.m_build_rhs.v + k * topodim*topodim,
	    &topodim);      
    }

#if 0
  std::cout << "build rhs" <<std::endl;
  std::cout << data_.m_build_rhs << std::endl;
#endif
#if 0
  for (wmesh_int_t i=0;i<this->m_build.m;++i)
    {
      local_matrix_.v[i]=0;
    }
#endif

  //  std::cout << "build_rhs "  << std::endl;
  //  std::cout << data_.m_build_rhs << std::endl;
  if (local_matrix_.ld == local_matrix_.m)
    {

#if 0
      std::cout << "mv local_matrix, alpha " <<  alpha_ << std::endl;
      std::cout << "mv local_matrix, m " <<  local_matrix_.m << std::endl;
      std::cout << "mv local_matrix, n " <<  local_matrix_.n << std::endl;
      std::cout << "mv local_matrix, ld " <<  local_matrix_.ld << std::endl;
      std::cout << "mv build, m " <<  this->m_build.m << std::endl;
      std::cout << "mv build, n " <<  this->m_build.n << std::endl;
      std::cout << "mv build, ld " <<  this->m_build.ld << std::endl;
      std::cout << "mv build, m " <<  data_.m_build_rhs.m << std::endl;
      std::cout << "mv build, n " <<  data_.m_build_rhs.n << std::endl;
      std::cout << "mv build, ld " <<  data_.m_build_rhs.ld << std::endl;
      for (wmesh_int_t i=0;i<this->m_build.m;++i)
	{
	  //	  fprintf(stdout," %e\n",local_matrix_.v[i]);
	  //	  std::cout << local_matrix_.v[i] << std::endl;
	}
      //      std::cout << data_.m_build_rhs << std::endl;
      //      std::cout << this->m_build << std::endl;
      //      std::cout << r1 << std::endl;
      //      std::cout << n1 << std::endl;
#endif
      //      std::cout << this->m_build.m << " " << this->m_build.n << std::endl;
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
struct wmesh_template_integral_diffusion_t<float>;
template
struct wmesh_template_integral_diffusion_t<double>;
