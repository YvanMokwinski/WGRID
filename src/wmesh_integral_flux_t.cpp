#include "bms.hpp"
#include "wmesh_shape_t.hpp"
#include "wmesh_cubature_t.hpp"
#include "wmesh_shape_eval_t.hpp"
#include "wmesh_cubature_info_t.hpp"
#include "wmesh_shape_info_t.hpp"
#include "wmesh_integral_flux_t.hpp"
#include "wmesh_shape_eval_factory_t.hpp"
#include "wmesh_cubature_factory_t.hpp"
#include "wmesh-blas.hpp"
#include "bms_templates.hpp"
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

template<typename T>
wmesh_template_integral_flux_t<T>::wmesh_template_integral_flux_t(wmesh_int_t 				element_,
								  const wmesh_cubature_info_t& 		cubature_info_,
								  const wmesh_shape_info_t& 		shape_info_element_,
								  const wmesh_shape_info_t& 		shape_info_velocity_,
								  const wmesh_shape_info_t& 		shape_info_trial_,
								  const wmesh_shape_info_t& 		shape_info_test_)
{
  wmesh_status_t status;
  
  status = wmesh_shape_def(&this->m_shape_element,
			   element_,
			   shape_info_element_.m_family,
			   shape_info_element_.m_degree);
  WMESH_STATUS_CHECK_FAIL(status);
  
  status = wmesh_shape_def(&this->m_shape_velocity,
			   element_,
			   shape_info_velocity_.m_family,
			   shape_info_velocity_.m_degree);
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
  this->m_build_residual_storage	= WMESH_STORAGE_INTERLEAVE;
  
  const wmesh_int_t q_n 	= this->m_cubature->m_w.n;
  const wmesh_int_t ni 		= this->m_shape_test.m_ndofs;
  const wmesh_int_t nj 		= this->m_shape_trial.m_ndofs;
  const wmesh_int_t n1 		= 1;

  wmesh_mat_t<T>::alloc(&this->m_build_residual,
			ni,
			q_n);

  wmesh_mat_t<T>::zero(this->m_build_residual);

  wmesh_mat_t<T>::alloc(&this->m_build,
			ni * nj,
			q_n);

  wmesh_mat_t<T>::zero(this->m_build);

  const wmesh_int_t inc_i = (this->m_shape_eval_test->m_diff_storage == WMESH_STORAGE_INTERLEAVE) ? 1 : this->m_shape_eval_test->m_diff[0].ld;
  const wmesh_int_t inc_j = (this->m_shape_eval_trial->m_diff_storage == WMESH_STORAGE_INTERLEAVE) ? 1 : this->m_shape_eval_trial->m_diff[0].ld;
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      const T wk = this->m_cubature->m_w.v[this->m_cubature->m_w.ld*k];
      const T * __restrict__  test_i  = this->m_shape_eval_test->m_f.v + this->m_shape_eval_test->m_f.ld * k;
      const T * __restrict__  trial_j = this->m_shape_eval_trial->m_f.v + this->m_shape_eval_trial->m_f.ld * k;
      xger(&ni,
	   &nj,
	   &wk,
	   test_i,
	   &inc_i,
	   trial_j,
	   &inc_j,
	   this->m_build.v + this->m_build.ld * k,
	   &ni);
    }      

  for (wmesh_int_t k=0;k<q_n;++k)
    {
      const T wk = this->m_cubature->m_w.v[this->m_cubature->m_w.ld*k];
      const T * __restrict__  test_i  = this->m_shape_eval_test->m_f.v + this->m_shape_eval_test->m_f.ld * k;
      xcopy(&ni,test_i,&n1,this->m_build_residual.v + this->m_build_residual.ld * k,&n1);
      xscal(&ni,&wk,this->m_build_residual.v + this->m_build_residual.ld * k,&n1);
    }      
};

template
struct wmesh_template_integral_flux_t<float>;
template
struct wmesh_template_integral_flux_t<double>;



template<typename T>
wmesh_template_integral_flux_t<T>::data_t::data_t(const wmesh_template_integral_flux_t<T>&  parent_)
{

  
  const wmesh_int_t q_n 		= parent_.m_cubature->m_w.n;
  const wmesh_int_t ndofs_velocity 	= parent_.m_shape_velocity.m_ndofs;
  const wmesh_int_t ndofs_element 	= parent_.m_shape_element.m_ndofs;
  const wmesh_int_t ndofs_test 		= parent_.m_shape_test.m_ndofs;
  const wmesh_int_t ndofs_trial 	= parent_.m_shape_trial.m_ndofs;
  const wmesh_int_t topodim 		= parent_.m_topodim;
  
  this->m_dofs_element_storage  = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&this->m_dofs_element, topodim + 1, ndofs_element);

  this->m_dofs_velocity_storage  = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&this->m_dofs_velocity, topodim + 1, ndofs_velocity);
  
  this->m_q_velocity_storage 		= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>::alloc(&this->m_q_velocity, topodim + 1, q_n);

  this->m_q_element_diff_storage 	= WMESH_STORAGE_INTERLEAVE;
  for (wmesh_int_t i=0;i<topodim;++i)
    {
      wmesh_mat_t<T>::alloc(&this->m_q_element_diff[i], topodim + 1, q_n);
    }
  
  wmesh_mat_t<T>::alloc(&this->m_build_rhs,1, q_n);
  wmesh_mat_t<T>::alloc(&this->m_local_matrix,ndofs_test,ndofs_trial);
  wmesh_mat_t<T>::alloc(&this->m_local_rhs,ndofs_test,1);
  wmesh_mat_t<T>::alloc(&this->m_q_a,q_n,1);
  wmesh_mat_t<T>::alloc(&this->m_q_coo,topodim+1,q_n);
}



template<typename T>
wmesh_status_t wmesh_template_integral_flux_t<T>::eval(wmesh_template_integral_flux_t<T>::data_t&  		data_,
						       const T 							alpha_,
						       wmesh_mat_t<T>& 						local_matrix_) const
{
  const wmesh_int_t n1 		= static_cast<wmesh_int_t>(1);
  const T r1 			= static_cast<T>(1);

  const wmesh_int_t q_n 	= this->m_cubature->m_w.n;
  const wmesh_int_t topodim 	= this->m_topodim;

  //
  //  Evaluate the velocity (velocity is in dimension (topodim+1) )
  //
  wmesh_mat_gemm(static_cast<T>(1),
		 data_.m_dofs_velocity,
		 this->m_shape_eval_velocity->m_f,
		 static_cast<T>(0),
		 data_.m_q_velocity);


  //
  // Evaluate the gradient of the reference coordinates.
  //
  for (wmesh_int_t i=0;i<topodim;++i)
    {
      wmesh_mat_gemm(static_cast<T>(1),
		     data_.m_dofs_element,
		     this->m_shape_eval_element->m_diff[i],
		     static_cast<T>(0),
		     data_.m_q_element_diff[i]);
#if 0
      std::cout << "aaaaaaaaaaaaaaaaaaaaaaaaaaaa " << std::endl;
      std::cout << data_.m_dofs_element << std::endl;
      std::cout << "bbbbbbbbbbbbbbbbbbbbbbbbbbbb " << std::endl;
      std::cout << this->m_shape_eval_element->m_diff[i] << std::endl;
      std::cout << "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhh " << std::endl;
      std::cout << data_.m_q_element_diff[i] << std::endl;
#endif      
    }
  
  T tmpu[3]{};
  T normal[3]{};
  T dx[3]{};
  T dy[3]{};
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      if (data_.m_q_velocity_storage==WMESH_STORAGE_BLOCK)
	{
	  for (wmesh_int_t i=0;i<=topodim;++i)
	    {
	      tmpu[i] = *(data_.m_q_velocity.v + data_.m_q_velocity.ld*i+k);
	    }
	}
      else
	{
	  for (wmesh_int_t i=0;i<=topodim;++i)
	    {
	      tmpu[i] = *(data_.m_q_velocity.v + data_.m_q_velocity.ld*k+i);
	    }
	}
      
      //
      // (3xn1) * (n1 * qn) = 3xqn = dx/dr dy/dr dz/dr
      //
      for (wmesh_int_t i=0;i<=topodim;++i)
	{
	  dx[i] = *(data_.m_q_element_diff[0].v + data_.m_q_element_diff[0].ld*k+i);
	}
      
      if (topodim==2)
	{
	  for (wmesh_int_t i=0;i<=topodim;++i)
	    {
	      dy[i] = *(data_.m_q_element_diff[1].v + data_.m_q_element_diff[1].ld*k+i);
	    }
	  //
	  //  0  1  
	  // -1  0 
	  //
	  normal[0] = dx[1] * dy[2] - dx[2] * dy[1];
	  normal[1] = dx[2] * dy[0] - dx[0] * dy[2];
	  normal[2] = dx[0] * dy[1] - dx[1] * dy[0];
	  T udotn = tmpu[0] * normal[0] + tmpu[1] * normal[1] + tmpu[2] * normal[2];
	  data_.m_build_rhs.v[k] = (udotn < static_cast<T>(0) ) ? -udotn : static_cast<T>(0);
	}
      else
	{
	  //
	  //  0  1  dx[0]
	  // -1  0  dx[1]
	  //
	  
	  normal[0] = dx[1];
	  normal[1] = -dx[0];
#if 0
	  std::cout << "######################## NORMAL" << std::endl;
	  std::cout << " " << normal[0] << " " << normal[1] << std::endl;
	  std::cout << "######################## TMPU" << std::endl;
	  std::cout << " " << tmpu[0] << " " << tmpu[1] << std::endl;
#endif
	  T udotn = tmpu[0] * normal[0] + tmpu[1] * normal[1];
#if 0
	  std::cout << "######################## UDOTN" << std::endl;
	  std::cout << " " <<  udotn << std::endl;
#endif
	  data_.m_build_rhs.v[k] = (udotn < static_cast<T>(0) ) ? -udotn : static_cast<T>(0);
	}
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






template<typename T>
wmesh_status_t wmesh_template_integral_flux_t<T>::eval_residual(wmesh_template_integral_flux_t<T>::data_t&  		data_,
								const T 						alpha_,
								wmesh_mat_t<T>& 					local_rhs_,
								T * 							q_a_,
								wmesh_int_t						q_a_inc_) const
{
  const wmesh_int_t n1 		= static_cast<wmesh_int_t>(1);
  const T r1 			= static_cast<T>(1);

  const wmesh_int_t q_n 	= this->m_cubature->m_w.n;
  const wmesh_int_t topodim 	= this->m_topodim;
  
  //
  //  Evaluate the velocity (velocity is in dimension (topodim+1) )
  //
  wmesh_mat_gemm(static_cast<T>(1),
		 data_.m_dofs_velocity,
		 this->m_shape_eval_velocity->m_f,
		 static_cast<T>(0),
		 data_.m_q_velocity);


  //
  // Evaluate the gradient of the reference coordinates.
  //
  std::cout <<data_.m_dofs_element << std::endl;
  for (wmesh_int_t i=0;i<topodim;++i)
    {
      wmesh_mat_gemm(static_cast<T>(1),
		     data_.m_dofs_element,
		     this->m_shape_eval_element->m_diff[i],
		     static_cast<T>(0),
		     data_.m_q_element_diff[i]);
    }
  
  T tmpu[3]{};
  T normal[3]{};
  T dx[3]{};
  T dy[3]{};
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      if (data_.m_q_velocity_storage==WMESH_STORAGE_BLOCK)
	{
	  for (wmesh_int_t i=0;i<=topodim;++i)
	    {
	      tmpu[i] = *(data_.m_q_velocity.v + data_.m_q_velocity.ld*i+k);
	    }
	}
      else
	{
	  for (wmesh_int_t i=0;i<=topodim;++i)
	    {
	      tmpu[i] = *(data_.m_q_velocity.v + data_.m_q_velocity.ld*k+i);
	    }
	}
      
      //
      // (3xn1) * (n1 * qn) = 3xqn = dx/dr dy/dr dz/dr
      //
      for (wmesh_int_t i=0;i<=topodim;++i)
	{
	  dx[i] = *(data_.m_q_element_diff[0].v + data_.m_q_element_diff[0].ld*k+i);
	}
      
      if (topodim==2)
	{
	  for (wmesh_int_t i=0;i<=topodim;++i)
	    {
	      dy[i] = *(data_.m_q_element_diff[1].v + data_.m_q_element_diff[1].ld*k+i);
	    }
	  //
	  //  0  1  
	  // -1  0 
	  //
	  normal[0] = dx[1] * dy[2] - dx[2] * dy[1];
	  normal[1] = dx[2] * dy[0] - dx[0] * dy[2];
	  normal[2] = dx[0] * dy[1] - dx[1] * dy[0];
	  T udotn = tmpu[0] * normal[0] + tmpu[1] * normal[1] + tmpu[2] * normal[2];
	  data_.m_build_rhs.v[k] = (udotn < static_cast<T>(0) ) ? -udotn : static_cast<T>(0);
	  data_.m_build_rhs.v[k] *= q_a_[q_a_inc_*k];
	}
      else
	{
	  //
	  //  0  1  dx[0]
	  // -1  0  dx[1]
	  //
	  
	  normal[0] = dx[1];
	  normal[1] = -dx[0];
#if 1
	  std::cout << "######################## NORMAL" << std::endl;
	  std::cout << " " << normal[0] << " " << normal[1] << std::endl;
	  std::cout << "######################## TMPU" << std::endl;
	  std::cout << " " << tmpu[0] << " " << tmpu[1] << std::endl;
#endif
	  T udotn = tmpu[0] * normal[0] + tmpu[1] * normal[1];
#if 1
	  std::cout << "######################## UDOTN" << std::endl;
	  std::cout << " " <<  udotn << " " << q_a_[q_a_inc_*k] << std::endl;
#endif
	  data_.m_build_rhs.v[k] = (udotn < static_cast<T>(0) ) ? -udotn : static_cast<T>(0);
	  data_.m_build_rhs.v[k] *= q_a_[q_a_inc_*k];
	}
    }
  
  //
  // Build the matrix.
  //
  //  if (local_rhs_.ld == 1)
    {
      xgemv("N",
	    &this->m_build_residual.m,
	    &this->m_build_residual.n,
	    &r1,
	    this->m_build_residual.v,
	    &this->m_build_residual.ld,
	    data_.m_build_rhs.v,
	    &n1,
	    &alpha_,
	    local_rhs_.v,
	    &n1);
    }
#if 0
  else
    {
      std::cerr << "need a correction." << std::endl;
      exit(1);
    }
#endif
  return WMESH_STATUS_SUCCESS;
};
