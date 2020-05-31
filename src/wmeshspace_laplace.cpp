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

template<typename T>
struct wmesh_cubature_t
{
  wmesh_int_t 		m_element;
  wmesh_int_t 		m_family;
  wmesh_int_t 		m_degree;
  wmesh_int_t 		m_c_storage;
  wmesh_mat_t<T> 	m_c;
  wmesh_mat_t<T> 	m_w;
};

template<typename T>
wmesh_status_t wmesh_cubature_def(wmesh_cubature_t<T>*__restrict__ self_,
				  wmesh_int_t 		element_,
				  wmesh_int_t 		family_,
				  wmesh_int_t 		degree_)
{
  wmesh_int_t rw_n;
  T * __restrict__ rw;
  
  memset(self_,0,sizeof(wmesh_cubature_t<T>));
  self_->m_element 	= element_;
  self_->m_family 	= family_;
  self_->m_degree 	= degree_;
  wmesh_status_t status;
  
  wmesh_int_t topodim;
  status = bms_element2topodim(element_,
			       &topodim);

  wmesh_int_t num_nodes1d;
  status = bms_cubature_num_nodes(WMESH_ELEMENT_EDGE,
				  family_,
				  degree_,
				  &num_nodes1d);
  WMESH_STATUS_CHECK(status);

  wmesh_int_t num_nodes;    
  status = bms_cubature_num_nodes(element_,
				  family_,
				  degree_,
				  &num_nodes);
  WMESH_STATUS_CHECK(status);

  self_->m_c_storage = WMESH_STORAGE_INTERLEAVE;	      
  
  wmesh_mat_t<T>::define(&self_->m_c,
			 topodim,
			 num_nodes,
			 (T*)malloc(sizeof(T) * topodim * num_nodes),
			 topodim);
  
  wmesh_mat_t<T>::define(&self_->m_w,
			 1,
			 num_nodes,
			 (T*)malloc(sizeof(T) * 1 * num_nodes),
			 1);
  
  status = bms_cubature_buffer_size(element_,
				    family_,
				    num_nodes1d,			
				    &rw_n);
  WMESH_STATUS_CHECK(status);
  
  rw = (rw_n > 0) ? (T*)malloc(sizeof(T)*rw_n) : nullptr;
  if (rw_n > 0 && !rw)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
    }
  
  status = bms_cubature(element_,
			family_,
			num_nodes1d,
			
			self_->m_c_storage,
			self_->m_c.m,
			self_->m_c.n,
			self_->m_c.v,
			self_->m_c.ld,
			    
			self_->m_w.n,
			self_->m_w.v,
			self_->m_w.ld,
			
			rw_n,
			rw);
  
  WMESH_STATUS_CHECK(status);
  if (rw)
    {
      free(rw);
    }
  rw = nullptr;  
  return WMESH_STATUS_SUCCESS;
}

struct wmesh_shape_t
{
  wmesh_int_t		m_element;
  wmesh_int_t		m_family;
  wmesh_int_t 		m_degree;
};

wmesh_status_t wmesh_shape_def(wmesh_shape_t*__restrict__ self_,
			       wmesh_int_t 		element_,
			       wmesh_int_t 		family_,
			       wmesh_int_t 		degree_)
{
  WMESH_CHECK_POINTER(self_);
  memset(self_,0,sizeof(wmesh_shape_t));
  self_->m_element 	= element_;
  self_->m_family 	= family_;
  self_->m_degree 	= degree_;
  return WMESH_STATUS_SUCCESS;
}

template<typename T>
struct wmesh_shape_eval_t
{
  wmesh_shape_t 	m_shape;


  wmesh_int_t 		m_f_storage;
  wmesh_mat_t<T> 	m_f;
  wmesh_int_t 		m_diff_storage;
  wmesh_mat_t<T> 	m_diff[3];

  wmesh_int_t 		m_wf_storage;
  wmesh_mat_t<T> 	m_wf;
  wmesh_int_t 		m_wdiff_storage;
  wmesh_mat_t<T> 	m_wdiff[3];

  
};



template<typename T>
wmesh_status_t wmesh_shape_eval_def_init(wmesh_shape_eval_t<T>*__restrict__ 	self_,
					 wmesh_int_t 				nodes_storage_,
					 const wmesh_mat_t<T> * 		nodes_)
{
  WMESH_CHECK_POINTER(self_);
  WMESH_CHECK_POINTER(nodes_);

  
  wmesh_int_t iw_n,rw_n;
  wmesh_status_t status;

  const wmesh_int_t element = self_->m_shape.m_element;
  const wmesh_int_t family  = self_->m_shape.m_family;
  const wmesh_int_t degree  = self_->m_shape.m_degree;

  wmesh_int_t topodim;
  status = bms_element2topodim(element,
			       &topodim);  
  WMESH_STATUS_CHECK(status);

  status = bms_shape_buffer_size(element,
				 self_->m_shape.m_family,
				 self_->m_shape.m_degree,
				 &iw_n,
				 &rw_n);
  WMESH_STATUS_CHECK(status);
  
  wmesh_int_p iw = (iw_n > 0) ? (wmesh_int_p)malloc(sizeof(T)*iw_n) : nullptr;
  if (iw_n > 0 && !iw)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
    }
  
  T * rw = (rw_n > 0) ? (T*)malloc(sizeof(T)*rw_n) : nullptr;
  if (rw_n > 0 && !rw)
    {
      if (iw)
	{
	  free(iw);	  
	}
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
    }
  
  wmesh_int_t diff[3] = {0,0,0};  
  status = bms_template_shape(element,
			      family,
			      degree,
			    
			      diff,
			    
			      nodes_storage_,			      
			      nodes_->m,
			      nodes_->n,
			      nodes_->v,
			      nodes_->ld,
			      
			      WMESH_STORAGE_INTERLEAVE,			    
			      self_->m_f.m,
			      self_->m_f.n,
			      self_->m_f.v,
			      self_->m_f.ld,
			    
			      iw_n,
			      iw,
			      rw_n,
			      rw);
  if (WMESH_STATUS_SUCCESS != status)
    {
      if (rw)
	{
	  free(rw);
	}
      if (iw)
	{
	  free(iw);
	}
      WMESH_STATUS_CHECK(status);
    }
	
  for (wmesh_int_t j=0;j<topodim;++j)
    {
      diff[j] = 1;
      status = bms_template_shape(element,
				  family,
				  degree,
			     
				  diff,
			     
				  nodes_storage_,
				  nodes_->m,
				  nodes_->n,
				  nodes_->v,
				  nodes_->ld,
				  
				  WMESH_STORAGE_INTERLEAVE,			    
				  self_->m_diff[j].m,
				  self_->m_diff[j].n,
				  self_->m_diff[j].v,
				  self_->m_diff[j].ld,
			      
				  iw_n,
				  iw,
				  rw_n,
				  rw);
      if (WMESH_STATUS_SUCCESS != status)
	{
	  if (rw)
	    {
	      free(rw);
	    }
	  if (iw)
	    {
	      free(iw);
	    }
	  WMESH_STATUS_CHECK(status);
	}
      diff[j] = 0;
    }
	
  if (rw)
    {
      free(rw);
    }
  if (iw)
    {
      free(iw);
    }

  
  return WMESH_STATUS_SUCCESS;  
}

template<typename T>
wmesh_status_t wmesh_shape_eval_def(wmesh_shape_eval_t<T>*__restrict__ 	self_,
				    wmesh_int_t 			element_,				
				    wmesh_int_t 			shape_family_,
				    wmesh_int_t 			shape_degree_,				
				    wmesh_int_t 			nodes_storage_,
				    const wmesh_mat_t<T> * 		nodes_,
				    const wmesh_mat_t<T> * 		weights_ = nullptr)
{
  WMESH_CHECK_POINTER(self_);
  WMESH_CHECK_POINTER(nodes_);

  wmesh_status_t status;

  memset(self_,0,sizeof(wmesh_shape_t));
  wmesh_shape_def(&self_->m_shape, element_, shape_family_,shape_degree_);
  
  const wmesh_int_t
    element 		= element_,
    degree		= shape_degree_,
    num_nodes 		= (nodes_storage_ == WMESH_STORAGE_INTERLEAVE) ? nodes_->n : nodes_->m;
  
  wmesh_int_t
    topodim,
    num_dofs_per_element;
  
  status = bms_element2topodim(element,&topodim);
  WMESH_STATUS_CHECK(status);

  status = bms_ndofs(element,
		     degree,
		     &num_dofs_per_element);
  WMESH_STATUS_CHECK(status);
  
  //
  //
  //
  {
    T*__restrict__ f_ptr = (T*)malloc(sizeof(T)*num_dofs_per_element*num_nodes);
  
    if (!f_ptr)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }

    T*__restrict__ diff_ptr[3];
    for (wmesh_int_t j=0;j<topodim;++j)
      {
	diff_ptr[j] = (T*)malloc(sizeof(T)*num_dofs_per_element*num_nodes);
      }

    { wmesh_int_t j;
      for (j=0;j<topodim;++j)
	{
	  if (!diff_ptr[j])
	    {
	      free(f_ptr);
	      f_ptr = nullptr;
	      for (wmesh_int_t k=0;k<j;++k)
		{
		  free(diff_ptr[k]);
		  diff_ptr[k] = nullptr;
		}
	      break;
	    }      
	}
      if (j < topodim)
	{
	  for (wmesh_int_t k=0;k<topodim;++k)
	    {
	      if (diff_ptr[k])
		{
		  free(diff_ptr[k]);
		  diff_ptr[k] = nullptr;
		}
	    }
	  WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	}
    }





    
    wmesh_mat_t<T>::define(&self_->m_f,
			   num_dofs_per_element,
			   num_nodes,
			   f_ptr,
			   num_dofs_per_element);
  
    for (wmesh_int_t j=0;j<topodim;++j)
      {
	wmesh_mat_t<T>::define(&self_->m_diff[j],
			       num_dofs_per_element,
			       num_nodes,
			       diff_ptr[j],
			       num_dofs_per_element);
      }
  }




  

  //
  //
  //

  status = wmesh_shape_eval_def_init(self_,
				     nodes_storage_,
				     nodes_);
  WMESH_STATUS_CHECK(status);


  if (weights_)
    {

      //
      // Allocate.
      //
      {
	T*__restrict__ f_ptr = (T*)malloc(sizeof(T)*num_dofs_per_element*num_nodes);
  
	if (!f_ptr)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }

	T*__restrict__ diff_ptr[3];
	for (wmesh_int_t j=0;j<topodim;++j)
	  {
	    diff_ptr[j] = (T*)malloc(sizeof(T)*num_dofs_per_element*num_nodes);
	  }

	{ wmesh_int_t j;
	  for (j=0;j<topodim;++j)
	    {
	      if (!diff_ptr[j])
		{
		  free(f_ptr);
		  f_ptr = nullptr;
		  for (wmesh_int_t k=0;k<j;++k)
		    {
		      free(diff_ptr[k]);
		      diff_ptr[k] = nullptr;
		    }
		  break;
		}      
	    }
	  if (j < topodim)
	    {
	      for (wmesh_int_t k=0;k<topodim;++k)
		{
		  if (diff_ptr[k])
		    {
		      free(diff_ptr[k]);
		      diff_ptr[k] = nullptr;
		    }
		}
	      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	    }
	}

	wmesh_mat_t<T>::define(&self_->m_wf,
			       num_dofs_per_element,
			       num_nodes,
			       f_ptr,
			       num_dofs_per_element);
    
	for (wmesh_int_t j=0;j<topodim;++j)
	  {
	    wmesh_mat_t<T>::define(&self_->m_wdiff[j],
				   num_dofs_per_element,
				   num_nodes,
				   diff_ptr[j],
				   num_dofs_per_element);
	  }
      }
      
      for (wmesh_int_t j=0;j<weights_->n;++j)
	{
	  T scal = weights_->v[weights_->ld*j];
	  for (wmesh_int_t i=0;i<self_->m_f.m;++i)
	    {
	      self_->m_wf.v[self_->m_wf.ld * j + i] = self_->m_f.v[self_->m_f.ld * j + i]  * scal;	      
	    }
	}
      
      for (wmesh_int_t idim=0;idim<topodim;++idim)
	{
	  for (wmesh_int_t j=0;j<weights_->n;++j)
	    {
	      T scal = weights_->v[weights_->ld*j];
	      for (wmesh_int_t i=0;i<self_->m_diff[idim].m;++i)
		{
		  self_->m_wdiff[idim].v[self_->m_wdiff[idim].ld * j + i] = self_->m_diff[idim].v[self_->m_diff[idim].ld * j + i]  * scal; 
		}
	    }
	}      
    }

  return WMESH_STATUS_SUCCESS;
}


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
struct wmesh_fem_laplace_t
{
  wmesh_cubature_t<T>  	m_cubature[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_t<T>	m_shape_eval_element[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_t<T>	m_shape_eval_f[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_t<T>	m_shape_eval_element_weighted[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_t<T>	m_shape_eval_f_weighted[WMESH_ELEMENT_ALL];
};

template <typename T>
wmesh_status_t wmesh_fem_laplace_def(wmesh_fem_laplace_t<T> *	self_,
				     wmesh_int_t 		element_,
				     wmesh_int_t 		shape_element_family_,
				     wmesh_int_t 		shape_element_degree_,
				     wmesh_int_t 		shape_f_family_,
				     wmesh_int_t 		shape_f_degree_,
				     wmesh_int_t 		cubature_family_,
				     wmesh_int_t 		cubature_degree_)
{
  WMESH_CHECK_POINTER(self_);  
  wmesh_status_t status;
  status = wmesh_cubature_def(&self_->m_cubature[element_],
			      element_,
			      cubature_family_,
			      cubature_degree_);
  WMESH_STATUS_CHECK(status);

  status = wmesh_shape_eval_def(&self_->m_shape_eval_element[element_],
				element_,				
				shape_element_family_,
				shape_element_degree_,
				self_->m_cubature[element_].m_c_storage,
				&self_->m_cubature[element_].m_c,
				&self_->m_cubature[element_].m_w);
  WMESH_STATUS_CHECK(status);

  status = wmesh_shape_eval_def(&self_->m_shape_eval_f[element_],
				element_,				
				shape_element_family_,
				shape_element_degree_,
				self_->m_cubature[element_].m_c_storage,
				&self_->m_cubature[element_].m_c,
				&self_->m_cubature[element_].m_w);
  WMESH_STATUS_CHECK(status);
  
  return WMESH_STATUS_SUCCESS;
}


template <typename T>
wmesh_status_t wmesh_fem_laplace_local_system(const wmesh_fem_laplace_t<T> *	self_,
					      wmesh_int_t 			element_,
					      wmesh_int_t 			cooelm_m_,
					      wmesh_int_t 			cooelm_n_,
					      const T*__restrict__		cooelm_,
					      wmesh_int_t 			cooelm_ld_,
					      wmesh_int_t 			lmat_m_,
					      wmesh_int_t 			lmat_n_,
					      T*__restrict__			lmat_,
					      wmesh_int_t 			lmat_ld_,
					      T*__restrict__			lrhs_,
					      wmesh_int_t 			lrhs_inc_,
					      wmesh_int_t 			rw_n_,
					      T*__restrict__			rw_)
{
  wmesh_shape_eval_t<T>	* shape_eval_element = self_->m_shape_eval_element[element_];
  wmesh_shape_eval_t<T>	* shape_eval_f = self_->m_shape_eval_f[element_];

  wmesh_int_t 	q_n = self_->m_cubature[element_].m_w.n;
  const T* 	q_w = self_->m_cubature[element_].m_w;
  wmesh_int_t 	topodim  	= cooelm_m_;
  wmesh_int_t  	topodimXtopodim = topodim*topodim;;
  T * 		jacobians 	= rw_;
  T * 		dets 		= topodimXtopodim * q_n;

  //
  // cooelm (3 X nn)
  // mat    nn x q_n
  //
  // cooelm * diff[0] = [ dx/dr(R0) dy/dr(R0) dz/dr(R0) ... ]  (3 X q_n)
  //
  //
  T r1=1, r0=0,B[9];
  for (wmesh_int_t i=0;i<topodim;++i)
    {
      BLAS_dgemm("N",
		 "N",
		 &topodim,
		 &shape_eval_element->m_diff[i].n,
		 &shape_eval_element->m_diff[i].m,
		 &r1,
		 cooelm_,
		 &cooelm_ld_,
		 shape_eval_element->m_diff[i].v,
		 &shape_eval_element->m_diff[i].ld,
		 &r0,
		 jacobians + topodim*i,
		 &topodimXtopodim);
    }

  for (wmesh_int_t j=0;j<q_n;++j)
    {
      T*jacobian = jacobians + topodimXtopodim * j;
      T det;
      if (topodim==2)
	{
	  T a00 = jacobian[0];
	  T a10 = jacobian[1];
	  T a01 = jacobian[2];
	  T a11 = jacobian[3];
	  det = a00 * a11 - a01 * a10;
	  jacobian[0] = a11 / det;
	  jacobian[1] = -a10 / det;
	  jacobian[2] = -a01 / det;
	  jacobian[3] = a00 / det;
	  if (det < 0.0)
	    det = -det;
	}
      else
	{
	  T a00 = jacobian[0];
	  T a10 = jacobian[1];
	  T a20 = jacobian[2];
	  T a01 = jacobian[3];
	  T a11 = jacobian[4];
	  T a21 = jacobian[5];
	  T a02 = jacobian[6];
	  T a12 = jacobian[7];
	  T a22 = jacobian[8];
	  T det0 = a11 * a22 - a12 * a21;
	  T det1 = a01 * a22 - a02 * a21;
	  T det2 = a01 * a12 - a02 * a11;
	  T det  = a00 * det0 - a10 * det1 + a20 * det2;
	  if (det < 0.0)
	    det = -det;			
	}
      dets[j] = det;

      //
      // Now inverse the jacobian.
      //
      
      for (wmesh_int_t i=0;i<topodim*topodim;++i) B[i]=0;
      for (wmesh_int_t i=0;i<topodim;++i) B[i]=1.0;
      wmesh_int_t perm[3],info_lapack;
      LAPACK_dgesv((wmesh_int_p)&topodim,
		   (wmesh_int_p)&topodim,
		   jacobian,
		   (wmesh_int_p)&topodim,
		   perm,
		   B,
		   (wmesh_int_p)&topodim,
		   (wmesh_int_p)&info_lapack);
		    
      for (wmesh_int_t i=0;i<topodim*topodim;++i)
	{
	  jacobians[topodim*topodim * j + i] = B[i];
	}
    }


  T lnabla_phi_i[3];
  T lnabla_phi_j[3];
  T nabla_phi_i[3];
  T nabla_phi_j[3];
  for (wmesh_int_t j=0;j<lmat_n_;++j)
    {
      for (wmesh_int_t i=0;i<lmat_m_;++i)
	{
	  lmat_[lmat_ld_*j+i] = 0;	      
	}
    }
      
  static constexpr wmesh_int_t n1 = static_cast<wmesh_int_t>(1);
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      for (wmesh_int_t j=0;j<lmat_n_;++j)
	{
	  T*jacobian = jacobians + topodimXtopodim * j;
	  for (wmesh_int_t idim=0;idim<topodim;++idim)
	    {
	      lnabla_phi_j[idim] = shape_eval_f->m_diff[idim].v[shape_eval_f->m_diff[idim].ld * k + j];
	    }
	  
	  BLAS_dgemm("N",
		     "N",
		     &topodim,
		     &n1,
		     &topodim,
		     &r1,
		     jacobian,
		     &topodim,
		     lnabla_phi_j,
		     &topodim,
		     &r0,
		     nabla_phi_j,
		     &topodim);

	  for (wmesh_int_t i=0;i<lmat_m_;++i)
	    {
	      for (wmesh_int_t idim=0;idim<topodim;++idim)
		{
		  lnabla_phi_i[idim] = shape_eval_f->m_diff[idim].v[shape_eval_f->m_diff[idim].ld * k + i];
		}

	      BLAS_dgemm("N",
			 "N",
			 &topodim,
			 &n1,
			 &topodim,
			 &r1,
			 jacobian,
			 &topodim,
			 lnabla_phi_i,
			 &topodim,
			 &r0,
			 nabla_phi_i,
			 &topodim);

	      T dot = static_cast<T>(0);
	      for (wmesh_int_t idim=0;idim<topodim;++idim)
		{
		  dot += nabla_phi_i[idim]*nabla_phi_j[idim];
		}
	      lmat_[lmat_ld_*j+i] += dot * dets[k] * q_w[k];
	    }
	  
	}      
    }



  for (wmesh_int_t i=0;i<lmat_m_;++i)
    {
      lrhs_[lrhs_inc_*i] = 0;	      
    }


  for (wmesh_int_t k=0;k<q_n;++k)
    {
      for (wmesh_int_t i=0;i<lmat_m_;++i)
	{
	  lrhs_[lrhs_inc_*i] += dets[k] * q_w[k] * shape_eval_f->m_f[k*shape_eval_f->m_ld + i] * 1;
	}      
    }
  
  //
  // Calculate the jacobian.
  //  
  return WMESH_STATUS_SUCCESS;
}




extern "C"
{
  
  wmesh_status_t wmeshspace_laplace(const wmeshspace_t*__restrict__ 	self_,
				    wmesh_int_t				csr_size_,
				    const_wmesh_int_p			csr_ptr_,
				    const_wmesh_int_p			csr_ind_,
				    double * 				csr_val_,
				    double * 				rhs_)
  {

    wmesh_int_t shape_family = WMESH_SHAPE_FAMILY_LAGRANGE;
    wmesh_int_t shape_degree = self_->m_degree;    
    wmesh_int_t topodim = self_->m_mesh->m_topology_dimension;

    const wmesh_int_t shape_element_family = WMESH_SHAPE_FAMILY_LAGRANGE;
    const wmesh_int_t shape_element_degree = 1;    

    wmesh_int_t cubature_family = WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE;
    wmesh_int_t cubature_degree = 2 * shape_degree;
    
    //
    // Initialize data.
    //
    wmesh_fem_laplace_t<double> fem{};
    wmesh_status_t status;

    
    wmesh_int_t ntypes;
    wmesh_int_t elements[4];
    status =  bms_topodim2elements(topodim,
				   &ntypes,
				   elements);
    WMESH_STATUS_CHECK(status);
    
    for (wmesh_int_t itype=0;itype<self_->m_c2d.m_size;++itype)
      {
	if (self_->m_c2d.m_n[itype] > 0)
	  {
	    status = wmesh_fem_laplace_def(&fem,
					   elements[itype],
					   shape_element_family,
					   shape_element_degree,				   
					   shape_family,
					   shape_degree,
					   cubature_family,
					   cubature_degree);
	    WMESH_STATUS_CHECK(status);
	  }
      }


    
    
#if 0
    WMESH_CHECK_POINTER(self_);
    WMESH_CHECK_POINTER(csr_size_);
    WMESH_CHECK_POINTER(csr_ptr_);
    WMESH_CHECK_POINTER(csr_ind_);
    WMESH_CHECK_POINTER(csr_val_);
    WMESH_CHECK_POINTER(rhs_);

    wmesh_status_t status;

    wmesh_mat_t<double> cubature_coordinates[4];
    wmesh_mat_t<double> cubature_weights[4];

    wmesh_mat_t<double> shape_f[4];
    wmesh_mat_t<double> shape_f_diff[4][3];

    wmesh_mat_t<double> shape_element[4];
    wmesh_mat_t<double> shape_element_diff[4][3];


    double lrhs[32];
    double lmat[32*32];


    //
    // Define cubatures.
    //
    wmesh_int_t n1d;
    wmesh_int_t cubature_num_nodes;

    //
    // Now loop over each zone.
    //
    double r1 = 1.0;
    double r0 = 0.0;
    double cooelm[32];
    for (wmesh_int_t itype=0;itype<ntypes;++itype)
      {
	if (self_->m_c2d.m_n[itype]>0)
	  {
	    for (wmesh_int_t ielm=0;ielm<self_->m_c2d.m_n[itype];++ielm)
	      {
		//
		// Reset.
		//
		wmesh_int_t num_dofs_per_elem = self_->m_c2d.m_m[itype];		
		for (wmesh_int_t j=0;j<num_dofs_per_elem;++j)
		  {
		    for (wmesh_int_t i=0;i<num_dofs_per_elem;++i)
		      {
			lmat[j * num_dofs_per_elem + i] = (double)0.0;
		      }
		  }
		for (wmesh_int_t i=0;i<num_dofs_per_elem;++i)
		  {
		    lrhs[i] = (double)0.0;
		  }


		//
		// Get coordinates of the elements.
		//
		wmesh_int_t cooelm_m = topodim;
		wmesh_int_t cooelm_ld = cooelm_ld;
		wmesh_int_t num_nodes_per_elem = self_->m_mesh->m_c2n.m_m[itype];
		for (wmesh_int_t i=0;i<num_nodes_per_elem;++i)
		  {
		    for (wmesh_int_t k=0;k<topodim;++k)
		      {			
			cooelm[i*cooelm_ld + k] =
			  self_->m_mesh->m_coo[(self_->m_mesh->m_c2n.m_data[self_->m_mesh->m_c2n.m_ptr[itype] +self_->m_mesh->m_c2n.m_ld[itype]    *ielm+i]-1)*self_->m_mesh->m_coo_ld + k];
		      }
		  }
		wmesh_mat_t<double> * mat_element_diff 	= &shape_element_diff[itype][0];		
		wmesh_mat_t<double> * mat_f_diff 	= &shape_f_diff[itype][0];		
		wmesh_mat_t<double> * eval_element_diff	= &shape_eval_element_diff[itype][0];		
		wmesh_mat_t<double> * eval_f_diff 	= &shape_eval_f_diff[itype][0];		
		wmesh_mat_t<double> * q_coo 		= &cubature_coordinates[itype];
		wmesh_int_t q_m = q_coo->m;
		wmesh_int_t q_n = q_coo->n;
		wmesh_int_t q_ld = q_coo->ld;
		double * q_v = q_coo->v;
		double det;
		double jacobian[9],g_gradient[9];
		double jacobians[1024*9];
		double dets[1024];


		double  diff_wj[3];
		double  diff_wi[3];
		double  gradient_wj[3];
		double  gradient_wi[3];
		wmesh_int_t diff_wi_inc[3],n1=1;
		wmesh_int_t diff_wj_inc[3];

		    wmesh_int_t perm[3], info_lapack;
		    double B[9];

		for (wmesh_int_t iq=0;iq<q_n;++iq)
		  {
		    for(wmesh_int_t l=0;l<topodim;++l)
		      {
			//
			// Compute the local gradient.
			//
			BLAS_dgemm("N",
				   "N",
				   &cooelm_m,
				   &q_n,
				   &num_nodes_per_elem,
				   &r1,
				   cooelm,
				   &cooelm_ld,
				   mat_element_diff[l].v,
				   &mat_element_diff[l].ld,
				   &r0,
				   &eval_element_diff[l].v,
				   &eval_element_diff[l].ld);
		      }
		  }

		//
		// Now dphi_j
		// 
		for (wmesh_int_t iq=0;iq<q_n;++iq)
		  {
		    //
		    // Compute jacobian
		    //
		    for(wmesh_int_t j=0;j<topodim;++j)
		      {
			for(wmesh_int_t i=0;i<topodim;++i)
			  {
			    jacobian[j*topodim+i] = element_diff[j].v[iq * element_diff[j].ld + i];
			  }
		      }

		    //
		    // Inverse the jacobian and compute the absolute determinant.
		    //
		    
		    if (topodim==2)
		      {
			double a00 = jacobian[0];
			double a10 = jacobian[1];
			double a01 = jacobian[2];
			double a11 = jacobian[3];
			det = a00 * a11 - a01 * a10;
			jacobian[0] = a11 / det;
			jacobian[1] = -a10 / det;
			jacobian[2] = -a01 / det;
			jacobian[3] = a00 / det;
			if (det < 0.0)
			  det = -det;
		      }
		    else
		      {
			double a00 = jacobian[0];
			double a10 = jacobian[1];
			double a20 = jacobian[2];
			double a01 = jacobian[3];
			double a11 = jacobian[4];
			double a21 = jacobian[5];
			double a02 = jacobian[6];
			double a12 = jacobian[7];
			double a22 = jacobian[8];
			double det0 = a11 * a22 - a12 * a21;
			double det1 = a01 * a22 - a02 * a21;
			double det2 = a01 * a12 - a02 * a11;
			double det  = a00 * det0 - a10 * det1 + a20 * det2;
			if (det < 0.0)
			  det = -det;			
		      }
		    
		    for (wmesh_int_t i=0;i<topodim*topodim;++i) B[i]=0;
		    for (wmesh_int_t i=0;i<topodim;++i) B[i]=1.0;
		    
		    LAPACK_dgesv((wmesh_int_p)&topodim,
				 (wmesh_int_p)&topodim,
				 jacobian,
				 (wmesh_int_p)&topodim,
				 perm,
				 B,
				 (wmesh_int_p)&topodim,
				 (wmesh_int_p)&info_lapack);
		    
		    for (wmesh_int_t i=0;i<topodim*topodim;++i)
		      {
			jacobians[topodim*topodim * iq + i] = B[i];
		      }		    
		    dets[iq] = det;
		  }


		for (wmesh_int_t iq=0;iq<q_n;++iq)
		  {
		    wmesh_int_t topodimXtopodim = topodim*topodim;
		    		  
		    for(wmesh_int_t j=0;j<num_dofs_per_elem;++j)
		      {
			
			for (wmesh_int_t l=0;l<topodim;++l)
			  {
			    diff_wj[l] = &mat_f_diff[l].v[mat_f_diff[l].ld * iq + j];
			    diff_wj_inc[l] = mat_f_diff[l].ld;
			  }

			BLAS_dgemm("N",
				   "N",
				   &topodim,
				   &n1,
				   &topodim,
				   &r1,
				   &jacobians[topodim*topodim*iq],
				   &topodimXtopodim,
				   diff_wj,
				   &topodim,
				   &r0,
				   gradient_wj,
				   &topodim);

			for(wmesh_int_t i=0;i<num_dofs_per_elem;++i)
			  {			    
			    
			    for (wmesh_int_t l=0;l<topodim;++l)
			      {
				diff_wi[l] = &mat_f_diff[l].v[mat_f_diff[l].ld * iq + i];
				diff_wi_inc[l] = mat_f_diff[l].ld;
			      }

			    BLAS_dgemm("N",
				       "N",
				       &topodim,
				       &n1,
				       &topodim,
				       &r1,
				       &jacobians[topodim*topodim*iq],
				       &topodimXtopodim,
				       diff_wi,
				       &topodim,
				       &r0,
				       gradient_wi,
				       &topodim);
			    
			    if (topodim==2)
			      {
				lmat[j*num_dofs_per_elem + i] += dets[iq] * qw[iq] * (  -(gradient_wi[0]*gradient_wj[0]+gradient_wi[1]*gradient_wi[1])   );
			      }
			    else if (topodim==3)
			      {
				lmat[j*num_dofs_per_elem + i] += dets[iq] * qw[iq] * (  -(gradient_wi[0]*gradient_wj[0]+gradient_wi[1]*gradient_wi[1]+gradient_wi[2]*gradient_wi[2])   );
			      }
			    
			  }			    
		      }
		  }

		
		//
		// Compute the right hand side.
		//
		for (wmesh_int_t iq=0;iq<q_n;++iq)
		  {
		    f[iq] = 1 * q_w[iq] * dets[iq];
		  }
		
		BLAS_dgemm("N",
			   "N",
			   &num_dofs_per_elem,
			   &n1,
			   &q_n,
			   &r1,
			   mat_f[itype].v,
			   &mat_f[itype].ld,
			   f,
			   &q_n,
			   &r1,
			   lrhs,
			   &num_dofs_per_elem);					    
		

		//
		// Assembly.
		//
		for (wmesh_int_t i=0;i<num_dofs_per_elem;++i)
		  {
		    wmesh_int_t idof = self_->m_c2d.m_data[self_->m_c2d.m_ptr[itype] +self_->m_c2d.m_ld[itype] * ielm + i]-1;
		    for (wmesh_int_t j=0;j<num_dofs_per_elem;++j)
		      {
			wmesh_int_t jdof = self_->m_c2d.m_data[self_->m_c2d.m_ptr[itype] + self_->m_c2d.m_ld[itype] * ielm + j]-1;
			//
			// Search in the matrix.
			//
			for (wmesh_int_t at = csr_ptr_[idof];at < csr_ptr_[idof+1];++at)
			  {
			    //  std::cout << csr_ind_[at] << std::endl;
			    if (csr_ind_[at]-1 == jdof)
			      {
				csr_val_[at] += lmat[j*num_dofs_per_elem + i];
			      }
			  }
		      }
		  }
		
		for (wmesh_int_t i=0;i<num_dofs_per_elem;++i)
		  {
		    rhs_[self_->m_c2d.m_data[self_->m_c2d.m_ptr[itype] + self_->m_c2d.m_ld[itype] *ielm+i]-1] += lrhs[i];
		  }
		
		//
		//
		//
	      }	      
	  }
      }
#endif
    return WMESH_STATUS_SUCCESS;
  }
  
  
}

#if 0
    //
    // Prepare data.
    //
    wmesh_int_t cubature_n[4]{};
    double * 	cubature_p[4]{};
    double * 	cubature_w[4]{};
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++i)
      {
	if (self_->m_c2d.m_n[l] > 0)
	  {
	    
	  }
      }
    
    //
    // Element shape.
    //
    double * 	eval_shape_element[4]{};
    double * 	eval_shape_element_dr[4]{};
    double * 	eval_shape_element_ds[4]{};
    double * 	eval_shape_element_dt[4]{};
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++i)
      {
	if (self_->m_c2d.m_n[l] > 0)
	  {
	    
	  }
      }

    //
    // Variable shape.
    //
    double * 	eval_shape_f[4]{};
    double * 	eval_shape_f_dr[4]{};
    double * 	eval_shape_f_ds[4]{};
    double * 	eval_shape_f_dt[4]{};
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++i)
      {
	if (self_->m_c2d.m_n[l] > 0)
	  {
	    
	  }
      }
    
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++i)
      {

	if (self_->m_c2d.m_n[l] > 0)
	  {
	    
	    //
	    // Need integrations points.
	    //

	    //
	    // Need to evaluate shape of the variable F on integrations points.
	    //

	    //
	    // Need to evaluate shape of the element on integrations points.
	    //

	    //
	    // Now loop over the element of this type.
	    //
	    for (wmesh_int_t j=0;j<self_->m_c2d.m_n[l];++j)
	      {

		//
		// Get the coordinates of the elements.
		//
		cooelm[0] = 0;
		
		//
		// Evaluate the derivatives of the transformation on  the cubature.
		//

		
		//
		// Element jacobian and determinant on the cubature.
		//


		//
		// Element jacobian and determinant on the cubature.
		//

		
		//
		// Load the dof indices.
		//
		for (wmesh_int_t i=0;i<self_->m_c2d.m_m[l];++i)
		  {
		    c2d[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[l] + j * self_->m_c2d.m_ld[l] + i] - 1;
		  }

		
		//
		// Assembly.
		//

	      }
	  }
      }

    
    wmesh_int_t num_dofs 		= self_->m_ndofs;
    wmesh_int_t iw_n, num_table_coeffs;

    status = bms_sparse_buffer_size	(self_->m_ndofs,
					 self_->m_c2d.m_size,
					 self_->m_c2d.m_m,
					 self_->m_c2d.m_n,
					 &iw_n,
					 &num_table_coeffs);
    WMESH_STATUS_CHECK(status);
    
    
    wmesh_int_p iw  = (iw_n > 0) ? (wmesh_int_p)malloc(iw_n * sizeof(wmesh_int_t)) : nullptr;
    if (iw_n > 0 && !iw)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    

    wmesh_int_p csr_ptr = (wmesh_int_p)malloc( (num_dofs + 1) * sizeof(wmesh_int_t));
    if (!csr_ptr)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }

    std::cout << "ptr" << std::endl;
    status = bms_sparse_ptr(num_dofs,
			    WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2d),			 
			    csr_ptr,
			    iw_n,
			    iw);

    WMESH_STATUS_CHECK(status);
    std::cout << "ptr done" << std::endl;

    wmesh_int_p csr_ind = (wmesh_int_p)malloc( csr_ptr[num_dofs] * sizeof(wmesh_int_t));
    if (!csr_ind)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }
    std::cout << "sp .." << std::endl;
    
    status = bms_sparse	(num_dofs,
			 WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2d),			 
			 csr_ptr,
			 csr_ind,
			 iw_n,
			 iw);
    WMESH_STATUS_CHECK(status);
    std::cout << "sp done." << std::endl;

    //
    // Free iw.
    //
    free(iw);
    iw = nullptr;

    csr_size_[0] 	= num_dofs;
    csr_ptr_[0] 	= csr_ptr;
    csr_ind_[0] 	= csr_ind;
    
    return WMESH_STATUS_SUCCESS;
  }
  
};
#endif
