
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
#if 0
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
struct wmesh_shape_eval_boundary_t
{
  //
  // [edge|triangle|quad][0:pos,1:neg][rot:1,2]
  //
  wmesh_shape_t 	m_shape;

  wmesh_cubature_t<T>   m_cubature_facets[2][2][4];
  wmesh_shape_eval_t<T> m_shape_eval_facets[2][2][4];
};


template<typename T>
wmesh_status_t wmesh_shape_eval_boundary_def(wmesh_shape_eval_boundary_t<T>*__restrict__ self_,
					     wmesh_int_t 			element_,				
					     wmesh_int_t 			shape_family_,
					     wmesh_int_t 			shape_degree_,				
					     wmesh_int_t 			cubature_family_,				
					     wmesh_int_t 			cubature_degree_)
{
  WMESH_CHECK_POINTER(self_);
  
  wmesh_status_t status;

  //
  // Reset.
  //
  memset(self_, 0, sizeof(wmesh_shape_eval_boundary_t<T>));

  //
  // Define shape.
  //
  status = wmesh_shape_def(&self_->m_shape,
			   element_,
			   shape_family_,
			   shape_degree_);

  WMESH_STATUS_CHECK(status);


  wmesh_int_t
    topodim,
    topodim_boundary;
  
  status = bms_element2topodim(element_,&topodim);
  WMESH_STATUS_CHECK(status);
  
  topodim_boundary = topodim - 1;

  //
  // Now we need to know the number of facets.
  //
  wmesh_int_t num_facets;
  wmesh_int_t facets[6];
  wmesh_int_t facets_num_nodes[6];
  status = bms_element_facets(element_,
				  &num_facets,
				  facets);
  WMESH_STATUS_CHECK(status);


  //
  // Defines cubatures.
  //
  status = bms_elements_num_entities(num_facets,
				     facets,
				     WMESH_ELEMENT_NODE,
				     facets_num_nodes);
  WMESH_STATUS_CHECK(status);

  //
  // Prepare quadratures.
  // [edge|triangle|quad][0:pos,1:neg][rot:1,2]
  //
  //
  // Prepare quadratures.
  // [edge|triangle|quad][0:pos,1:neg][rot:1,2]
  //
  for (wmesh_int_t ifacet = 0;ifacet < num_facets;++ifacet)
    {      
      const wmesh_int_t element_facet 	= facets[ifacet];
      const wmesh_int_t facet_num_nodes	= facets_num_nodes[ifacet];
      
      {
	const wmesh_int_t signed_rotation = 1;
	wmesh_int_t i,j,k;
	i = element_facet - ( (topodim_boundary==2) ? 2 : 1);
	j = (signed_rotation > 0) ? 0 : 1;
	k = signed_rotation;
	
	status =  wmesh_cubature_def(&self_->m_cubature_facets[i][j][k],
				     element_facet,
				     cubature_family_,
				     cubature_degree_);
	WMESH_STATUS_CHECK(status);
      }
    }
  
  for (wmesh_int_t ifacet = 0;ifacet < num_facets;++ifacet)
    {
      
      const wmesh_int_t element_facet 	= facets[ifacet];
      const wmesh_int_t facet_num_nodes	= facets_num_nodes[ifacet];
      
      for (wmesh_int_t 	signed_rotation = -facet_num_nodes;signed_rotation<=facet_num_nodes;++signed_rotation)
	{
	  if (signed_rotation != 0 && signed_rotation!=1)
	    {
	      wmesh_int_t i,j,k;
	      i = element_facet - ( (topodim_boundary==2) ? 2 : 1);
	      j = (signed_rotation > 0) ? 0 : 1;
	      k = signed_rotation;
	      
	      // for (wmesh_int_t signed_rotation=-
	      status =  wmesh_cubature_def(&self_->m_cubature_facets[i][j][k],
					   element_facet,
					   cubature_family_,
					   cubature_degree_);
	      WMESH_STATUS_CHECK(status);

	      status = bms_mirrored_local_coordinates(element_facet,
						      signed_rotation,
						      
						      self_->m_cubature_facets[i][0][1].m_c_storage,
						      self_->m_cubature_facets[i][0][1].m_c.m,
						      self_->m_cubature_facets[i][0][1].m_c.n,
						      self_->m_cubature_facets[i][0][1].m_c.v,
						      self_->m_cubature_facets[i][0][1].m_c.ld,
						      
						      self_->m_cubature_facets[i][j][k].m_c_storage,
						      self_->m_cubature_facets[i][j][k].m_c.m,
						      self_->m_cubature_facets[i][j][k].m_c.n,
						      self_->m_cubature_facets[i][j][k].m_c.v,
						      self_->m_cubature_facets[i][j][k].m_c.ld);
	      WMESH_STATUS_CHECK(status);

	    }
	}
    }


  //
  // Now eval shape over quadratures.
  // [edge|triangle|quad][0:pos,1:neg][rot:1,2]
  //
  for (wmesh_int_t ifacet = 0;ifacet < num_facets;++ifacet)
    {
      const wmesh_int_t element_facet 	= facets[ifacet];
      const wmesh_int_t facet_num_nodes	= facets_num_nodes[ifacet];
      for (wmesh_int_t 	signed_rotation = -facet_num_nodes;signed_rotation<=facet_num_nodes;++signed_rotation)
	{
	  if (signed_rotation != 0)
	    {
	      wmesh_int_t i,j,k;
	      i = element_facet - ( (topodim_boundary==2) ? 2 : 1);
	      j = (signed_rotation > 0) ? 0 : 1;
	      k = signed_rotation;
	      status = wmesh_shape_eval_def(&self_->m_shape_eval_facets[i][j][k],
					    element_facet,				
					    shape_family_,
					    shape_degree_,
					    self_->m_cubature_facets[i][j][k].m_c_storage,
					    &self_->m_cubature_facets[i][j][k].m_c,
					    &self_->m_cubature_facets[i][j][k].m_w);
	      
	      WMESH_STATUS_CHECK(status);
	    }
	}
    }

  return WMESH_STATUS_SUCCESS;
}






template<typename T>
struct wmeshspacedg_fem_advection_t
{
  wmesh_cubature_t<T>  			m_cubature[WMESH_ELEMENT_ALL];
  
  //
  // 
  //
  wmesh_shape_eval_t<T>			m_shape_eval_element[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_t<T>			m_shape_eval_f[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_t<T>			m_shape_eval_u[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_t<T>			m_shape_eval_test[WMESH_ELEMENT_ALL];

  //
  // 
  //
  wmesh_shape_eval_boundary_t<T>	m_shape_eval_boundary_element[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_boundary_t<T>	m_shape_eval_boundary_f[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_boundary_t<T> 	m_shape_eval_boundary_test[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_boundary_t<T> 	m_shape_eval_boundary_u[WMESH_ELEMENT_ALL];


#if 0
  //
  // Buffer.
  //
  wmesh_shape_eval_t<T>			m_shape_eval_fglobal[WMESH_ELEMENT_ALL];

  
  status = wmesh_shape_eval_def(&self_->m_shape_eval_fglobal[element_],
				element_,				
				shape_f_family_,
				shape_f_degree_,
				self_->m_cubature[element_].m_c_storage,
				&self_->m_cubature[element_].m_c,
				&self_->m_cubature[element_].m_w);
  WMESH_STATUS_CHECK(status);
#endif

};



template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_add(wmeshspacedg_fem_advection_t<T> *	self_,
					      wmesh_int_t 			element_,
					      wmesh_int_t 			cubature_family_,
					      wmesh_int_t 			cubature_degree_,
					      wmesh_int_t 			shape_element_family_,
					      wmesh_int_t 			shape_element_degree_,
					      wmesh_int_t 			shape_u_family_,
					      wmesh_int_t 			shape_u_degree_,
					      wmesh_int_t 			shape_f_family_,
					      wmesh_int_t 			shape_f_degree_,
					      wmesh_int_t 			shape_test_family_,
					      wmesh_int_t 			shape_test_degree_)
{
  wmesh_status_t status;
  status = wmesh_shape_eval_boundary_def(&self_->m_shape_eval_boundary_element[element_],
					 element_,
					 shape_element_family_,
					 shape_element_degree_,
					 cubature_family_,				
					 cubature_degree_);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_eval_boundary_def(&self_->m_shape_eval_boundary_u[element_],
					 element_,
					 shape_u_family_,
					 shape_u_degree_,
					 cubature_family_,				
					 cubature_degree_);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_eval_boundary_def(&self_->m_shape_eval_boundary_f[element_],
					 element_,
					 shape_f_family_,
					 shape_f_degree_,
					 cubature_family_,				
					 cubature_degree_);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_eval_boundary_def(&self_->m_shape_eval_boundary_test[element_],
					 element_,
					 shape_test_family_,
					 shape_test_degree_,
					 cubature_family_,				
					 cubature_degree_);
  WMESH_STATUS_CHECK(status);      




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
				shape_f_family_,
				shape_f_degree_,
				self_->m_cubature[element_].m_c_storage,
				&self_->m_cubature[element_].m_c,
				&self_->m_cubature[element_].m_w);
  WMESH_STATUS_CHECK(status);

  status = wmesh_shape_eval_def(&self_->m_shape_eval_u[element_],
				element_,				
				shape_u_family_,
				shape_u_degree_,
				self_->m_cubature[element_].m_c_storage,
				&self_->m_cubature[element_].m_c,
				&self_->m_cubature[element_].m_w);
  WMESH_STATUS_CHECK(status);
  
  status = wmesh_shape_eval_def(&self_->m_shape_eval_test[element_],
				element_,				
				shape_test_family_,
				shape_test_degree_,
				self_->m_cubature[element_].m_c_storage,
				&self_->m_cubature[element_].m_c,
				&self_->m_cubature[element_].m_w);
  WMESH_STATUS_CHECK(status);

  return WMESH_STATUS_SUCCESS;
}


template<typename T>
static inline T determinant1x1(const T * __restrict__ 	jacobian_,
			       wmesh_int_t 		jacobian_ld_)
{
  return jacobian_[jacobian_ld_*0+0];
}

template<typename T>
static inline T determinant2x2(const T * __restrict__ 	jacobian_,
			       wmesh_int_t 		jacobian_ld_)
{
  return jacobian_[jacobian_ld_*0+0] * jacobian_[jacobian_ld_*1+1] - jacobian_[jacobian_ld_*0+1] *jacobian_[jacobian_ld_*1+0];
}

template<typename T>
static inline T determinant3x3(const T * __restrict__ 	jacobian_,
			       wmesh_int_t 		jacobian_ld_)
{
  const T a00 = jacobian_[jacobian_ld_ * 0 + 0];
  const T a10 = jacobian_[jacobian_ld_ * 0 + 1];
  const T a20 = jacobian_[jacobian_ld_ * 0 + 2];
  const T a01 = jacobian_[jacobian_ld_ * 1 + 0];
  const T a11 = jacobian_[jacobian_ld_ * 1 + 1];
  const T a21 = jacobian_[jacobian_ld_ * 1 + 2];
  const T a02 = jacobian_[jacobian_ld_ * 2 + 0];
  const T a12 = jacobian_[jacobian_ld_ * 2 + 1];
  const T a22 = jacobian_[jacobian_ld_ * 2 + 2];
  const T det0 = a11 * a22 - a12 * a21;
  const T det1 = a01 * a22 - a02 * a21;
  const T det2 = a01 * a12 - a02 * a11;
  return a00 * det0 - a10 * det1 + a20 * det2;
}







template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_local_system_buffer_size(const wmeshspacedg_fem_advection_t<T> *__restrict__	self_,
							  wmesh_int_t 					element_,							
							  wmesh_int_p 					rw_n_)
{
  const wmesh_cubature_t<T>*__restrict__ 	cubature		= &self_->m_cubature[element_];
  const wmesh_int_t 				q_n 			= cubature->m_w.n;
  wmesh_status_t status;
  wmesh_int_t topodim;
  status = bms_element2topodim	(element_,
				 &topodim);
  WMESH_STATUS_CHECK(status);
  rw_n_[0] = topodim * topodim * q_n + q_n;
  return WMESH_STATUS_SUCCESS;
}

template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_local_system(const wmeshspacedg_fem_advection_t<T> *__restrict__	self_,
						       wmesh_int_t 					element_,
						       wmesh_int_t 					cooelm_storage_,
						       wmesh_int_t 					cooelm_m_,
						       wmesh_int_t 					cooelm_n_,
						       const T*__restrict__				cooelm_,
						       wmesh_int_t 					cooelm_ld_,
						       wmesh_int_t 					uelm_storage_,
						       wmesh_int_t 					uelm_m_,
						       wmesh_int_t 					uelm_n_,
						       const T*__restrict__				uelm_,
						       wmesh_int_t 					uelm_ld_,
						       wmesh_int_t 					lmat_m_,
						       wmesh_int_t 					lmat_n_,
						       T*__restrict__					lmat_,
						       wmesh_int_t 					lmat_ld_,
						       T*__restrict__					lrhs_,
						       wmesh_int_t 					lrhs_inc_,
						       wmesh_int_t 					rw_n_,
						       T*__restrict__					rw_)
{
  static constexpr T r0 = static_cast<T>(0);
  static constexpr T r1 = static_cast<T>(1);

  WMESH_CHECK_POINTER(self_);
  
  wmesh_status_t status;
  wmesh_int_t required_rw_n;
  status = wmeshspacedg_fem_advection_local_system_buffer_size(self_,
							       element_,
							       &required_rw_n);
  WMESH_STATUS_CHECK(status);
  
  if (rw_n_ < required_rw_n)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);       
    }
  
  const wmesh_shape_eval_t<T> * __restrict__ 	shape_eval_element 	= &self_->m_shape_eval_element[element_];
  const wmesh_shape_eval_t<T> * __restrict__	shape_eval_f 		= &self_->m_shape_eval_f[element_];
  const wmesh_shape_eval_t<T> * __restrict__	shape_eval_u 		= &self_->m_shape_eval_u[element_];
  const wmesh_shape_eval_t<T> * __restrict__	shape_eval_test 	= &self_->m_shape_eval_test[element_];  
  const wmesh_cubature_t<T>*__restrict__ 	cubature		= &self_->m_cubature[element_];

  const wmesh_int_t 				u_n 			= (uelm_storage_ == WMESH_STORAGE_INTERLEAVE) ? uelm_n_ : uelm_m_;
  const wmesh_int_t 				q_n 			= cubature->m_w.n;
  const T* __restrict__				q_w 			= cubature->m_w.v;  
  const wmesh_int_t 	topodim  					= cooelm_m_;
  const wmesh_int_t  	topodimXtopodim 				= topodim*topodim;;
  
  //
  // Compute Jacobians and determinants.
  //
  T *__restrict__ jacobians 	= rw_;
  T *__restrict__ dets 		= rw_ + topodimXtopodim * q_n;
  rw_ += required_rw_n;
  rw_n_ -= required_rw_n;


  


  
  for (wmesh_int_t idim=0;idim<topodim;++idim)
    {
      BLAS_dgemm("N",
		 "N",
		 &topodim,
		 &shape_eval_element->m_diff[idim].n,
		 &shape_eval_element->m_diff[idim].m,
		 &r1,
		 cooelm_,
		 &cooelm_ld_,
		 shape_eval_element->m_diff[idim].v,
		 &shape_eval_element->m_diff[idim].ld,
		 &r0,
		 jacobians + topodim*idim,
		 &topodimXtopodim);
    }
  
  T B[9];
  for (wmesh_int_t j=0;j<q_n;++j)
    {
      T*jacobian = jacobians + topodimXtopodim * j;
      T det;
      if (topodim==2)
	{
	  det = determinant2x2(jacobian,
			       topodim);
	}
      else if (topodim==3)
	{
	  det = determinant3x3(jacobian,
			       topodim);
	}
      else
	{
	  det = determinant1x1(jacobian,
			       topodim);
	}
          
      if (det < 0.0)
	det = -det;
      dets[j] = det;
      
      //
      // Reset identity
      //
      for (wmesh_int_t i=0;i<topodim*topodim;++i)
	{
	  B[i] = r0;
	}
      
      for (wmesh_int_t i=0;i<topodim;++i)
	{
	  B[i*topodim+i] = r1;
	}

      //
      // Inverse the jacobian.
      //
      
      {	wmesh_int_t info_lapack,perm[3];
	LAPACK_dgesv((wmesh_int_p)&topodim,
		     (wmesh_int_p)&topodim,
		     jacobian,
		     (wmesh_int_p)&topodim,
		     perm,
		     B,
		     (wmesh_int_p)&topodim,
		     (wmesh_int_p)&info_lapack);
	if (info_lapack)
	  {
	    fprintf(stderr,"error lapack\n");
	    exit(1);
	  } }

      //
      // Store the inverse.
      //
      for (wmesh_int_t i=0;i<topodimXtopodim;++i)
	{
	  jacobians[topodim*topodim * j + i] = B[i];
	}
    }


  //
  // Reset the local matrix and the local rhs.
  //
  for (wmesh_int_t j=0;j<lmat_n_;++j)
    {
      for (wmesh_int_t i=0;i<lmat_m_;++i)
	{
	  lmat_[lmat_ld_*j+i] = r0;	      
	}
    }

  for (wmesh_int_t i=0;i<lmat_m_;++i)
    {
      lrhs_[lrhs_inc_*i] = r0;	      
    }
#if 0



  //
  // Eval u.
  //

  
  {
    //
    // Over the cell: evaluate velocity field.
    //
    cell_eval_u = cell_build_eval_u * cell_u;

    //
    // Over the cell: evaluate global gradient phi_j.
    //
    
    //
    // Over the cell: evaluate global gradient phi_j.
    //
    
    //
    // Over the cell: evaluate global gradient phi_j.
    //
  }

  for (wmesh_int_t cell_idx=0;cell_idx<num_cells;++cell_idx)
    {
      //
      // Get velocity.
      //

      
      //
      // Now we need to integrate over the boundary.
      //
      for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
	{
	  wmesh_int_t facet = facets[ifacet];

	  //
	  // Get the faces.
	  //
	  facet2n[];


	  //
	  // Where is the minimum?
	  //
	  wmesh_int_t c=facet2n[0];
	  wmesh_int_t signed_rotation=1;
	  for (wmesh_int_t i=1;i<4;++i)
	    {
	      if (facet2n[i] < c)
		{
		  c = facet2n[i];
		  signed_rotation = i+1;
		}
	    }
	  if (facet2n[ (signed_rotation-1) ])
	  
	  
	  //
	  // Over the cell: evaluate velocity field.
	  //
	  cell_eval_u = cell_build_eval_u * cell_u;
	  
	  //
	  // Over the cell: evaluate global gradient phi_j.
	  //
	  
	  //
	  // Over the cell: evaluate global gradient phi_j.
	  //
	  
	  //
	  // Over the cell: evaluate global gradient phi_j.
	  //
	}      
    }
  

  //
  // Advection equation.
  //
  //
  // ( u.nabla(phi_j) ) test_i - \int_{Boundary} max(0,-u.n) phi_j test_i  - \int_{Boundary} max(0,-u.n) phi_j^{-} test_i = \int_{dirichlet} f test_i
  //
  // Diagonal block ( u.nabla(phi_j) ) test_i - \int_{Boundary} max(0,-u.n) phi_j test_i
  // Extra-diagonal block   - \int_{Boundary} max(0,-u.n) phi_j^{-} test_i
  //

  //
  // Compute global gradient.
  //
  
  //
  // Diagonal block over the cell.
  //
  {
    T nabla_phi_j[3],velocity[3];
    for (wmesh_int_t j=0;j<lmat_n_;++j)
      {
	for (wmesh_int_t i=0;i<lmat_n_;++i)
	  {
	    T mij =  static_cast<T>(0);
	    for (wmesh_int_t k=0;k<q_n;++k)
	      {
		//
		// Get velocity.
		//
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    velocity[l] = q_velocity[k*topodim+l];
		  }		
		
		//
		// Get nabla_phi_j
		//
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    nabla_phi_j[l] = shape_eval_fglobal->m_diff[l].v[shape_eval_fglobal->m_diff[l].ld*k+j];
		  }
		
		//
		// Get test_i
		//
		T test_i 	= shape_eval_test->m_f.v[shape_eval_test->m_f.ld*k+i];
		T dot 		= static_cast<T>(0);
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    dot += nabla_phi_j[l] * velocity[l];
		  }
		
		mij += dets[k] * q_w[k] * (  dot * test_i  );
	      }	  	  
	  }
      }
  }

  //
  // Diagonal block over the boundary faces.
  //

  //
  // Loop over the facets.
  //  
  {
    T nabla_phi_j[3],velocity[3];
    for (wmesh_int_t j=0;j<lmat_n_;++j)
      {
	for (wmesh_int_t i=0;i<lmat_n_;++i)
	  {
	    T mij =  static_cast<T>(0);
	    for (wmesh_int_t k=0;k<q_n;++k)
	      {
		//
		// Get velocity.
		//
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    velocity[l] = q_velocity[k*topodim+l];
		  }		
		
		//
		// Get nabla_phi_j
		//
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    nabla_phi_j[l] = shape_eval_fglobal->m_diff[l].v[shape_eval_fglobal->m_diff[l].ld*k+j];
		  }
		
		//
		// Get test_i
		//
		T test_i 	= shape_eval_test->m_f.v[shape_eval_test->m_f.ld*k+i];
		T dot 		= static_cast<T>(0);
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    dot += nabla_phi_j[l] * velocity[l];
		  }
		
		mij += dets[k] * q_w[k] * (  dot * test_i  );
	      }	  	  
	  }
      }
  }

  
  
  
  //
  // Evaluate velocity field over gauss points over the cell.
  //
  BLAS_dgemm("N",
	     "N",
	     &topodim,
	     &q_n,
	     &u_n,
	     &r1,
	     uelm_,
	     &uelm_ld_,
	     shape_eval_u->m_f.v,
	     shape_eval_u->m_f.ld,
	     &r0,
	     q_velocity,
	     &topodim);

  //
  // Compute the global gradient of phi_j
  //
  {
    T lnabla_phi_j[3];
    T nabla_phi_j[3];
    for (wmesh_int_t k=0;k<q_n;++k)
      {
	T*jacobian = jacobians + topodimXtopodim * k;
	for (wmesh_int_t j=0;j<lmat_n_;++j)
	  {
	    //
	    // Compute the global gradient of shape function j at point k.
	    //
	    for (wmesh_int_t idim=0;idim<topodim;++idim)
	      {
		lnabla_phi_j[idim] = shape_eval_f->m_diff[idim].v[shape_eval_f->m_diff[idim].ld * k + j];
	      }
	    
	    BLAS_dgemm("T",
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
	    
	    //
	    // Compute the global gradient of shape function j at point k.
	    //
	    for (wmesh_int_t idim=0;idim<topodim;++idim)
	      {
		shape_eval_fglobal->m_diff[idim].v[shape_eval_fglobal->m_diff[idim].ld * k + j] = lnabla_phi_j[idim];
	      }	  
	  }
      }
  }



  
  // 
  // diagonal block integral_K (u(x_k) dr_j + v(x_k) ds_j(x_k) + w(x_k)*dt_j(x_k)) * test_i(x_k)
  //
  for (wmesh_int_t j=0;j<lmat_n_;++j)
    {
      for (wmesh_int_t i=0;i<lmat_n_;++i)
	{
	  double val=0.0;
	  for (wmesh_int_t k=0;k<q_n;++k)
	    {
	      double dot = 0.0;
	      for (wmesh_int_t idim=0;idim<topodim;++idim)
		{
		  dot += q_velocity[topodim * k + idim] * shape_eval_fglobal->m_diff[idim].v[shape_eval_fglobal->m_diff[idim].ld * k + j];
		}
	      
	      val += dot * shape_eval_test->m_f.v[shape_eval_test->m_f.ld * k + idim] * q_w[k] * dets[k];
	    }
	  
	  lmat[j*lmat_n_ + i] += val;
	}
    }

  // 
  // extra diagonal block integral_K (u(x_k) dr_j + v(x_k) ds_j(x_k) + w(x_k)*dt_j(x_k)) * test_i(x_k)
  //
  
  // 
  // sum over faces integral_K -max(0, u(x_k) nx(x_k) + v(x_k) ny(x_k) + w(x_k)*nz(x_k) ) * test_i(x_k)
  //
  for (wmesh_int_t i=0;i<numfaces;++i)
    {
      //
      // Compute normals.
      //
      
      //
      // Compute velocity over the face.
      //
      BLAS_dgemm("N",
		 "N",
		 &topodim,
		 &q_n,
		 &u_n,
		 &r1,
		 uelm_,
		 &uelm_ld_,
		 shape_eval_u->m_f.v,
		 shape_eval_u->m_f.ld,
		 &r0,
		 q_velocity,
		 &topodim);

      //
      // 
      //
      
      BLAS_dgemm("N",
		 "N",
		 &topodim,
		 &q_n,
		 &u_n,
		 &r1,
		 uelm_,
		 &uelm_ld_,
		 shape_eval_u->m_f.v,
		 shape_eval_u->m_f.ld,
		 &r0,
		 q_velocity,
		 &topodim);
      
    }

  
  static constexpr wmesh_int_t n1 = static_cast<wmesh_int_t>(1);
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      T*jacobian = jacobians + topodimXtopodim * k;
      for (wmesh_int_t j=0;j<lmat_n_;++j)
	{
	  //
	  // Compute the global gradient of shape function j at point k.
	  //
	  for (wmesh_int_t idim=0;idim<topodim;++idim)
	    {
	      lnabla_phi_j[idim] = shape_eval_f->m_diff[idim].v[shape_eval_f->m_diff[idim].ld * k + j];
	    }
	  
	  BLAS_dgemm("T",
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
	      //
	      // Compute dot product of the global gradient of shape function i
	      // and the global gradient of shape function j.
	      //
	      T dot = r0;
	      for (wmesh_int_t idim=0;idim<topodim;++idim)
		{
		  dot += q_velocity[idim]*nabla_phi_j[idim];
		}
#if 0
	      if (i==2)
		{
		  std::cout << "lnabla phi_i = " << lnabla_phi_i[0] << " " << lnabla_phi_i[1] << std::endl;
		  std::cout << "lnabla phi_j = " << lnabla_phi_j[0] << " " << lnabla_phi_j[1] << std::endl;
		  std::cout << "nabla phi_i = " << nabla_phi_i[0] << " " << nabla_phi_i[1] << std::endl;
		  std::cout << "nabla phi_j = " << nabla_phi_j[0] << " " << nabla_phi_j[1] << std::endl;
		  std::cout << "j= " << j << std::endl;
		  std::cout << "j= " << j << std::endl;
		  std::cout << "j= " << j << std::endl;
		  std::cout << "dot " << dot<<std::endl;
		  std::cout << "q_w " << q_w[k]<<std::endl;
		  std::cout << "dets " << dets[k]<<std::endl;
		}
#endif
	      lmat_[lmat_ld_*j+i] += dot *  q_w[k]  *  dets[k];
	    }
	  
	}      
    }

  //
  // Reset the local rhs.
  //
  for (wmesh_int_t i=0;i<lmat_m_;++i)
    {
      lrhs_[lrhs_inc_*i] = r0;	      
    }

  for (wmesh_int_t k=0;k<q_n;++k)
    {
      for (wmesh_int_t i=0;i<lmat_m_;++i)
	{
	  lrhs_[lrhs_inc_*i] += shape_eval_f->m_f.v[k*shape_eval_f->m_f.ld + i] * dets[k] * q_w[k];
	}      
    }
#endif  
  //
  // Calculate the jacobian.
  //  
  return WMESH_STATUS_SUCCESS;
}



template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_global_system_zone(const wmesh_t *		self_,
							     wmesh_int_t			itype_,
							     const wmeshspacedg_fem_advection_t<T> *	fem_,
							     const wmeshspace_t *		velocity_space_,
							     wmesh_int_t 			velocity_storage_,
							     wmesh_int_t 			velocity_m_,
							     wmesh_int_t 			velocity_n_,
							     const double * 			velocity_,
							     wmesh_int_t 			velocity_ld_,
							     wmesh_int_t			csr_size_,
							     const_wmesh_int_p			csr_ptr_,
							     const_wmesh_int_p			csr_ind_,
							     T * 				csr_val_,
							     T * 				rhs_)
{
  const wmesh_int_t num_cells = self_->m_c2n.m_n[itype_];
  for (wmesh_int_t ielm=0;ielm < num_cells;++ielm)
    {
      
      //
      // 1/ Load element connectivity.
      //
      //
      // 2/ Load velocity connectivity.
      //
      
      //
      // 3/ Load element coordinates.
      //
      //
      // 4/ Load velocity dofs.
      //

      // 
      // 5/ Compute jacobians on cell cubature.
      //
      
      //
      // 6/ Eval velocity on cell cubature.
      //

      // 
      // 7/ Compute nabla_phi_j on  cell cubature.
      //

      //
      // Diagonal part of the flux (no need to rotate here. )
      // 


      
      // 
      // 7/ Compute normals on boundary.
      //
      
      // 
      // 8/ Compute velocity on boundary.
      //
      
      // 
      // 9/ Evaluate max(0,-u.n)
      //


      for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
	{
	  c2c_tmpi[ifacet] = c2c_v[c2c_ptr[itype] + c2c_ld[itype] * ielm  + ifacet];
	}
      
      //
      // Now compute extra-diagonal blocks
      //
      for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
	{

	  auto c2c_info_i = c2c_tmpi[ifacet]; 	  
	  if (c2c_info_i != 0)
	    {
	      wmesh_int_t isigned_rotation;
	      //
	      // Eextract the connectivity of ifacet-th facet in the ielm-th element of type itype.
	      //
	      
	      //
	      // Get the position of the minimum.
	      //

	      wmesh_int_t jfacet;
	      wmesh_int_t jsigned_rotation;
	      wmesh_int_t jelm;
	      wmesh_int_t jtype;
	      
	      status = bms_c2c_cindex(c2c_info_i,&jelm);
	      status = bms_c2c_ctype(c2c_info_i,&jtype);
	      
	      //
	      // Load c2c_j
	      //
	      for (wmesh_int_t l=0;l<num_facets;++l)
		{
		  c2c_tmpj[l] = c2c_v[c2c_ptr[jtype] + c2c_ld[jtype] * jelm  + l];
		}
	      
	      for (jfacet=0;jfacet<num_facets;++jfacet)
		{
		  auto c2c_info_j = c2c_tmpj[jfacet]; 		  
		  status = bms_c2c_cindex(c2c_info_j,&kelm);
		  status = bms_c2c_ctype(c2c_info_j,&ktype);
		  if ( (ielm == kelm) && (itype == ktype) )
		    {
		      break;
		    }
		}

	      //
	      // Now extract the connectivity of jfacet-th facet in the jelm-th element of type jtype.
	      //
	      
	      //
	      // Get the position of the minimum.
	      //

	      //
	      // Cubature.
	      //

	      //
	      // Assembly.
	      //

	    }
	}
    }

  
#if 0  
  const wmesh_int_t num_dofs_per_element 	= ;
  const wmesh_int_t num_elements 		= self_->m_c2n.m_n[itype_];
  WMESH_CHECK( num_elements > 0 );
  T cooelm[32];
  const wmesh_int_t cooelm_storage 	= WMESH_STORAGE_INTERLEAVE;
  const wmesh_int_t cooelm_m 		= self_->m_topology_dimension;
  const wmesh_int_t cooelm_n 		= self_->m_c2n.m_m[itype_];
  const wmesh_int_t cooelm_ld 		= self_->m_topology_dimension;
  
  wmesh_status_t status;
  const wmesh_int_t element = (cooelm_m == 1) ? 1 + itype_ : ((cooelm_m == 2) ? 2 + itype_ : 4 + itype_);
  

  const wmesh_int_t lmat_m  = num_dofs_per_element;
  const wmesh_int_t lmat_n  = num_dofs_per_element;
  const wmesh_int_t lmat_ld = num_dofs_per_element;
  const wmesh_int_t lrhs_inc = 1;

  
  T * lmat = (T*)malloc(sizeof(T)*lmat_ld*lmat_n);
  T * lrhs = (T*)malloc(sizeof(T)*lmat_m);
  const wmesh_int_t velocity_num_ldofs = velocity_space_->m_c2d.m_m[itype_];
  T * umat = (T*)malloc(sizeof(T)* cooelm_m * velocity_num_ldofs);




  
  wmesh_int_t rw_n;
  status = wmeshspacedg_fem_advection_local_system_buffer_size(fem_,
							       element,
							       &rw_n);
  WMESH_STATUS_CHECK(status);
  
  T * rw =  (T*)malloc(sizeof(T)*rw_n);
  wmesh_int_p c2d  = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*num_dofs_per_element);
  
  const wmesh_t * mesh = self_;

  for (wmesh_int_t ielm=0;ielm<num_elements;++ielm)
    {      
      //
      // get coordinates of element.
      //
      status = wmeshspace_get_cooelm(mesh,
				     itype_,
				     ielm,
				     cooelm_storage,
				     cooelm_m,
				     cooelm_n,
				     cooelm,
				     cooelm_ld);
      WMESH_STATUS_CHECK(status);

      //
      // get velocity.
      //
      status = wmeshspace_get_uelm(velocity_space_,
				   itype_,
				   ielm,
				   velocity_storage_,
				   velocity_m_,
				   velocity_n_,
				   velocity_,
				   velocity_ld_,
				   uelm_storage,
				   uelm_m,
				   uelm_n,
				   uelm,
				   uelm_ld);
      WMESH_STATUS_CHECK(status);

      
      status = wmeshspacedg_fem_advection_local_system(fem_,
						       element,
						       
						       cooelm_storage,
						       cooelm_m,
						       cooelm_n,
						       cooelm,
						       cooelm_ld,
						       
						       uelm_storage,
						       uelm_m,
						       uelm_n,
						       uelm,
						       uelm_ld,
						       
						       lmat_m,
						       lmat_n,
						       lmat,
						       lmat_ld,
						       lrhs,
						       lrhs_inc,
						       rw_n,
						       rw);
      WMESH_STATUS_CHECK(status);


      
#if 0
      

      
      for (wmesh_int_t i=0;i<lmat_n;++i)
	{
	  std::cout << "rhs before " << lrhs[lrhs_inc*i];	      
	  std::cout << std::endl;
	}
      std::cout << "lmat_n " << lmat_n << std::endl;
      for (int i=0;i<3+3*(self_->m_degree-1);++i)
	{
#if 1
	  for (int j=0;j<lmat_n;++j)
	    {
	      lmat[j*lmat_ld + i] = 0;
	    }
	  lmat[i*lmat_ld + i] = 1;
#endif
	  lrhs[i*lrhs_inc] = 0;
	}
      for (wmesh_int_t i=0;i<lmat_n;++i)
	{
	  std::cout << "rhs after["<<i <<"] " << lrhs[lrhs_inc*i];	      
	  std::cout << std::endl;
	}
     
      std::cout << "-----------------------------------------------------" << std::endl;
      std::cout << "-----------------------------------------------------" << std::endl;
      std::cout << "-----------------------------------------------------" << std::endl;
      for (wmesh_int_t i=0;i<lmat_n;++i)
	{
	  for (wmesh_int_t j=0;j<lmat_n;++j)
	    {	      
	      std::cout << " " << lmat[j*lmat_ld + i];	      
	    }
	  std::cout << std::endl;
	}
      std::cout << "-----------------------------------------------------" << std::endl;
      for (wmesh_int_t i=0;i<lmat_n;++i)
	{
	  std::cout << " " << lrhs[lrhs_inc*i];	      
	  std::cout << std::endl;
	}
      std::cout << "-----------------------------------------------------" << std::endl;
      std::cout << "-----------------------------------------------------" << std::endl;
      std::cout << "-----------------------------------------------------" << std::endl;
      wmesh_int_t perm[2048],n1=1,info_lapack;
      LAPACK_dgesv((wmesh_int_p)&lmat_n,
		   (wmesh_int_p)&n1,
		   lmat,
		   (wmesh_int_p)&lmat_n,
		   perm,
		   lrhs,
		   (wmesh_int_p)&lmat_n,
		   (wmesh_int_p)&info_lapack);
      std::cout << "-----------------------------------------------------" << std::endl;
      for (wmesh_int_t i=0;i<lmat_n;++i)
	{
	  std::cout << " " << lrhs[lrhs_inc*i];	      
	  std::cout << std::endl;
	}
      std::cout << "-----------------------------------------------------" << std::endl;
      std::cout << "-----------------------------------------------------" << std::endl;
      std::cout << "-----------------------------------------------------" << std::endl;
      
      exit(1);
#endif 



      
      
      //
      // Get dof indices.
      //
      for (wmesh_int_t i=0;i<self_->m_c2d.m_m[itype_];++i)
	{
	  c2d[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[itype_] + self_->m_c2d.m_ld[itype_] * ielm + i] - 1;
	}

      

      //
      // Assembly.
      //
      for (wmesh_int_t i=0;i<num_dofs_per_element;++i)
	{
	  const wmesh_int_t idof = c2d[i];
	  for (wmesh_int_t j=0;j<num_dofs_per_element;++j)
	    {
	      const wmesh_int_t jdof = c2d[j];
	      
	      //
	      // Search in the csr matrix.
	      //
	      bool found = false;
	      for (wmesh_int_t at = csr_ptr_[idof];at < csr_ptr_[idof+1];++at)
		{
		  if (csr_ind_[at]==0)
		    {
		  std::cerr << "pblm aaaaaaaassembly " << std::endl;
		  exit(1);

		    }
		  //  std::cout << csr_ind_[at] << std::endl;
		  if (csr_ind_[at] - 1 == jdof)
		    {
		      csr_val_[at] += lmat[j*num_dofs_per_element + i];
		      found = true;
		      break;
		    }
		}
	      if (!found)
		{
		  std::cerr << "pblm assembly " << std::endl;
		  exit(1);
		}
	    }
	}
      
      for (wmesh_int_t i=0;i<num_dofs_per_element;++i)
	{
	  rhs_[c2d[i]] += lrhs[i];
	}
      
    }


  //
  // Now apply bc.
  //
  for (wmesh_int_t i=0;i<self_->m_ndofs;++i)
    {
      if (self_->m_dof_codes[i] == 1)
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

	  rhs_[i] = 0.0;
	  
	}
    }
  free(rw);
#endif
  return WMESH_STATUS_SUCCESS;
}


template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_global_system(const wmesh_t *				self_,							
							const wmeshspace_t *			velocity_space_,
							wmesh_int_t 				velocity_storage_,
							wmesh_int_t 				velocity_m_,
							wmesh_int_t 				velocity_n_,
							const double * 				velocity_,
							wmesh_int_t 				velocity_ld_,
							const wmeshspacedg_fem_advection_t<T> *	fem_,
							wmesh_int_t				csr_size_,
							const_wmesh_int_p			csr_ptr_,
							const_wmesh_int_p			csr_ind_,
							T * 					csr_val_,
							T * 					rhs_)
{
  const wmesh_int_t ntypes = self_->m_c2d.m_size;
  for (wmesh_int_t itype=0;itype<ntypes;++itype)
    {
      const wmesh_int_t num_elements = self_->m_c2d.m_n[itype];
      if (num_elements > 0)
	{
	  wmesh_status_t status = wmeshspacedg_fem_advection_global_system_zone(self_,
										itype,
										fem_,
										velocity_space_,
										velocity_storage_,
										velocity_m_,
										velocity_n_,
										velocity_,
										velocity_ld_,
										csr_size_,
										csr_ptr_,
										csr_ind_,
										csr_val_,
										rhs_);
	  WMESH_STATUS_CHECK(status);
	}
    }
  return WMESH_STATUS_SUCCESS;
}


  


extern "C"
{
  
  wmesh_status_t wmeshspacedg_advection(const wmesh_t*				self_,
					const wmesh_shape_t* __restrict__	shape_f_,
					const wmesh_shape_t* __restrict__	shape_u_,
					const wmesh_shape_t* __restrict__	shape_test_,
					
					const wmeshspace_t *			velocity_space_,

					wmesh_int_t 				velocity_storage_,
					wmesh_int_t 				velocity_m_,
					wmesh_int_t 				velocity_n_,
					const double * 				velocity_,
					wmesh_int_t 				velocity_ld_,
					
					wmesh_int_t				csr_size_,
					const_wmesh_int_p			csr_ptr_,
					const_wmesh_int_p			csr_ind_,
					double * 				csr_val_,
					
					double * 				rhs_)
  {
    static constexpr double r0 = static_cast<double>(0);

    static constexpr  wmesh_int_t 	shape_element_family 	= WMESH_SHAPE_FAMILY_LAGRANGE;
    static constexpr  wmesh_int_t 	shape_element_degree 	= 1;    
    static constexpr  wmesh_int_t 	cubature_family 	= WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE;

    wmesh_status_t status;
    
    //
    // Set csr_val and rhs_ to zero.
    //
    const wmesh_int_t N = csr_ptr_[csr_size_];
    for (wmesh_int_t i=0;i<N;++i)
      csr_val_[i] = r0; 
    for (wmesh_int_t i=0;i<csr_size_;++i)
      rhs_[i] = r0;
    

    const wmesh_int_t topodim 	= self_->m_topology_dimension;
    const wmesh_int_t 	cubature_degree 	= (shape_u_->m_degree + (shape_f_->m_degree-1) + shape_test_->m_degree) * 2;
    wmesh_int_t ntypes;
    wmesh_int_t elements[4];
    status =  bms_topodim2elements(topodim,
				   &ntypes,
				   elements);
    WMESH_STATUS_CHECK(status);

    
    //
    // Initialize data.
    //
    wmeshspacedg_fem_advection_t<double> fem{};
        
    std::cout << "dgadvection def .. " << std::endl;
    for (wmesh_int_t itype=0;itype<self_->m_c2n.m_size;++itype)
      {
	if (self_->m_c2n.m_n[itype] > 0)
	  {
	    status = wmeshspacedg_fem_advection_add(&fem,
						    elements[itype],
						    cubature_family,
						    cubature_degree,
						    shape_element_family,
						    shape_element_degree,				   
						    shape_f_->m_family,
						    shape_f_->m_degree,
						    shape_u_->m_family,
						    shape_u_->m_degree,
						    shape_f_->m_family,
						    shape_f_->m_degree);
	    WMESH_STATUS_CHECK(status);
	  }
      }

    status = wmeshspacedg_fem_advection_global_system(self_,
						      velocity_space_,
						      velocity_storage_,
						      velocity_m_,
						      velocity_n_,
						      velocity_,
						      velocity_ld_,
						      &fem,
						      csr_size_,
						      csr_ptr_,
						      csr_ind_,
						      csr_val_,
						      rhs_);
    WMESH_STATUS_CHECK(status);
    
    return WMESH_STATUS_SUCCESS;
  }
}

#if 0

template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_local_system_buffer_size(const wmeshspacedg_fem_advection_t<T> *__restrict__	self_,
							  wmesh_int_t 					element_,							
							  wmesh_int_p 					rw_n_)
{
  const wmesh_cubature_t<T>*__restrict__ 	cubature		= &self_->m_cubature[element_];
  const wmesh_int_t 				q_n 			= cubature->m_w.n;
  wmesh_status_t status;
  wmesh_int_t topodim;
  status = bms_element2topodim	(element_,
				 &topodim);
  WMESH_STATUS_CHECK(status);
  rw_n_[0] = topodim * topodim * q_n + q_n;
  return WMESH_STATUS_SUCCESS;
}

template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_local_system(const wmeshspacedg_fem_advection_t<T> *__restrict__	self_,
						       wmesh_int_t 					element_,
						       wmesh_int_t 					cooelm_storage_,
						       wmesh_int_t 					cooelm_m_,
						       wmesh_int_t 					cooelm_n_,
						       const T*__restrict__				cooelm_,
						       wmesh_int_t 					cooelm_ld_,
						       wmesh_int_t 					uelm_storage_,
						       wmesh_int_t 					uelm_m_,
						       wmesh_int_t 					uelm_n_,
						       const T*__restrict__				uelm_,
						       wmesh_int_t 					uelm_ld_,
						       wmesh_int_t 					lmat_m_,
						       wmesh_int_t 					lmat_n_,
						       T*__restrict__					lmat_,
						       wmesh_int_t 					lmat_ld_,
						       T*__restrict__					lrhs_,
						       wmesh_int_t 					lrhs_inc_,
						       wmesh_int_t 					rw_n_,
						       T*__restrict__					rw_)
{
  static constexpr T r0 = static_cast<T>(0);
  static constexpr T r1 = static_cast<T>(1);

  WMESH_CHECK_POINTER(self_);
  
  wmesh_status_t status;
  wmesh_int_t required_rw_n;
  status = wmeshspacedg_fem_advection_local_system_buffer_size(self_,
							       element_,
							       &required_rw_n);
  WMESH_STATUS_CHECK(status);
  
  if (rw_n_ < required_rw_n)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);       
    }
  
  const wmesh_shape_eval_t<T> * __restrict__ 	shape_eval_element 	= &self_->m_shape_eval_element[element_];
  const wmesh_shape_eval_t<T> * __restrict__	shape_eval_f 		= &self_->m_shape_eval_f[element_];
  const wmesh_shape_eval_t<T> * __restrict__	shape_eval_u 		= &self_->m_shape_eval_u[element_];
  const wmesh_shape_eval_t<T> * __restrict__	shape_eval_test 	= &self_->m_shape_eval_test[element_];  
  const wmesh_cubature_t<T>*__restrict__ 	cubature		= &self_->m_cubature[element_];

  const wmesh_int_t 				u_n 			= (uelm_storage_ == WMESH_STORAGE_INTERLEAVE) ? uelm_n_ : uelm_m_;
  const wmesh_int_t 				q_n 			= cubature->m_w.n;
  const T* __restrict__				q_w 			= cubature->m_w.v;  
  const wmesh_int_t 	topodim  					= cooelm_m_;
  const wmesh_int_t  	topodimXtopodim 				= topodim*topodim;;
  
  //
  // Compute Jacobians and determinants.
  //
  T *__restrict__ jacobians 	= rw_;
  T *__restrict__ dets 		= rw_ + topodimXtopodim * q_n;
  rw_ += required_rw_n;
  rw_n_ -= required_rw_n;

  
  for (wmesh_int_t idim=0;idim<topodim;++idim)
    {
      BLAS_dgemm("N",
		 "N",
		 &topodim,
		 &shape_eval_element->m_diff[idim].n,
		 &shape_eval_element->m_diff[idim].m,
		 &r1,
		 cooelm_,
		 &cooelm_ld_,
		 shape_eval_element->m_diff[idim].v,
		 &shape_eval_element->m_diff[idim].ld,
		 &r0,
		 jacobians + topodim*idim,
		 &topodimXtopodim);
    }
  
  T B[9];
  for (wmesh_int_t j=0;j<q_n;++j)
    {
      T*jacobian = jacobians + topodimXtopodim * j;
      T det;
      if (topodim==2)
	{
	  det = determinant2x2(jacobian,
			       topodim);
	}
      else if (topodim==3)
	{
	  det = determinant3x3(jacobian,
			       topodim);
	}
      else
	{
	  det = determinant1x1(jacobian,
			       topodim);
	}
          
      if (det < 0.0)
	det = -det;
      dets[j] = det;
      
      //
      // Reset identity
      //
      for (wmesh_int_t i=0;i<topodim*topodim;++i)
	{
	  B[i] = r0;
	}
      
      for (wmesh_int_t i=0;i<topodim;++i)
	{
	  B[i*topodim+i] = r1;
	}

      //
      // Inverse the jacobian.
      //
      
      {	wmesh_int_t info_lapack,perm[3];
	LAPACK_dgesv((wmesh_int_p)&topodim,
		     (wmesh_int_p)&topodim,
		     jacobian,
		     (wmesh_int_p)&topodim,
		     perm,
		     B,
		     (wmesh_int_p)&topodim,
		     (wmesh_int_p)&info_lapack);
	if (info_lapack)
	  {
	    fprintf(stderr,"error lapack\n");
	    exit(1);
	  } }

      //
      // Store the inverse.
      //
      for (wmesh_int_t i=0;i<topodimXtopodim;++i)
	{
	  jacobians[topodim*topodim * j + i] = B[i];
	}
    }


  //
  // Reset the local matrix and the local rhs.
  //
  for (wmesh_int_t j=0;j<lmat_n_;++j)
    {
      for (wmesh_int_t i=0;i<lmat_m_;++i)
	{
	  lmat_[lmat_ld_*j+i] = r0;	      
	}
    }

  for (wmesh_int_t i=0;i<lmat_m_;++i)
    {
      lrhs_[lrhs_inc_*i] = r0;	      
    }



  //
  // Advection equation.
  //
  //
  // ( u.nabla(phi_j) ) test_i - \int_{Boundary} max(0,-u.n) phi_j test_i  - \int_{Boundary} max(0,-u.n) phi_j^{-} test_i = \int_{dirichlet} f test_i
  //
  // Diagonal block ( u.nabla(phi_j) ) test_i - \int_{Boundary} max(0,-u.n) phi_j test_i
  // Extra-diagonal block   - \int_{Boundary} max(0,-u.n) phi_j^{-} test_i
  //

  //
  // Compute global gradient.
  //
  
  //
  // Diagonal block over the cell.
  //
  {
    T nabla_phi_j[3],velocity[3];
    for (wmesh_int_t j=0;j<lmat_n_;++j)
      {
	for (wmesh_int_t i=0;i<lmat_n_;++i)
	  {
	    T mij =  static_cast<T>(0);
	    for (wmesh_int_t k=0;k<q_n;++k)
	      {
		//
		// Get velocity.
		//
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    velocity[l] = q_velocity[k*topodim+l];
		  }		
		
		//
		// Get nabla_phi_j
		//
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    nabla_phi_j[l] = shape_eval_fglobal->m_diff[l].v[shape_eval_fglobal->m_diff[l].ld*k+j];
		  }
		
		//
		// Get test_i
		//
		T test_i 	= shape_eval_test->m_f.v[shape_eval_test->m_f.ld*k+i];
		T dot 		= static_cast<T>(0);
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    dot += nabla_phi_j[l] * velocity[l];
		  }
		
		mij += dets[k] * q_w[k] * (  dot * test_i  );
	      }	  	  
	  }
      }
  }

  //
  // Diagonal block over the boundary faces.
  //

  //
  // Loop over the facets.
  //  
  {
    T nabla_phi_j[3],velocity[3];
    for (wmesh_int_t j=0;j<lmat_n_;++j)
      {
	for (wmesh_int_t i=0;i<lmat_n_;++i)
	  {
	    T mij =  static_cast<T>(0);
	    for (wmesh_int_t k=0;k<q_n;++k)
	      {
		//
		// Get velocity.
		//
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    velocity[l] = q_velocity[k*topodim+l];
		  }		
		
		//
		// Get nabla_phi_j
		//
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    nabla_phi_j[l] = shape_eval_fglobal->m_diff[l].v[shape_eval_fglobal->m_diff[l].ld*k+j];
		  }
		
		//
		// Get test_i
		//
		T test_i 	= shape_eval_test->m_f.v[shape_eval_test->m_f.ld*k+i];
		T dot 		= static_cast<T>(0);
		for (wmesh_int_t l=0;l<topodim;++l)
		  {
		    dot += nabla_phi_j[l] * velocity[l];
		  }
		
		mij += dets[k] * q_w[k] * (  dot * test_i  );
	      }	  	  
	  }
      }
  }

  
  
  
  //
  // Evaluate velocity field over gauss points over the cell.
  //
  BLAS_dgemm("N",
	     "N",
	     &topodim,
	     &q_n,
	     &u_n,
	     &r1,
	     uelm_,
	     &uelm_ld_,
	     shape_eval_u->m_f.v,
	     shape_eval_u->m_f.ld,
	     &r0,
	     q_velocity,
	     &topodim);

  //
  // Compute the global gradient of phi_j
  //
  {
    T lnabla_phi_j[3];
    T nabla_phi_j[3];
    for (wmesh_int_t k=0;k<q_n;++k)
      {
	T*jacobian = jacobians + topodimXtopodim * k;
	for (wmesh_int_t j=0;j<lmat_n_;++j)
	  {
	    //
	    // Compute the global gradient of shape function j at point k.
	    //
	    for (wmesh_int_t idim=0;idim<topodim;++idim)
	      {
		lnabla_phi_j[idim] = shape_eval_f->m_diff[idim].v[shape_eval_f->m_diff[idim].ld * k + j];
	      }
	    
	    BLAS_dgemm("T",
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
	    
	    //
	    // Compute the global gradient of shape function j at point k.
	    //
	    for (wmesh_int_t idim=0;idim<topodim;++idim)
	      {
		shape_eval_fglobal->m_diff[idim].v[shape_eval_fglobal->m_diff[idim].ld * k + j] = lnabla_phi_j[idim];
	      }	  
	  }
      }
  }



  
  // 
  // diagonal block integral_K (u(x_k) dr_j + v(x_k) ds_j(x_k) + w(x_k)*dt_j(x_k)) * test_i(x_k)
  //
  for (wmesh_int_t j=0;j<lmat_n_;++j)
    {
      for (wmesh_int_t i=0;i<lmat_n_;++i)
	{
	  double val=0.0;
	  for (wmesh_int_t k=0;k<q_n;++k)
	    {
	      double dot = 0.0;
	      for (wmesh_int_t idim=0;idim<topodim;++idim)
		{
		  dot += q_velocity[topodim * k + idim] * shape_eval_fglobal->m_diff[idim].v[shape_eval_fglobal->m_diff[idim].ld * k + j];
		}
	      
	      val += dot * shape_eval_test->m_f.v[shape_eval_test->m_f.ld * k + idim] * q_w[k] * dets[k];
	    }
	  
	  lmat[j*lmat_n_ + i] += val;
	}
    }

  // 
  // extra diagonal block integral_K (u(x_k) dr_j + v(x_k) ds_j(x_k) + w(x_k)*dt_j(x_k)) * test_i(x_k)
  //
  
  // 
  // sum over faces integral_K -max(0, u(x_k) nx(x_k) + v(x_k) ny(x_k) + w(x_k)*nz(x_k) ) * test_i(x_k)
  //
  for (wmesh_int_t i=0;i<numfaces;++i)
    {
      //
      // Compute normals.
      //
      
      //
      // Compute velocity over the face.
      //
      BLAS_dgemm("N",
		 "N",
		 &topodim,
		 &q_n,
		 &u_n,
		 &r1,
		 uelm_,
		 &uelm_ld_,
		 shape_eval_u->m_f.v,
		 shape_eval_u->m_f.ld,
		 &r0,
		 q_velocity,
		 &topodim);

      //
      // 
      //
      
      BLAS_dgemm("N",
		 "N",
		 &topodim,
		 &q_n,
		 &u_n,
		 &r1,
		 uelm_,
		 &uelm_ld_,
		 shape_eval_u->m_f.v,
		 shape_eval_u->m_f.ld,
		 &r0,
		 q_velocity,
		 &topodim);
      
    }

  
  static constexpr wmesh_int_t n1 = static_cast<wmesh_int_t>(1);
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      T*jacobian = jacobians + topodimXtopodim * k;
      for (wmesh_int_t j=0;j<lmat_n_;++j)
	{
	  //
	  // Compute the global gradient of shape function j at point k.
	  //
	  for (wmesh_int_t idim=0;idim<topodim;++idim)
	    {
	      lnabla_phi_j[idim] = shape_eval_f->m_diff[idim].v[shape_eval_f->m_diff[idim].ld * k + j];
	    }
	  
	  BLAS_dgemm("T",
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
	      //
	      // Compute dot product of the global gradient of shape function i
	      // and the global gradient of shape function j.
	      //
	      T dot = r0;
	      for (wmesh_int_t idim=0;idim<topodim;++idim)
		{
		  dot += q_velocity[idim]*nabla_phi_j[idim];
		}
#if 0
	      if (i==2)
		{
		  std::cout << "lnabla phi_i = " << lnabla_phi_i[0] << " " << lnabla_phi_i[1] << std::endl;
		  std::cout << "lnabla phi_j = " << lnabla_phi_j[0] << " " << lnabla_phi_j[1] << std::endl;
		  std::cout << "nabla phi_i = " << nabla_phi_i[0] << " " << nabla_phi_i[1] << std::endl;
		  std::cout << "nabla phi_j = " << nabla_phi_j[0] << " " << nabla_phi_j[1] << std::endl;
		  std::cout << "j= " << j << std::endl;
		  std::cout << "j= " << j << std::endl;
		  std::cout << "j= " << j << std::endl;
		  std::cout << "dot " << dot<<std::endl;
		  std::cout << "q_w " << q_w[k]<<std::endl;
		  std::cout << "dets " << dets[k]<<std::endl;
		}
#endif
	      lmat_[lmat_ld_*j+i] += dot *  q_w[k]  *  dets[k];
	    }
	  
	}      
    }

  //
  // Reset the local rhs.
  //
  for (wmesh_int_t i=0;i<lmat_m_;++i)
    {
      lrhs_[lrhs_inc_*i] = r0;	      
    }

  for (wmesh_int_t k=0;k<q_n;++k)
    {
      for (wmesh_int_t i=0;i<lmat_m_;++i)
	{
	  lrhs_[lrhs_inc_*i] += shape_eval_f->m_f.v[k*shape_eval_f->m_f.ld + i] * dets[k] * q_w[k];
	}      
    }
  
  //
  // Calculate the jacobian.
  //  
  return WMESH_STATUS_SUCCESS;
}


template <typename T>
wmesh_status_t wmeshspace_get_cooelm(const wmesh_t * __restrict__	self_,
				     wmesh_int_t			itype_,
				     wmesh_int_t			ielm_,
				     wmesh_int_t			cooelm_storage_,
				     wmesh_int_t			cooelm_m_,
				     wmesh_int_t			cooelm_n_,
				     T * __restrict__			cooelm_,
				     wmesh_int_t			cooelm_ld_)
{
  WMESH_CHECK_POINTER(self_);
  WMESH_CHECK_POINTER(cooelm_);
  switch(cooelm_storage_)
    {
    case WMESH_STORAGE_INTERLEAVE:
      {
	for (wmesh_int_t j=0;j<cooelm_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<cooelm_m_;++i)
	      {
		cooelm_[cooelm_ld_*j+i]
		  = self_->m_coo[self_->m_coo_ld * ( self_->m_c2n.m_data[self_->m_c2n.m_ptr[itype_] + self_->m_c2n.m_ld[itype_] * ielm_ + j] - 1) + i];
	      }
	  }
	return WMESH_STATUS_SUCCESS;  
      }

    case WMESH_STORAGE_BLOCK:
      {
	for (wmesh_int_t j=0;j<cooelm_n_;++j)
	  {
	    for (wmesh_int_t i=0;i<cooelm_m_;++i)
	      {
		cooelm_[cooelm_ld_*i+j]
		  = self_->m_coo[self_->m_coo_ld * ( self_->m_c2n.m_data[self_->m_c2n.m_ptr[itype_] + self_->m_c2n.m_ld[itype_] * ielm_ + j] - 1) + i];
	      }
	  }
	return WMESH_STATUS_SUCCESS;  
      }

    }
  return WMESH_STATUS_INVALID_ENUM;  
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
#endif
#endif
