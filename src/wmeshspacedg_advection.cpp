/*
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
#include "wmesh-blas.hpp"
#include "bms.hpp"

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
inline void wmesh_mat_gemm(T 				alpha_,
			   const wmesh_mat_t<T>&	a_,
			   const wmesh_mat_t<T>&	b_,
			   T 				beta_,
			   wmesh_mat_t<T>&		c_)
{
  xgemm("N",
	"N",
	&c_.m,
	&c_.n,
	&a_.n,
	&alpha_,
	a_.v,
	&a_.ld,
	b_.v,
	&b_.ld,
	&beta_,
	c_.v,
	&c_.ld);  
}

template<typename T>
inline wmesh_status_t wmesh_mat_gemm(const char * 		transa_,
			   const char * 		transb_,
			   T 				alpha_,
			   const wmesh_mat_t<T>&	a_,
			   const wmesh_mat_t<T>&	b_,
			   T 				beta_,
			   wmesh_mat_t<T>&		c_)
{
#ifndef NDEBUG
  wmesh_int_t am = (transa_[0] == 'N') ? &a_.m : &a_.n;
  wmesh_int_t an = (transa_[0] == 'N') ? &a_.n : &a_.m;
  wmesh_int_t bm = (transb_[0] == 'N') ? &b_.m : &b_.n;
  wmesh_int_t bn = (transb_[0] == 'N') ? &b_.n : &b_.m;
  wmesh_int_t cm = b_.m;
  wmesh_int_t cn = b_.n;
  WMESH_CHECK(an == bm);
  WMESH_CHECK(am == cm);
  WMESH_CHECK(bn == cn);
#endif
  
  xgemm(transa_,
	transb_,
	&c_.m,
	&c_.n,
	(transa_[0] == 'N') ? &a_.n : &a_.m,
	&alpha_,
	a_.v,
	&a_.ld,
	b_.v,
	&b_.ld,
	&beta_,
	c_.v,
	&c_.ld);
  return WMESH_STATUS_SUCCESS;
}


template <typename T>
wmesh_status_t bms_element_jacobians(wmesh_int_t 					element_,
				     wmesh_int_t 					cooelm_storage_,
				     const wmesh_mat_t<T>&				cooelm_,
				     wmesh_int_t 					element_shape_eval_diff_storage_,
				     const wmesh_mat_t<T>*				element_shape_eval_diff_,
				     wmesh_mat_t<T>&					jacobians_,
				     wmesh_mat_t<T>&					jacobians_det_)
{
  static constexpr T r0 = static_cast<T>(0);
  static constexpr T r1 = static_cast<T>(1);
  
  wmesh_status_t 	status;
  const wmesh_int_t 	topodim 	= (cooelm_storage_ == WMESH_STORAGE_INTERLEAVE) ? cooelm_.m : cooelm_.n;
  const wmesh_int_t  	topodimXtopodim	= topodim*topodim;
  const wmesh_int_t	num_eval_nodes  = (element_shape_eval_diff_storage_ == WMESH_STORAGE_INTERLEAVE) ? element_shape_eval_diff_[0].n : element_shape_eval_diff_[0].m;
  wmesh_mat_t<T> 	dim_jacobians;
  T B[9];
  for (wmesh_int_t idim=0;idim<topodim;++idim)
    {
      wmesh_mat_t<T>::define(&dim_jacobians,
			     topodim,
			     num_eval_nodes,
			     jacobians_.v + topodim*idim,
			     topodimXtopodim);
      
      wmesh_mat_gemm((cooelm_storage_ == WMESH_STORAGE_INTERLEAVE) ? "N" : "T",
		     (element_shape_eval_diff_storage_ == WMESH_STORAGE_INTERLEAVE) ? "N" : "T",
		     r1,		     
		     cooelm_,
		     element_shape_eval_diff_[idim],
		     r0,
		     dim_jacobians);      
    }
  
  for (wmesh_int_t j=0;j<num_eval_nodes;++j)
    {
      T*jacobian = jacobians_.v + jacobians_.ld * j;
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
	{
	  std::cout << "// WARNING " << "NEGATIVE DETERMINANT" << std::endl;
	  det = -det;
	}
      
      jacobians_det_.v[jacobians_det_.ld * j + 0] = det;
      
      //
      // Reset identity
      //
      for (wmesh_int_t i=0;i<topodimXtopodim;++i)
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
      {
	wmesh_int_t info_lapack, perm[3];
	LAPACK_dgesv((wmesh_int_p)&topodim,
		     (wmesh_int_p)&topodim,
		     jacobian,
		     (wmesh_int_p)&topodim,
		     perm,
		     B,
		     (wmesh_int_p)&topodim,
		     (wmesh_int_p)&info_lapack);
	WMESH_CHECK(info_lapack == 0);
      }
      
      for (wmesh_int_t i=0;i<topodimXtopodim;++i)
	{
	  jacobians_.v[jacobians_.ld * j + i] = B[i];
	}
    }
  return WMESH_STATUS_SUCCESS;
}


#if 0
template<typename T>
struct wmesh_shape_eval_boundary_restriction_t
{
  //
  //       -2 -1  1 2
  //    -3 -2 -1  1 2 3 
  // -4 -3 -2 -1  1 2 3 4  
  wmesh_shape_t 		m_cell_shape;
  wmesh_shape_eval_t<T>   	m_cell_shape_eval_facets[2][8];
  wmesh_mat_t<T>   		m_facets_nodes_facets[2][8];
};


template<typename T>
void wmesh_shape_eval_boundary_restriction_def(wmesh_shape_eval_boundary_restriction_t*self_,
					       wmesh_int_t element_,
					       wmesh_int_t family_,
					       wmesh_int_t degree_)
{
  memset(self_,0,sizeof(wmesh_shape_eval_boundary_restriction_t<T>));
  wmesh_status status = wmesh_shape_def(&self_->m_cell_shape,
					element_,
					family_,
					degree_);
  WMESH_STATUS_CHECK(status);
  
}
#endif


template<typename T>
struct wmesh_shape_eval_boundary_t
{
  wmesh_shape_t 	m_shape;
  wmesh_int_t 		m_facets_nodes_storage;								  
  wmesh_mat_t<T>   	m_facets_nodes[2][8];
  wmesh_shape_eval_t<T> m_facets_shape_eval[2][8];
};

template<typename T>
wmesh_status_t wmesh_shape_eval_boundary_def(wmesh_shape_eval_boundary_t<T>*__restrict__ 	self_,
					     wmesh_int_t 					element_,				
					     wmesh_int_t 					shape_family_,
					     wmesh_int_t 					shape_degree_,
					     wmesh_int_t 					nodes_storage_,
					     const wmesh_mat_t<T>				nodes_)
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
  T ref_coordinates[32];
  status =  bms_ordering_vertices(element_,
				  ref_coordinates);  
  WMESH_STATUS_CHECK(status);
  
  for (wmesh_int_t ifacet = 0;ifacet < num_facets;++ifacet)
    {      
      const wmesh_int_t element_facet 	= facets[ifacet];
      const wmesh_int_t facet_num_nodes	= facets_num_nodes[ifacet];

#if 0
      //
      // Transform the nodes into the world of the cell coordinates.
      // 3xn = 3x2 2xn (r,s) -> (u,v,w)
      wmesh_mat_gemm();
#endif      
      //
      // Define the nodes.
      //
      {
	wmesh_int_t i = element_facet - ( (topodim_boundary==2) ? 2 : 1);
	status =  wmesh_mat_t<T>::alloc(&self_->m_facets_nodes[i][0],
					topodim,
					facet_num_nodes);
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
	      const wmesh_int_t i 	= element_facet - ( (topodim_boundary==2) ? 2 : 1);
	      const wmesh_int_t lidx 	= (signed_rotation < 0) ? (facet_num_nodes + (-signed_rotation - 1)) : (signed_rotation - 1);
	      
	      status =  wmesh_mat_t<T>::alloc(&self_->m_facets_nodes[i][lidx],
					      topodim,
					      facet_num_nodes);
	      WMESH_STATUS_CHECK(status);	      
	      status = bms_mirrored_local_coordinates(element_facet,
						      signed_rotation,						      
						      self_->m_facets_nodes_storage,
						      WMESH_MAT_FORWARD(self_->m_facets_nodes[i][0]),
						      self_->m_facets_nodes_storage,
						      WMESH_MAT_FORWARD(self_->m_facets_nodes[i][lidx]) );
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
	      k = (signed_rotation > 0) ? signed_rotation-1 : -signed_rotation-1;
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
wmesh_status_t wmeshspacedg_fem_advection_diagonal_block_flux_buffer_size(const wmeshspacedg_fem_advection_t<T> *__restrict__	self_,
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



wmesh_status_t wmesh_get_neighbor_info(const wmesh_t*self_)
{
  const wmesh_int_t num_cells = mesh->m_c2n.m_n[itype_];
  wmesh_int_t c2nelm[8],c2nelm_inc=1;
  wmesh_int_t c2celm[6],c2celm_inc=1;
  wmesh_int_t c2celm_types[6],c2celm_types_inc=1;
  for (wmesh_int_t ielm=0; ielm < num_cells;++ielm)
    {      
      //
      // Get connectivity of element.
      //
      status = wmesh_get_c2nelm(mesh,
				itype_,
				ielm,
				c2nelm,
				c2nelm_inc);
      WMESH_STATUS_CHECK(status);

      //
      // Get neighbors information.
      //
      wmesh_int_t num_facets = 0;
      status = wmesh_get_c2celm(mesh,
				itype_,
				ielm,
				&num_facets,
				c2celm,
				c2celm_inc,
				c2celm_types,
				c2celm_types_inc);
      WMESH_STATUS_CHECK(status);
      
      for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
	{
	  get_q2n(c2nelm,
		  lquadrilateral_idx,
		  q2n,
		  s_q2n_m_,
		  s_q2n_n_,
		  s_q2n_v_ + s_q2n_ptr_,
		  s_q2n_ld_);
	  
	  wmesh_int_t ifacet_type = 1;
	  	  
	  wmesh_int_t neicell_idx 	= c2celm[ifacet];
	  wmesh_int_t neicell_type	= c2celm_types[ifacet];
	  wmesh_int_t neicell_facet_idx	= -1;
	  
	  //
	  // Need to know jfacet.
	  //	  
	  wmesh_int_t neicell_num_facets = 0;
	  status = wmesh_get_c2celm(self_,
				    neicell_type,
				    neicell_idx,
				    &neicell_num_facets,
				    neicell_nei_ids,
				    1,
				    neicell_nei_types,
				    1);
	  WMESH_STATUS_CHECK(status);
	  for (wmesh_int_t l=0;l<neicell_num_facets;++l)
	    {
	      if (neicell_nei_types[l]>=0)
		{
		  if (neicell_nei_ids[l] == ielm)
		    {
		      neicell_facet_idx = l;
		      break;
		    }
		}
	    }
	  if (neicell_facet_idx < 0)
	    {
	      std::cout << "inconsistent c2c"  << std::endl;
	    }
	  
	  
	  //
	  // Get connectivity of the face.
	  //
	  wmesh_int_t c2n_facet[4],c2n_facet_inc=1;
	  status = wmesh_get_c2nfacet(self_,
				      itype_,
				      ielm,
				      ifacet,
				      c2n_facet,
				      c2n_facet_inc);
	  WMESH_STATUS_CHECK(status);
	  

	  wmesh_int_t nei_c2n_facet[4],nei_c2n_facet_inc=1;
	  status = wmesh_get_c2nfacet(self_,
				      jtype,
				      jelm,
				      jfacet,
				      nei_c2n_facet,
				      nei_c2n_facet_inc);
	  WMESH_STATUS_CHECK(status);


	  //
	  //
	  //
	  wmesh_int_t signed_rotation = 0;
	  {
	    wmesh_int_t min = c2n_facet[0];
	    for (wmesh_int_t l=1;l<;++l)
	      {
		if (min > c2n_facet[l])
		  {
		    min = c2n_facet[l];
		    signed_rotation = l;
		  }
	      }
	  }
	  
	  wmesh_int_t nei_signed_rotation = 0;
	  {
	    wmesh_int_t min = nei_c2n_facet[0];
	    for (wmesh_int_t l=1;l<;++l)
	      {
		if (min > nei_c2n_facet[l])
		  {
		    min = nei_c2n_facet[l];
		    nei_signed_rotation = l;
		  }
	      }
	  }
	}
}

//
// Only diagonal part.
//
template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_diagonal_block(const wmeshspacedg_fem_advection_t<T> *__restrict__	self_,
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
						       wmesh_int_t 					rw_n_,
						       T*__restrict__					rw_)
{
  static constexpr T r0 = static_cast<T>(0);
  static constexpr T r1 = static_cast<T>(1);

  WMESH_CHECK_POINTER(self_);
  
  wmesh_status_t status;
  wmesh_int_t required_rw_n;
  status = wmeshspacedg_fem_advection_diagonal_block_flux_buffer_size(self_,
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

  //  const wmesh_int_t 				u_n 			= (uelm_storage_ == WMESH_STORAGE_INTERLEAVE) ? uelm_n_ : uelm_m_;
  const wmesh_int_t 				q_n 			= cubature->m_w.n;
  const T* __restrict__				q_w 			= cubature->m_w.v;  
  const wmesh_int_t 	topodim  					= cooelm_m_;
  const wmesh_int_t  	topodimXtopodim 				= topodim*topodim;;
  
  //
  // Compute Jacobians and determinants.
  //
  std::cout << "    compute jacobians  " << std::endl;
  T *__restrict__ jacobians 	= rw_;
  T *__restrict__ dets 		= rw_ + topodimXtopodim * q_n;
  rw_ += required_rw_n;
  rw_n_ -= required_rw_n;

  
  for (wmesh_int_t idim=0;idim<topodim;++idim)
    {
      //      std::cout << "? " << shape_eval_element->m_diff[idim].m << std::endl;
      //      std::cout << "? " << shape_eval_element->m_diff[idim].ld << std::endl;
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

  std::cout << "    compute jacobians done.  " << std::endl;

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


  const wmesh_int_t num_cells = mesh->m_c2n.m_n[itype_];
  wmesh_int_t c2nelm[8],c2nelm_inc=1;
  wmesh_int_t c2celm[6],c2celm_inc=1;
  wmesh_int_t c2celm_types[6],c2celm_types_inc=1;
  for (wmesh_int_t ielm=0; ielm < num_cells;++ielm)
    {      
      //
      // Get connectivity of element.
      //
      status = wmesh_get_c2nelm(mesh,
				itype_,
				ielm,
				c2nelm,
				c2nelm_inc);
      WMESH_STATUS_CHECK(status);

      //
      // Get neighbors information.
      //
      wmesh_int_t num_facets = 0;
      status = wmesh_get_c2celm(mesh,
				itype_,
				ielm,
				&num_facets,
				c2celm,
				c2celm_inc,
				c2celm_types,
				c2celm_types_inc);
      WMESH_STATUS_CHECK(status);
      
      for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
	{
	  get_q2n(c2nelm,
		  lquadrilateral_idx,
		  q2n,
		  s_q2n_m_,
		  s_q2n_n_,
		  s_q2n_v_ + s_q2n_ptr_,
		  s_q2n_ld_);
	  
	  wmesh_int_t ifacet_type = 1;
	  	  
	  wmesh_int_t neicell_idx 	= c2celm[ifacet];
	  wmesh_int_t neicell_type	= c2celm_types[ifacet];
	  wmesh_int_t neicell_facet_idx	= -1;
	  
	  //
	  // Need to know jfacet.
	  //	  
	  wmesh_int_t neicell_num_facets = 0;
	  status = wmesh_get_c2celm(self_,
				    neicell_type,
				    neicell_idx,
				    &neicell_num_facets,
				    neicell_nei_ids,
				    1,
				    neicell_nei_types,
				    1);
	  WMESH_STATUS_CHECK(status);
	  for (wmesh_int_t l=0;l<neicell_num_facets;++l)
	    {
	      if (neicell_nei_types[l]>=0)
		{
		  if (neicell_nei_ids[l] == ielm)
		    {
		      neicell_facet_idx = l;
		      break;
		    }
		}
	    }
	  if (neicell_facet_idx < 0)
	    {
	      std::cout << "inconsistent c2c"  << std::endl;
	    }
	  
	  
	  //
	  // Get connectivity of the face.
	  //
	  wmesh_int_t c2n_facet[4],c2n_facet_inc=1;
	  status = wmesh_get_c2nfacet(self_,
				      itype_,
				      ielm,
				      ifacet,
				      c2n_facet,
				      c2n_facet_inc);
	  WMESH_STATUS_CHECK(status);
	  

	  wmesh_int_t nei_c2n_facet[4],nei_c2n_facet_inc=1;
	  status = wmesh_get_c2nfacet(self_,
				      jtype,
				      jelm,
				      jfacet,
				      nei_c2n_facet,
				      nei_c2n_facet_inc);
	  WMESH_STATUS_CHECK(status);


	  //
	  //
	  //
	  wmesh_int_t signed_rotation = 0;
	  {
	    wmesh_int_t min = c2n_facet[0];
	    for (wmesh_int_t l=1;l<;++l)
	      {
		if (min > c2n_facet[l])
		  {
		    min = c2n_facet[l];
		    signed_rotation = l;
		  }
	      }
	  }
	  
	  wmesh_int_t nei_signed_rotation = 0;
	  {
	    wmesh_int_t min = nei_c2n_facet[0];
	    for (wmesh_int_t l=1;l<;++l)
	      {
		if (min > nei_c2n_facet[l])
		  {
		    min = nei_c2n_facet[l];
		    nei_signed_rotation = l;
		  }
	      }
	  }


	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_element 	= &fem_->m_shape_eval_boundary_element[element];
	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_test 	= &fem_->m_shape_eval_boundary_test[element];
	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_f 	= &fem_->m_shape_eval_boundary_f[element];
	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_f_nei 	= &fem_->m_shape_eval_boundary_f[nei_element];
	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_u 	= &fem_->m_shape_eval_boundary_f[element];

	  const wmesh_cubature_t<T> * 	cubature_element 	= &eval_boundary_element->m_cubature_facets[facet_type][0][signed_rotation];
	  const wmesh_shape_eval_t<T> * eval_element 		= &eval_boundary_element->m_shape_eval_facets[facet_type][0][signed_rotation];
	  const wmesh_shape_eval_t<T> * eval_test		= &eval_boundary_test->m_shape_eval_facets[facet_type][0][signed_rotation];
	  const wmesh_shape_eval_t<T> * eval_f			= &eval_boundary_f->m_shape_eval_facets[facet_type][0][signed_rotation];
	  const wmesh_shape_eval_t<T> * eval_f_nei		= &eval_boundary_f->m_shape_eval_facets[nei_facet_type][1][nei_signed_rotation];
	  
	  //
	  // Compute the normals.
	  //
	  for (wmesh_int_t idim=0;idim<topodim;++idim)
	    {
	      BLAS_dgemm("N",
			 "N",
			 &topodim,
			 &eval_element->m_diff[idim].n,
			 &eval_element->m_diff[idim].m,
			 &r1,
			 cooelm_,
			 &cooelm_ld_,
			 eval_element->m_diff[idim].v,
			 &eval_element->m_diff[idim].ld,
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
	  // Compute the reference normal.
	  //	  
	  for (wmesh_int_t j=0;j<q_n;++j)
	    {
	      BLAS_dgemm("T",
			 "N",
			 &topodim,
			 &n1,
			 &topodim,
			 &r1,
			 jacobian,
			 &topodim,
			 reference_normal,
			 &topodim,
			 &r0,
			 physical_normal,
			 &topodim);
	    }


	  //
	  // Evaluate u.
	  //

	  //
	  // u.dot(n).
	  //	  
	  for (wmesh_int_t j=0;j<q_n;++j)
	    {
	      T dot = static_cast<T>(0);
	      for (wmesh_int_t i=0;i<topodim;++i)
		{
		  dot += n[i]*u[i];
		}
	      dot = ( dot < 0.0 ) ? -dot : 0.0;
	      dot * qfacet_w[j] * dets[j];
	    }

  

  //  std::cout << "    evaluate velocity ....  " << std::endl;
  
  //
  // Evaluate velocity.
  //
  T * ucubature = (T *)malloc(sizeof(T) *  uelm_m_ * q_n);
  dgemm("N",
	"N",
	&uelm_m_,
	&q_n,
	&uelm_n_,
	&r1,
	uelm_,
	&uelm_ld_,
	shape_eval_u->m_f.v,
	&shape_eval_u->m_f.ld,	
	&r0,
	ucubature,
	&uelm_m_);
  
  //  std::cout << "    evaluate velocity done.  " << std::endl;
  //
  // Compute the local matrix.
  //
  T lnabla_phi_j[3];
  T nabla_phi_j[3];
  static constexpr wmesh_int_t n1 = static_cast<wmesh_int_t>(1);
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      //  std::cout << "    integration point  " << k << "/" << q_n << std::endl;
      T*jacobian = jacobians + topodimXtopodim * k;
      for (wmesh_int_t j=0;j<lmat_n_;++j)
	{
	  //
	  // Compute the global gradient of shape function j at point k.
	  //
	  for (wmesh_int_t idim=0;idim<topodim;++idim)
	    {
	      const wmesh_int_t i 	= element_facet - ( (topodim_boundary==2) ? 2 : 1);
	      const wmesh_int_t lidx 	= (signed_rotation < 0) ? (facet_num_nodes + (-signed_rotation - 1)) : (signed_rotation - 1);
	      status = wmesh_shape_eval_def(&self_->m_facets_shape_eval[i][lidx],
					    element_facet,				
					    shape_family_,
					    shape_degree_,
					    self_->m_facets_nodes_storage,
					    &self_->m_facets_nodes[i][lidx]);	      
	      WMESH_STATUS_CHECK(status);
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
		  dot += ucubature[k*topodim + idim] * nabla_phi_j[idim];
		}
	      lmat_[lmat_ld_*j+i] += ( dot  * shape_eval_test->m_f.v[shape_eval_test->m_f.ld * k + i]  ) *  q_w[k]  *  dets[k];
	    }
	  
	}      
    }

  return WMESH_STATUS_SUCCESS;
}


template<typename T>
struct wmesh_facet_restrict_t
{
  wmesh_int_t 		m_shape_family;
  wmesh_int_t 		m_shape_degree;
  wmesh_int_t 		m_nodes_family;
  
  wmesh_int_t 		m_shape_family_over_facet;
  wmesh_int_t 		m_shape_degree_over_facet;
  wmesh_int_t 		m_nodes_family_over_facet;

  wmesh_mat_t<T>   	m_nodes_over_facet[WMESH_ELEMENT_ALL][2];
  wmesh_mat_t<T> 	m_restrict[WMESH_ELEMENT_ALL][6];
};

  //
  //
  //

  wmesh_shape_eval_boundary_t<T>	m_shape_eval_boundary_element[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_boundary_t<T>	m_shape_eval_boundary_f[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_boundary_t<T> 	m_shape_eval_boundary_test[WMESH_ELEMENT_ALL];
  wmesh_shape_eval_boundary_t<T> 	m_shape_eval_boundary_u[WMESH_ELEMENT_ALL];

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
#if 0
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
  status = wmesh_shape_eval_boundary_def(&self_->m_shape_eval_boundary_test[element_],
					 element_,
					 shape_test_family_,
					 shape_test_degree_,
					 cubature_family_,				
					 cubature_degree_);
  WMESH_STATUS_CHECK(status);      
#endif
#if 0
  
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



#if 0
template<typename T>
wmesh_status_t wmesh_facet_restrict_def(wmesh_facet_restrict_t * self_,
					wmesh_int_t 		shape_family_,
					wmesh_int_t 		shape_degree_,
					wmesh_int_t 		nodes_family_,

					wmesh_int_t 		shape_family_,
					wmesh_int_t 		shape_degree_,
					wmesh_int_t 		nodes_family_)
{
  
  //
  // Build nodes over facets.
  //
  
  //
  // Build nodal evaluation of the shape over the element.
  //  
}
#endif






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
wmesh_status_t wmeshspacedg_fem_advection_diagonal_block_buffer_size(const wmeshspacedg_fem_advection_t<T> *__restrict__	self_,
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


//
// Only diagonal part.
//
template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_diagonal_block(const wmeshspacedg_fem_advection_t<T> *__restrict__	self_,
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
						       wmesh_int_t 					rw_n_,
						       T*__restrict__					rw_)
{
  static constexpr T r0 = static_cast<T>(0);
  static constexpr T r1 = static_cast<T>(1);

  WMESH_CHECK_POINTER(self_);
  
  wmesh_status_t status;
  wmesh_int_t required_rw_n;
  status = wmeshspacedg_fem_advection_diagonal_block_buffer_size(self_,
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

  //  const wmesh_int_t 				u_n 			= (uelm_storage_ == WMESH_STORAGE_INTERLEAVE) ? uelm_n_ : uelm_m_;
  const wmesh_int_t 				q_n 			= cubature->m_w.n;
  const T* __restrict__				q_w 			= cubature->m_w.v;  
  const wmesh_int_t 	topodim  					= cooelm_m_;
  const wmesh_int_t  	topodimXtopodim 				= topodim*topodim;;
  
  //
  // Compute Jacobians and determinants.
  //
  std::cout << "    compute jacobians  " << std::endl;




  
  T *__restrict__ jacobians 	= rw_;
  T *__restrict__ dets 		= rw_ + topodimXtopodim * q_n;
  rw_ += required_rw_n;
  rw_n_ -= required_rw_n;

  
  for (wmesh_int_t idim=0;idim<topodim;++idim)
    {
      //      std::cout << "? " << shape_eval_element->m_diff[idim].m << std::endl;
      //      std::cout << "? " << shape_eval_element->m_diff[idim].ld << std::endl;
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

  std::cout << "    compute jacobians done.  " << std::endl;

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


  //  std::cout << "    evaluate velocity ....  " << std::endl;
  
  //
  // Evaluate velocity.
  //
  T * ucubature = (T *)malloc(sizeof(T) *  uelm_m_ * q_n);
  dgemm("N",
	"N",
	&uelm_m_,
	&q_n,
	&uelm_n_,
	&r1,
	uelm_,
	&uelm_ld_,
	shape_eval_u->m_f.v,
	&shape_eval_u->m_f.ld,	
	&r0,
	ucubature,
	&uelm_m_);
  
  //  std::cout << "    evaluate velocity done.  " << std::endl;
  //
  // Compute the local matrix.
  //
  T lnabla_phi_j[3];
  T nabla_phi_j[3];
  static constexpr wmesh_int_t n1 = static_cast<wmesh_int_t>(1);
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      //  std::cout << "    integration point  " << k << "/" << q_n << std::endl;
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
		  dot += ucubature[k*topodim + idim] * nabla_phi_j[idim];
		}
	      lmat_[lmat_ld_*j+i] += ( dot  * shape_eval_test->m_f.v[shape_eval_test->m_f.ld * k + i]  ) *  q_w[k]  *  dets[k];
	    }
	  
	}      
    }


  //
  //
  //
  
  

  
#if 0
  for (wmesh_int_t i=0;i<lmat_m_;++i)
    {
      lrhs_[lrhs_inc_*i] = r0;	      
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
wmesh_status_t bms_restrict_operator(wmesh_int_t 		element_,
				     wmesh_int_t 		shape_family_,
				     wmesh_int_t 		shape_degree_,
				     wmesh_int_t 		facet_local_idx_,
				     wmesh_int_t 		facet_,
				     wmesh_int_t 		facet_type_,
				     wmesh_int_t 		facet_nodes_storage_,
				     const wmesh_mat_t<T>& 	facet_nodes_,
				     const wmesh_mat_t<T>& 	facet_eval_to_element_,
				     wmesh_int_t 		mat_storage_,
				     wmesh_mat_t<T>& 		mat_,
				     wmesh_mat_t<T>& 		facet_nodes_into_element_)
{
  

  
  wmesh_int_t iw_n,rw_n;
  wmesh_status_t status;

  //
  // Facet nodes is has a topological dimension equal to (topodim-1).
  //

  //
  // Get reference coordinates.
  //
  T ref_element_coordinates[32];
  T ref_facet_coordinates[32];
  status = bms_ordering_vertices(element_,
				 ref_element_coordinates);
  WMESH_STATUS_CHECK(status);

  //
  // Get facet coordinates.
  //
  
  //
  // Transformation (1d-2d) (2d-3d)
  //
  wmesh_mat_gemm(static_cast<T>(1.0),
		 ref_facet_coordinates,
		 facet_eval_to_element_,
		 static_cast<T>(0.0),
		 facet_nodes_into_element_);


  

  
  //
  // Now eval the shape over this.
  //
  wmesh_shape_eval_t<T> shape_eval;
  wmesh_shape_eval_def(&shape_eval,
		       element_,
		       shape_family_,
		       shape_degree_,
		       facet_nodes_storage_,
		       facet_nodes_into_element_);

  wmesh_mat_copy(mat_,
		 shape_eval.m_f);

  return WMESH_STATUS_SUCCESS;  
    
}

template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_extra_diagonal_block2(const wmeshspacedg_fem_advection_t<T> *__restrict__	self_,
								wmesh_int_t 						facet_type_,
								
								wmesh_int_t 						cooelm_storage_,
								const wmesh_mat_t<T>&					cooelm_,

								wmesh_int_t 						cooelm_restrict_op_storage_,
								const wmesh_mat_t<T>&					cooelm_restrict_op_,

								wmesh_int_t 						coofacet_eval_storage_,
								const wmesh_mat_t<T>&					coofacet_eval_,
								
								wmesh_int_t 						uelm_storage_,
								const wmesh_mat_t<T>&					uelm_,
								wmesh_int_t 						uelm_restrict_op_storage_,
								const wmesh_mat_t<T>&					uelm_restrict_op_,
								
								wmesh_int_t 						ufacet_eval_storage_,
								const wmesh_mat_t<T>&					ufacet_eval_,
								
								wmesh_int_t 						phij_restrict_op_storage_,
								const wmesh_mat_t<T>&					phij_restrict_op_,

								wmesh_int_t 						phij_eval_storage_,
								const wmesh_mat_t<T>&					phij_eval_,

								wmesh_int_t 						phii_restrict_op_storage_,
								const wmesh_mat_t<T>&					phii_restrict_op_,

								wmesh_int_t 						phii_eval_storage_,
								const wmesh_mat_t<T>&					phii_eval_,

								wmesh_mat_t<T>&						lmat_,
								
								wmesh_int_t 						ufacet_storage_,
								wmesh_mat_t<T>&						ufacet_,
								
								wmesh_int_t 						coofacet_storage_,
								wmesh_mat_t<T>&						coofacet_,
								
								const wmesh_cubature_t<T>& 				cubature_facet_,
								
								wmesh_int_t  						cubature_velocity_storage_,
								wmesh_mat_t<T>& 					cubature_velocity_,
								
								wmesh_int_t 						rw_n_,
								T*__restrict__						rw_)
{

  //
  // BLOCK nx6         BLOCK 6xdim = n * dim
  // uelm_restrict_op_ uelm_  
  //
  // ILV 6xn           ILV
  // uelm_             uelm_restrict_op_  
  //
  static constexpr T r0 = static_cast<T>(0);
  static constexpr T r1 = static_cast<T>(1);

  
  //
  // Extract facet coordinates.
  //
  wmesh_mat_gemm(r1,
		 cooelm_,
		 cooelm_restrict_op_,
		 r0,
		 coofacet_);

  
  //
  // Extract velocity on facets.
  //
  wmesh_mat_gemm(r1,
		 uelm_,
		 uelm_restrict_op_,
		 r0,
		 ufacet_);

  
  //
  // Interpolate the normal.
  //
  
  
  //
  // Extract normal on facets.
  //
  wmesh_mat_gemm(r1,
		 uelm_,
		 uelm_restrict_op_,
		 r0,
		 ufacet_);


  
  //
  // Extract coordinates, the result has a block storage.
  //
  // coofacet = R_c * cooelm 
  //
  {
    const char * transa = (WMESH_STORAGE_BLOCK == cooelm_restrict_op_storage_) ? "N" : "T";
    const char * transb = (WMESH_STORAGE_BLOCK == cooelm_storage_) ? "N" : "T";
    const wmesh_int_t m = (WMESH_STORAGE_INTERLEAVE == cooelm_restrict_op_storage_) ? cooelm_restrict_op_.n : cooelm_restrict_op_.m;
    const wmesh_int_t n = (WMESH_STORAGE_INTERLEAVE == cooelm_storage_) ? cooelm_.m : cooelm_.n;
    const wmesh_int_t k = (WMESH_STORAGE_INTERLEAVE == cooelm_restrict_op_storage_) ? cooelm_restrict_op_.m : cooelm_restrict_op_.n;

    
    xgemm(transa,
	  transb,
	  &m,
	  &n,
	  &k,
	  &r1,
	  cooelm_restrict_op_.v,
	  &cooelm_restrict_op_.ld,
	  cooelm_.v,
	  &cooelm_.ld,
	  &r0,
	  coofacet_.v,
	  &coofacet_.ld);
    
  }
  
  //
  // Extract velocity, the result has a block storage.
  //
  // ufacet = R_u * uelm 
  //
  {
    const char * transa = (WMESH_STORAGE_BLOCK == uelm_restrict_op_storage_) ? "N" : "T";
    const char * transb = (WMESH_STORAGE_BLOCK == uelm_storage_) ? "N" : "T";
    const wmesh_int_t m = (WMESH_STORAGE_INTERLEAVE == uelm_restrict_op_storage_) ? uelm_restrict_op_.n : uelm_restrict_op_.m;
    const wmesh_int_t n = (WMESH_STORAGE_INTERLEAVE == uelm_storage_) ? uelm_.m : uelm_.n;
    const wmesh_int_t k = (WMESH_STORAGE_INTERLEAVE == uelm_restrict_op_storage_) ? uelm_restrict_op_.m : uelm_restrict_op_.n;
    
    xgemm(transa,
	  transb,
	  &m,
	  &n,
	  &k,
	  &r1,
	  uelm_restrict_op_.v,
	  &uelm_restrict_op_.ld,
	  uelm_.v,
	  &uelm_.ld,
	  &r0,
	  ufacet_.v,
	  &ufacet_.ld);
  }

  //
  // dx/dr dx/dr * dr / ds
  // dy/dr dy/ds * ds / dr
  //
  // u -> (r= u,s = 0)
  // u -> (r= 1-u,s = u)

  //
  // (u,v) -> (r= 1-u-v,s=u,t=v)
  //
  // dr/du = -1 -1     ( 2 -1  )/3  (-1 1 0)
  //          1  0       -1 2
  //          0  1
  // u -> (r= 1-u,s = u)
  // 
  // dx/dr = dx/dt * dt/dr
  // 
  // X = (1-t)P0 + t * P1   
  // dx/dt = P1-P0 (1;0) 2x1 (1,0) * (1;0) = 1
  // dx/dt = P1-P0 (1;-1) 2x1 (1,-1) * (1;-1) = 2 
  // (1;-1) (1,-1) (1,-1
  //                -1,1)
  // dx/dt dt/dr
  // Evaluate the normals
  //
  // 
  //
  //
  //
  //
  
  {
    const char * transa = (WMESH_STORAGE_INTERLEAVE == coofacet_eval_storage_) ? "T" : "N";
    const char * transb = "N";
    const wmesh_int_t m = (WMESH_STORAGE_INTERLEAVE == coofacet_eval_storage_) ? coofacet_eval_.n : coofacet_eval_.m;
    const wmesh_int_t n = (WMESH_STORAGE_INTERLEAVE == coofacet_storage_) ? coofacet_.m : coofacet_.n;
    const wmesh_int_t k = (WMESH_STORAGE_INTERLEAVE == coofacet_eval_storage_) ? coofacet_eval_.m : coofacet_eval_.n;
    
    xgemm(transa,
	  transb,
	  &m,
	  &n,
	  &k,
	  &r1,
	  coofacet_eval_.v,
	  &coofacet_eval_.ld,
	  coofacet_.v,
	  &coofacet_.ld,
	  &r0,
	  cubature_facet_.v,
	  &cubature_facet_.ld);    
  }
  
  

  //
  // Evaluate physical coordinates at the cubature points.
  // Maybe not needed.
  //
  {
    const char * transa = (WMESH_STORAGE_INTERLEAVE == coofacet_eval_storage_) ? "T" : "N";
    const char * transb = "N";
    const wmesh_int_t m = (WMESH_STORAGE_INTERLEAVE == coofacet_eval_storage_) ? coofacet_eval_.n : coofacet_eval_.m;
    const wmesh_int_t n = (WMESH_STORAGE_INTERLEAVE == coofacet_storage_) ? coofacet_.m : coofacet_.n;
    const wmesh_int_t k = (WMESH_STORAGE_INTERLEAVE == coofacet_eval_storage_) ? coofacet_eval_.m : coofacet_eval_.n;
    
    xgemm(transa,
	  transb,
	  &m,
	  &n,
	  &k,
	  &r1,
	  coofacet_eval_.v,
	  &coofacet_eval_.ld,
	  coofacet_.v,
	  &coofacet_.ld,
	  &r0,
	  cubature_facet_.v,
	  &cubature_facet_.ld);    
  }
  
  

  //
  // Evaluate velocity at the cubature points.
  //
  {
    const char * transa = (WMESH_STORAGE_INTERLEAVE == ufacet_eval_storage_) ? "T" : "N";
    const char * transb = "N";
    const wmesh_int_t m = (WMESH_STORAGE_INTERLEAVE == ufacet_eval_storage_) ? ufacet_eval_.n : ufacet_eval_.m;
    const wmesh_int_t n = (WMESH_STORAGE_INTERLEAVE == ufacet_storage_) ? ufacet_.m : ufacet_.n;
    const wmesh_int_t k = (WMESH_STORAGE_INTERLEAVE == ufacet_eval_storage_) ? ufacet_eval_.m : ufacet_eval_.n;
    
    xgemm(transa,
	  transb,
	  &m,
	  &n,
	  &k,
	  &r1,
	  ufacet_eval_.v,
	  &ufacet_eval_.ld,
	  ufacet_.v,
	  &ufacet_.ld,
	  &r0,
	  cubature_velocity_.v,
	  &cubature_velocity_.ld);    
  }

  
  

  //
  // Now we operate in 2d.
  //

  

  //
  //
  //
  // phi_j_restrict_op_ (9x27)  or (16x27)
  // phi_i_restrict_op_ (9x27) 
  // G_j = (9x100) 
  // 27 x 27 = phi_j^T * phi_i * (u.n) * q_ 
  // G_j^{ILV ? T} phi_j^{ILV?T} 100x9 9x27 100x27
  // G_j^{ILV ? T} phi_j^{ILV?T} |u.n|
  //

  

  //
  // Compute coordinates.
  //
  
  //
  // Interpolate the normal.
  //
  {
    const char * transa = (WMESH_STORAGE_BLOCK == cooelm_restrict_op_storage_) ? "N" : "T";
    const char * transb = (WMESH_STORAGE_BLOCK == cooelm_storage_) ? "N" : "T";
    const wmesh_int_t m = (WMESH_STORAGE_INTERLEAVE == cooelm_restrict_op_storage_) ? cooelm_restrict_op_.n : cooelm_restrict_op_.m;
    const wmesh_int_t n = (WMESH_STORAGE_INTERLEAVE == cooelm_storage_) ? cooelm_.m : cooelm_.n;
    const wmesh_int_t k = (WMESH_STORAGE_INTERLEAVE == cooelm_restrict_op_storage_) ? cooelm_restrict_op_.m : cooelm_restrict_op_.n;    
    xgemm(transa,
	  transb,
	  &m,
	  &n,
	  &k,
	  &r1,
	  cooelm_restrict_op_.v,
	  &cooelm_restrict_op_.ld,
	  cooelm_.v,
	  &cooelm_.ld,
	  &r0,
	  coofacet_.v,
	  &coofacet_.ld);
    
  }

  
  
#if 0
  //
  // Get the restriction operator for uelm_
  //
  const wmesh_mat_t<T>& h = 


      //
      // Restrict the coordinates of the element.
      //
      dgemm("N",
	    "N",
	    &cooelm_m_,
	    &facet_restrict_element->m_restrict[element_][ifacet].n,
	    &cooelm_n_,
	    &r1,
	    cooelm_,
	    &cooelm_ld_,
	    facet_restrict_element->m_restrict[element_][ifacet].v,
	    &facet_restrict_element->m_restrict[element_][ifacet].ld,	
	    &r0,
	    efacet,
	    &uelm_m_);

  
  
  static constexpr T r0 = static_cast<T>(0);
  static constexpr T r1 = static_cast<T>(1);

  WMESH_CHECK_POINTER(self_);
  
  wmesh_status_t status;
  wmesh_int_t required_rw_n,element_=ielement_;
  status = wmeshspacedg_fem_advection_diagonal_block_buffer_size(self_,
								 element_,
								 &required_rw_n);
  WMESH_STATUS_CHECK(status);
  
  if (rw_n_ < required_rw_n)
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_WORKSPACE);       
    }


  

  //
  // Get the face coordinates.
  //

  //
  // Get the facet nodes
  //
  
  //
  // Get the restriction of the basis over the face.
  //

  const wmesh_shape_eval_boundary_t<T> * eval_boundary_element 	= &self_->m_shape_eval_boundary_element[element_];
  const wmesh_shape_eval_boundary_t<T> * eval_boundary_test 	= &self_->m_shape_eval_boundary_test[element_];
  const wmesh_shape_eval_boundary_t<T> * eval_boundary_f 	= &self_->m_shape_eval_boundary_f[element_];
  const wmesh_shape_eval_boundary_t<T> * eval_boundary_u 	= &self_->m_shape_eval_boundary_f[element_];
  //  const wmesh_shape_eval_boundary_t<T> * eval_boundary_f_nei 	= &self_->m_shape_eval_boundary_f[jelement_];
  wmesh_int_t signed_rotation = 1;
  const wmesh_cubature_t<T> * 	cubature_facet 			= &eval_boundary_element->m_cubature_facets[facet_type_][0][signed_rotation];
  const wmesh_shape_eval_t<T> * eval_facet 			= &eval_boundary_element->m_shape_eval_facets[facet_type_][0][signed_rotation];
  const wmesh_shape_eval_t<T> * eval_test			= &eval_boundary_test->m_shape_eval_facets[facet_type_][0][signed_rotation];
  const wmesh_shape_eval_t<T> * eval_f				= &eval_boundary_f->m_shape_eval_facets[facet_type_][0][signed_rotation];
  //  const wmesh_shape_eval_t<T> * eval_f_nei			= &eval_boundary_f->m_shape_eval_facets[facet_type][1][nei_signed_rotation];

#if 0
  const wmesh_shape_eval_t<T> * __restrict__ 	shape_eval_element 	= &self_->m_shape_eval_element[element_];
  const wmesh_shape_eval_t<T> * __restrict__	shape_eval_f 		= &self_->m_shape_eval_f[element_];
  const wmesh_shape_eval_t<T> * __restrict__	shape_eval_u 		= &self_->m_shape_eval_u[element_];
  const wmesh_shape_eval_t<T> * __restrict__	shape_eval_test 	= &self_->m_shape_eval_test[element_];  
#endif
  const wmesh_cubature_t<T>*__restrict__ 	cubature		= &self_->m_cubature[element_];
  //  const wmesh_int_t 				u_n 			= (uelm_storage_ == WMESH_STORAGE_INTERLEAVE) ? uelm_n_ : uelm_m_;
  const wmesh_int_t 				q_n 			= cubature->m_w.n;
  const T* __restrict__				q_w 			= cubature->m_w.v;  
  wmesh_int_t 	topodim;
  status = bms_element2topodim(element_,&topodim);
  const wmesh_int_t  	topodimXtopodim 				= topodim*topodim;;
  
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
  
  //  wmesh_facet_restrict_t<T> facet_restrict;
  wmesh_mat_t<T> 	m_restrict[WMESH_ELEMENT_ALL][6];

  T ufacet[];
  for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
    {

      //
      // Restrict the coordinates of the element.
      //
      dgemm("N",
	    "N",
	    &cooelm_m_,
	    &facet_restrict_element->m_restrict[element_][ifacet].n,
	    &cooelm_n_,
	    &r1,
	    cooelm_,
	    &cooelm_ld_,
	    facet_restrict_element->m_restrict[element_][ifacet].v,
	    &facet_restrict_element->m_restrict[element_][ifacet].ld,	
	    &r0,
	    efacet,
	    &uelm_m_);
      
      //
      // Restrict the velocity.
      // u = dim x u_ndofs
      // restrict = u_ndofs x facet_ndofs
      //
      dgemm("N",
	    "N",
	    &uelm_m_,
	    &facet_restrict_u->m_restrict[element_][ifacet].n,
	    &uelm_n_,
	    &r1,
	    uelm_,
	    &uelm_ld_,
	    facet_restrict_u->m_restrict[element_][ifacet].v,
	    &facet_restrict_u->m_restrict[element_][ifacet].ld,	
	    &r0,
	    ufacet,
	    &uelm_m_);
      
      //
      // Restrict the test functions: transpose of the current restrict.
      // restrict: (ngauss_facet x facet_ndofs)_{facet} ^T *  (test_ndofs x facet_ndofs)_element^T
      //
      //
      // Restrict the f functions: transpose of the current restrict.
      // restrict facet_ndofs x test_ndofs
      //
      for (wmesh_int_t k=0;k<q_n;++k)
	{	  
	  xger(&lmat_m_,&lmat_n_,
	       alpha,
	       phi_j + phi_j_ld * k,
	       1,
	       phi_i + phi_i_ld * k,
	       1,
	       lmat_,
	       lmat_ld_);
	}
      
      //
      // Evaluate r-velocity.
      //

      //
      // Evaluate normal.
      //

      //
      // Evaluate determinant.
      //
      
      //
      // Evaluate r-test-i.
      //

      //
      // Evaluate r-test-i.
      //
      dgemm("T",
	    "N",
	    &facet_restrict_test->m_restrict[element_][ifacet].m,
	    &uelm_m_,
	    &facet_restrict_test->m_restrict[element_][ifacet].n,
	    &uelm_n_,
	    &r1,
	    uelm_,
	    &uelm_ld_,
	    facet_restrict_f->m_restrict[element_][ifacet].v,
	    &facet_restrict_f->m_restrict[element_][ifacet].ld,	
	    facet_restrict_test->m_restrict[element_][ifacet].v,
	    &facet_restrict_test->m_restrict[element_][ifacet].ld,	
	    &r0,
	    fufacet,
	    &uelm_m_);


      //
      // Now evaluate the velocity over the facet from its restriction.
      //
      
      //
      // Now evaluate the velocity over the facet from its restriction.
      //

      //
      // Interpolate the non-normalized normal.
      //
      
      
    }
  
  //
  // Restrict the test functions.
  //
  
  
  //  std::cout << "    evaluate velocity ....  " << std::endl;
  
  //
  // Evaluate velocity.
  //
  T * ucubature = (T *)malloc(sizeof(T) *  uelm_m_ * q_n);
  dgemm("N",
	"N",
	&uelm_m_,
	&q_n,
	&uelm_n_,
	&r1,
	uelm_,
	&uelm_ld_,
	eval_u->m_f.v,
	&eval_u->m_f.ld,	
	&r0,
	ucubature,
	&uelm_m_);
  
  //  std::cout << "    evaluate velocity done.  " << std::endl;
  //
  // Compute the local matrix.
  //
  T lnabla_phi_j[3];
  T nabla_phi_j[3];
  static constexpr wmesh_int_t n1 = static_cast<wmesh_int_t>(1);
  for (wmesh_int_t k=0;k<q_n;++k)
    {
      //  std::cout << "    integration point  " << k << "/" << q_n << std::endl;
      T*jacobian = jacobians + topodimXtopodim * k;
      for (wmesh_int_t j=0;j<lmat_n_;++j)
	{
	  //
	  // Compute the global gradient of shape function j at point k.
	  //
	  for (wmesh_int_t idim=0;idim<topodim;++idim)
	    {
	      lnabla_phi_j[idim] = eval_f->m_diff[idim].v[eval_f->m_diff[idim].ld * k + j];
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
		  dot += ucubature[k*topodim + idim] * nabla_phi_j[idim];
		}
	      lmat_[lmat_ld_*j+i] += ( dot  * eval_test->m_f.v[eval_test->m_f.ld * k + i]  ) *  q_w[k]  *  dets[k];
	    }
	  
	}      
    }


  //
  //
  //
  
  

  
#if 0
  for (wmesh_int_t i=0;i<lmat_m_;++i)
    {
      lrhs_[lrhs_inc_*i] = r0;	      
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
#endif
  return WMESH_STATUS_SUCCESS;
}
#endif





template <typename T>
wmesh_status_t wmeshspace_get_elmdofs(const wmeshspace_t * __restrict__	self_,
				      wmesh_int_t			itype_,
				      wmesh_int_t			ielm_,
				      wmesh_int_t			dofs_storage_,
				      wmesh_int_t			dofs_m_,
				      wmesh_int_t			dofs_n_,
				      const T * 			dofs_,		
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

template <typename T>
wmesh_status_t wmesh_get_cooelm(const wmesh_t * __restrict__	self_,
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





wmesh_status_t wmesh_get_c2nelm(const wmesh_t * __restrict__	self_,
				wmesh_int_t			itype_,
				wmesh_int_t			ielm_,
				wmesh_int_p			c2n_,
				wmesh_int_t			c2n_inc_)
{
  WMESH_CHECK_POINTER(self_);
  WMESH_CHECK_POINTER(c2n_);
  for (wmesh_int_t i=0;i<self_->m_c2n.m_m[itype_];++i)
    {
      c2n_[c2n_inc_*i] = self_->m_c2n.m_data[self_->m_c2n.m_ptr[itype_] + self_->m_c2n.m_ld[itype_] * ielm_ + i] - 1;
    }
  return WMESH_STATUS_SUCCESS;  
}

wmesh_status_t wmesh_get_c2celm(const wmesh_t * __restrict__	self_,
				wmesh_int_t			itype_,
				wmesh_int_t			ielm_,
				wmesh_int_p			num_facets_,
				wmesh_int_p			c2c_,
				wmesh_int_t			c2c_inc_,
				wmesh_int_p			c2c_types_,
				wmesh_int_t			c2c_types_inc_)
{
  WMESH_CHECK_POINTER(self_);
  WMESH_CHECK_POINTER(c2c_);
  WMESH_CHECK_POINTER(c2c_types_);
  const wmesh_int_t num_facets = self_->m_c2c.m_m[itype_];
  for (wmesh_int_t i=0;i<num_facets;++i)
    {
      auto c2c_info = self_->m_c2c.m_data[self_->m_c2c.m_ptr[itype_] + self_->m_c2c.m_ld[itype_] * ielm_ + i];
      if (c2c_info != 0)
	{
	  wmesh_status_t status = bms_c2c_cindex(c2c_info,&c2c_[c2c_inc_*i]);
	  WMESH_STATUS_CHECK(status);
	  status = bms_c2c_ctype(c2c_info,&c2c_types_[c2c_types_inc_*i]);
	  WMESH_STATUS_CHECK(status);
	}
      else
	{
	  c2c_[c2c_inc_*i] = 0;
	  c2c_types_[c2c_types_inc_*i] = -1;
	}
    }
  
  num_facets_[0] = num_facets;
  return WMESH_STATUS_SUCCESS;  
}


wmesh_status_t bms_dg_size_blocks(wmesh_int_t 			degree_,
				  wmesh_int_p 			size_blocks_)
{
  wmesh_status_t status;
  for (wmesh_int_t element=WMESH_ELEMENT_EDGE;element<WMESH_ELEMENT_ALL;++element)
    {
      status = bms_ndofs(element,
			 degree_,
			 &size_blocks_[element]);
      WMESH_STATUS_CHECK(status);
    }
  return WMESH_STATUS_SUCCESS;
}



template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_global_system_zone(const wmeshspacedg_t *			self_,
							     wmesh_int_t				itype_,
							     const wmeshspacedg_fem_advection_t<T> *	fem_,
							     const wmeshspace_t *			velocity_space_,
							     wmesh_int_t 				velocity_storage_,
							     wmesh_int_t 				velocity_m_,
							     wmesh_int_t 				velocity_n_,
							     const double * 				velocity_,
							     wmesh_int_t 				velocity_ld_,
							     wmesh_int_t				csr_size_,
							     const_wmesh_int_p				csr_ptr_,
							     const_wmesh_int_p				csr_ind_,
							     T * 					csr_val_,
							     T * 					rhs_)
{
  const wmesh_t *				mesh = self_->m_mesh;

  wmesh_status_t status;
  const wmesh_int_t topodim 			= mesh->m_topology_dimension;
  const wmesh_int_t element 			= itype_ + ( (topodim==3) ? 4 : ((topodim==2)?2:1) );

  const wmesh_int_t dg_ndofs 		= self_->m_dofs_m[itype_];
  const wmesh_int_t test_ndofs 		= dg_ndofs;

  //
  // Alloc cooelm.
  //
  const wmesh_int_t cooelm_storage = WMESH_STORAGE_INTERLEAVE;
  const wmesh_int_t uelm_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>
    lrhs,
    lmat,
    cooelm,
    uelm;

  const wmesh_int_t elm_ndofs 	= mesh->m_c2n.m_m[itype_];
  const wmesh_int_t u_ndofs 	= velocity_space_->m_c2d.m_m[itype_];
  
  wmesh_mat_t<T>::alloc(&cooelm,topodim,elm_ndofs);
  wmesh_mat_t<T>::alloc(&uelm,topodim,u_ndofs);
  wmesh_mat_t<T>::alloc(&lmat,dg_ndofs,dg_ndofs);
  wmesh_mat_t<T>::alloc(&lrhs,1,dg_ndofs);
  
  const wmesh_int_t nei_cooelm_storage 		= WMESH_STORAGE_INTERLEAVE;
  const wmesh_int_t nei_cooelm_m 			= topodim;
  const wmesh_int_t nei_cooelm_n 			= mesh->m_c2n.m_m[itype_];
  const wmesh_int_t nei_cooelm_ld 			= nei_cooelm_m;
  T* 	nei_cooelm = (T*)malloc(sizeof(T)* nei_cooelm_m * nei_cooelm_n );
  
  wmesh_int_t rw_n;
  status = wmeshspacedg_fem_advection_diagonal_block_buffer_size(fem_,
							       element,
							       &rw_n);
  WMESH_STATUS_CHECK(status);
  T * rw =  (T*)malloc(sizeof(T)*rw_n);
  
  wmesh_int_t max_size_block = 0;
  wmesh_int_t size_blocks[4];
  wmesh_int_t num_types = mesh->m_c2n.m_size;
  for (wmesh_int_t i=0;i<num_types;++i)
    {
      wmesh_int_t element 	= ( (topodim==3) ? 4 : (topodim==2) ? 2 : 1 ) + i;
      size_blocks[i] 		= fem_->m_shape_eval_f[element].m_shape.m_ndofs;
      max_size_block 		= (max_size_block < size_blocks[i]) ? size_blocks[i] : max_size_block;
    }

  wmesh_int_t dofs_ptr[4 + 1];
  
  dofs_ptr[0] = 0;
  for (wmesh_int_t i=0;i<num_types;++i)
    {
      dofs_ptr[i+1] = dofs_ptr[i] + size_blocks[i] * mesh->m_c2n.m_n[i];
    }


  wmesh_int_p idofs = (wmesh_int_p)malloc(max_size_block * sizeof(wmesh_int_t));
  wmesh_int_p jdofs = (wmesh_int_p)malloc(max_size_block * sizeof(wmesh_int_t));
  
  const wmesh_int_t num_cells = mesh->m_c2n.m_n[itype_];
  wmesh_int_t c2nelm[8],c2nelm_inc=1;
  wmesh_int_t c2celm[6],c2celm_inc=1;
  wmesh_int_t c2celm_types[6],c2celm_types_inc=1;


#if 0
  template <typename T>
    wmesh_status_t bms_element_jacobians(wmesh_int_t 					element_,
					 wmesh_int_t 					cooelm_storage_,
					 const wmesh_mat_t<T>&				cooelm_,
					 wmesh_int_t 					element_shape_eval_diff_storage_,
					 const wmesh_mat_t<T>*				element_shape_eval_diff_,
					 wmesh_mat_t<T>&				jacobians_,
					 wmesh_mat_t<T>&				jacobians_det_);
#endif
  
  std::cout << " Calculation and Assembly over  " << num_cells << " cells." << std::endl;
  for (wmesh_int_t ielm=0; ielm < num_cells;++ielm)
    {
      std::cout << "loop  " << ielm << std::endl;
      std::cout << "    get_c2nelm  " << std::endl;
      
      //
      // Get connectivity of element.
      //
      status = wmesh_get_c2nelm(mesh,
				itype_,
				ielm,
				c2nelm,
				c2nelm_inc);
      WMESH_STATUS_CHECK(status);

      std::cout << "    get_c2celm  " << std::endl;
      
      //
      // Get neighbors information.
      //
      wmesh_int_t num_facets = 0;
      status = wmesh_get_c2celm(mesh,
				itype_,
				ielm,
				&num_facets,
				c2celm,
				c2celm_inc,
				c2celm_types,
				c2celm_types_inc);
      WMESH_STATUS_CHECK(status);
      std::cout << "    get_cooelm  " << std::endl;
      
      //
      // Get coordinates of element.
      //
      status = wmesh_get_cooelm(mesh,
				itype_,
				ielm,
				cooelm_storage,
				WMESH_MAT_FORWARD(cooelm));
      WMESH_STATUS_CHECK(status);

      std::cout << "    get_velocity dofs  " << std::endl;

      //
      // Get velocity dofs.
      //
      status = wmeshspace_get_elmdofs(velocity_space_,
				      itype_,
				      ielm,
				      velocity_storage_,
				      velocity_m_,
				      velocity_n_,
				      velocity_,
				      velocity_ld_,
				      uelm_storage,
				      WMESH_MAT_FORWARD(uelm));
      WMESH_STATUS_CHECK(status);
      std::cout << "    diagonal block  " << std::endl;
      
      //
      // Only diagonal.
      //
      status = wmeshspacedg_fem_advection_diagonal_block(fem_,
							 element,
							 
							 cooelm_storage,
							 WMESH_MAT_FORWARD(cooelm),
						       
							 uelm_storage,
							 WMESH_MAT_FORWARD(uelm),
							 
							 WMESH_MAT_FORWARD(lmat),

							 rw_n,
							 rw);
      WMESH_STATUS_CHECK(status);

      status = wmeshspacedg_fem_advection_diagonal_block_flux(fem_,
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


      status = wmeshspacedg_get_dofs(self_,itype_,ielm,idofs,1);
      WMESH_STATUS_CHECK(status);

      memcpy(jdofs,idofs,sizeof(wmesh_int_t)*size_blocks[itype_]);
      wmesh_int_t ni = size_blocks[itype_];
      wmesh_int_t nj = ni;

      WMESH_STATUS_CHECK(status);
      status =  bms_sparse_add(ni,
			       idofs,
			       1,
			       
			       nj,
			       jdofs,
			       1,
			       
			       lmat.v,
			       lmat.ld,
			       
			       csr_size_,
			       csr_ptr_,
			       csr_ind_,
			       csr_val_);
      
    }



  //
  // Initialize normals.
  //
  wmesh_int_t refnormals_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T> refnormals;
  wmesh_mat_t<T>::alloc(&refnormals,topodim,mesh->m_c2c.m_m[itype_]);
  memset(refnormals.v,0,sizeof(T)*refnormals.m*refnormals.n);
  
  std::cout << "refnormals " << std::endl;
  std::cout << refnormals << std::endl;
  status = bms_ordering_normals(element,
				refnormals_storage,
				WMESH_MAT_FORWARD(refnormals) );
  WMESH_STATUS_CHECK(status);
  std::cout << "refnormals " << std::endl;
  std::cout << refnormals << std::endl;


  wmesh_int_t num_facets;
  wmesh_int_t facets[6];
  wmesh_int_t facets_num_nodes[6];
  status = bms_element_facets(element,
			      &num_facets,
			      facets);
  WMESH_STATUS_CHECK(status);

  std::cout << "num_facets " << num_facets << std::endl;
  for (wmesh_int_t i=0; i < num_facets;++i)
    {
      std::cout << " facet " << facets[i] << std::endl;
    }
  
  std::cout << "num_cells " << num_cells << std::endl;
  
  wmesh_int_t facet_elm_ndofs[2];
  wmesh_int_t facet_u_ndofs[2];
  wmesh_int_t facet_dg_ndofs[2];
  wmesh_int_t facet_test_ndofs[2];

  const wmesh_int_t 	facet_cubature_family = WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE;
  const wmesh_int_t 	facet_cubature_degree = (velocity_space_->m_degree + 1)  + self_->m_degree + self_->m_degree + 2; // * 2 : degree jacobian
  wmesh_cubature_t<T> 	facet_cubatures[2];
  wmesh_int_t qn[2];
  if (topodim==2)
    {
      status = wmesh_cubature_def(&facet_cubatures[0],
				  WMESH_ELEMENT_EDGE,
				  facet_cubature_family,
				  facet_cubature_degree);
      WMESH_STATUS_CHECK(status);
      qn[0] = facet_cubatures[0].m_c.n;
    }
  else if (topodim==3)
    {
      status = wmesh_cubature_def(&facet_cubatures[0],
				  WMESH_ELEMENT_TRIANGLE,
				  facet_cubature_family,
				  facet_cubature_degree);
      WMESH_STATUS_CHECK(status);
      qn[0] = facet_cubatures[0].m_c.n;
      
      status = wmesh_cubature_def(&facet_cubatures[1],
				  WMESH_ELEMENT_QUADRILATERAL,
				  facet_cubature_family,
				  facet_cubature_degree);
      
      WMESH_STATUS_CHECK(status);
      qn[1] = facet_cubatures[1].m_c.n;
    }
  else if (topodim==1)
    {
    }

  
  
  if (topodim==2)
    {
      status = bms_ndofs(WMESH_ELEMENT_EDGE,
			 1,
			 &facet_elm_ndofs[0]);      
      WMESH_STATUS_CHECK(status);

      status = bms_ndofs(WMESH_ELEMENT_EDGE,
			 velocity_space_->m_degree,
			 &facet_u_ndofs[0]);
      
      WMESH_STATUS_CHECK(status);

      status = bms_ndofs(WMESH_ELEMENT_EDGE,
			 self_->m_degree,
			 &facet_dg_ndofs[0]);
      
      WMESH_STATUS_CHECK(status);
      
      status = bms_ndofs(WMESH_ELEMENT_EDGE,
			 self_->m_degree,
			 &facet_test_ndofs[0]);
      
      WMESH_STATUS_CHECK(status);
    }
  else
    {
      status = bms_ndofs(WMESH_ELEMENT_TRIANGLE,
			 1,
			 &facet_elm_ndofs[0]);      
      WMESH_STATUS_CHECK(status);

      status = bms_ndofs(WMESH_ELEMENT_TRIANGLE,
			 velocity_space_->m_degree,
			 &facet_u_ndofs[0]);
      
      WMESH_STATUS_CHECK(status);

      status = bms_ndofs(WMESH_ELEMENT_TRIANGLE,
			 self_->m_degree,
			 &facet_dg_ndofs[0]);
      
      WMESH_STATUS_CHECK(status);
      
      status = bms_ndofs(WMESH_ELEMENT_TRIANGLE,
			 self_->m_degree,
			 &facet_test_ndofs[0]);
      
      WMESH_STATUS_CHECK(status);


      status = bms_ndofs(WMESH_ELEMENT_QUADRILATERAL,
			 1,
			 &facet_elm_ndofs[1]);      
      WMESH_STATUS_CHECK(status);

      status = bms_ndofs(WMESH_ELEMENT_QUADRILATERAL,
			 velocity_space_->m_degree,
			 &facet_u_ndofs[1]);
      
      WMESH_STATUS_CHECK(status);

      status = bms_ndofs(WMESH_ELEMENT_QUADRILATERAL,
			 self_->m_degree,
			 &facet_dg_ndofs[1]);
      
      WMESH_STATUS_CHECK(status);
      
      status = bms_ndofs(WMESH_ELEMENT_QUADRILATERAL,
			 self_->m_degree,
			 &facet_test_ndofs[1]);
      
      WMESH_STATUS_CHECK(status);

    }
  
  wmesh_int_t 				coofacet_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			coofacet[2];
  wmesh_int_t facet_num_types = (topodim==3) ? 2 : ( (topodim==2) ? 1 : 1);
  for (wmesh_int_t jtype=0;jtype<facet_num_types;++jtype)
    {
      wmesh_mat_t<T>::alloc(&coofacet[jtype], topodim, facet_elm_ndofs[jtype]);
    }

  wmesh_int_t 				ufacet_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			ufacet[2];
  for (wmesh_int_t jtype=0;jtype<facet_num_types;++jtype)
    {
      wmesh_mat_t<T>::alloc(&ufacet[jtype], topodim, facet_u_ndofs[jtype]);
    }

  //
  // EVALS
  //
  wmesh_int_t 				coofacet_eval_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			coofacet_eval[2];
  for (wmesh_int_t i=0;i<facet_num_types;++i)
    {
      wmesh_mat_t<T>::alloc(&coofacet_eval[i], facet_elm_ndofs[i],qn[i]);
    }
  
  wmesh_int_t 				ufacet_eval_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			ufacet_eval[2];  
  for (wmesh_int_t i=0;i<facet_num_types;++i)
    {
      wmesh_mat_t<T>::alloc(&ufacet_eval[i], facet_u_ndofs[i], qn[i]);
    }
  
  wmesh_int_t 				phij_eval_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			phij_eval[2];
  for (wmesh_int_t i=0;i<facet_num_types;++i)
    {
      wmesh_mat_t<T>::alloc(&phij_eval[i], facet_dg_ndofs[i], qn[i]);
    }
  
  wmesh_int_t 				phii_eval_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			phii_eval[2];
  for (wmesh_int_t i=0;i<facet_num_types;++i)
    {
      wmesh_mat_t<T>::alloc(&phii_eval[i], facet_test_ndofs[i], qn[i]);
    }


  
  //
  // RESTRICTS
  //
  wmesh_int_t 				cooelm_restrict_op_storage 	= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			cooelm_restrict_op[6];
  for (wmesh_int_t i=0;i<num_facets;++i)
    {
      wmesh_mat_t<T>::alloc(&cooelm_restrict_op[i],facet_elm_ndofs[i],elm_ndofs);
    }
  
  wmesh_int_t 				uelm_restrict_op_storage 	= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			uelm_restrict_op[6];
  for (wmesh_int_t i=0;i<num_facets;++i)
    {
      wmesh_mat_t<T>::alloc(&uelm_restrict_op[i],facet_u_ndofs[i],u_ndofs);
    }

  wmesh_int_t 				phij_restrict_op_storage 	= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			phij_restrict_op[6];
  for (wmesh_int_t i=0;i<num_facets;++i)
    {
      wmesh_mat_t<T>::alloc(&phij_restrict_op[i],facet_dg_ndofs[i],dg_ndofs);
    }
  
  wmesh_int_t 				phii_restrict_op_storage 	= WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T>			phii_restrict_op[6];
  for (wmesh_int_t i=0;i<num_facets;++i)
    {
      wmesh_mat_t<T>::alloc(&phii_restrict_op[i], facet_test_ndofs[i], test_ndofs);
    }

  
  
#if 0  
  bms_restrict_op(WMESH_ELEMENT_TRIANGLE,
		  WMESH_MAT_FORWARD(uelm_restrict_op[0]));
  
  
  wmesh_int_t  				cubature_facet_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T> 			cubature_facet[2];
  for (wmesh_int_t jtype=0;jtype<facet_num_types;++jtype)
    {
      wmesh_mat_t<T>::alloc(&cubature_facet[jtype], topodim, qn[jtype]);
    }
  
  wmesh_int_t  				cubature_velocity_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T> 			cubature_velocity[2];

#endif
  
  for (wmesh_int_t ielm=0; ielm < num_cells;++ielm)
    {
      
      wmesh_int_t verif_num_facets = 0;
      status = wmesh_get_c2celm(mesh,
				itype_,
				ielm,
				&verif_num_facets,
				c2celm,
				c2celm_inc,
				c2celm_types,
				c2celm_types_inc);
      WMESH_STATUS_CHECK(status);


      for (wmesh_int_t j=0;j<num_facets;++j)
	{
	  wmesh_int_t nei_idx  = c2celm[c2celm_inc*j];
	  wmesh_int_t nei_type = c2celm_types[c2celm_types_inc*j];
	  std::cout << "(" << nei_idx << "," << nei_type << ")" << std::endl; 
#if 0
	  ststus = wmeshspacedg_fem_advection_extra_diagonal_block2(fem_,
								    nei_type,
								    
								    cooelm_storage_,
								    cooelm_,
								    
								    cooelm_restrict_op_storage_,
								    cooelm_restrict_op_,
								    
								    coofacet_eval_storage_,
								    coofacet_eval_,
								
								    uelm_storage_,
								    uelm_,
								    uelm_restrict_op_storage_,
								    uelm_restrict_op_,
								    
								    ufacet_eval_storage_,
								    ufacet_eval_,
								
								    phij_restrict_op_storage_,
								    phij_restrict_op_,

								    phij_eval_storage_,
								    phij_eval_,

								    phii_restrict_op_storage_,
								    phii_restrict_op_,

								    phii_eval_storage_,
								    phii_eval_,
								    
								    lmat_,
								    
								    ufacet_storage_,
								    ufacet_,
								    
								    coofacet_storage_,
								    coofacet_,
								    
								    cubature_facet_storage_,
								    cubature_facet_,
								    
								    cubature_velocity_storage_,
								    cubature_velocity_,
								    
								    rw_n_,
								    rw_);
#endif	    
	}
      
      
    }

  
#if 0
  
  for (wmesh_int_t ielm=0; ielm < num_cells;++ielm)
    {
      
      //
      // Get connectivity of element.
      //
      status = wmesh_get_c2nelm(self_,
				itype_,
				ielm,
				c2nelm,
				c2nelm_inc);
      WMESH_STATUS_CHECK(status);

      //
      // Get neighbors information.
      //
      wmesh_int_t num_facets = 0;
      status = wmesh_get_c2celm(self_,
				itype_,
				ielm,
				&num_facets,
				c2celm,
				c2celm_inc,
				c2celm_types,
				c2celm_types_inc);
      WMESH_STATUS_CHECK(status);
      
      
      //
      // Get coordinates of element.
      //
      status = wmesh_get_cooelm(self_,
				itype_,
				ielm,
				cooelm_storage,
				cooelm_m,
				cooelm_n,
				cooelm,
				cooelm_ld);
      WMESH_STATUS_CHECK(status);


      //
      // Get velocity dofs.
      //
      status = wmeshspace_get_elmdofs(velocity_space_,
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

      
      //
      // Only diagonal.
      //
      status = wmeshspacedg_fem_advection_diagonal_block(fem_,
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

							 rw_n,
							 rw);
      WMESH_STATUS_CHECK(status);
      
      //
      // Now extra-diagonal.
      //
      
      //
      // Get dof indices within the zone.
      //
      for (wmesh_int_t l=0;l<size_blocks[itype_];++l)
	{
	  idofs[l] = dofs_ptr[itype_] + ielm * size_blocks[itype_] + l;
	}
      for (wmesh_int_t ifacet=0;ifacet<num_facets;++ifacet)
	{

	  bool has_nei = (c2celm_types[ifacet] >= 0);
	  if (has_nei)
	    continue;
	  
	  wmesh_int_t jelm 	= c2celm[ifacet];
	  wmesh_int_t jtype 	= c2celm_types[ifacet];
	  wmesh_int_t jfacet    = -1;
	  //
	  // Need to know jfacet.
	  //	  
	  wmesh_int_t nei_num_facets = 0;
	  status = wmesh_get_c2celm(self_,
				    jtype,
				    jelm,
				    &nei_num_facets,
				    nei_c2celm,
				    1,
				    nei_c2celm_types,
				    1);
	  WMESH_STATUS_CHECK(status);
	  for (wmesh_int_t l=0;l<nei_num_facets;++l)
	    {
	      if (nei_c2celm_types[l]>=0)
		{
		  if (nei_c2celm[l] == ielm)
		    {
		      jfacet = l;
		      break;
		    }
		}
	    }
	  if (jfacet < 0)
	    {
	      std::cout << "inconsistent c2c"  << std::endl;
	    }


	  //
	  // Get matrix interpolation on the face (k*k)x(k*k*k)
	  //
	  
	  //
	  // Get restriction.
	  //
	  
	  
	  //
	  // Get connectivity of the face.
	  //
	  wmesh_int_t c2n_facet[4],c2n_facet_inc=1;
	  status = wmesh_get_c2nfacet(self_,
				      itype_,
				      ielm,
				      ifacet,
				      c2n_facet,
				      c2n_facet_inc);
	  WMESH_STATUS_CHECK(status);
	  

	  wmesh_int_t nei_c2n_facet[4],nei_c2n_facet_inc=1;
	  status = wmesh_get_c2nfacet(self_,
				      jtype,
				      jelm,
				      jfacet,
				      nei_c2n_facet,
				      nei_c2n_facet_inc);
	  WMESH_STATUS_CHECK(status);


	  //
	  //
	  //
	  wmesh_int_t signed_rotation = 0;
	  {
	    wmesh_int_t min = c2n_facet[0];
	    for (wmesh_int_t l=1;l<;++l)
	      {
		if (min > c2n_facet[l])
		  {
		    min = c2n_facet[l];
		    signed_rotation = l;
		  }
	      }
	  }
	  
	  wmesh_int_t nei_signed_rotation = 0;
	  {
	    wmesh_int_t min = nei_c2n_facet[0];
	    for (wmesh_int_t l=1;l<;++l)
	      {
		if (min > nei_c2n_facet[l])
		  {
		    min = nei_c2n_facet[l];
		    nei_signed_rotation = l;
		  }
	      }
	  }


	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_element 	= &fem_->m_shape_eval_boundary_element[element];
	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_test 	= &fem_->m_shape_eval_boundary_test[element];
	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_f 	= &fem_->m_shape_eval_boundary_f[element];
	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_f_nei 	= &fem_->m_shape_eval_boundary_f[nei_element];
	  const wmesh_shape_eval_boundary_t<T> * eval_boundary_u 	= &fem_->m_shape_eval_boundary_f[element];

	  const wmesh_cubature_t<T> * 	cubature_element 	= &eval_boundary_element->m_cubature_facets[facet_type][0][signed_rotation];
	  const wmesh_shape_eval_t<T> * eval_element 		= &eval_boundary_element->m_shape_eval_facets[facet_type][0][signed_rotation];
	  const wmesh_shape_eval_t<T> * eval_test		= &eval_boundary_test->m_shape_eval_facets[facet_type][0][signed_rotation];
	  const wmesh_shape_eval_t<T> * eval_f			= &eval_boundary_f->m_shape_eval_facets[facet_type][0][signed_rotation];
	  const wmesh_shape_eval_t<T> * eval_f_nei		= &eval_boundary_f->m_shape_eval_facets[nei_facet_type][1][nei_signed_rotation];
	  
	  //
	  // Compute the normals.
	  //
	  for (wmesh_int_t idim=0;idim<topodim;++idim)
	    {
	      BLAS_dgemm("N",
			 "N",
			 &topodim,
			 &eval_element->m_diff[idim].n,
			 &eval_element->m_diff[idim].m,
			 &r1,
			 cooelm_,
			 &cooelm_ld_,
			 eval_element->m_diff[idim].v,
			 &eval_element->m_diff[idim].ld,
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
	  // Compute the reference normal.
	  //	  
	  for (wmesh_int_t j=0;j<q_n;++j)
	    {
	      BLAS_dgemm("T",
			 "N",
			 &topodim,
			 &n1,
			 &topodim,
			 &r1,
			 jacobian,
			 &topodim,
			 reference_normal,
			 &topodim,
			 &r0,
			 physical_normal,
			 &topodim);
	    }


	  //
	  // Evaluate u.
	  //

	  //
	  // u.dot(n).
	  //	  
	  for (wmesh_int_t j=0;j<q_n;++j)
	    {
	      T dot = static_cast<T>(0);
	      for (wmesh_int_t i=0;i<topodim;++i)
		{
		  dot += n[i]*u[i];
		}
	      dot = ( dot < 0.0 ) ? -dot : 0.0;
	      dot * qfacet_w[j] * dets[j];
	    }

	  
	  
	}
	    
    }


      
      //
      // Assembly.
      //
      for (wmesh_int_t i=0;i<size_blocks[itype_];++i)
	{
	  const wmesh_int_t idof = idofs[i];
	  for (wmesh_int_t j=0;j<size_blocks[itype_];++j)
	    {
	      const wmesh_int_t jdof = idofs[j];
	      
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
		      csr_val_[at] += lmat[j*lmat_ld + i];
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

#if 0
      for (wmesh_int_t i=0;i<num_local_dofs;++i)
	{
	  rhs_[c2d[i]] += lrhs[i];
	}
#endif
      
    }

  //
  // No force ...
  //

  
  //
  // Apply weak boundary conditions.
  //
  for (wmesh_int_t ielm =0;ielm<num_cells;++ielm)
    {

      for (wmesh_int_t ifacet=0;ifacet<self_->m_c2c.m_m[itype_];++ifacet)
	{
	  auto nei_info = self_->m_c2c.m_data[self_->m_c2c.m_ptr[itype_]  + self_->m_c2c.m_ld[itype_] * ielm + ifacet ];
	  if (nei_info == 0)
	    {
	      
	      //
	      // This is a candidate for boundary condition.
	      //

	      //
	      // int_ (u.n) bc(x,y) test_i 
	      //	      
	      wmesh_int_t jelm,jtype;
	      status = bms_c2c_cindex(nei_info,&jelm);
	      WMESH_STATUS_CHECK(status);
	      status = bms_c2c_ctype(nei_info,&jtype);
	      WMESH_STATUS_CHECK(status);	  
	    }
	}
    }
  
  
#if 0
  
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

#endif  
  return WMESH_STATUS_SUCCESS;
}


template <typename T>
wmesh_status_t wmeshspacedg_fem_advection_global_system(const wmeshspacedg_t *				self_,							
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
  const wmesh_t * 	mesh 	= self_->m_mesh;  
  const wmesh_int_t 	ntypes 	= mesh->m_c2n.m_size;
  for (wmesh_int_t itype=0;itype<ntypes;++itype)
    {
      const wmesh_int_t num_elements = mesh->m_c2n.m_n[itype];
      if (num_elements > 0)
	{
	  std::cout << "Compute global system for cells of type " << itype << std::endl;
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
  
  wmesh_status_t wmeshspacedg_advection(  const wmeshspacedg_t *			self_,
					const wmesh_shape_info_t* __restrict__	shape_info_f_,
					const wmesh_shape_info_t* __restrict__	shape_info_u_,
					const wmesh_shape_info_t* __restrict__	shape_info_test_,
					
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
    const wmesh_t *				mesh = self_->m_mesh;
    std::cout << "wmeshspacedg_advection ... " << std::endl;
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
    
    
    //
    // Initialize data.
    //

    const wmesh_int_t topodim 	= mesh->m_topology_dimension;
    const wmesh_int_t 	cubature_degree 	= (shape_info_u_->m_degree + (shape_info_f_->m_degree-1) + shape_info_test_->m_degree) * 2;
    wmesh_int_t ntypes;
    wmesh_int_t elements[4];
    status =  bms_topodim2elements(topodim,
				   &ntypes,
				   elements);
    WMESH_STATUS_CHECK(status);
    WMESH_CHECK(ntypes == mesh->m_c2n.m_size);


    std::cout << "Prepare data. " << std::endl;    
    wmeshspacedg_fem_advection_t<double> fem;
    memset(&fem,0,sizeof(wmeshspacedg_fem_advection_t<double>));
    for (wmesh_int_t itype=0;itype< ntypes;++itype)
      {
	if (mesh->m_c2n.m_n[itype] > 0)
	  {
	    status = wmeshspacedg_fem_advection_add(&fem,
						    elements[itype],
						    cubature_family,
						    cubature_degree,
						    shape_element_family,
						    shape_element_degree,				   
						    shape_info_u_->m_family,
						    shape_info_u_->m_degree,
						    shape_info_f_->m_family,
						    shape_info_f_->m_degree,
						    shape_info_test_->m_family,
						    shape_info_test_->m_degree);
	    WMESH_STATUS_CHECK(status);
	  }
      }
    std::cout << "Prepare data done." << std::endl;    
    std::cout << "Compute global system. " << std::endl;
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
    std::cout << "Compute global system done. " << std::endl;
    
    return WMESH_STATUS_SUCCESS;
  }
}
#endif
*/
