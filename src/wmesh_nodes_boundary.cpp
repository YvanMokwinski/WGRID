
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh.hpp"
#include <chrono>

#include "bms.h"
#include "wmesh-utils.hpp"
#include "wmesh-blas.hpp"
#include "bms_templates.hpp"
#include "bms.hpp"
#include "wmesh_nodes_boundary_t.hpp"
#include "wmesh_shape_eval_factory_t.hpp"
#include "wmesh_nodes_factory_t.hpp"

template<typename T>
static std::ostream& operator<<(std::ostream&out_,
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

//
// 
//

template<typename T>
wmesh_status_t wmesh_nodes_boundary_def(wmesh_nodes_boundary_t<T>*__restrict__ 	self_,
					wmesh_int_t 				element_,
					wmesh_int_t 				nodes_family_,
					wmesh_int_t 				nodes_degree_)
{
  //
  // Initialize the data.
  //
  memset(self_,0,sizeof(wmesh_nodes_boundary_t<T>));

  wmesh_int_t topodim;
  wmesh_status_t status = bms_element2topodim(element_,
			       &topodim);
  WMESH_STATUS_CHECK(status);

  wmesh_int_t element_num_nodes;
  status = bms_elements_num_nodes(1,
				  &element_,
				  &element_num_nodes);
  WMESH_STATUS_CHECK(status);
  
  wmesh_mat_t<T> ref_cooelm;
  wmesh_mat_t<T>::alloc(&ref_cooelm,topodim,element_num_nodes);
  WMESH_STATUS_CHECK(status);
  status = bms_element_geometry(element_,ref_cooelm.v);
  WMESH_STATUS_CHECK(status);
  std::cout << ref_cooelm << std::endl;;

  wmesh_int_t num_facets;
  status = bms_element_facets(element_,			 
			      &num_facets,
			      self_->m_facets);
  WMESH_STATUS_CHECK(status);
  self_->m_num_facets = num_facets;

  self_->m_facets_nodes_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_int_t s_m[2],s_n[2],s_ld[2],s_v[2][64];
  const wmesh_int_t element_type = (topodim==3) ? (element_ -4):  ( (topodim==2) ? (element_-2) : (element_-1) );
  const wmesh_shape_eval_t<T>*facet_shape_eval_coodofs[2];

  if (topodim==2)
    {
      const wmesh_nodes_t<T> * nodes = wmesh_nodes_factory_t<T>::nodes_instance(WMESH_ELEMENT_EDGE,
										nodes_family_,
										nodes_degree_);
      
      facet_shape_eval_coodofs[0] = wmesh_shape_eval_factory_t<T>::shape_eval_instance(WMESH_ELEMENT_EDGE,
										       WMESH_SHAPE_FAMILY_LAGRANGE,
										       1,
										       nodes->m_c_storage,
										       &nodes->m_c);
      status =  bms_s_e2n_type(element_type,
			       topodim,
			       &s_m[0],
			       &s_n[0],
			       &s_v[0][0],
			       &s_ld[0]);
    }
  else   if (topodim==3)
    {
      {
	const wmesh_nodes_t<T> * nodes = wmesh_nodes_factory_t<T>::nodes_instance(WMESH_ELEMENT_TRIANGLE,
										  nodes_family_,
										  nodes_degree_);
	
	facet_shape_eval_coodofs[0] = wmesh_shape_eval_factory_t<T>::shape_eval_instance(WMESH_ELEMENT_TRIANGLE,
											 WMESH_SHAPE_FAMILY_LAGRANGE,
											 1,
											 nodes->m_c_storage,
											 &nodes->m_c);
      }
      
      
      {
	const wmesh_nodes_t<T> * nodes = wmesh_nodes_factory_t<T>::nodes_instance(WMESH_ELEMENT_QUADRILATERAL,
										  nodes_family_,
										  nodes_degree_);
	
	facet_shape_eval_coodofs[1] = wmesh_shape_eval_factory_t<T>::shape_eval_instance(WMESH_ELEMENT_QUADRILATERAL,
											 WMESH_SHAPE_FAMILY_LAGRANGE,
											 1,
											 nodes->m_c_storage,
											 &nodes->m_c);
      }

      status =  bms_s_t2n_type(element_type,
			       &s_m[0],
			       &s_n[0],
			       &s_v[0][0],
			       &s_ld[0]);
      status =  bms_s_q2n_type(element_type,
			       &s_m[1],
			       &s_n[1],
			       &s_v[1][0],
			       &s_ld[1]);
    }
  else
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_CONFIG);
    }

  wmesh_mat_t<T> coofacet;
  T coofacet_data[32];
  for (wmesh_int_t ifacet=0;ifacet < num_facets;++ifacet)
    {      
      const wmesh_int_t 	facet 		= self_->m_facets[ifacet];
      const wmesh_int_t 	facet_type	= (topodim==3) ? (facet - 2) : (facet-1);
      const wmesh_int_t 	facet_num_nodes = (topodim==3) ? (facet_type + 3) : ( (topodim==2) ? 2 : 1 );// be careful, trick.

      wmesh_mat_t<T>::define(&coofacet,
			     topodim,
			     facet_num_nodes,
			     coofacet_data,
			     topodim);

      
      const wmesh_int_t s_f2n_m 	= s_m[facet_type];
      const wmesh_int_t s_f2n_n 	= s_n[facet_type];
      const_wmesh_int_p s_f2n_v 	= &s_v[facet_type][0];
      const wmesh_int_t s_f2n_ld 	= s_ld[facet_type];
      for (wmesh_int_t i=0;i<s_f2n_m;++i)
	{
	  const wmesh_int_t idx = s_f2n_v[s_f2n_ld*ifacet+i];
	  for (wmesh_int_t k=0;k<topodim;++k)
	    {
	      coofacet.v[coofacet.ld*i+k] = ref_cooelm.v[ref_cooelm.ld*idx + k];
	    }
	}

      for (wmesh_int_t signed_rotation_idx=0;signed_rotation_idx<facet_num_nodes;++signed_rotation_idx)
	{

	  wmesh_mat_t<T>::alloc(&self_->m_facets_nodes[ifacet][signed_rotation_idx],coofacet.m,facet_shape_eval_coodofs[facet_type]->m_f.n);
	  //
	  // Generate the coordinates.
	  //
	  T r0 = static_cast<T>(0);
	  T r1 = static_cast<T>(1);
	  
	  wmesh_mat_gemm(r1,
			 coofacet,
			 facet_shape_eval_coodofs[facet_type]->m_f,
			 r0,
			 self_->m_facets_nodes[ifacet][signed_rotation_idx]);
	  
	  //
	  // transform for the next rotation, signed_rotation_idx is the position of the first node.
	  //
	  //	  if (signed_rotation_idx<facet_num_nodes-1)
	  {
	    T tmp[topodim];
	    for (wmesh_int_t i=0;i<topodim;++i)
	      {
		tmp[i] = coofacet.v[coofacet.ld*(facet_num_nodes-1)+i];
	      }
	    for (wmesh_int_t j=facet_num_nodes-1;j>0;--j)
	      {
		for (wmesh_int_t i=0;i<topodim;++i)
		  {
		    coofacet.v[coofacet.ld*j+i] = coofacet.v[coofacet.ld*(j-1)+i];
		  }
	      }
	    for (wmesh_int_t i=0;i<topodim;++i)
	      {
		coofacet.v[coofacet.ld*0+i] = tmp[i];
	      }
	  }
	}

      //
      // now reverse.
      //      
      {
	T tmp[topodim];
	for (wmesh_int_t j=1;j<=facet_num_nodes/2;++j)
	  {
	    for (wmesh_int_t i=0;i<topodim;++i)
	      {
		tmp[i] = coofacet.v[coofacet.ld*j+i];
	      }
	    for (wmesh_int_t i=0;i<topodim;++i)
	      {
		coofacet.v[coofacet.ld*j+i] = coofacet.v[coofacet.ld * ( facet_num_nodes-j ) + i];
	      }
	    for (wmesh_int_t i=0;i<topodim;++i)
	      {
		coofacet.v[coofacet.ld * ( facet_num_nodes-j ) + i] = tmp[i];
	      }
	  }
      }
      
      for (wmesh_int_t signed_rotation_idx=0;signed_rotation_idx<facet_num_nodes;++signed_rotation_idx)
	{
	  wmesh_mat_t<T>::alloc(&self_->m_facets_nodes[ifacet][facet_num_nodes+signed_rotation_idx],coofacet.m,facet_shape_eval_coodofs[facet_type]->m_f.n);
	  //
	  // Generate the coordinates.
	  //
	  T r0 = static_cast<T>(0);
	  T r1 = static_cast<T>(1);
	  wmesh_mat_gemm(r1,
			 coofacet,
			 facet_shape_eval_coodofs[facet_type]->m_f,
			 r0,
			 self_->m_facets_nodes[ifacet][facet_num_nodes + signed_rotation_idx]);
	  
	  //
	  // transform for the next rotation, signed_rotation_idx is the position of the first node.
	  //
	  if (signed_rotation_idx<facet_num_nodes-1)
	  {
	    T tmp[topodim];
	    for (wmesh_int_t i=0;i<topodim-1;++i)
	      {
		tmp[i] = coofacet.v[coofacet.ld*(facet_num_nodes-1)+i];
	      }
	    for (wmesh_int_t j=facet_num_nodes-1;j>0;--j)
	      {
		for (wmesh_int_t i=0;i<topodim-1;++i)
		  {
		    coofacet.v[coofacet.ld*j+i] = coofacet.v[coofacet.ld*(j-1)+i];
		  }
	      }
	    for (wmesh_int_t i=0;i<topodim-1;++i)
	      {
		coofacet.v[coofacet.ld*0+i] = tmp[i];
	      }
	  }
	}

      
    }      
  
  return WMESH_STATUS_SUCCESS;
}

template
wmesh_status_t wmesh_nodes_boundary_def<float>(wmesh_nodes_boundary_t<float>*__restrict__ 	self_,
					       wmesh_int_t 					element_,
					       wmesh_int_t 					nodes_family_,
					       wmesh_int_t 					nodes_degree_);
template
wmesh_status_t wmesh_nodes_boundary_def<double>(wmesh_nodes_boundary_t<double>*__restrict__ 	self_,
						wmesh_int_t 					element_,
						wmesh_int_t 					nodes_family_,
						wmesh_int_t 					nodes_degree_);
