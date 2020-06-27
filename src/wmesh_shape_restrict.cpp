
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
wmesh_status_t wmesh_shape_restrict_def(wmesh_shape_restrict_t<T>*__restrict__ 	self_,
					wmesh_int_t 				element_,
					wmesh_int_t 				shape_family_,
					wmesh_int_t 				shape_degree_)
{
  //
  // Initialize the data.
  //
  memset(self_,0,sizeof(wmesh_shape_restrict_t<T>));

  //
  // Define shape over element.
  // 
  wmesh_status_t status;
  status = wmesh_shape_def(&self_->m_shape,
			   element_,
			   shape_family_,
			   shape_degree_);
  WMESH_STATUS_CHECK(status);

  
  wmesh_int_t element_num_nodes;
  status = bms_elements_num_nodes(1,
				  &element_,
				  &element_num_nodes);
  WMESH_STATUS_CHECK(status);

  
  //
  // topodim.
  //
  wmesh_int_t topodim;
  status = bms_element2topodim(element_,
			       &topodim);
  WMESH_STATUS_CHECK(status);
  const wmesh_int_t facet_topodim = topodim - 1;
  
  //
  // Register facets.
  //
  wmesh_int_t num_facets;
  status = bms_element_facets(element_,			 
			      &num_facets,
			      self_->m_facets);
  WMESH_STATUS_CHECK(status);
  self_->m_num_facets = num_facets;

  
  for (wmesh_int_t ifacet=0;ifacet < num_facets;++ifacet)
    {
      const wmesh_int_t facet = self_->m_facets[ifacet];
      self_->m_facet_types[ifacet] = (topodim==3) ? (facet-2) : ( (topodim==2) ? (facet-1) : (facet-0) );
    }
  
  for (wmesh_int_t ifacet=0;ifacet < num_facets;++ifacet)
    {
      const wmesh_int_t facet = self_->m_facets[ifacet];
      status = bms_ndofs(facet,
			 shape_degree_,
			 &self_->m_facets_num_dofs[ifacet]);
      WMESH_STATUS_CHECK(status);
      
      status = bms_ndofs(facet,
			 1,
			 &self_->m_facets_num_nodes[ifacet]);
      WMESH_STATUS_CHECK(status);
    }
  
  for (wmesh_int_t ifacet=0;ifacet < num_facets;++ifacet)
    {
      for (wmesh_int_t signed_rotation_idx=0;signed_rotation_idx < self_->m_facets_num_nodes[ifacet];++signed_rotation_idx)
	{
	  wmesh_mat_t<T>::alloc(&self_->m_restrict[ifacet][signed_rotation_idx],
				self_->m_shape.m_ndofs,
				self_->m_facets_num_dofs[ifacet]);
	}
    }

  //
  // Get the reference geometry of the cell.
  //
  wmesh_mat_t<T> ref_cooelm;
  wmesh_mat_t<T>::alloc(&ref_cooelm,topodim,element_num_nodes);
  WMESH_STATUS_CHECK(status);
  status = bms_element_geometry(element_,ref_cooelm.v);
  WMESH_STATUS_CHECK(status);
  std::cout << ref_cooelm << std::endl;;

    
  //
  // Traverse the facets.
  //
  wmesh_int_t s_m[2],s_n[2],s_ld[2],s_v[2][64];
  const wmesh_int_t element_type = (topodim==3) ? (element_ -4):  ( (topodim==2) ? (element_-2) : (element_-1) );
  if (topodim==2)
    {
      status =  bms_s_e2n_type(element_type,
			       topodim,
			       &s_m[0],
			       &s_n[0],
			       &s_v[0][0],
			       &s_ld[0]);
    }
  else   if (topodim==3)
    {
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
  
  
  T facet_nodes[32];
  wmesh_int_t f2n[4];
  T coofacet_data[32];
  wmesh_int_t coofacet_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T> coofacet;
  wmesh_int_t coodofs_storage = WMESH_STORAGE_INTERLEAVE;
  wmesh_mat_t<T> coodofs;

  wmesh_mat_t<T> facet_eval_nodes[2];
  
  {
  bool has[2]{};
  for (wmesh_int_t ifacet=0;ifacet < num_facets;++ifacet)
    {
      const wmesh_int_t 	facet 		= self_->m_facets[ifacet];
      const wmesh_int_t 	facet_type	= self_->m_facet_types[ifacet];
      const wmesh_int_t 	facet_num_dofs 	= self_->m_facets_num_dofs[ifacet];
      const wmesh_int_t 	facet_num_nodes = (topodim==3) ? (facet_type + 3) : ( (topodim==2) ? 2 : 1 );// be careful, trick.
      if (!has[facet_type])
	{
  has[facet_type]=true;
  wmesh_shape_t shape_facet;
  wmesh_shape_def(&shape_facet,facet,shape_family_,shape_degree_);
  status = wmesh_shape_calculate_eval<T>(shape_facet,coodofs_storage,coodofs,facet_eval_nodes[facet_type]);
	  WMESH_STATUS_CHECK(status);
}
    }
 }

#if 0
  wmesh_shape_t shape_element;
  wmesh_shape_def(&shape_element,element,shape_family_,shape_degree_);
#endif

  //  wmesh_int_t  = ( (topodim==3) ?  ( (self_->facet_num_dofs[0]<self_->facet_num_dofs[1]) ? self_->facet_num_dofs[1] : self_->facet_num_dofs[0] ) : self_->facet_num_dofs[0] );
  //  wmesh_mat_t<T> coodofs;
  //  wmesh_mat_t<T>::alloc(&coodofs,topodim,element_num_nodes);
  //
  //
  // trick be careful !!!!!!!!!!!!!!!!!!! (take the last one, since triangle are before quadrilateral if (any).
  T * coodofs_data = (T*)malloc(sizeof(wmesh_int_t) * self_->m_facets_num_dofs[self_->m_num_facets-1]);
  
  for (wmesh_int_t ifacet=0;ifacet < num_facets;++ifacet)
    {

      const wmesh_int_t 	facet 		= self_->m_facets[ifacet];
      const wmesh_int_t 	facet_type	= self_->m_facet_types[ifacet];
      const wmesh_int_t 	facet_num_dofs 	= self_->m_facets_num_dofs[ifacet];
      const wmesh_int_t 	facet_num_nodes = (topodim==3) ? (facet_type + 3) : ( (topodim==2) ? 2 : 1 );// be careful, trick.

      wmesh_mat_t<T>::define(&coofacet,topodim,facet_num_nodes,coofacet_data,topodim);
      wmesh_mat_t<T>::define(&coodofs,topodim,facet_num_dofs,coodofs_data,topodim);
      
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
      std::cout << "ifacet " << ifacet << std::endl;
      std::cout << "     ifacet " << ifacet << std::endl;
      std::cout << "     coo " << std::endl << coofacet  << std::endl;
      
      for (wmesh_int_t signed_rotation_idx=0;signed_rotation_idx<facet_num_nodes;++signed_rotation_idx)
	{
	  const wmesh_int_t signed_rotation = signed_rotation_idx+1;// (j < facet_num_nodes) ?  (j+1) : (facet_num_nodes - j - 1);
	  wmesh_mat_t<T>& restrict_mat = self_->m_restrict[ifacet][signed_rotation_idx];


	  // 
	  // coofacet_cell = coofacet * facet_eval
	  // 
	  //
	  // Transform (n-1)d to nd coordinates.
	  //
	  T r1 = static_cast<T>(1);
	  T r0 = static_cast<T>(0);
	  wmesh_mat_gemm(r1, coofacet, facet_eval_nodes[facet_type], r0, coodofs);
	  

	  status = wmesh_shape_calculate_eval(self_->m_shape,
					      coodofs_storage,
					      coodofs,
					      self_->m_restrict[ifacet][signed_rotation_idx]);
	  //

	  // transform for the next rotation, signed_rotation_idx is the position of the first node.
	  //
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
  
#if 0

  
  //
  // Allocate facet nodes coordinates.
  //
  
      
  //
  // Traverse the facets.
  //
  for (wmesh_int_t ifacet=0;ifacet < num_facets;++ifacet)
    {

      //
      // Get facet informations.
      //
      const wmesh_int_t facet 		= self_->m_facets[ifacet];
      const wmesh_int_t	facet_type	= self_->m_facet_types[ifacet];
      const wmesh_int_t	facet_num_dofs 	= self_->m_facet_num_dofs[ifacet];      
      const wmesh_int_t facet_num_nodes = self_->m_facet_num_nodes[ifacet];

      //
      // Allocate facet coordinates.
      //
      wmesh_mat_t<T> facet_coonodes;
      status = wmesh_mat_t<T>::alloc(&facet_coonodes,topodim-1,facet_num_nodes);
      WMESH_STATUS_CHECK(status);

      //
      // Extract facet coordinates
      //

      
      //
      // Allocate facet dofs coordinates.
      //
      const wmesh_int_t facet_coodofs_storage = WMESH_STORAGE_INTERLEAVE;
      wmesh_mat_t<T> facet_coodofs;
      status = wmesh_mat_t<T>::alloc(&facet_coodofs,topodim,facet_num_dofs);
      WMESH_STATUS_CHECK(status);
      
      //
      // For each rotation.
      //
      for (wmesh_int_t signed_rotation=1;signed_rotation <= facet_num_nodes;++signed_rotation)
	{

	  //
	  // Evaluate facet dofs coordinates.
	  //
	  wmesh_mat_gemm(r1,
			 facet_coonodes,
			 facet_eval_coodofs,
			 r0,
			 facet_coodofs);

	  //
	  // Evaluate cell shape. over facet_coodofs
	  //
	  
	  
	  //
	  // Transform for next rotation by shifting to the right.
	  //
	  {
	    T tmp[topodim - 1];
	    for (wmesh_int_t i=0;i<topodim-1;++i)
	      {
		tmp[i] = facet_coonodes.v[facet_coonodes.ld*(facet_num_nodes-1)+i];
	      }
	    for (wmesh_int_t j=facet_num_nodes-1;j>0;--j)
	      {
		for (wmesh_int_t i=0;i<topodim-1;++i)
		  {
		    facet_coonodes.v[facet_coonodes.ld*j+i] = facet_coonodes.v[facet_coonodes.ld*(j-1)+i];
		  }
	      }
	    for (wmesh_int_t i=0;i<topodim-1;++i)
	      {
		facet_coonodes.v[facet_coonodes.ld*0+i] = tmp[i];
	      }
	  }
	  

	}
      
    }
  //
  // Get the facets.
  //
  wmesh_int_p facets 		=  &self_->m_facets[0];
  //
  // Get reference coordinates of the dofs.
  //
  wmesh_mat_t<T> ref_facet_coodofs;
  wmesh_mat_t<T>::alloc(&ref_facet_coodofs,
			ref_facet_coodofs,
			facet_num_dofs);
  
  status = wmesh_shape_calculate_nodes(facet,
				       shape_family_,
				       shape_degree_,
				       ref_facet_coodofs);
  
  //
  // Get the number of 1d nodes.
  //
  wmesh_int_t num_nodes1d;
  status = bms_cubature_num_nodes(WMESH_ELEMENT_EDGE,
				  cubature_family_,
				  cubature_degree_,
				  &num_nodes1d);
  
  //
  // Init cubatures.
  //
  wmesh_cubature_t<T> ref_facets_cubature[2];
  if (topodim==3)
    {
      status = wmesh_cubature_def(&ref_facets_cubature[0],
				  WMESH_ELEMENT_TRIANGLE,
				  cubature_family_,
				  cubature_degree_);
      WMESH_STATUS_CHECK(status);
      status = wmesh_cubature_def(&ref_facets_cubature[1],
				  WMESH_ELEMENT_QUADRILATERAL,
				  cubature_family_,
				  cubature_degree_);
      WMESH_STATUS_CHECK(status);
    }
  else if (topodim==2)
    {
      status = wmesh_cubature_def(&ref_facets_cubature[0],
				  WMESH_ELEMENT_EDGE,
				  cubature_family_,
				  cubature_degree_);
      WMESH_STATUS_CHECK(status);
    }
  else
    {
      WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
    }


  //
  //
  //
  T ref_nodes[32];
  status = bms_element_geometry(element_,ref_nodes);
  WMESH_STATUS_CHECK(status);

  
  //
  // Traverse the facets.
  //
  wmesh_int_t s_m[2],s_n[2],s_ld[2],s_v[2][64];
  wmesh_int_t element_type = (topodim==3) ? (element_ -4):  ( (topodim==2) ? (element_-2) : (element_-1) );
  if (topodim==2)
    {
      status =  bms_s_e2n_type(element_type,
			       topodim,
			       &s_m[0],
			       &s_n[0],
			       &s_v[0][0],
			       &s_ld[0]);
    }
  else   if (topodim==3)
    {
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
  
  T facet_nodes[32];
  wmesh_int_t f2n[4];
  for (wmesh_int_t ifacet=0;ifacet < num_facets;++ifacet)
    {
      const wmesh_int_t facet 			= self_->m_facets[ifacet];
      wmesh_int_t 	facet_type		= self_->m_facet_types[ifacet];
      wmesh_int_t 	facet_num_dofs 	= self_->m_facet_num_dofs[ifacet];
      wmesh_int_t 	facet_num_nodes = (topodim==3) ? (facet_type + 3) : ( (topodim==2) ? 2 : 1 );// be careful, trick.
      for (wmesh_int_t j=0;j<facet_num_nodes;++j)
	{
	  const wmesh_int_t signed_rotation = j;// (j < facet_num_nodes) ?  (j+1) : (facet_num_nodes - j - 1);

	  
	  wmesh_mat_t<T>::alloc(&self_->m_restrict[ifacet][j],
				facet_num_dofs,
				cell_num_dofs);

	  //
	  // Rotate (n-1)d coordinates.
	  //
	  status =  bms_mirrored_local_coordinates(facet,
						   signed_rotation,
						   ref_facet_nodes_storage,
						   WMESH_MAT_FORWARD(ref_facet_nodes[facet_type]),
						   self_->m_restrict_storage,
						   WMESH_MAT_FORWARD(ref_rotated_facet_nodes) );
	  WMESH_STATUS_CHECK(status);
	  
	  //
	  // Transform (n-1)d to nd coordinates.
	  //
	  xgemm(ref_coofacet,eval_element,ref_rotated_facet_nodes_into_cell);
	  
	  //
	  // Now evaluate the cell shape basis over these points.
	  //
	  wmesh_shape_eval_t<T>::build_eval(ref_rotated_facet_nodes_into_cell,
					    self_->m_restrict_storage,
					    self_->m_restrict[ifacet][j]);


	 
	  //
	  // Now transform the coordinates to cell coordinates.
	  //
	  //
	  // We need the coordinates of the facet on the reference element.
	  //
	  for (wmesh_int_t i=0;i<facet_num_nodes;++i)
	    {
	      if (facet_type==0) // edge or triangle.
		{
		  f2n[i] = s_v[facet_type][s_ld[facet_type] * ifacet + i];
		  //		  std::cout << "#### " << s_v[facet_type][s_ld[facet_type] * ifacet + i] << std::endl;
		}
	      else
		{
		  //
		  // quad.
		  //
		  wmesh_int_t kfacet = 0;
		  if (element_ == WMESH_ELEMENT_PYRAMID)
		    {
		      kfacet = ifacet-4;
		    }
		  else if (element_ == WMESH_ELEMENT_HEXAHEDRON)
		    {
		      kfacet = ifacet;
		    }
		  else if (element_ == WMESH_ELEMENT_WEDGE)
		    {
		      kfacet = ifacet-2;
		    }
		  //
		  // For quad.
		  //
		  f2n[i] = s_v[facet_type][s_ld[facet_type] * kfacet + i];
		}
	    }
	  
	  for (wmesh_int_t i=0;i<facet_num_nodes;++i)
	    {
	      //	      std::cout << "??????????? " << f2n[i] << std::endl;
	      for (wmesh_int_t k=0;k<topodim;++k)
		{  
		  facet_nodes[topodim * i+k] = ref_nodes[topodim * f2n[i]+k];
		}
	    }
	  
	  T cellcoo[3];
	  T facetcoo[3]{};
	  
	  if (facet_num_nodes==2)
	    {
	      for (wmesh_int_t i=0;i<cubature_num_nodes;++i)
		{
		  for (wmesh_int_t k=0;k<topodim-1;++k)
		    {
		      facetcoo[k] = self_->m_facets_cubature[ifacet][j].m_c.v[self_->m_facets_cubature[ifacet][j].m_c.ld * i + k];
		    }
		  
		  T r  = facetcoo[0];
		  T l0 = ( static_cast<T>(1) - (r) ) / static_cast<T>(2);
		  T l1 = ( static_cast<T>(1) + (r) ) / static_cast<T>(2);
		  for (wmesh_int_t k=0;k<topodim;++k)
		    {
		      cellcoo[k] = l0 * facet_nodes[topodim * 0 + k] + l1 * facet_nodes[topodim * 1 + k];
		    }
		  //
		  // copy
		  //
		  for (wmesh_int_t k=0;k<topodim;++k)
		    {
		      self_->m_facets_cubature[ifacet][j].m_c.v[self_->m_facets_cubature[ifacet][j].m_c.ld * i + k] = cellcoo[k];
		    }		  
		}
	    }
	  else if (facet_num_nodes==3)
	    {
	      for (wmesh_int_t i=0;i<cubature_num_nodes;++i)
		{
		  
		  for (wmesh_int_t k=0;k<topodim-1;++k)
		    {
		      facetcoo[k] = self_->m_facets_cubature[ifacet][j].m_c.v[self_->m_facets_cubature[ifacet][j].m_c.ld * i + k];
		    }
		  
		  T r = facetcoo[0];
		  T s = facetcoo[1];
		  T l0 = ( static_cast<T>(1) - (r+s) );
		  T l1 = r;
		  T l2 = s;
#if 0
		  for (wmesh_int_t k=0;k<topodim;++k)
		    {
		      std::cout << "coo0 " << facet_nodes[topodim * 0 + k] << std::endl;
		    }
		  for (wmesh_int_t k=0;k<topodim;++k)
		    {
		      std::cout << "coo1 " << facet_nodes[topodim * 1 + k] << std::endl;
		    }
		  for (wmesh_int_t k=0;k<topodim;++k)
		    {
		      std::cout << "coo2 " << facet_nodes[topodim * 2 + k] << std::endl;
		    }
#endif
		  for (wmesh_int_t k=0;k<topodim;++k)
		    {
		      cellcoo[k] = l0 * facet_nodes[topodim * 0 + k] + l1 * facet_nodes[topodim * 1 + k] + l2 * facet_nodes[topodim * 2 + k];
		    }

		  //
		  // copy
		  //
		  for (wmesh_int_t k=0;k<topodim;++k)
		    {
		      self_->m_facets_cubature[ifacet][j].m_c.v[self_->m_facets_cubature[ifacet][j].m_c.ld * i + k] = cellcoo[k];
		    }		  
		}
	    }
	  else if (facet_num_nodes==4)
	    {
	      
	      for (wmesh_int_t i=0;i<cubature_num_nodes;++i)
		{

		  for (wmesh_int_t k=0;k<topodim-1;++k)
		    {
		      facetcoo[k] = self_->m_facets_cubature[ifacet][j].m_c.v[self_->m_facets_cubature[ifacet][j].m_c.ld * i + k];
		    }
	      
		  T r  = facetcoo[0];
		  T s  = facetcoo[1];
		  T l0 = ( static_cast<T>(1) - r)*( static_cast<T>(1) - s)/static_cast<T>(4);
		  T l1 = ( static_cast<T>(1) + r)*( static_cast<T>(1) - s)/static_cast<T>(4);
		  T l2 = ( static_cast<T>(1) + r)*( static_cast<T>(1) + s)/static_cast<T>(4);
		  T l3 = ( static_cast<T>(1) - r)*( static_cast<T>(1) + s)/static_cast<T>(4);
		  for (wmesh_int_t k=0;k<topodim;++k)
		    {
		      cellcoo[k] = l0 * facet_nodes[topodim * 0 + k] + l1 * facet_nodes[topodim * 1 + k] + l2 * facet_nodes[topodim * 2 + k]  + l3 * facet_nodes[topodim * 3 + k];
		    }

		  //
		  // sum_k p(x_k) w_k det(x_k) //  (1-r)/2 (0,0)  (1+r)/2 (1,0) (-1/2 1/2)
		  // 
		  //
		  //
		  
		  //
		  // copy 
		  //
		  for (wmesh_int_t k=0;k<topodim;++k)
		    {
		      self_->m_facets_cubature[ifacet][j].m_c.v[self_->m_facets_cubature[ifacet][j].m_c.ld * i + k] = cellcoo[k];
		    }		  
		}
	    }	  	  
	}
    }
#endif  
  return WMESH_STATUS_SUCCESS;
}


template
wmesh_status_t wmesh_shape_restrict_def<float>(wmesh_shape_restrict_t<float>*__restrict__ 	self_,
					       wmesh_int_t 					element_,
					       wmesh_int_t 					shape_family_,
					       wmesh_int_t 					shape_degree_);

template
wmesh_status_t wmesh_shape_restrict_def<double>(wmesh_shape_restrict_t<double>*__restrict__ 	self_,
						wmesh_int_t 					element_,
						wmesh_int_t 					shape_family_,
						wmesh_int_t 					shape_degree_);
