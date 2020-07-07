#if 0
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
#include "bms.hpp"

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
wmesh_status_t wmesh_cubature_boundary_def(wmesh_cubature_boundary_t<T>*__restrict__ 	self_,
					   wmesh_int_t 				element_,
					   wmesh_int_t 				cubature_family_,
					   wmesh_int_t 				cubature_degree_)
{

  memset(self_,0,sizeof(wmesh_cubature_boundary_t<T>));

  self_->m_element 		= element_;
  self_->m_cubature_family 	= cubature_family_;
  self_->m_cubature_degree 	= cubature_degree_;
  wmesh_status_t status;

  wmesh_int_t topodim;
  status = bms_element2topodim(element_,
			       &topodim);
  WMESH_STATUS_CHECK(status);

  //
  // Get the facets.
  //
  wmesh_int_t num_facets;
  status = bms_element_facets(element_,			 
			      &num_facets,
			      self_->m_facets);
  WMESH_STATUS_CHECK(status);
  self_->m_num_facets = num_facets;
  wmesh_int_t topodim_boundary 	= topodim - 1;
  wmesh_int_p facets 		=  &self_->m_facets[0];
  
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
      wmesh_int_t facet = facets[ifacet];
      wmesh_int_t facet_type = (topodim==3) ? (facet - 2) : 0;
      wmesh_int_t cubature_num_nodes = ref_facets_cubature[facet_type].m_c.n;
      wmesh_int_t facet_num_nodes = (topodim==3) ? (facet_type + 3) : ( (topodim==2) ? 2 : 1 );// be careful, trick.
      for (wmesh_int_t j=0;j<facet_num_nodes * 2;++j)
	{
	  const wmesh_int_t signed_rotation = (j < facet_num_nodes) ?  (j+1) : (facet_num_nodes - j - 1);
	  
	  self_->m_facets_cubature[ifacet][j].m_element   = facet;
	  self_->m_facets_cubature[ifacet][j].m_family    = cubature_family_;
	  self_->m_facets_cubature[ifacet][j].m_degree    = cubature_degree_;
	  self_->m_facets_cubature[ifacet][j].m_c_storage = WMESH_STORAGE_INTERLEAVE;	      
	  
	  wmesh_mat_t<T>::alloc(&self_->m_facets_cubature[ifacet][j].m_c,
				topodim,
				cubature_num_nodes);
	  
	  memset(self_->m_facets_cubature[ifacet][j].m_c.v,0,
		 sizeof(T)*topodim*
				cubature_num_nodes);
	  
	  wmesh_mat_t<T>::alloc(&self_->m_facets_cubature[ifacet][j].m_w,
				1,
				cubature_num_nodes);
#if 0
	  std::cout << "ref facet  " << signed_rotation << std::endl;
	  std::cout << ref_facets_cubature[facet_type].m_c << std::endl;
	  std::cout << "signed  " << signed_rotation << std::endl;
#endif
	  status =  bms_mirrored_local_coordinates(facet,
						   signed_rotation,
						   ref_facets_cubature[facet_type].m_c_storage,
						   WMESH_MAT_FORWARD(ref_facets_cubature[facet_type].m_c),
						   self_->m_facets_cubature[ifacet][j].m_c_storage,
						   WMESH_MAT_FORWARD(self_->m_facets_cubature[ifacet][j].m_c) );
	  WMESH_STATUS_CHECK(status);
#if 0
	  std::cout << "ref cell " << std::endl;
	  std::cout << self_->m_facets_cubature[ifacet][j].m_c << std::endl;
#endif
	  //
	  // memcpy the weights.
	  //
	  memcpy(self_->m_facets_cubature[ifacet][j].m_w.v,
		 ref_facets_cubature[facet_type].m_w.v,
		 sizeof(T)*cubature_num_nodes);
	  T fd = 0;
	  for (wmesh_int_t i=0;i<ref_facets_cubature[facet_type].m_w.n;++i)
	    {
	      fd += ref_facets_cubature[facet_type].m_w.v[i];
	    }
	  std::cout << "sum weights " << fd << std::endl;
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
  
  return WMESH_STATUS_SUCCESS;
}


template
wmesh_status_t wmesh_cubature_boundary_def<float>(wmesh_cubature_boundary_t<float>*__restrict__ 	self_,
						  wmesh_int_t 						element_,
						  wmesh_int_t 						cubature_family_,
						  wmesh_int_t 						cubature_degree_);


template
wmesh_status_t wmesh_cubature_boundary_def<double>(wmesh_cubature_boundary_t<double>*__restrict__ 	self_,
						   wmesh_int_t 						element_,
						   wmesh_int_t 						cubature_family_,
						   wmesh_int_t 						cubature_degree_);
#endif
