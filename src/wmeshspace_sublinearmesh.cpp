//
// treat 2d
//
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-math.hpp"
#include "wmesh-status.h"
#include "wmesh_t.hpp"
#include "wmesh_utils.hpp"
#include "wmesh-blas.h"
#include "bms.h"

#include <chrono>
#include <iostream>
#include <math.h>
#include "bms_templates.hpp"
#include "wmeshspace_t.hpp"

using namespace std::chrono;
#if 0
template<typename T>
wmesh_status_t wmeshspace_generate_element_coodofs(const wmeshspace_t * __restrict__ 	self_,
						   wmesh_int_t 				element_type_,
						   wmesh_int_t 				element_idx_,
						   wmesh_int_t 				cooelm_storage_,
						   const wmesh_mat_t<T>& 		cooelm_,
						   wmesh_int_t 				coodofs_storage_,
						   wmesh_mat_t<T>& 			coodofs_)
{
  static constexpr T r0 = static_cast<T>(0);
  static constexpr T r1 = static_cast<T>(1);
  
  wmesh_mat_gemm(shape_eval_element->m_storage,
		 cooelm_storage_,
		 coodofs_storage_,
		 r1,
		 shape_eval_element->m_f,
		 cooelm_,
		 r0,
		 coodofs_);

  
  wmesh_status_t 	status;
  const wmesh_t * 	mesh 	= self_->get_mesh();
  wmesh_int_t         topodim = mesh->m_topology_dimension;
  wmesh_int_t 	coo_m  	= mesh->m_coo_m;
  
  wmesh_int_t 	num_types;
  wmesh_int_t 	elements[4];
  T 		cell_xyz[32];
  wmesh_int_t 	cell_xyz_ld = coo_m;
  wmesh_mat_t<T> eval_basis[4];

  status = bms_topodim2elements(topodim,
				&num_types,
				elements);
  WMESH_STATUS_CHECK(status);        
  for (wmesh_int_t l=0;l<num_types;++l)
    {
      if (mesh->m_c2n.m_n[l]>0)
	{
	  wmesh_int_t element = (topodim==3) ? (4+l) : ((topodim==2)? (2+l) : 1+l);
	  wmesh_int_t element_num_nodes;

	  status = bms_elements_num_nodes(1,&element,&element_num_nodes);
	  WMESH_STATUS_CHECK(status);        
	    
	  const wmesh_t* 	rmacro			= self_->get_refinement_pattern(l);
	  const T * 		rmacro_coo 		= rmacro->m_coo;
	  wmesh_int_t 		rmacro_coo_ld 		= rmacro->m_coo_ld;
	  wmesh_int_t 		rmacro_num_nodes 	= rmacro->m_num_nodes;
	  const wmesh_int_t  	mat_rmacro_coo_storage = WMESH_STORAGE_INTERLEAVE;
	  wmesh_mat_t<T> mat_rmacro_coo;
	  wmesh_mat_t<T>::define(&mat_rmacro_coo,topodim,rmacro_num_nodes,(T*)rmacro_coo,rmacro_coo_ld);

	  const wmesh_int_t  	eval_basis_storage = WMESH_STORAGE_INTERLEAVE;
	  wmesh_mat_t<T>::alloc(&eval_basis[l], element_num_nodes, rmacro_num_nodes);

	  wmesh_shape_t shape_element(element, WMESH_SHAPE_FAMILY_LAGRANGE, 1);
	    
	  wmesh_shape_calculate_eval(shape_element,
				     mat_rmacro_coo_storage,
				     mat_rmacro_coo,
				     eval_basis_storage,
				     eval_basis[l]);
	}

    }
    
  wmesh_int_t rw_n = 0;
  for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
    {
      wmesh_int_t k = self_->m_c2d.m_m[l]*topodim;
      rw_n = (rw_n < k) ? k : rw_n;
    }
  T * rw = (T*)malloc(sizeof(T)*rw_n);
  for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
    {
      //	double * 	refeval = refevals[l];
	
      //
      // Local c2d.
      //	    
      // wmesh_int_t c2d_n	= self_->m_c2d.m_n[l];
      wmesh_int_t c2d_m 		= self_->m_c2d.m_m[l];
      wmesh_int_t c2d_ld 		= self_->m_c2d.m_ld[l];
      wmesh_int_p c2d_v 		= self_->m_c2d.m_data + self_->m_c2d.m_ptr[l];
	
      //
      // Local c2n.
      //
      wmesh_int_t c2n_n 		= mesh->m_c2n.m_n[l];
      wmesh_int_t c2n_m 		= mesh->m_c2n.m_m[l];
      wmesh_int_t c2n_ld		= mesh->m_c2n.m_ld[l];
      wmesh_int_p c2n_v 		= mesh->m_c2n.m_data + mesh->m_c2n.m_ptr[l];
	
      for (wmesh_int_t j=0;j<c2n_n;++j)
	{
	    
	  //
	  // Get the coordinates of the cell.
	  //		
	  for (wmesh_int_t i=0;i<c2n_m;++i)
	    {
	      wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
	      for (wmesh_int_t k=0;k<mesh->m_coo_m;++k)
		{
		  cell_xyz[cell_xyz_ld * i + k] = mesh->m_coo[mesh->m_coo_ld * idx + k];
		}
	    }

	  //	    wmesh_mat_gemm(static_cast<T>(1),cell_xyz,eval_basis[l],static_cast<T>(0),physical_coordinates);
	    
	  xgemm("N",
		"N",
		&coo_m ,
		&c2d_m,
		&c2n_m ,
		&r1,
		cell_xyz,
		&cell_xyz_ld,
		eval_basis[l].v,
		&eval_basis[l].ld,
		&r0,
		rw,
		&coo_m);

	  //
	  // Copy.
	  //
	  if (WMESH_STORAGE_INTERLEAVE == coo_storage_)
	    {
	      for (wmesh_int_t i=0;i<c2d_m;++i)
		{
		  wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		  for (wmesh_int_t k = 0;k<coo_m;++k)
		    {
		      coo_[coo_ld_ * idx + k] = rw[coo_m * i + k];
		    }
		}
	    }
	  else
	    {
	      for (wmesh_int_t i=0;i<c2d_m;++i)
		{
		  wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		  for (wmesh_int_t k = 0;k<coo_m;++k)
		    {
		      coo_[coo_ld_ * k + idx] = rw[coo_m * i + k];
		    }
		}
	    }	    
	}
    }
  return WMESH_STATUS_SUCCESS;
};

#endif
template<typename T>
wmesh_status_t wmeshspace_generate_coodofs(const wmeshspace_t * __restrict__ 	self_,
					   wmesh_int_t 				coo_storage_,
					   wmesh_int_t 				coo_m_,
					   wmesh_int_t 				coo_n_,
					   T *  __restrict__			coo_,
					   wmesh_int_t 				coo_ld_)
{    
  static constexpr T r0 = static_cast<T>(0);
  static constexpr T r1 = static_cast<T>(1);
  wmesh_status_t 	status;
  const wmesh_t * 	mesh 	= self_->get_mesh();
  wmesh_int_t         topodim = mesh->m_topology_dimension;
  wmesh_int_t 	coo_m  	= mesh->m_coo_m;
  
    wmesh_int_t 	num_types;
    wmesh_int_t 	elements[4];
    T 		cell_xyz[32];
    wmesh_int_t 	cell_xyz_ld = coo_m;
    wmesh_mat_t<T> eval_basis[4];

    status = bms_topodim2elements(topodim,
				  &num_types,
				  elements);
    WMESH_STATUS_CHECK(status);        
    for (wmesh_int_t l=0;l<num_types;++l)
      {
	if (mesh->m_c2n.m_n[l]>0)
	  {
	    wmesh_int_t element = (topodim==3) ? (4+l) : ((topodim==2)? (2+l) : 1+l);
	    wmesh_int_t element_num_nodes;

	    status = bms_elements_num_nodes(1,&element,&element_num_nodes);
	    WMESH_STATUS_CHECK(status);        
	    
	    const wmesh_t* 	rmacro			= self_->get_refinement_pattern(l);
	    const T * 		rmacro_coo 		= rmacro->m_coo;
	    wmesh_int_t 	rmacro_coo_ld 		= rmacro->m_coo_ld;
	    wmesh_int_t 	rmacro_num_nodes 	= rmacro->m_num_nodes;
	    const wmesh_int_t  	mat_rmacro_coo_storage = WMESH_STORAGE_INTERLEAVE;
	    wmesh_mat_t<T> mat_rmacro_coo;
	    wmesh_mat_t<T>::define(&mat_rmacro_coo,topodim,rmacro_num_nodes,(T*)rmacro_coo,rmacro_coo_ld);

	    const wmesh_int_t  	eval_basis_storage = WMESH_STORAGE_INTERLEAVE;
	    wmesh_mat_t<T>::alloc(&eval_basis[l], element_num_nodes, rmacro_num_nodes);

	    wmesh_shape_t shape_element(element, WMESH_SHAPE_FAMILY_LAGRANGE, 1);
	    
	    wmesh_shape_calculate_eval(shape_element,
				       mat_rmacro_coo_storage,
				       mat_rmacro_coo,
				       eval_basis_storage,
				       eval_basis[l]);
	  }

      }
    
    wmesh_int_t rw_n = 0;
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
      {
	wmesh_int_t k = self_->m_c2d.m_m[l]*topodim;
	rw_n = (rw_n < k) ? k : rw_n;
      }
    T * rw = (T*)malloc(sizeof(T)*rw_n);
    for (wmesh_int_t l=0;l<self_->m_c2d.m_size;++l)
      {
	//	double * 	refeval = refevals[l];
	
	//
	// Local c2d.
	//	    
	// wmesh_int_t c2d_n	= self_->m_c2d.m_n[l];
	wmesh_int_t c2d_m 		= self_->m_c2d.m_m[l];
	wmesh_int_t c2d_ld 		= self_->m_c2d.m_ld[l];
	wmesh_int_p c2d_v 		= self_->m_c2d.m_data + self_->m_c2d.m_ptr[l];
	
	//
	// Local c2n.
	//
	wmesh_int_t c2n_n 		= mesh->m_c2n.m_n[l];
	wmesh_int_t c2n_m 		= mesh->m_c2n.m_m[l];
	wmesh_int_t c2n_ld		= mesh->m_c2n.m_ld[l];
	wmesh_int_p c2n_v 		= mesh->m_c2n.m_data + mesh->m_c2n.m_ptr[l];
	
	for (wmesh_int_t j=0;j<c2n_n;++j)
	  {
	    
	    //
	    // Get the coordinates of the cell.
	    //		
	    for (wmesh_int_t i=0;i<c2n_m;++i)
	      {
		wmesh_int_t idx = c2n_v[c2n_ld * j + i] - 1;
		for (wmesh_int_t k=0;k<mesh->m_coo_m;++k)
		  {
		    cell_xyz[cell_xyz_ld * i + k] = mesh->m_coo[mesh->m_coo_ld * idx + k];
		  }
	      }

	    //	    wmesh_mat_gemm(static_cast<T>(1),cell_xyz,eval_basis[l],static_cast<T>(0),physical_coordinates);
	    
	    xgemm("N",
		  "N",
		  &coo_m ,
		  &c2d_m,
		  &c2n_m ,
		  &r1,
		  cell_xyz,
		  &cell_xyz_ld,
		  eval_basis[l].v,
		  &eval_basis[l].ld,
		  &r0,
		  rw,
		  &coo_m);

	    //
	    // Copy.
	    //
	    if (WMESH_STORAGE_INTERLEAVE == coo_storage_)
	      {
		for (wmesh_int_t i=0;i<c2d_m;++i)
		  {
		    wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		    for (wmesh_int_t k = 0;k<coo_m;++k)
		      {
			coo_[coo_ld_ * idx + k] = rw[coo_m * i + k];
		      }
		  }
	      }
	    else
	      {
		for (wmesh_int_t i=0;i<c2d_m;++i)
		  {
		    wmesh_int_t idx = c2d_v[c2d_ld * j + i] - 1;
		    for (wmesh_int_t k = 0;k<coo_m;++k)
		      {
			coo_[coo_ld_ * k + idx] = rw[coo_m * i + k];
		      }
		  }
	      }	    
	  }
      }
    return WMESH_STATUS_SUCCESS;
  };



template
wmesh_status_t wmeshspace_generate_coodofs<double>(const wmeshspace_t * __restrict__ 	self_,
						   wmesh_int_t 				coo_storage_,
						   wmesh_int_t 				coo_m_,
						   wmesh_int_t 				coo_n_,
						   double *  __restrict__			coo_,
						   wmesh_int_t 				coo_ld_);
#if 0
template
wmesh_status_t wmeshspace_generate_coodofs<float>(const wmeshspace_t * __restrict__ 	self_,
						  wmesh_int_t 				coo_storage_,
						  wmesh_int_t 				coo_m_,
						  wmesh_int_t 				coo_n_,
						  float *  __restrict__			coo_,
						  wmesh_int_t 				coo_ld_);

  wmesh_status_t wmeshspace_generate_scoodofs(const wmeshspace_t * __restrict__ 	self_,
					      wmesh_int_t 				coo_storage_,
					      wmesh_int_t 				coo_m_,
					      wmesh_int_t 				coo_n_,
					      float *  __restrict__			coo_,
					      wmesh_int_t 				coo_ld_)
  {
    return wmeshspace_generate_coodofs(self_,
				       coo_storage_,
				       coo_m_,
				       coo_n_,
				       coo_,
				       coo_ld_);
  }
#endif

extern "C"
{

  wmesh_status_t bms_djacobip(wmesh_int_t 	alpha_,
			      wmesh_int_t 	beta_,
			      wmesh_int_t 	N_,
			      wmesh_int_t 	x_n_,
			      const double * 	x_,
			      wmesh_int_t  	x_ld_,
			      double * 		y_,
			      wmesh_int_t  	y_ld_,
			      wmesh_int_t 	work_n_,			   
			      double * 		work_)
  
  {
    return bms_jacobip(alpha_,
		       beta_,
		       N_,
		       x_n_,
		       x_,
		       x_ld_,
		       y_,
		       y_ld_,
		       work_n_,			   
		       work_);
  }

  wmesh_status_t wmeshspace_generate_dcoodofs(const wmeshspace_t * __restrict__ 	self_,
					      wmesh_int_t 				coo_storage_,
					      wmesh_int_t 				coo_m_,
					      wmesh_int_t 				coo_n_,
					      double *  __restrict__			coo_,
					      wmesh_int_t 				coo_ld_)
  {
    return wmeshspace_generate_coodofs(self_,
				       coo_storage_,
				       coo_m_,
				       coo_n_,
				       coo_,
				       coo_ld_);
  }

  wmesh_status_t wmeshspace_sublinearmesh(wmeshspace_t * 	self_,
					  wmesh_t ** 		mesh__)
  {    
    wmesh_status_t 	status;
    const wmesh_t * 	mesh 	= self_->get_mesh();
    wmesh_int_t         topodim = mesh->m_topology_dimension;
    wmesh_int_t 	coo_m  	= mesh->m_coo_m;
    
    wmesh_int_t 	num_types;
    wmesh_int_t 	elements[4];


    status = bms_topodim2elements(topodim,
				  &num_types,
				  elements);
    WMESH_STATUS_CHECK(status);    
    
    wmesh_int_t coo_dofs_m  	= coo_m;
    wmesh_int_t coo_dofs_n  	= self_->get_ndofs();
    wmesh_int_t coo_dofs_ld 	= coo_dofs_m;
    double * 	coo_dofs 	= (double*)malloc(sizeof(double) * coo_dofs_n * coo_dofs_ld);
    
    status =  wmeshspace_generate_coodofs(self_,
					  WMESH_STORAGE_INTERLEAVE,
					  coo_dofs_m,
					  coo_dofs_n,
					  coo_dofs,
					  coo_dofs_ld);
    WMESH_STATUS_CHECK(status);

    //
    // GENERATE SUBLINEAR CONNECTIVITY
    //
    {
      
      wmesh_int_t c2n_size = num_types;
      wmesh_int_t c2n_ptr[5];
      wmesh_int_t c2n_m[4]{};
      wmesh_int_t c2n_n[4]{};	
      wmesh_int_t c2n_ld[4]{};
      
      status = bms_elements_num_nodes(num_types,
				      elements,
				      c2n_m);
      WMESH_STATUS_CHECK(status);
      
      for (int i=0;i<num_types;++i)
	{
	  const wmesh_t * pattern = self_->get_refinement_pattern(i);
	  if (pattern!=nullptr)
	    {
	      if (i==1 && topodim == 3)
		{
		  c2n_n[0] += mesh->m_c2n.m_n[i] * pattern->m_c2n.m_n[0];
		  c2n_n[1] += mesh->m_c2n.m_n[i] * pattern->m_c2n.m_n[1];
		}
	      else
		{
		  c2n_n[i] = mesh->m_c2n.m_n[i] * pattern->m_num_cells;
		}
	    }
	}

      status = wmesh_int_sparsemat_init(c2n_size,
					c2n_ptr,
					c2n_m,
					c2n_n,
					c2n_ld);
      WMESH_STATUS_CHECK(status);      
      wmesh_int_p c2n_v = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*c2n_ptr[c2n_size]);
      if (!c2n_v)
	{
	  WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);      
	}

      wmesh_int_t mxdofs = 0;
      for (int i=0;i<num_types;++i)
	{
	  const wmesh_t * pattern = self_->get_refinement_pattern(i);
	  //	  if (c2n_n[i] > 0) mxdofs = (self_->m_patterns[i]->m_num_nodes > mxdofs) ? self_->m_patterns[i]->m_num_nodes : mxdofs;
	  if (nullptr != pattern)
	    {
	      mxdofs = (pattern->m_num_nodes > mxdofs) ? pattern->m_num_nodes : mxdofs;
	    }
	}
      
      wmesh_int_p dofs = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);
      wmesh_int_p lidx = (wmesh_int_p)malloc(sizeof(wmesh_int_t)*mxdofs);
      
      printf("generate connectivity.\n");
      wmesh_int_t idx[4]{0,0,0,0};
      for (wmesh_int_t l=0;l<num_types;++l)
	{
	  auto ref_c2n = &self_->get_refinement_pattern(l)->m_c2n;
	  
	  wmesh_int_t ncells = self_->m_c2d.m_n[l];
	  for (wmesh_int_t j=0;j<ncells;++j)
	    {
	      
	      //
	      // extract dofs.
	      //
	      for (wmesh_int_t i=0;i<self_->m_c2d.m_m[l];++i)
		{
		  dofs[i] = self_->m_c2d.m_data[self_->m_c2d.m_ptr[l] + j * self_->m_c2d.m_ld[l] + i];
		}

		
	      //
	      // For each.
	      //
	      for (wmesh_int_t ref_cell_type=0;ref_cell_type<ref_c2n->m_size;++ref_cell_type)
		{
		  wmesh_int_t ref_ncells = ref_c2n->m_n[ref_cell_type];
		  wmesh_int_t ref_nnodes_per_cell = ref_c2n->m_m[ref_cell_type];
		  for (wmesh_int_t sj=0;sj<ref_ncells;++sj)
		    {
			
		      for (wmesh_int_t si=0;si<ref_nnodes_per_cell;++si)
			{
			  lidx[si] = ref_c2n->m_data[ref_c2n->m_ptr[ref_cell_type] + sj * ref_c2n->m_ld[ref_cell_type]  + si ] - 1;
			}
			
		      for (wmesh_int_t si=0;si<ref_c2n->m_m[ref_cell_type];++si)
			{
			  c2n_v[ c2n_ptr[ref_cell_type] + c2n_ld[ref_cell_type] * idx[ref_cell_type] + si] = dofs[lidx[si]];
			}
		      ++idx[ref_cell_type];		       
		    }
		}
	    }
	}
      
      //
      // Define the mesh.
      //
      status =  wmesh_def(mesh__,
			  mesh->m_topology_dimension,				 
			  c2n_size,
			  c2n_ptr,
			  c2n_m,
			  c2n_n,
			  c2n_v,
			  c2n_ld,
			  coo_dofs_m,
			  coo_dofs_n,
			  coo_dofs,
			  coo_dofs_ld);

      //
      // Copy the dof codes.
      // It's a bit tricky, need to be changed.
      //
      const wmesh_int_t ndofs = self_->get_ndofs();
      for (wmesh_int_t i=0;i<ndofs;++i)
	{
	  mesh__[0]->m_n_c.v[i] = self_->m_dof_codes[i];
	}
      
      WMESH_STATUS_CHECK(status);
    }
    return WMESH_STATUS_SUCCESS;
  }

};




