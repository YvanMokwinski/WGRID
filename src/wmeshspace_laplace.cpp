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

extern "C"
{
  
  wmesh_status_t wmeshspace_laplace	(const wmeshspace_t*__restrict__ 	self_,
					 wmesh_int_t 				csr_size_,
					 const_wmesh_int_p			csr_ptr_,
					 const_wmesh_int_p			csr_ind_,
					 double * 				csr_val_,
					 double * 				rhs_)
  {
    WMESH_CHECK_POINTER(self_);
    WMESH_CHECK_POINTER(csr_size_);
    WMESH_CHECK_POINTER(csr_ptr_);
    WMESH_CHECK_POINTER(csr_ind_);
    WMESH_CHECK_POINTER(csr_val_);
    WMESH_CHECK_POINTER(rhs_);


#if 0
    wmesh_t * 			m_mesh;
    wmesh_int_t 		m_degree;
    wmesh_t * 			m_patterns[4];

    wmesh_int_sparsemat_t 	m_c2d;
    wmesh_int_sparsemat_t 	m_c2d_n;
    wmesh_int_sparsemat_t 	m_c2d_e;
    wmesh_int_sparsemat_t 	m_c2d_t;
    wmesh_int_sparsemat_t 	m_c2d_q;
    wmesh_int_sparsemat_t 	m_c2d_i;
    wmesh_int_t 		m_ndofs;
#endif

    wmesh_int_t shape_family = WMESH_SHAPE_FAMILY_LAGRANGE;
    wmesh_int_t topodim = self_->m_mesh->m_topology_dimension;
    wmesh_int_t degree = self_->m_degree;    
    wmesh_status_t status;
    wmesh_int_t ntypes;

    wmesh_int_t cubature_family = WMESH_CUBATURE_FAMILY_GAUSSLEGENDRE;
    wmesh_int_t cubature_degree = 2 * degree;


    wmesh_mat_t<double> cubature_coordinates[4];
    wmesh_mat_t<double> cubature_weights[4];

    wmesh_mat_t<double> shape_f[4];
    wmesh_mat_t<double> shape_f_diff[4][3];

    wmesh_mat_t<double> shape_element[4];
    wmesh_mat_t<double> shape_element_diff[4][3];

    
    wmesh_int_t elements[4];
    status =  bms_topodim2elements(topodim,
				   &ntypes,
				   elements);
    WMESH_STATUS_CHECK(status);


    //
    // Define cubatures.
    //
    wmesh_int_t n1d;
    wmesh_int_t cubature_num_nodes;

    status = bms_cubature_num_nodes(WMESH_ELEMENT_EDGE,
				    cubature_family,
				    cubature_degree,
				    &n1d);
    WMESH_STATUS_CHECK(status);
    
    for (wmesh_int_t i=0;i<ntypes;++i)
      {
	status = bms_cubature_num_nodes(elements[i],
					cubature_family,
					cubature_degree,
					&cubature_num_nodes);
	WMESH_STATUS_CHECK(status);
	std::cout << "define cubature "  << cubature_num_nodes << std::endl;

	wmesh_mat_t<double>::define(&cubature_coordinates[i],
				    topodim,
				    cubature_num_nodes,
				    (double*)malloc(sizeof(double)*topodim*cubature_num_nodes),
				    topodim);
	
	wmesh_mat_t<double>::define(&cubature_weights[i],
				    1,
				    cubature_num_nodes,
				    (double*)malloc(sizeof(double)*1*cubature_num_nodes),
				    1);

	
	wmesh_int_t rw_n;
	status = bms_cubature_buffer_size(elements[i],
					  cubature_family,
					  n1d,			
					  &rw_n);
	WMESH_STATUS_CHECK(status);


	double * rw = (rw_n > 0) ? (double*)malloc(sizeof(double)*rw_n) : nullptr;
	if (rw_n > 0 && !rw)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }
	
	status = bms_cubature(elements[i],
			      cubature_family,
			      n1d,

			      WMESH_STORAGE_INTERLEAVE,
			      cubature_coordinates[i].m,
			      cubature_coordinates[i].n,
			      cubature_coordinates[i].v,
			      cubature_coordinates[i].ld,

			      cubature_weights[i].n,
			      cubature_weights[i].v,
			      cubature_weights[i].ld,
			      
			      rw_n,
			      rw);  
	WMESH_STATUS_CHECK(status);
	if (rw)
	  {
	    free(rw);
	  }
	rw = nullptr;
      }


    wmesh_int_t degree_element = 1;
    for (wmesh_int_t i=0;i<ntypes;++i)
      {
       	
	wmesh_mat_t<double>::define(&shape_element[i],
				    self_->m_mesh->m_c2n.m_m[i],
				    cubature_num_nodes,
				    (double*)malloc(sizeof(double)*self_->m_mesh->m_c2n.m_m[i]*cubature_num_nodes),
				    self_->m_mesh->m_c2n.m_m[i]);

	for (wmesh_int_t j=0;j<topodim;++j)
	  {
	    wmesh_mat_t<double>::define(&shape_element_diff[i][j],
					self_->m_mesh->m_c2n.m_m[i],
					cubature_num_nodes,
					(double*)malloc(sizeof(double)*self_->m_mesh->m_c2n.m_m[i]*cubature_num_nodes),
					self_->m_mesh->m_c2n.m_m[i]);
	  }
	
	wmesh_int_t iw_n,rw_n;
	status = bms_shape_buffer_size(elements[i],
				       WMESH_SHAPE_FAMILY_LAGRANGE,
				       degree_element,
				       &iw_n,
				       &rw_n);
	WMESH_STATUS_CHECK(status);

	wmesh_int_p iw = (iw_n > 0) ? (wmesh_int_p)malloc(sizeof(double)*iw_n) : nullptr;
	if (iw_n > 0 && !iw)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }

	double * rw = (rw_n > 0) ? (double*)malloc(sizeof(double)*rw_n) : nullptr;
	if (rw_n > 0 && !rw)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }

	wmesh_int_t diff[3] = {0,0,0};
	status = bms_dshape(elements[i],
			    shape_family,
			    degree_element,
			    
			    diff,
			    
			    WMESH_STORAGE_INTERLEAVE,
			    cubature_coordinates[i].m,
			    cubature_coordinates[i].n,
			    cubature_coordinates[i].v,
			    cubature_coordinates[i].ld,
			    
			    WMESH_STORAGE_INTERLEAVE,			    
			    shape_element[i].m,
			    shape_element[i].n,
			    shape_element[i].v,
			    shape_element[i].ld,
			    
			    iw_n,
			    iw,
			    rw_n,
			    rw);
	
	WMESH_STATUS_CHECK(status);
	
	for (wmesh_int_t j=0;j<topodim;++j)
	{
	  diff[j] = 1;
	  status = bms_dshape(elements[i],
			      shape_family,
			      degree_element,
			    
			      diff,
			    
			      WMESH_STORAGE_INTERLEAVE,
			      cubature_coordinates[i].m,
			      cubature_coordinates[i].n,
			      cubature_coordinates[i].v,
			      cubature_coordinates[i].ld,
			    
			      WMESH_STORAGE_INTERLEAVE,			    
			      shape_element_diff[i][j].m,
			      shape_element_diff[i][j].n,
			      shape_element_diff[i][j].v,
			      shape_element_diff[i][j].ld,
			      
			      iw_n,
			      iw,
			      rw_n,
			      rw);

	  WMESH_STATUS_CHECK(status);
	  diff[j] = 0;
	}
	
	if (rw)
	  {
	    free(rw);
	  }
	rw = nullptr;
	if (iw)
	  {
	    free(iw);
	  }
	iw = nullptr;
      }




    
    for (wmesh_int_t i=0;i<ntypes;++i)
      {
	wmesh_int_t ndofs_per_element = self_->m_c2d.m_m[i];
	wmesh_mat_t<double>::define(&shape_f[i],
				    ndofs_per_element,
				    cubature_num_nodes,
				    (double*)malloc(sizeof(double)*ndofs_per_element*cubature_num_nodes),
				    ndofs_per_element);

	for (wmesh_int_t j=0;j<topodim;++j)
	  {
	    wmesh_mat_t<double>::define(&shape_f_diff[i][j],
					ndofs_per_element,
					cubature_num_nodes,
					(double*)malloc(sizeof(double)*ndofs_per_element*cubature_num_nodes),
					ndofs_per_element);
	  }
	
	wmesh_int_t iw_n,rw_n;
	status = bms_shape_buffer_size(elements[i],
				       shape_family,
				       degree,
				       &iw_n,
				       &rw_n);
	WMESH_STATUS_CHECK(status);

	wmesh_int_p iw = (iw_n > 0) ? (wmesh_int_p)malloc(sizeof(double)*iw_n) : nullptr;
	if (iw_n > 0 && !iw)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }

	double * rw = (rw_n > 0) ? (double*)malloc(sizeof(double)*rw_n) : nullptr;
	if (rw_n > 0 && !rw)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }

	wmesh_int_t diff[3] = {0,0,0};
	status = bms_dshape(elements[i],
			    shape_family,
			    degree,
			    
			    diff,
			    
			    WMESH_STORAGE_INTERLEAVE,
			    cubature_coordinates[i].m,
			    cubature_coordinates[i].n,
			    cubature_coordinates[i].v,
			    cubature_coordinates[i].ld,
			    
			    WMESH_STORAGE_INTERLEAVE,			    
			    shape_f[i].m,
			    shape_f[i].n,
			    shape_f[i].v,
			    shape_f[i].ld,
			    
			    iw_n,
			    iw,
			    rw_n,
			    rw);
	
	WMESH_STATUS_CHECK(status);
	
	for (wmesh_int_t j=0;j<topodim;++j)
	{
	  diff[j] = 1;
	  status = bms_dshape(elements[i],
			      shape_family,
			      degree,
			    
			      diff,
			    
			      WMESH_STORAGE_INTERLEAVE,
			      cubature_coordinates[i].m,
			      cubature_coordinates[i].n,
			      cubature_coordinates[i].v,
			      cubature_coordinates[i].ld,
			    
			      WMESH_STORAGE_INTERLEAVE,			    
			      shape_f_diff[i][j].m,
			      shape_f_diff[i][j].n,
			      shape_f_diff[i][j].v,
			      shape_f_diff[i][j].ld,
			      
			      iw_n,
			      iw,
			      rw_n,
			      rw);

	  WMESH_STATUS_CHECK(status);
	  diff[j] = 0;
	}
	
	if (rw)
	  {
	    free(rw);
	  }
	rw = nullptr;
	if (iw)
	  {
	    free(iw);
	  }
	iw = nullptr;
      }

    //
    // Now loop over each zone.
    //
    double cooelm[32];
    for (wmesh_int_t itype=0;itype<ntypes;++itype)
      {
	if (self_->m_c2d.m_n[itype]>0)
	  {
	    for (wmesh_int_t ielm=0;ielm<self_->m_c2d.m_n[itype];++ielm)
	      {

		//
		// Get coordinates of the elements.
		//
		wmesh_int_t num_nodes_per_elem = self_->m_mesh->m_c2n.m_m[itype];
		for (wmesh_int_t i=0;i<num_nodes_per_elem;++i)
		  {
		    for (wmesh_int_t k=0;k<topodim;++k)
		      {			
			cooelm[i*topodim + k] =
			  self_->m_mesh->m_coo[(self_->m_mesh->m_c2n.m_data[self_->m_mesh->m_c2n.m_ptr[itype] +self_->m_mesh->m_c2n.m_ld[itype]    *ielm+i]-1)*self_->m_mesh->m_coo_ld + k];
		      }
		  }

		
		//
		// Evaluate the inverse of the Jacobian for each gauss points.
		//
		for (wmesh_int_t k=0;k<cubature_coordinates[itype].n;++k)
		  {
		    double ref_nabla_element[9];
		    //
		    // cooelm 3xNumNodesPerElm
		    //
		    
		    for (wmesh_int_t l=0;l<topodim;++l)
		      {
			double r1 = 1.0;
			double r0 = 0.0;
			BLAS_dgemm("N",
				   "N",
				   &topodim,
				   &cubature_coordinates[itype].n,
				   &num_nodes_per_elem,
				   &r1,
				   cooelm,
				   &topodim,
				   shape_element_diff[itype][l].v,
				   &shape_element_diff[itype][l].ld,
				   &r0,
				   &ref_nabla_element[topodim*l],
				   &topodim);
		      }
		    
		  }
		
		
		//
		// Evaluate the absolute determinant of the Jacobian.
		//

		//
		// Evaluate the global gradient.
		//
#if 0
		wmesh_int_t num_dofs_per_elem = self_->m_c2d.m_m[itype];
		for (wmesh_int_t j=0;j<num_dofs_per_elem;++j)
		  {
		    for (wmesh_int_t i=0;i<num_dofs_per_elem;++i)
		      {			
			for (wmesh_int_t k=0;k<cubature_num_nodes;++k)
			  {
			    dger(cubature_weights[itype].n,
				 cubature_weights[itype].v[k],
				 shape_diff_g[itype][k].v + shape_diff_g[itype][k].ld * i,
				 shape_diff_g[itype][k].v + shape_diff_g[itype][k].ld * j);
			  }
		      }
		  }
#endif
		
	      }	      
	  }
      }
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
