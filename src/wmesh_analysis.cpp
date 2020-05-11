
#include <stdlib.h>
#include <string.h>
#include "wmesh-types.hpp"
#include "wmesh-status.h"
#include "wmesh_medit.hpp"
#include "bms.h"
#include "wmesh.hpp"
#include <chrono>
#include <iostream>

using namespace std::chrono;
extern "C"
{

  wmesh_status_t wmesh_shape_eval_basis(wmesh_int_t     cell_type_,
					wmesh_int_t 	coo_m_,
					wmesh_int_t 	coo_n_,
					double * 	coo_x_,
					wmesh_int_t 	coo_ld_,
					double * 	x_,
					wmesh_int_t 	x_ld_)
  {
    
    static const double s_zero = ((double)0.0);
    static const double s_one = ((double)1.0);
    static constexpr double one = 1.0;
    if (cell_type_==0)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	    
	    x_[k*x_ld_+0] = s_one - (r+s+t);
	    x_[k*x_ld_+1] = r;
	    x_[k*x_ld_+2] = s;
	    x_[k*x_ld_+3] = t;
	  }		  
      }
    else if (cell_type_==1)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	    
	    const auto
	      tscale = (t < s_one) ? s_one / (s_one - t) : s_zero,
	      lts = s_one - t - s,
	      ltr = s_one - t - r;
	    
	    x_[k*x_ld_+0] = ltr  * lts * tscale;
	    x_[k*x_ld_+1] = r * lts * tscale;
	    x_[k*x_ld_+2] = r * s * tscale;
	    x_[k*x_ld_+3] = ltr * s  * tscale;
	    x_[k*x_ld_+4] = t;
	  }		  
      }
    else if (cell_type_==2)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];
	   
	    const auto
	      l = one - (r+s),
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = l * lt;
	    x_[k*x_ld_+1] = r * lt;
	    x_[k*x_ld_+2] = s * lt;
	    x_[k*x_ld_+3] = l * t;
	    x_[k*x_ld_+4] = r * t;
	    x_[k*x_ld_+5] = s * t;
	  }		  
      }
    else if (cell_type_==3)
      {
	for (wmesh_int_t k=0;k<coo_n_;++k)
	  {
	    const auto
	      r = coo_x_[coo_ld_*k+0],
	      s = coo_x_[coo_ld_*k+1],
	      t = coo_x_[coo_ld_*k+2];

	    const auto
	      lr = one - r,
	      ls = one - s,
	      lt = one - t;
	    
	    x_[k*x_ld_+0] = lr* ls * lt;
	    x_[k*x_ld_+1] = r * ls * lt;
	    x_[k*x_ld_+2] = r *  s * lt;
	    x_[k*x_ld_+3] = lr*  s * lt;
	    x_[k*x_ld_+4] = lr* ls * t;
	    x_[k*x_ld_+5] =  r* ls * t;
	    x_[k*x_ld_+6] =  r*  s * t;
	    x_[k*x_ld_+7] = lr*  s * t;
	  }
      }
    else
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_INVALID_ARGUMENT);
      }
    return WMESH_STATUS_SUCCESS;
  }

  




  
  
  
  wmesh_status_t wmesh_analysis_edges(wmesh_t* 		self_,
				      wmesh_int_p 	num_edges_)
  {
    wmesh_int_t edge_idx = 0;
    
    wmesh_int_t work_n;
    wmesh_status_t status =  wmesh_indexing_edges_buffer_size(4,
							      self_->m_c2n.m_n,
							      &work_n);      
    WMESH_STATUS_CHECK(status);
    wmesh_int_t * work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      } 
    status =  wmesh_indexing_edges(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
				   WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2e),
				   WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_e2n),				     
				   &edge_idx,
				   work_n,
				   work);
    free(work);
    WMESH_STATUS_CHECK(status);
    num_edges_[0] = edge_idx;
    return WMESH_STATUS_SUCCESS;
  }
  

  wmesh_status_t wmesh_analysis_faces(wmesh_t* 		self_)
  {
    wmesh_status_t status;
    
    wmesh_int_t work_n;
    status =  wmesh_indexing_triangles_buffer_size(self_->m_c2n.m_size,
						   self_->m_c2n.m_n,
						   &work_n);
    WMESH_STATUS_CHECK(status);
    wmesh_int_p work = (wmesh_int_p)malloc(work_n*sizeof(wmesh_int_t));
    if (!work)
      {
	WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
      }      
    
    //
    // Enumerates triangles faces.
    //
    wmesh_int_t face_idx = 0;
    status =  wmesh_indexing_triangles(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
				       WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2f_t),
				       WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_t2n),				     
				       &face_idx,
				       work_n,
				       work);      
    
    WMESH_STATUS_CHECK(status);      
    face_idx = 0;
    status = wmesh_indexing_quadrilaterals(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
					   WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2f_q),
					   WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_q2n),				     
					   &face_idx,
					   work_n,
					   work);      
    WMESH_STATUS_CHECK(status);
    return WMESH_STATUS_SUCCESS;
  }

  
  wmesh_status_t wmesh_c2c(wmesh_t* 		self_,
			   wmesh_int_p 		num_hyperfaces_,
			   wmesh_int_p 		num_boundary_hyperfaces_)
  {
    const wmesh_int_t topodim = self_->m_topology_dimension;
    //
    // Initialize c2c.
    //
    wmesh_status_t status = wmesh_init_c2c	(topodim,
						 &self_->m_c2n,
						 &self_->m_c2c);
    WMESH_STATUS_CHECK(status);	

    if (topodim == 3)
      {
	//
	// Initialize sub c2c for triangles.
	//
	status = wmesh_init_c2f_t(&self_->m_c2c,
				  &self_->m_c2c_t);
	WMESH_STATUS_CHECK(status);
	
	//
	// Initialize sub c2c for quadrilaterals.
	//
	status = wmesh_init_c2f_q(&self_->m_c2c,
				  &self_->m_c2c_q);
	WMESH_STATUS_CHECK(status);      
      }

    wmesh_int_t work_n;
    wmesh_int_t * work  = nullptr;

    if (self_->m_topology_dimension==3)
      {	
	status =  bms_c2c_buffer_size(self_->m_c2n.m_size,
				      self_->m_c2n.m_n,
				      &work_n);    
	WMESH_STATUS_CHECK(status);    
	work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
	if (!work)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }      
	
	num_boundary_hyperfaces_[0] = 0;
	num_boundary_hyperfaces_[1] = 0;
	num_hyperfaces_[0] 		= 0;
	num_hyperfaces_[1] 		= 0;
	
	status = bms_c2c(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
			 WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2c_t),
			 WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2c_q),
			 WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_t2n),
			 WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_q2n),
			 work_n,
			 work,
			 num_hyperfaces_,
			 num_boundary_hyperfaces_);
	
	free(work);      
	WMESH_STATUS_CHECK(status);    
	return WMESH_STATUS_SUCCESS;
      }
    else
      {
	status =  bms_c2c_e_buffer_size(self_->m_c2n.m_size,
					self_->m_c2n.m_n,
					&work_n);    
	WMESH_STATUS_CHECK(status);    
	work = (wmesh_int_t*)malloc(work_n*sizeof(wmesh_int_t));
	if (!work)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_ERROR_MEMORY);
	  }      
	
	num_boundary_hyperfaces_[0] = 0;
		
	status = bms_c2c_e(WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2n),
			   WMESH_INT_SPARSEMAT_FORWARD(self_->m_c2c),
			   WMESH_INT_SPARSEMAT_FORWARD(self_->m_s_e2n),
			   num_boundary_hyperfaces_,
			   work_n,
			   work);
	
	num_hyperfaces_[0] = (num_boundary_hyperfaces_[0] + self_->m_c2n.m_n[0]*3+self_->m_c2n.m_n[1]*4)/2;
	
	free(work);      
	WMESH_STATUS_CHECK(status);    
	return WMESH_STATUS_SUCCESS;

      }
   
  }

  
  using timing_t = high_resolution_clock::time_point;
  inline timing_t timing_stop(){ return high_resolution_clock::now(); }
  inline timing_t timing_start(){ return high_resolution_clock::now(); }
  inline double timing_seconds(timing_t&t,timing_t&t1){ return duration_cast<duration<double>>(t1-t).count(); }


  
  wmesh_status_t wmesh_analysis(wmesh_t*    self_)
  {
    const wmesh_int_t topology_dimension = self_->m_topology_dimension;
    wmesh_status_t 	status;
    wmesh_int_t 	num_facets[2];
    wmesh_int_t 	num_bfacets[2];
    
#ifndef NDEBUG
    std::cerr << "wmesh_analysis: topology_dimension " << topology_dimension << std::endl;
#endif

    //
    // 3/ Create the cells to edges.
    //    
    {
      wmesh_status_t status;
      status = wmesh_init_c2e(&self_->m_c2n,
			      &self_->m_c2e,
			      self_->m_topology_dimension);
      WMESH_STATUS_CHECK(status);
    }
    
    //
    // 4/ Create the cells to faces.
    //
    if (topology_dimension == 3)
    {
      wmesh_status_t status;
      status = wmesh_init_c2f(&self_->m_c2n,
			      &self_->m_c2f);
      WMESH_STATUS_CHECK(status);
      
      //
      // 4.a/ Define the part of triangles.
      //      
      status = wmesh_init_c2f_t	(&self_->m_c2f,
				 &self_->m_c2f_t);
      WMESH_STATUS_CHECK(status);
      
      //
      // 4.b/ Define the part of quadrilaterals.
      //      
      status = wmesh_init_c2f_q	(&self_->m_c2f,
				 &self_->m_c2f_q);
      WMESH_STATUS_CHECK(status);      
    }
    
    //
    // Initialize the reference shape edges to nodes.
    //


    

      auto start = timing_start();      
	status    = wmesh_c2c(self_,
			      num_facets,
			      num_bfacets);

	if (self_->m_topology_dimension==3)
	  {
	    self_->m_num_triangles      = num_facets[0];
	    self_->m_num_quadrilaterals = num_facets[1];
	    self_->m_num_bfaces[0]      = num_bfacets[0];
	    self_->m_num_bfaces[1]      = num_bfacets[1];
	  }
	else if (self_->m_topology_dimension==2)
	  {
	    self_->m_num_edges      = num_facets[0];
	    self_->m_num_bfaces[0]  = num_bfacets[0];
	  }
	else if (self_->m_topology_dimension==1)
	  {
	    WMESH_STATUS_CHECK(WMESH_STATUS_NOT_IMPLEMENTED);
	  }
	
	auto stop = timing_stop();
	std::cout << "neighbors detection, elapsed " << timing_seconds(start,stop) << std::endl;

	//      WMESH_STATUS_CHECK(status);
	
	if (self_->m_topology_dimension==3)
	  {
	    start = timing_start();      
	    status = wmesh_analysis_edges(self_, &self_->m_num_edges);
	    stop = timing_stop();
	    std::cout << "edge analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	    WMESH_STATUS_CHECK(status);

	    start = timing_start();      
	    status    = wmesh_analysis_faces(self_);
	    stop = timing_stop();
	    std::cout << "faces analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	    WMESH_STATUS_CHECK(status);
	  }	
	else if (self_->m_topology_dimension==2)
	  {
	    start = timing_start();      
	    status = wmesh_analysis_edges(self_, &self_->m_num_edges);
	    stop = timing_stop();
	    std::cout << "edge analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	    WMESH_STATUS_CHECK(status);

#if 0
	    start = timing_start();      
	    status    = wmesh_analysis_faces(self_);
	    stop = timing_stop();
	    std::cout << "faces analysis, elapsed " << timing_seconds(start,stop) << std::endl;
	    WMESH_STATUS_CHECK(status);
#endif
	  }
	

	WMESH_STATUS_CHECK(status);
    
    if (self_->m_topology_dimension==3)
      {

      std::cout << "boundary faces (triangles, quadrilaterals):  "
		<< num_bfacets[0]
		<< ", "
		<< num_bfacets[1]
		<< std::endl;
    
      std::cout << "total num faces (triangles, quadrilaterals):  "
		<< num_facets[0]
		<< ", "
		<< num_facets[1] << std::endl;
    
      std::cout << "total num edges:  " << self_->m_num_edges  << std::endl;
    }
    else if (self_->m_topology_dimension==2)
      {

      std::cout << "boundary edges:  "
		<< num_bfacets[0]
		<< std::endl;
    
      std::cout << "total num edges:  "
		<< num_facets[0]
		<< std::endl;
    
      std::cout << "total num edges:  " << self_->m_num_edges  << std::endl;
    }
    return WMESH_STATUS_SUCCESS;
  }
  


};
